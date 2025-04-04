# Horizontal Plane Solution Python-Fortran Model
# Author G. Stevens 04/05/25
# Licensed under GPU version 3.0
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import seaborn as sns; sns.set()
from tkinter import*
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import tkinter.font as tkFont
from tkinter import filedialog
from PIL import ImageTk, Image
from fpdf import FPDF
import datetime
import csv
import os
from scipy.interpolate import griddata
import time
import hpsmodel
from osgeo import gdal
root = Tk()

# HPS Function
def HPS(pjt_name,k,grad,por,aq_th,disp_x, disp_y, disp_z, R, Lambda,flow_dir,nsrs,\
        model_t,model_x,model_y,model_z,nint,ctr_cut): 
    
    # store starting time. Used as refernce for time duration (see end of program)
    begin = time.time() 
    # Obtain values from Entry widgets and validate
    try:
      k = float(k.get())
    except ValueError:
        messagebox.showerror("Invalid Input","k must be a number!")
        return False    
    try:
      grad = float(grad.get())
    except ValueError:
         messagebox.showerror("Invalid Input","Gradient must be a number!")
         return False
    try: 
      por = float(por.get()) # Porosity
    except ValueError:
        messagebox.showerror("Invalid Input","Porosity must be a number!")
        return False
    try: 
      aq_th = float(aq_th.get()) # Aquifer Thickness (meters)
    except ValueError:
         messagebox.showerror("Invalid Input","Aquifer Thickness must be a number!") 
         return False
    try:
      disp_x = float(disp_x.get()) # Dispersivity in X-Direction (meters)
    except ValueError:
         messagebox.showerror("Invalid Input","Dispersivity-x must be a number!") 
         return False
    try:
      disp_y = float(disp_y.get())  #Dispersivity in Y-Direction (meters)
    except ValueError:
         messagebox.showerror("Invalid Input","Dispersivity-y must be a number!") 
         return False
    try:
      disp_z = float(disp_z.get()) # Dispersivity in Z-Direction (meters)
    except ValueError:
         messagebox.showerror("Invalid Input","Dispersivity-z must be a number!") 
         return False
    try:
      R = float(R.get())  # Retardation
    except ValueError:
         messagebox.showerror("Invalid Input","Retardation must be a number!") 
         return False
    try:
      Lambda = float(Lambda.get()) # Contaminant Decay rate (1/Years)
    except ValueError:
         messagebox.showerror("Invalid Input","Contaminant decay must be a number!") 
         return False
    try:   
      model_t = float(model_t.get())  # Time of calculation (years)
    except ValueError:
         messagebox.showerror("Invalid Input","Time of calulation must be a number!") 
         return False
    try:
      model_x = int(model_x.get()) + 1 # Number of calulation points in the X-direction
      if isinstance(model_x, int):
          pass
    except ValueError:
         messagebox.showerror("Invalid Input","Model horizontal length"
                             " must be an integer number!") 
         return False
    try:
      model_y = int(model_y.get()) + 1 # Number of calulation points in the Y-direction
      if isinstance(model_y, int):
        pass
    except ValueError:
         messagebox.showerror("Invalid Input","Model horizontal length"
                           " must be an integer number!")
         return False
    try:
      model_z = int(model_z.get()) # z-dimension (depth) of calculation points (meters) 
    except ValueError:
            messagebox.showerror("Invalid Input, Time of calulation must be a number!") 
            return False  
    try:  
       nint = int(nint.get())  # Number of Time Steps for Integration
    except ValueError:
             messagebox.showerror("Invalid Input, Number of time steps must be a number!")                                 
             return False  
    try:
      ctr_cut = float(ctr_cut.get()) # Contour cutoff
    except ValueError:
         messagebox.showerror("Invalid Input","Contour cutoff must be a number!") 
         return False
    try:
         if ctr_cut < 0 or ctr_cut >360:
             raise Exception
    except Exception:
        messagebox.showerror("Input Error", "Contour cutoff must be between 0 and 360")
    try:
      flow_dir = float(flow_dir.get()) # Contour cutoff
    except ValueError:
         messagebox.showerror("Invalid Input","Contour cutoff must be a number!") 
         return False
    zobs = 0.0 # depth of observation points
    zs = 0.0
    C = 0.0 # initialize C
    try:
      nsrs = int(nsrs.get()) # Number of Sources
    except ValueError:
         messagebox.showerror("Invalid Input","Number of sources must be a number!") 
         return False
    pjt_name = pjt_name.get()
    #Obtain source data and create arrays
    results = []
    for c in range(1,7):
      column = []
      results.append(column)
      for r in range(3, nsrs+3):
            slaves = frame4.grid_slaves(row = r, column = c)
            entry = slaves[0]
            if c < 3:    # Source location needs to be int
              value = int(entry.get())
            else:
              value = float(entry.get()) # Source dimension, flux & conc float
            column.append(value)
    # Create arrays with source information
    srs_x = results[0]
    srs_y= results[1]
    srs_w= results[2]
    srs_l= results[3]
    srs_conc= results[4]
    srs_flx= results[5] 
       
    #Obtain obs point data and create arrays
    num_obs = int(nobs.get()) # Get number of obs points from tkinter Frame 5
    results_obs = []
    
    for c_obs in range(1,3):
      column_obs = []
      model_z_list = [model_z]
      model_z_list.append(model_z)
      results_obs.append(column_obs)
      for r_obs in range(3, num_obs+3):
            slaves_obs = frame5.grid_slaves(row = r_obs, column = c_obs)
            entry_obs = slaves_obs[0]
            value_obs = int(entry_obs.get())
            column_obs.append(value_obs)
            
    results_obs.append(model_z_list)           
    x_obs = results_obs[0]
    y_obs = results_obs[1]
    z_obs = results_obs[2]
    
    # Check if obs pts inside of model bounds
    for element in x_obs:
      try:  
         if element > model_x:
            raise Exception
      except Exception:
       messagebox.showerror("Input Error", "Obs points x-location must be within model boundaries") 
       return False
   
    for element in y_obs:
      try:  
         if element > model_y:
            raise Exception
      except Exception:
       messagebox.showerror("Input Error", "Obs points y-location must be within model boundaries") 
       return False
            
    # Insert "Model Running" text widget. Text insert located in HPS function.
    # Idle update inserted in HPS Function While loop so text will appear
    run_widget = tk.Text(frame1, width = 20, height =1,borderwidth=0, \
                 relief="flat",background="#F0F0F0",foreground ="blue")
    run_widget.grid(row = 10, column = 7,sticky = N)
    run_widget.insert("1.0", " Model Running ......")
    
    C_ARRAY = np.zeros((nsrs, model_y, model_x))
    C_ARRAY = np.asfortranarray(C_ARRAY)
                                           
    frame1.update()# Insert idle update so "Model Running" text will appear 
    
    #Run HPS Fortran code   
    hpsmodel.hps(srs_x, srs_y, srs_conc, srs_flx, srs_w, srs_l, C_ARRAY, 
            k, grad, por, aq_th, disp_x, disp_y, disp_z, R, Lambda, model_t, 
            nint, model_z, zs, nsrs, model_x, model_y)
   
    sum_conc = C_ARRAY.sum(axis=0) # New array of z levels sums of conc array
                                # to add conc plumes from multiple sources
    df_sum = pd.DataFrame(sum_conc) # Create dataframe of sum_conc 
    
    obs_pts = np.zeros((len(x_obs),3))
    # Obtain concentrations at observation points
    for j in range (len(x_obs)):
          obs_pts[j] = [x_obs[j],y_obs[j],df_sum.iat[y_obs[j],x_obs[j]]]

    # Add model_z value to obs_pts list for window display and pdf creation  
    obs_pts_xyzc = np.insert(obs_pts, 2, model_z, axis=1)
    df_obs = pd.DataFrame(obs_pts)
    df_obs.columns =['x','y','c']
    df_obs.insert(2, 'z', model_z) # add model depth to obs pt data
       
    # Create xyc column data for contour map and and export as csv file
    df_xyz = df_sum.stack().reset_index(name='c').rename(columns={'level_0':'x',\
                    'level_1': 'y'})
     
    df_csv = df_xyz[['y','x','c']].copy().rename(columns={'y':'y_new',\
                      'x':'x_new'}) # create df for contour map  
    
    df_csv.to_csv('hps_xyz.csv',encoding='utf-8', index=False) # create csv for export
   
    #Create contour map
    Z = pd.pivot(df_xyz, index='x', columns='y', values='c').T.values
    X_unique = np.sort(df_xyz.x.unique())
    Y_unique = np.sort(df_xyz.y.unique())
    X, Y = np.meshgrid(X_unique, Y_unique)
    Zmax = df_xyz.max(axis=0)['c']
    Zmin = df_xyz.min(axis=0)['c']
    Zmin = ctr_cut
  
    cont_int = []
    for m in range(11):
        interval = Zmin + (((Zmax-Zmin)/10) * m)
        cont_int.append(interval)
      
    # Create contour map for output tab2
    global fig
    fig = plt.figure()
    ax = fig.add_subplot(111)
    w = 6.5  # set width of figure
    h = 4.5 # set height of figure
    fig.set_figwidth(w)
    fig.set_figheight(h)
    fig.set_facecolor("#F0F0F0")
       
    global plot_size # Declare global to be used in get_pdf function
    plot_size = fig.get_size_inches()*fig.dpi # Determine plot size
    contour = ax.contourf(Y,X,Z,levels = cont_int,cmap='Spectral_r')#;scatter(Y,X)
    cont = plt.contour(Y, X, Z, levels = cont_int,colors='k', linewidths=0.5)
    
    today = datetime.date.today()
    today = today.strftime("%m/%d/%y")
    project_title = str(pjt_name + "     " + today)
    plt.title(pjt_name,x=0.10, y=1.0, fontsize=9)
    plt.tick_params(axis='both', labelsize=9)
    plt.xlabel('x (meters)', size = 9)
    plt.ylabel('y (meters)', size = 9)
    cbar = plt.colorbar(contour)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label('mg/l', size = 9)
    plt.clabel(cont, inline=1, fontsize=9)
   
    fig.tight_layout()
   
    # Create scatter plot of source locatiosn
    for num_srs in range(nsrs):
        ax.scatter(srs_x[num_srs],srs_y[num_srs],color='black', s=10)
        ax.annotate("   srs #{}".format(num_srs+1), (srs_x[num_srs], srs_y[num_srs]),fontsize=9)
    for no_obs in range(num_obs):
        ax.scatter(obs_pts[no_obs][0],obs_pts[no_obs][1],color='black', s=10)
        ax.annotate("   obs #{}".format(no_obs+1),(obs_pts[no_obs][0],obs_pts[no_obs][1]),fontsize=9)
       
    run_widget.grid_forget() # Remove "model Running..." text
    notebook.select(tab2) # change to output tab when execute model
         
    # Insert aquifer, model and source information
    # Delete and current data from tab2
    for widget in tab2.winfo_children():
      widget.destroy()  # deleting widget
      
    #Create matlibplot of contaminant contours
    canvas = FigureCanvasTkAgg(fig, master=tab2)
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=0, columnspan = 21, rowspan = 19, padx = 0, pady = 10)
    
    #Crete entry widgets and text for data results
    title_entry = Entry(tab2, width=32,font=('helvetica',10), bg =("#F0F0F0"),justify='left',bd=0)
    title_entry.grid(row=1, column=22, columnspan = 6,sticky ="NW", pady = 10)
    title_entry.insert(0, project_title) 
    title = title_entry.get()
        
    aq_info_title_entry = Entry(tab2, width=32,font=('helvetica',10,"underline"),bg =("#F0F0F0"), justify='left',bd=0)
    aq_info_title_entry.grid(row=2, column=22, columnspan = 6, sticky = 'W')
    aq_info_title_entry.insert(0, "Aquifer/Model Information: ")
    aq_info_title = aq_info_title_entry.get()
      
    hyd_cond_entry = Entry(tab2, width=30,font=('helvetica',10),bg =("#F0F0F0"), justify='left',bd=0)
    hyd_cond_entry.grid(row=3, column=22, columnspan = 6, sticky = 'W')
    hyd_cond_entry.insert(0,  "Hydraulic Conductivity (m/d) = {}".format(k))
    model_info = hyd_cond_entry.get()
       
    grad_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    grad_entry.grid(row=4, column=22, columnspan = 7, sticky = 'W')
    grad_entry.insert(0, "Gradient (m/m) = {}".format(grad))
    grad_info = grad_entry.get()
    
    por_entry = Entry(tab2, width=30,font=('Arial',10,),bg =("#F0F0F0"), justify='left',bd=0)
    por_entry.grid(row=5, column=22, columnspan = 6, sticky = 'W')
    por_entry.insert(0, "Porosity  = {}".format(por))
    por_info = por_entry.get()
    
    aq_th_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    aq_th_entry.grid(row=6, column=22, columnspan = 6, sticky = 'W')
    aq_th_entry.insert(0, "Aquifer Thicknes (m)  = {}\n".format(aq_th))
    aq_th_info = aq_th_entry.get()
    
    disp_x_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    disp_x_entry.grid(row=7, column=22, columnspan = 6, sticky = 'W')
    disp_x_entry.insert(0, "Dispersivity-x (m)  = {}\n".format(disp_x))
    disp_x_info = disp_x_entry.get()
    
    disp_y_entry = Entry(tab2, width=30,font=('Arial',10,),bg =("#F0F0F0"), justify='left',bd=0)
    disp_y_entry.grid(row=8, column=22, columnspan = 6, sticky = 'W')
    disp_y_entry.insert(0, "Dispersivity-y (m) = {}\n".format(disp_y))
    disp_y_info = disp_y_entry.get()
    
    disp_z_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    disp_z_entry.grid(row=9, column=22, columnspan = 6, sticky = 'W')
    disp_z_entry.insert(0, "Dispersivity-z (m) = {}\n".format(disp_z))
    disp_z_info =  disp_z_entry.get()
    
    R_entry = Entry(tab2, width=30,font=('Arial',10,),bg =("#F0F0F0"), justify='left',bd=0)
    R_entry.grid(row=10, column=22, columnspan = 6, sticky = 'W')
    R_entry.insert(0, "Retardation  = {}\n".format(R))
    R_info = R_entry.get()
    
    Lambda_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    Lambda_entry.grid(row=11, column=22, columnspan = 6, sticky = 'W')
    Lambda_entry.insert(0, "Decay Rate (1/years) = {}\n".format(Lambda))
    Lambda_info = Lambda_entry.get()
    
    model_t_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    model_t_entry.grid(row=12, column=22, columnspan = 6, sticky = 'W')
    model_t_entry.insert(0, "Time of Calculation (years) = {}\n".format(model_t))
    model_t_info =  model_t_entry.get()
    
    model_z_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    model_z_entry.grid(row=13, column=22, columnspan = 6, sticky = 'W')
    model_z_entry.insert(0, "Depth of Model (m) = {}\n".format(model_z))
    model_z_info =  model_z_entry.get()
    
    flow_dir_entry = Entry(tab2, width=30,font=('Arial',10),bg =("#F0F0F0"), justify='left',bd=0)
    flow_dir_entry.grid(row=14, column=22, columnspan = 6, sticky = 'W')
    flow_dir_entry.insert(0, "Flow Direction (deg) = {}\n".format(flow_dir))
    flow_dir_info =  flow_dir_entry.get()
    
    srs_header_entry = Entry(tab2, width=30,font=('Arial',10,"underline"),bg =("#F0F0F0"), justify='left',bd=0)
    srs_header_entry.grid(row=15, column=22, columnspan = 8, sticky = 'W')
    srs_header_entry.insert(0, "Contaminant Source Information")
    srs_header_info =  srs_header_entry.get()
    
    col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
    col_header_entry.grid(row=16, column=23)
    col_header_entry.insert(0, "x(m)")
    col_header_x =  col_header_entry.get()
    
    col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
    col_header_entry.grid(row=16, column=24)
    col_header_entry.insert(0, "y(m)")
    col_header_y=  col_header_entry.get()
    
    col_header_entry = Entry(tab2, width=8,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
    col_header_entry.grid(row=16, column=25)
    col_header_entry.insert(0, "wid(m)")
    col_header_w =  col_header_entry.get()
    
    col_header_entry = Entry(tab2, width=8,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
    col_header_entry.grid(row=16, column=26)
    col_header_entry.insert(0, "len(m)")
    col_header_l =  col_header_entry.get()
    
    col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
    col_header_entry.grid(row=16, column=27)
    col_header_entry.insert(0, "conc(mg/l)")
    col_header_c =  col_header_entry.get()
    
    col_header_entry = Entry(tab2, width=8,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
    col_header_entry.grid(row=16, column=28)
    col_header_entry.insert(0, "flux(gpd)")
    col_header_f =  col_header_entry.get()
    
    # Create list with source information    
    srs_data = []
    srs_data.append(srs_x)
    srs_data.append(srs_y)
    srs_data.append(srs_w)
    srs_data.append(srs_l)
    srs_data.append(srs_conc)
    srs_data.append(srs_flx)
   
    obs_data = []
    for i in range(num_obs):
        inner_list = []
        for j in range(4): # changed from 3 to 4
          inner_list.append(obs_pts_xyzc[i][j])
        obs_data.append(inner_list)
        
    # Create grid and print source information
    srs_total_rows = len(srs_flx) # Determine no. of columns
        
    for j in range(srs_total_rows):
       for i in range(6): # total # of columns
               
                srs_no_new = Entry(tab2, width=5,font=('Arial',10),\
                                justify ='center',bg =("#F0F0F0"),bd=0)
                srs_no_new.grid(row=j+17, column= 22, sticky = 'W')
                srs_no_new.insert(0, "srs #{}".format(j+1))    
                
                row_new = Entry(tab2, width=10,font=('Arial',10),\
                                justify='center',bg =("#F0F0F0"),bd=0)
                row_new.grid(row=j+17, column=i+23, sticky = 'W')
                row_new.insert(0, srs_data[i][j])   
                srs_entry = row_new.get()
    
    # Create grid and print obs points information
    obs_total_rows = srs_total_rows + 18
   
    if num_obs >= 1:    
         obs_header_entry = Entry(tab2, width=25,font=('Arial',10,"underline"),bg =("#F0F0F0"), justify='left',bd=0)
         obs_header_entry.grid(row=obs_total_rows, column=22, columnspan = 6, sticky = 'W')
         obs_header_entry.insert(0, "Observation Point Information")
         obs_Title =  obs_header_entry.get()
         
         col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
         col_header_entry.grid(row=obs_total_rows + 1, column=23)
         col_header_entry.insert(0, "x(m)")
         col_x_header =  col_header_entry.get()
         
         col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
         col_header_entry.grid(row=obs_total_rows + 1, column=24)
         col_header_entry.insert(0, "y(m)")
         col_y_header =  col_header_entry.get()
         
         col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
         col_header_entry.grid(row=obs_total_rows + 1, column=25)
         col_header_entry.insert(0, "z(m)")
         col_y_header =  col_header_entry.get()
         
         col_header_entry = Entry(tab2, width=10,font=('Arial',10, "bold"),justify='center',bg =("#F0F0F0"),bd=0)
         col_header_entry.grid(row=obs_total_rows + 1, column=26)
         col_header_entry.insert(0, "conc (mg/l")
         col_c_header =  col_header_entry.get()
         
         for y in range(num_obs):
             for x in range(4): # total # of columns
                   obs_no_new = Entry(tab2, width=8,font=('Arial',10),justify ='center',bg =("#F0F0F0"),bd=0)
                   obs_no_new.grid(row=y+obs_total_rows+2, column= 22, rowspan = 1,sticky = 'W')
                   obs_no_new.insert(0, "obs #{}".format(y+1))    
                   obs_new = Entry(tab2, width=10,font=('Arial',10),justify='center',bg =("#F0F0F0"),bd=0)
                   obs_new.grid(row=y+obs_total_rows+2, column=x+23, rowspan = 1, columnspan = 4,sticky = 'W')
                   if x <3:
                       data = int(obs_data[y][x])
                       obs_new.insert(0, data)   
                   else:
                       data = format(obs_data[y][x], '.4e')
                       obs_new.insert(0, data)   
                    
                   obs_entry = obs_new.get()
     #Create entry widgets and text  for data results
    frame12 = Frame(tab2)
    frame12.grid(row=20, column=8, rowspan = 5)
    
    label1 = Label(frame12, text='Source 1 Location') 
    label1.grid(row=21, column=8, columnspan = 2, sticky = 'W')
    label2 = Label(frame12, text='easting',pady=2) 
    label2.grid(row=22, column=8)
    label3 = Label(frame12, text='northing',pady=2) 
    label3.grid(row=22, column=9)
    label4 = Label(frame12, text='zone',pady=2) 
    label4.grid(row=22, column=10)
    label5 = Label(frame12, text='hemisphere',pady=2) 
    label5.grid(row=22, column=11)
    # Create Entry frames
    utm_x = Entry(frame12,width=10)
    utm_y = Entry(frame12,width=10)
    utm_zone = Entry(frame12,width=10)
    utm_hemi = Entry(frame12,width=10)
    #Prefill source 1 utm lovation values
    utm_x.insert(0, float("535983"))
    utm_y.insert(0, float("4408771"))
    utm_zone.insert(0, int("14"))
    utm_hemi.insert(0, str("N"))
    
    #Frame12 grid locations
    utm_x.grid(row=23,column = 8, columnspan = 1) 
    utm_y.grid(row=23, column = 9, columnspan = 1) 
    utm_zone.grid(row=23, column = 10, columnspan = 1) 
    utm_hemi.grid(row=23, column = 11, columnspan = 1) 
     
    button = tk.Button(tab2, text="Export pdf", command = \
             lambda:export_to_pdf(title, aq_info_title, model_info, grad_info,\
                    aq_th_info,disp_x_info, disp_y_info, disp_z_info,  R_info,\
                        Lambda_info, model_t_info, model_z_info, flow_dir_info,\
                            srs_data, obs_data, nsrs, num_obs))
    button.grid(row=20, column = 3, columnspan = 5, rowspan=3, padx=10, pady = 10, sticky = "N")
    
    button = tk.Button(frame12, text="Export geotif", command = \
             lambda:export_to_tif(utm_x, utm_y, utm_zone, utm_hemi,flow_dir,ctr_cut, srs_data))
    button.grid(row=20, column = 8, rowspan=1, columnspan=2, padx=0, pady = 10, sticky = "NW")
   
    time.sleep(1) 
    # store end time 
    end = time.time() 
     
    # total time taken 
    print(f"Total runtime of the program is {end - begin}") 
##############################################################################
#Create entries and labels for additional source rows
def get_num_rows(nsrs,srs_x,srs_y,srs_w,srs_l,srs_conc, srs_flx):
   
    try:
      nsrs = int(nsrs.get())  # Number of calulation points in the X-direction
      if isinstance(nsrs, int):
          pass
    except ValueError:
         messagebox.showerror("Invalid Input","No. of sources"
                             " must be an integer number!") 
         return False
     
    try:
           if nsrs >5:
               raise Exception
    except Exception:
          messagebox.showerror("Input Error", "No. of sources must be less than 5")
          return False
   
    try:
      srs_x = int(srs_x.get()) + 1 # Number of calulation points in the X-direction
      if isinstance(model_x, int):
          pass
    except ValueError:
         messagebox.showerror("Invalid Input","Source x location"
                             " must be an integer number!") 
         return False
   
    try:
      srs_y = int(srs_y.get()) + 1 # Number of calulation points in the X-direction
      if isinstance(model_y, int):
          pass
    except ValueError:
         messagebox.showerror("Invalid Input","Source y location"
                             " must be an integer number!") 
         return False
   
     
    global srs_rows_old
    srs_rows_new = nsrs
    #srs_rows_new = int(nsrs.get())
    if srs_rows_new - srs_rows_old >= 1: 
       for i in range(srs_rows_old,srs_rows_new):
           # Label row No
           label1 = Label(frame4, text=f'# {i+1}',pady=2) 
           label1.grid(row= 3+i, column=0,)
               
           #Enter Data Frame 4
           srs_x = Entry(frame4,width=10)
           srs_y = Entry(frame4,width=10)
           srs_w = Entry(frame4,width=10)
           srs_l = Entry(frame4,width=10)
           srs_conc = Entry(frame4,width=10)
           srs_flx = Entry(frame4,width=10)
          
           #Grid Location Frame 4
           srs_x.grid(row=i+3, column = 1)
           srs_y.grid(row=i+3, column = 2)
           srs_w.grid(row=i+3, column = 3)
           srs_l.grid(row=i+3, column = 4)
           srs_conc.grid(row=i+3, column = 5)
           srs_flx.grid(row=i+3, column = 6)
                  
    else:   
          # Remove all Frame 4 rows including column headers
          #max_row,max_col = frame4.grid_size()
          max_col,max_row = frame4.grid_size()
          #max_row = grid_size[0]
          for row_count in range(srs_rows_new, max_row):
                  for widget in frame4.grid_slaves(row =row_count+3):
                      widget.destroy() 
                                                     
    srs_rows_old = srs_rows_new # redfine new row num to old
    
##############################################################################      
def get_obs_points(nobs): # Make array of Value columns

       try:
         nobs = int(nobs.get())  # Number of calulation points in the X-direction
         if isinstance(nobs, int):
             pass
       except ValueError:
            messagebox.showerror("Invalid Input","No. of obs points"
                                " must be an integer number!") 
            return False
       try:
               if nobs >5:
                   raise Exception
       except Exception:
              messagebox.showerror("Input Error", "No. of sources must be less than 5")
              return False
          
       global obs_rows_old
       obs_rows_new = nobs
       #obs_rows_new = int(nobs.get())
       if obs_rows_new - obs_rows_old > 0: 
          for i in range(obs_rows_old,obs_rows_new):
             # Add Frame 1 rows 
             #Frame 1 labels 
             label1 = Label(frame5, text='x location',pady=2) 
             label1.grid(row= 2, column=1)
             label2 = Label(frame5, text='y location',pady=2) 
             label2.grid(row= 2, column=2)
             #Frame 1 row number label  
             label3 = Label(frame5, text=f'# {i+1}',pady=2) 
             label3.grid(row= 3+i, column=0)
             #Frame 1 entry widgets
             value1 = Entry(frame5,width=10)
             value2 = Entry(frame5,width=10)
             #Frame 1 grid locations
             value1= value1.grid(row=i+3, column = 1) 
             value2 = value2.grid(row=i+3, column = 2)
             
       else:   
             # Remove all Frame 1 rows including column headers
             #max_row,max_col = frame5.grid_size()
             max_col,max_row = frame5.grid_size()
             if obs_rows_new < 1:
                 for row_count in range(obs_rows_new, max_row):
                     for widget in frame5.grid_slaves(row =row_count+2):
                         widget.destroy() 
             else:
                # Remove unused Frame 1 rows
                for row_count in range(obs_rows_new, max_row):
                     for widget in frame5.grid_slaves(row=row_count+3):  
                         widget.destroy() 
                                                  
       obs_rows_old = obs_rows_new # redfine new row num to old
 #############################################################################
def export_to_pdf(title, aq_info_title, model_info, grad_info,\
       aq_th_info,disp_x_info, disp_y_info, disp_z_info,  R_info,\
           Lambda_info, model_t_info, model_z_info, flow_dir_info, srs_data,\
           obs_data, nsrs, num_obs):
    
       #Capture contaminant map
      fig.savefig("contaminant_plume.png")
      cont_map = Image.open("contaminant_plume.png")
      cont_map = cont_map.convert('RGBA')
      map_data = np.array(cont_map)   # "map_data" is a height x width x 4 numpy array
      red, green, blue, alpha = map_data.T # Temporarily unpack the bands for readability
      white_areas_map = (red == 240) & (blue == 240) & (green == 240)
      map_data[..., :-1][white_areas_map.T] = (255, 255, 255) # Transpose back needed
      map_image = Image.fromarray(map_data)
      map_image.save("contaminant_map.png")
           
      # Generate the PDF
      # Add header to pdf pages
      class PDF(FPDF):
       def header(self):
        # Set font
        self.set_font('Arial', 'B', 12)
        # Title
        self.image('HPS.ico', 5, 6, 15)
        self.cell(5, 7, '          HPS Model Results', ln = 2)
      # Create pdf   
      pdf = PDF()
      pdf.add_page()
      pdf.set_stretching(110)
      pdf.set_font("helvetica",size=12)
      pdf.cell(110, 6, txt = "  ",align = 'L',ln = 2)
      pdf.cell(110, 8, txt = title, align = 'L',ln = 2)
      pdf.set_font("helvetica", 'U',size=12)
      pdf.cell(110, 6, txt = aq_info_title,align = 'L',ln = 2)
      pdf.set_font("helvetica",size=12)
      pdf.cell(110, 6, txt = model_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = grad_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = aq_th_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = disp_x_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = disp_y_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = disp_z_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = R_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = Lambda_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = model_t_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = model_z_info,align = 'L',ln = 2)
      pdf.cell(110, 6, txt = flow_dir_info,align = 'L',ln = 2)
      
      pdf.set_font("helvetica", 'U',size=12)
      pdf.cell(110, 8, txt = "Contaminant Source Information",align = 'L',ln = 2)
      pdf.set_font("helvetica",size=12)
      # transpose srs_data list
      srs_data_trans = [[row[i] for row in srs_data] for i in range(len(srs_data[0]))]
      list_length = len(srs_data_trans) # Determine number of lists in srs_data_trans
      #add srs no to front of each list
      for i in range(list_length):
          srs_no = "srs #{}".format(i+1) 
         # Add an element to the second sublist
          srs_data_trans[i].insert(0, srs_no)
           
      # Create tuple of column header 
      tuple_srs_title =(("  ", "x(m)", "y(m)", "wid(m)", "len(m)", "conc(mg/l)",\
                      "flux (gpd)"))
      #convert list of lists to tuple of tuples
      tuple_srs_data = tuple(tuple(i) for x,i in enumerate(srs_data_trans))  
      tuple_srs_total = tuple([tuple(str(x) for x in tup) for tup in tuple_srs_data])
      #add title tupple to data tupple
      tuple_srs = (tuple_srs_title,) + tuple_srs_total
      
      #Works up to here     
      with pdf.table(width=130, col_widths=(12, 10, 10, 15, 15, 22, 20),align='L',\
                     text_align="C",first_row_as_headings=True,\
                     borders_layout="SINGLE_TOP_LINE") as table:
          for data_row in tuple_srs:
              row = table.row()
              for datum in data_row:
                  row.cell(datum)
      
      # Change concentration to scientific notation
      if num_obs >= 1:  
          for sublist in obs_data:
             sublist[3] = format(sublist[3],'4e')
        
          #add srs no to front of each list
          pdf.cell(110, 5, txt = "  ",align = 'L',ln = 2)
          pdf.set_font("helvetica", 'U',size=12)
          pdf.cell(110, 5, txt = "Observation Point Information",align = 'L',ln = 2)
          pdf.set_font("helvetica",size=12)
          obs_list_length = len(obs_data) # Determine number of lists in srs_data_trans
          for i in range(obs_list_length):
              obs_no = "obs #{}".format(i+1) 
              obs_data[i].insert(0, obs_no)
                  
          # Create tuple of column header 
          tuple_obs_title =(("  ", "x(m)", "y(m)", "z(m)", "conc(mg/l)"))
          #convert list of lists to tuple of tuples
          tuple_obs_data = tuple(tuple(i) for x,i in enumerate(obs_data))  
          tuple_obs_total = tuple([tuple(str(x) for x in tup) for tup in tuple_obs_data])
          #add title tupple to data tupple
          tuple_obs = (tuple_obs_title,) + tuple_obs_total
                
          # Create Table of obs values
          with pdf.table(width=140, col_widths=(12, 15, 15, 15, 22),align='LEFT',\
                         text_align="CENTER",first_row_as_headings=True,\
                         borders_layout="SINGLE_TOP_LINE") as table:
              for data_row in tuple_obs:
                  row = table.row()
                  for datum in data_row:
                      row.cell(datum)
            
      # Print contaminant map
      # increase y with incresae number of srs and obs. If row spacing exceeds
      # 150 then will trigger page break and force image to next page
      row_spacing = ((nsrs + num_obs)*8) +132
      if row_spacing < 150:
          pdf.image("contaminant_map.png", x = 10, y = row_spacing, w = 200)
      else:
          pdf.add_page()
          pdf.image("contaminant_map.png", x = 5, y = 25, w = 200)
     
      file_path = filedialog.asksaveasfilename(defaultextension=".pdf",\
                filetypes=[("pdf files", "*.pdf"),("All files", "*.*")])
      pdf.output(file_path, "F")
      os.remove("contaminant_plume.png")
#############################################################################
def export_to_tif (utm_x, utm_y, utm_zone, utm_hemi,flow_dir, ctr_cut, srs_data):
    try:
      utm_x = float(utm_x.get())
    except ValueError:
        messagebox.showerror("Invalid Input","utm x must be a number!")
        return False    
    try:
      utm_y = float(utm_y.get())
    except ValueError:
         messagebox.showerror("Invalid Input","utm y must be a number!")
         return False
    try: 
      utm_zone = float(utm_zone.get()) 
      utm_zone >= 1 or utm_zone <= 60
    except ValueError:
        messagebox.showerror("Invalid Input","Utm zone must be a number between 1 and 60!")
        return False
    try: 
      utm_hemi = str(utm_hemi.get()) 
      utm_hemi == "N" or utm_hemi == "n" or utm_hemi == "S" or utm_hemi == "s"
    except ValueError:
         messagebox.showerror("Invalid Input","Hemispeher must be north or south!") 
         return False
   
    #Open and read csv file created by Def HPS (one meter square grid)
    # with unrotated cartesian coordinates & concentration values
    with open('hps_xyz.csv', 'r') as file:
        csv_reader = csv.reader(file)
        list_csv_string = (list(csv_reader))
        list_csv_string.pop(0) #Remove header row

    #convert string values to float
    list_csv_num = []
    for sublist in list_csv_string:
        float_sublist = []
        for item in sublist:
            float_sublist.append(float(item))
        list_csv_num.append(float_sublist)

    # Define cartesian coordinates of pivot point 
    x_pivot =  srs_data[0][0]
    y_pivot = srs_data[1][0]

    # Define corresponding utm coordinates/refernce of pivot point
    utm_x_pivot = utm_x
    utm_y_pivot  = utm_y
    #utm_zone = 14
   # utm_hemi = 'n'

    # Define UTM ESPG number
    if utm_hemi =='n' or 'N':
        epsg = int(32600 + utm_zone)
    else:
        epsg = int(32700 + utm_zone)
    epsg = f'EPSG:{epsg}'

    # Reference unrotated cartesian to orgin (0,0)
    list_diff_num = [([None] * 3) for _ in range(len(list_csv_num))] 
    for i in range(len(list_csv_num)):
        list_diff_num[i][0] = list_csv_num[i][0] - x_pivot 
        list_diff_num[i][1] = list_csv_num[i][1] - y_pivot 
        list_diff_num[i][2] = list_csv_num[i][2]

    #Define angle of rotation
    deg = flow_dir
    rad_deg = math.radians(-deg + 90)

    # Create list for rotated cartesian points
    list_utm_rot = [([None] * 3) for _ in range(len(list_csv_num))] 

    #Rotate all the points
    for i in range(len(list_csv_num)):
        # Rotate X coordinates
        list_utm_rot[i][0] = (utm_x_pivot + ((list_diff_num[i][0]) *  
         math.cos(rad_deg)) - (list_diff_num[i][1]) * math.sin(rad_deg))
        
        # Rotate Y coordinates
        list_utm_rot[i][1] = (utm_y_pivot + ((list_diff_num[i][0]) * 
          math.sin(rad_deg)) + (list_diff_num[i][1]) * math.cos(rad_deg))  
        
        # Z values
        list_utm_rot[i][2] = list_diff_num[i][2]  

    #Convert rotated list into numpy array         
    list_utm_array = np.array(list_utm_rot, dtype = float)

    # Save/Export data as text file
    np.savetxt("hps_utm.csv", list_utm_array, delimiter=",")

    # Remove values below cutoff
    mask = (list_utm_array[:, 2] >= ctr_cut)
    list_utm_array = list_utm_array[mask]

    #Create equal distance grid spacing for gdal.Translate    
    # Create x,y,z arrays for gridding
    x = list_utm_array[:,0]
    y = list_utm_array[:,1]
    z = list_utm_array[:,2]

    # Determine min and max x,y,z values
    min_x = int(min(x))
    max_x = int(max(x))
    min_y = int(min(y))
    max_y = int(max(y))
    #min_z = int(min(z))
    #max_z = int(max(z)) * 1.2 # multiply by 1.2 for greater interpolated z values

    #Calculate variable line spacing in meshgrid
    x_diff = max_x - min_x
    y_diff = max_y - min_y

    if x_diff >= y_diff:
        lnsp = int(x_diff * 10 )
    else:
        y_diff > x_diff
        lnsp = int(y_diff * 10) 

    # Create empty x and y axes for meshgrid
    xi = np.linspace(min_x, max_x, lnsp)
    yi = np.linspace(min_y, max_y, lnsp)

    # Create a empty meshgrid 
    X, Y = np.meshgrid(xi, yi)

    # Interpolate rotated utm values onto the grid
    # from scipy.interpolate import griddata
    Z = griddata((x, y), z, (X, Y), method='cubic')

    #Create Numpy Array of gridded data
    xyz_arr = np.array(Z)

    # Caluate width, height and band No.
    no_of_bands = 1
    height = xyz_arr.shape[0]
    width = xyz_arr.shape[1]

    #sAssign X,Y max & min and calculate resolution
    x1 = min_x
    y1 = min_y
    x2 = max_x
    y2 = max_y
    x_res = (x2 - x1)/width
    y_res = (y2 - y1)/width

    file_path = filedialog.asksaveasfilename(defaultextension=".tif", 
       filetypes=[("GeoTIFF", "*.tif"), ("All Files", "*.*")])

    # Load driver, create dataset, transform and assign projection
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(file_path, width, height, no_of_bands, gdal.GDT_Float32)
    out_ds.SetGeoTransform((x1, x_res, 0, y1, 0, y_res))
    out_ds.SetProjection(epsg)  
    print(epsg)
    # Write the array to the dataset's band
    band = out_ds.GetRasterBand(1)
    band.WriteArray(xyz_arr)
    band.FlushCache()
    out_ds = None
     
##############################################################################       
# HPS Tkinter GUI

#setting tkinter window size
width= root.winfo_screenwidth() 
height= root.winfo_screenheight()
root.geometry("%dx%d" % (width, height))

# Create a notebook widget
notebook = ttk.Notebook(root)
notebook.grid(row =0, column =0)

# Create frames for each tab
tab1 = tk.Frame(notebook)
tab2 = tk.Frame(notebook)

# Add tabs to the notebook
notebook.add(tab1, text="    Input    ")
notebook.add(tab2, text="    Output    ")

# Add a frame to tab1
frame1 =  Frame(tab1)
# Add a frame to tab2
frame2 =  Frame(tab2)

normal_font = tkFont.Font(family="Arial", size=10, weight=tkFont.NORMAL)
root.title("HPS Model")
root.iconbitmap("HPS.ico")
root.state('zoomed')
width= root.winfo_screenwidth()               
height= root.winfo_screenheight()               
root.geometry("%dx%d" % (width, height))# Set GUI to full screen
frame1 = Frame(tab1, relief="groove")
frame1.grid(row=0, column=0)
frame2=LabelFrame(frame1, text='Aquifer Info',padx=18, pady=10)
frame2.grid(row=1, column=1,sticky = N)
frame3=LabelFrame(frame1, text='Model Info',padx=52, pady=10)
frame3.grid(row=1, column=5,sticky = NW)
frame3.columnconfigure(0, weight =2, pad = 20)
frame4=LabelFrame(frame1, text='Source Info',padx=5, pady=10)
frame4.grid(sticky = N,row=8, column=1)
frame5=LabelFrame(frame1, text='Observation Points',padx=10, pady=10)
frame5.grid(sticky = N,row=8, column=5)
frame6=LabelFrame(frame1, width=20,borderwidth=0) # add spcae between frames
frame6.grid(row=1, column=4)
frame7=LabelFrame(frame1,borderwidth=0,padx=10, pady=10) #Project name box
frame7.grid(row=0, column=1,sticky = N)
frame8=LabelFrame(frame1,width=10,padx=50, pady=10) # add spcae between frames
frame8.grid(row=0, column=0,sticky = N)
#Create input boxes and labels
#Frame 2 Aquifer Information
label1 = Label(frame2, text='Hydraulic Conductivity (m/d) = ',pady=2) 
label1.grid(row=1, column=0)
label2 = Label(frame2, text='Porosity = ',pady=2) 
label2.grid(row=2, column=0)
label3 = Label(frame2, text='gradient (m/m) = ',pady=2) 
label3.grid(row=3, column=0)
label4 = Label(frame2, text='Aquifer Thickness (m) = ',pady=2) 
label4.grid(row=4, column=0)
label5 = Label(frame2, text='x',padx = 4) 
label5.grid(row=5,column = 1)
label6 = Label(frame2, text='y',padx = 4) 
label6.grid(row=5,column = 2)
label7 = Label(frame2, text='z',padx = 4) 
label7.grid(row=5,column = 3)
label8 = Label(frame2, text='Dispersivity (m) = ',pady=2) 
label8.grid(row=6, column=0)
label9 = Label(frame2, text='Retardation = ',pady=2) 
label9.grid(row=7, column=0)
label10 = Label(frame2, text='Decay Rate (1/years) = ',pady=2) 
label10.grid(row=8, column=0)
label11 = Label(frame2, text='Flow Direction (0 - 359 deg) = ',pady=2) 
label11.grid(row=9, column=0)
label12 = Label(frame2, text='        ',pady=2) 
label12.grid(row=10, column=4)

#Frame 3 Model Information
label1 = Label(frame3, text='Model Horizontal Length (m) = ',pady=2) 
label1.grid(row=1, column=4)
label2 = Label(frame3, text='Model Vertical Length(m)  = ',pady=2) 
label2.grid(row=2, column=4)
label3 = Label(frame3, text='Depth of Model (m) = ',pady=2) 
label3.grid(row=3, column=4)
label4 = Label(frame3, text='Model Time Duration (years)= ',pady=2) 
label4.grid(row=4, column=4)
label5 = Label(frame3, text='No. Model Iterations =',pady=2) 
label5.grid(row=5, column=4)
label6 = Label(frame3, text='Contour Cutoff (mg/l) =',pady=2) 
label6.grid(row=6, column=4)
label7 = Label(frame3, text='   ',pady=2) 
label7.grid(row=7, column=4)
label8 = Label(frame3, text='   ',pady=2) 
label8.grid(row=8, column=4)
label9 = Label(frame3, text='   ',pady=2) 
label9.grid(row=9, column=4)
label10 = Label(frame3, text='   ',pady=1) 
label10.grid(row=10, column=4)

#Frame 4 Source Information
label1 = Label(frame4, text='No. of Sources = ',pady=2) 
label1.grid(row=1, column=0, columnspan=2)
nsrs = Entry(frame4,width=10)
nsrs.grid(row=1, column = 2)
nsrs.insert(0, int("1")) # Prefill No. of Values
label1 = Label(frame4, text='#1',pady=2) 
label1.grid(row= 3, column=0)
label2 = Label(frame4, text='x location',pady=2) 
label2.grid(row=2, column=1)
label3 = Label(frame4, text='y location',pady=2) 
label3.grid(row=2, column=2)
label4 = Label(frame4, text='width (m)',pady=2) 
label4.grid(row=2, column=3)
label5 = Label(frame4, text='length (m)',pady=2) 
label5.grid(row=2, column=4)
label6 = Label(frame4, text='conc (mg/l)',pady=2) 
label6.grid(row=2, column=5)
label7 = Label(frame4, text='flux (gpd)',pady=2) 
label7.grid(row=2, column=6)

#Frame 5 Obs Pts Information
label1 = Label(frame5, text='No of Obs Points = ',pady=2) 
label1.grid(row=1, column=0, columnspan=2)

#Frame 7
label1 = Label(frame7, text='Project Name ',pady=2) 
label1.grid(row=1, column=0, sticky=W)

#Enter Data Frame 2
k = Entry(frame2,width=10)    
por = Entry(frame2,width=10)
grad = Entry(frame2,width=10)
aq_th = Entry(frame2,width=10)
disp_x = Entry(frame2,width=10)
disp_y = Entry(frame2,width=10)
disp_z = Entry(frame2,width=10)
R = Entry(frame2,width=10)
Lambda = Entry(frame2,width=10)
flow_dir = Entry(frame2,width=10)

#Prefill Frame 2 values
k.insert(0, float("2.75"))
por.insert(0, float("0.2"))
grad.insert(0, float("0.002"))
aq_th.insert(0, float("50"))
disp_x.insert(0, float("2"))
disp_y.insert(0, float("0.2"))
disp_z.insert(0, float("0.1"))
R.insert(0, float("1.0"))
Lambda.insert(0, float("0.000"))
flow_dir.insert(0, float("45"))

#Enter Data Frame 3
model_x = Entry(frame3,width=10)
model_y = Entry(frame3,width=10)
model_z = Entry(frame3,width=10)
nint = Entry(frame3,width=10)
model_t =  Entry(frame3,width=10)
ctr_cut =  Entry(frame3,width=10)

#Prefill Frame 3 values
model_x.insert(0, int("30"))
model_y.insert(0, int("12"))
model_z.insert(0, int("0"))
nint.insert(0, int("250"))
model_t.insert(0, float("100"))
ctr_cut.insert(0, float("0.00"))

#Enter Data Frame 4
srs_x = Entry(frame4,width=10)
srs_y = Entry(frame4,width=10)
srs_w = Entry(frame4,width=10)
srs_l = Entry(frame4,width=10)
srs_conc = Entry(frame4,width=10)
srs_flx = Entry(frame4,width=10)

# Prefill values in all columns of Row 1 of Data Frame 4
srs_x.insert(0, int("2"))
srs_y.insert(0, int("6"))
srs_w.insert(0, int("10"))
srs_l.insert(0, int("1"))
srs_conc.insert(0, float("45"))
srs_flx.insert(0, float("5"))

#Enter Data Frame 5
nobs = Entry(frame5,width=10)
nobs.insert(0, int("0")) # Prefill No. of Obs Points

#Enter Data Frame 7
pjt_name =  Entry(frame7,width=55)
#Prefill Frame 2 values
pjt_name.insert(0, ("Project Name"))

#Grid Location Frame 2
k.grid(row=1, column = 1)
por.grid(row=2, column = 1)
grad.grid(row=3, column = 1)
aq_th.grid(row=4, column = 1)
disp_x.grid(row=6, column = 1)
disp_y.grid(row=6, column = 2)
disp_z.grid(row=6, column = 3)
R.grid(row=7, column = 1)
Lambda.grid(row=8, column = 1)
flow_dir.grid(row=9, column = 1)

#Grid Location Frame 3
model_x.grid(row=1, column = 5)
model_y.grid(row=2, column = 5)
model_z.grid(row=3, column = 5)
model_t.grid(row=4, column = 5)
nint.grid(row=5, column = 5)
ctr_cut.grid(row=6, column = 5)

#Grid Location Frame 4
srs_x.grid(row=3, column = 1)
srs_y.grid(row=3, column = 2)
srs_w.grid(row=3, column = 3)
srs_l.grid(row=3, column = 4)
srs_conc.grid(row=3, column = 5)
srs_flx.grid(row=3, column = 6)

srs_rows_old = 1 # Define initial srs row num value for use in function

#Grid Location Frame 5
nobs.grid(row=1, column = 2, padx=10)

obs_rows_old = 0 # Define initial obs row num value for use in function

pjt_name.grid(row=1, column = 1)

#Add Figure to Input
# Convert the image to a Tkinter-compatible photo image
#image = Image.open(r"C:\Users\rpahy\OneDrive\Documents\HPS\HPS Python\HPS_Figure.png")
image = Image.open(r"HPS_Figure.png")
# Resize the image using the new resampling attribute
image = image.resize((370, 270), Image.Resampling.LANCZOS)

# Convert the image to a Tkinter-compatible photo image
photo_image = ImageTk.PhotoImage(image,master=root)

# Create a label widget to hold the image
image_label = Label(frame1, image=photo_image)

# Place the label in the frame
image_label.grid(row=1, column = 7, padx = 20, pady = 0)

#Create submit button and call HPS function
myButton = Button(frame1, text='Run Model', command = lambda:HPS(pjt_name,k,grad,por,\
           aq_th,disp_x, disp_y, disp_z, R, Lambda,flow_dir,nsrs,model_t,\
           model_x,model_y,model_z,nint,ctr_cut))
myButton.grid(row = 8, column = 7, sticky = 'N', padx = 10, pady=20)

# Button for adding/removing sources
add_fields_button = Button(frame4, text="Enter No.of Sources", \
           command=lambda:get_num_rows(nsrs,srs_x,srs_y,srs_w,srs_l,srs_conc,\
           srs_flx))
add_fields_button.grid(row=1,column=3, columnspan = 2)  

# Button for adding/removing obs points
add_fields_button = Button(frame5, text="Enter No.of obs points", \
            command=lambda:get_obs_points(nobs))
add_fields_button.grid(row=1,column=3,columnspan = 2)  

root.mainloop()


      