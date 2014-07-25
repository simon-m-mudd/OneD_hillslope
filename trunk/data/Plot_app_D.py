## Plot_app_D.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This plots the evolution in time of apparent diffusivity after a 
## step change in actual diffusivity
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 17/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import matplotlib.pyplot as pp
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rcParams
import matplotlib.lines as mpllines
import os
from glob import glob

def Plot_app_D():

  # These set the font size
  label_size = 20
  #title_size = 30
  axis_size = 28

  # set some variables for the fonts
  #rcParams['font.family'] = 'sans-serif'
  #rcParams['font.sans-serif'] = ['arial']
  rcParams['font.size'] = label_size  

  
  for FileName in glob("*.CRNpdata_v3"): # assigns a number to each iteration (i.e. for every .tree file in the directory)
  #if label_size == 20:
    #FileName = "CRN_oneD_hs.0_0_1.CRNpdata"
    print "fname: "+FileName

    # open the file
    f1 = open(FileName, 'r')
    header1 = f1.readline()
    
    line = header1.strip()
    columns = line.split()
    
    erate = columns[0];
    start_D = columns[1];
    D_ratio = columns[2];
    
    
    split_fname = FileName.split('.')    
    no_tree_levs = len(split_fname)
    fnumbers = split_fname[no_tree_levs-2]
    figfiletitle = "Plot_appD_vs_that_v2" + fnumbers + ".eps"
    print "figtitle is: "  + figfiletitle
    
    title = "E: " + erate + " $m$ $yr^{-1}$, initial D: " + start_D + " $m^2$ $yr^{-1}$, D ratio: "+D_ratio
    print "title is: "+title  
    
    # initialise the lists for the data
    t_imet = []
    t_ime_hatt = []
    zeta_rtt = []
    app_E_10Bet = []
    app_E_14Ct  = []
    app_D_10Bet = []
    app_D_14Ct = []
    Derr_10Bet = []
    Derr_14Ct = []
    
    # get the data from the lines
    for line in f1:
      line = line.strip()
      columns = line.split()
      t_imet.append( float(columns[0]) )   
      t_ime_hatt.append (float(columns[1]) )
      zeta_rtt.append (float(columns[2]) )
      app_E_10Bet.append( float(columns[3]) )
      app_E_14Ct.append( float(columns[4]) ) 
      app_D_10Bet.append( float(columns[5]) )
      app_D_14Ct.append( float(columns[6]) )  
      Derr_10Bet.append( float(columns[7]) )
      Derr_14Ct.append( float(columns[8]) )              
    f1.close()
  
    # now get the data as arrays
    t_ime = np.asarray(t_imet)
    t_ime_hat =  np.asarray(t_ime_hatt)
    app_E_10Be = np.asarray(app_E_10Bet)
    app_E_14C = np.asarray(app_E_14Ct)    
    app_D_10Be = np.asarray(app_D_10Bet)
    app_D_14C = np.asarray(app_D_14Ct)  
    Derr_10Be = np.asarray(Derr_10Bet)
    Derr_14C = np.asarray(Derr_14Ct)  
    
    #print "t_ime_hat is: "
    #print t_ime_hat
    #print " Derr_10Be is: "
    #print Derr_10Be
    
    # now plot the results    
    fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
    ax1 = fig.add_subplot(1,1,1)
    ax1.semilogy(t_ime_hat,Derr_10Be,linewidth=3)
    pp.axis([0, 0.5, 0.001, 5])
 
    ax1.spines['top'].set_linewidth(2.5)
    ax1.spines['left'].set_linewidth(2.5)
    ax1.spines['right'].set_linewidth(2.5)
    ax1.spines['bottom'].set_linewidth(2.5) 
    ax1.tick_params(axis='both', width=2.5)    
  
    for line in ax1.get_xticklines():
      line.set_marker(mpllines.TICKDOWN)
  
    for line in ax1.get_yticklines():
      line.set_marker(mpllines.TICKLEFT)
  
    pp.xlabel('$\hat{t}$ after change in $D$ $m^2$ $yr^{-1}$',fontsize = axis_size)
    pp.ylabel('(Error in $D$)/$D$ (dimensionless)',fontsize = axis_size) 
    pp.title(title, fontsize = label_size)   

    pp.tight_layout()
        
    pp.savefig(figfiletitle, format='eps')  
    
    pp.cla()
    pp.clf()
    pp.close(fig)

  
if __name__ == "__main__":
    Plot_app_D()  