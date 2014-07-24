## Plot_relax.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes a contour plot of the relaxation time (0.4053*L_H^2/D)
## as a function of L_H and D
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 23/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import matplotlib.pyplot as pp
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rcParams
import matplotlib.lines as mpllines
import os
from glob import glob

def Plot_relax():

  # These set the font size
  label_size = 20
  #title_size = 30
  axis_size = 28

  # set some variables for the fonts
  #rcParams['font.family'] = 'sans-serif'
  #rcParams['font.sans-serif'] = ['arial']
  
  rcParams['font.size'] = label_size  
  delta = 0.05
  logL_h = np.arange(1, 3, delta)
  logD = np.arange(-5, -1, delta)
  LL_H, L_D = np.meshgrid(logL_h, logD)
  
  L_H = np.power(10,LL_H)
  D = np.power(10,L_D)
  
  L_H2 = np.power(L_H,2)
  
  top_term = np.multiply(L_H2,0.4053)
  
  relax = np.divide(top_term,D)
  Log10_relax = np.log10(relax)
  
  #print LL_H
  #print L_D
  #print relax
  

    
  # now plot the results    
  fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
  ax1 = fig.add_subplot(1,1,1)  
  #countour_range = np.arange(3, 9, 0.2)
  #im = pp.imshow(Log10_relax, interpolation='bilinear', 
  #              cmap=cmx.gray, extent=(10, 500, 0.0001, 0.1))
  CS = pp.contour(L_H, D, Log10_relax)
  pp.axis([10, 500, 0.0001, 0.1])
  
  
  #pp.clabel(CS, inline=1, fontsize=10)
  
  pp.rcParams['xtick.direction'] = 'out'
  pp.rcParams['ytick.direction'] = 'out'
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  
  ax1.spines['top'].set_linewidth(2.5)
  ax1.spines['left'].set_linewidth(2.5)
  ax1.spines['right'].set_linewidth(2.5)
  ax1.spines['bottom'].set_linewidth(2.5) 
  ax1.tick_params(axis='both', width=2.5)    

  for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)

  for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

  pp.ylabel('$D$ ($m^2$ $yr^{-1}$)',fontsize = axis_size)
  pp.xlabel('$L_H$ ($m$)',fontsize = axis_size) 

  pp.tight_layout()
      
  pp.savefig("Relax_contour.eps", format='eps')  
  
  pp.cla()
  pp.clf()
  pp.close(fig)

  
if __name__ == "__main__":
    Plot_relax()  