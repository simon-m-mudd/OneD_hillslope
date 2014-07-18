## Plot_ER_perturb.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Function to plot the reposnse in E* R* space after a step change in
## erosion rates
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 17/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import matplotlib.pyplot as pp
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rcParams
import matplotlib.lines as mpllines

def read_Es_Rs_data():

  # These set the font size
  label_size = 20
  #title_size = 30
  axis_size = 28

  # set some variables for the fonts
  rcParams['font.family'] = 'sans-serif'
  rcParams['font.sans-serif'] = ['arial']
  rcParams['font.size'] = label_size  
  


  f = open('test_params_data1.data', 'r')
  header1 = f.readline()
  maxerr1 = []
  logr1 = []
  tmax1 = []
  
  i = 0
  for line in f:
    line = line.strip()
    columns = line.split()
    logr1[i] = float(columns[4])
    maxerr1[i]  = float(columns[4])
    tmax1[i]  = float(columns[4])
    
  
  fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
  ax1 = fig.add_subplot(2,1,1)
  ax1.plot(logr1,maxerr1)
  
  pp.xlabel('Log ratio',fontsize = axis_size)
  pp.ylabel('Max error',fontsize = axis_size)    
  
  pp.show
  
if __name__ == "__main__":
    read_Es_Rs_data()  