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
  #rcParams['font.family'] = 'sans-serif'
  #rcParams['font.sans-serif'] = ['arial']
  rcParams['font.size'] = label_size  
  
  
  # get the initial erosion rate
  # set these to zero. There are four initial rates
  initial_erosion = np.zeros(4)
  
  # read the data from the first file
  f1 = open('test_params_data1.data', 'r')
  header1 = f1.readline()
  maxerr1 = []
  logr1 = []
  ratio1 = []
  tmax1 = []
  
  for line in f1:
    line = line.strip()
    columns = line.split()
    initial_erosion[0] = float(columns[0])
    logr1.append( float(columns[2]) )   
    ratio1.append (float(columns[1]) )
    maxerr1.append( float(columns[4]) )
    tmax1.append( float(columns[6]) )   
  f1.close()

  lr1 = np.asarray(logr1)
  r1 =  np.asarray(ratio1)
  me1 = np.asarray(maxerr1)
  tm1 = np.asarray(tmax1)
  
  print r1

  f2 = open('test_params_data2.data', 'r')
  header2 = f2.readline()
  maxerr2 = []
  logr2 = []
  ratio2 = []
  tmax2 = []
  
  for line in f2:
    line = line.strip()
    columns = line.split()
    initial_erosion[1] = float(columns[0])
    logr2.append( float(columns[2]) )   
    ratio2.append (float(columns[1]) )
    maxerr2.append( float(columns[4]) )
    tmax2.append( float(columns[6]) )   
  f2.close()

  lr2 = np.asarray(logr2)
  r2 =  np.asarray(ratio2)
  me2 = np.asarray(maxerr2)
  tm2 = np.asarray(tmax2)    
  
  print r2

  f3 = open('test_params_data3.data', 'r')
  header3 = f3.readline()
  maxerr3 = []
  logr3 = []
  ratio3 = []
  tmax3 = []
  
  for line in f3:
    line = line.strip()
    columns = line.split()
    initial_erosion[2] = float(columns[0])
    logr3.append( float(columns[2]) )   
    ratio3.append (float(columns[1]) )
    maxerr3.append( float(columns[4]) )
    tmax3.append( float(columns[6]) )   
  f3.close()

  lr3 = np.asarray(logr3)
  r3 =  np.asarray(ratio3)
  me3 = np.asarray(maxerr3)
  tm3 = np.asarray(tmax3)   

  f4 = open('test_params_data4.data', 'r')
  header4 = f4.readline()
  maxerr4 = []
  logr4 = []
  ratio4 = []
  tmax4 = []
  
  for line in f4:
    line = line.strip()
    columns = line.split()
    initial_erosion[3] = float(columns[0])
    logr4.append( float(columns[2]) )   
    ratio4.append (float(columns[1]) )
    maxerr4.append( float(columns[4]) )
    tmax4.append( float(columns[6]) )   
  f4.close()

  lr4 = np.asarray(logr4)
  r4 =  np.asarray(ratio4)
  me4 = np.asarray(maxerr4)
  tm4 = np.asarray(tmax4)   

  
  # now get the labels
  label1 = str(initial_erosion[0])
  label2 = str(initial_erosion[1])
  label3 = str(initial_erosion[2])
  label4 = str(initial_erosion[3])
  
  fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
  ax1 = fig.add_subplot(2,1,1)
  ax1.semilogx(r1,me1,linewidth=3, label = label1)
  ax1.semilogx(r2,me2,linewidth=3, label = label2)
  ax1.semilogx(r3,me3,linewidth=3, label = label3)
  ax1.semilogx(r4,me4,linewidth=3, label = label4)
  #,r2,me2,linewidth=3,r3,me3,linewidth=3,r4,me4,linewidth=3
  ax1.legend()

  ax1.spines['top'].set_linewidth(2.5)
  ax1.spines['left'].set_linewidth(2.5)
  ax1.spines['right'].set_linewidth(2.5)
  ax1.spines['bottom'].set_linewidth(2.5) 
  ax1.tick_params(axis='both', width=2.5)    

  for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)

  for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

  pp.xlabel('Ratio of $\hat{U}_{final}/\hat{U}_{initial}$',fontsize = axis_size)
  pp.ylabel('Max error in $\hat{R}$',fontsize = axis_size)    



  # add the next subplot
  ax2 = fig.add_subplot(2,1,2)
  ax2.semilogx(r1,tm1,linewidth=3)
  ax2.semilogx(r2,tm2,linewidth=3)
  ax2.semilogx(r3,tm3,linewidth=3)
  ax2.semilogx(r4,tm4,linewidth=3)
  #,r2,tm2,linewidth=3,r3,tm3,linewidth=3,r4,tm4,linewidth=3

  ax2.spines['top'].set_linewidth(2.5)
  ax2.spines['left'].set_linewidth(2.5)
  ax2.spines['right'].set_linewidth(2.5)
  ax2.spines['bottom'].set_linewidth(2.5) 
  ax2.tick_params(axis='both', width=2.5)    

  for line in ax2.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)

  for line in ax2.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

  pp.xlabel('Ratio of $\hat{U}_{final}/\hat{U}_{initial}$',fontsize = axis_size)
  pp.ylabel('$\hat{t}$ to max error',fontsize = axis_size)    
  
  pp.tight_layout()
        
  pp.savefig('Perturbation_EsRs.png', format='png')
  
if __name__ == "__main__":
    read_Es_Rs_data()  