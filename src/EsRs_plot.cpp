// this program takes a file that contains the mean and
// standard deviation of the three hillslope parameters
// t_star_peak
// U_star_peak
// U_star_width
// and prints out results from the mean, minimum
// and maximum versions of the model runs
// run from command line with something like
// EsRs_plotting.exe fitted_U_star_params.param DBPR_Estar_Rstar_data_SD fit_EsRs_trajectories

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "OneDImplicitHillslope.hpp"
using namespace std;

string itoa(int num)
{
    stringstream converter;
    converter << num;
    return converter.str();
}

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: not enough inputs. The program needs 1)fitted param filename" << endl;
		cout << "2) the fname of the data and 3) the prefix of the outfile" << endl;
		exit(EXIT_SUCCESS);
	}

	string pathname = argv[1];
	cout << "the path is: " << pathname << endl;

	string param_name = argv[2];
	cout << "param_name is: " << param_name << endl;

	string param_suffix = argv[3];
	cout << "out param_suffix is: " << param_suffix << endl;

	string of_prefix = param_name;
	string dot = ".";
	string num;

  string pfname = pathname+param_name+dot+param_suffix;
  cout << "Parameter filename is: " << pfname << endl;
    
  // open the parameter file
  ifstream infile;
  infile.open(pfname.c_str());

  // read in the paramaters
  // first parameters for the starting values of Ustar
  double logUstart_low; 
  double logUstart_high;
  int N_Ustart;
  
  infile >> logUstart_low >> logUstart_high >> N_Ustart;
  infile.close();
  
  // create a hillslope
  OneDImplicitHillslope Hillslope;

  // first create a file for the theoretical E* vs R* curve
  string EsRs_fname =  pathname+"EstarvsRstar.data";
  cout << "EsRs filename is: " << EsRs_fname << endl;
  ofstream esRs_out;
  esRs_out.open(EsRs_fname.c_str());
  int NEstar_nodes = 201;
  double log_start_Es = -0.7;
  double dLog = 0.02;
  double log_Estar;
  double Estar;
  double Rstar;
  cout << "Printing the SS Estar vs Rstar curve" << endl;
  for (int i = 0; i< NEstar_nodes; i++)
  {
     log_Estar = log_start_Es+double(i)*dLog;
     Estar = pow(10,log_Estar);
     Rstar = Hillslope.analytical_R_star(Estar);
     
     esRs_out << Estar << "\t" << Rstar << endl;  
  }
  esRs_out.close();

  // now do some time series. 
  double start_estar;
  double end_estar;
  for(int i = 1; i<=4; i++)
  {
    string uscore = "_";
    string this_num = itoa(i);
    string data_ext = "EsRsdata";
    string this_outfile = pathname+param_name+uscore+this_num+dot+data_ext;
    cout << "Now doing timeseries, filename is: " << this_outfile << endl;
    
    ofstream EsRs_transient_out;
    EsRs_transient_out.open(this_outfile.c_str());
    
    if(i == 1)
    {
      start_estar = 0.2;
      end_estar = 30;
    }
    else if(i == 2)
    {
      start_estar = 1;
      end_estar = 10;
    }    
    else if(i == 3)
    {
      start_estar = 10;
      end_estar = 1;
    }     
    else if(i == 4)
    {
      start_estar = 30;
      end_estar = 0.2;
    }  
    
    // now initiate a hillslope          
    OneDImplicitHillslope thisHillslope;
    
    // set to steady state
    thisHillslope.set_analytical_steady(start_estar);
    double this_E_star = thisHillslope.calculate_E_star();
    EsRs_transient_out << "0\t" << thisHillslope.calculate_E_star()
                       << "\t" << thisHillslope.calculate_R_star() 
                       << "\t" << thisHillslope.analytical_R_star(this_E_star) << endl;
        
    // now run for 3 relaxation times
    double dt_hat = 0.0005;
    double t_ime_hat = 0;
    double end_time_hat = 0.5;
    double tolerance = 0.00001;
    double print_spacing = 0.005;
    double next_print = print_spacing;
    cout << "Starting hillslope loop" << endl;
    while (t_ime_hat < end_time_hat)
    {
      thisHillslope.hillslope_timestep(dt_hat, t_ime_hat, end_estar, tolerance);
      
      if(t_ime_hat >= next_print)
      {
        //cout << "Printing to file, time is: " << t_ime_hat << " and end time: " << end_time_hat << endl;
        double this_E_star = thisHillslope.calculate_E_star();
        EsRs_transient_out << t_ime_hat << "\t"
                           << thisHillslope.calculate_E_star() << "\t" 
                           << thisHillslope.calculate_R_star() << "\t"
                           << thisHillslope.analytical_R_star(this_E_star) << endl;
        next_print += print_spacing;                  
      }
    
    }
    
    EsRs_transient_out.close();
    
  }

  // now print the hillslopes for the max pulse of incision 
  start_estar = 1;
  end_estar = 10;
  OneDImplicitHillslope profHillslope(start_estar,start_estar,start_estar, 0.1);
  
  vector<double> x_hat;
  Array1D<double> zeta_init;
  Array1D<double> zeta_final;
  Array1D<double> zeta_intermediate;

  x_hat = profHillslope.get_x_hat();
  profHillslope.set_analytical_steady(start_estar);
  zeta_init = profHillslope.get_zeta_hat();

  // now run until the you get a big difference
  double dt_hatt = 0.0005;
  double t_ime_hatt = 0;
  double end_time_hatt = 0.05;
  double ttolerance = 0.000001;
  while (t_ime_hatt < end_time_hatt)
  {
    profHillslope.hillslope_timestep(dt_hatt, t_ime_hatt, end_estar, ttolerance);
  }
  
  // now get zeta
  zeta_intermediate =   profHillslope.get_zeta_hat();
  
  // and get the zeta with the same E* but with the proper R*
  double this_estar = profHillslope.calculate_E_star();
  profHillslope.set_analytical_steady(this_estar);
  zeta_final = profHillslope.get_zeta_hat();
  
  // now print to file
  string HS_prof_name = "HS_prof";
  string HS_profiles = pathname + param_name+"_"+HS_prof_name+".data";
  ofstream HS_prof_out;
  HS_prof_out.open(HS_profiles.c_str());
  
  int sz_prof = zeta_final.dim1();
  for (int i = 0; i<sz_prof; i++)
  {
    HS_prof_out << x_hat[i] << "\t" << zeta_init[i] << "\t" 
                << zeta_intermediate[i] << "\t" << zeta_final[i] << endl;
  }
  HS_prof_out.close();
}

