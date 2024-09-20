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
  OneDImplicitHillslope thisHillslope;

  double start_estar = 0.5;
  double end_estar = 10;


  // set to steady state
  thisHillslope.set_analytical_steady(start_estar);


  string uscore = "_";
  string data_ext = "EsRsdata";
  string this_outfile = pathname+param_name+uscore+dot+data_ext;
  cout << "Now doing timeseries, filename is: " << this_outfile << endl;
    
  ofstream EsRs_transient_out;
  EsRs_transient_out.open(this_outfile.c_str());

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


}

