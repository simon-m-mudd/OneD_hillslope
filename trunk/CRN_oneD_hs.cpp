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
#include <list>
#include "OneDImplicitHillslope.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDParticle.hpp"
#include "CRN_funcs.hpp"
using namespace std;

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

  // read data elements from infile
  string stringin;
  double start_erosion;
  double rho_r;
  double S_c;
  double start_D;
  double particle_spacing;
  double start_depth;
  double L_H;

  // read in data from parameter file
  infile >> stringin >> start_erosion >> stringin >> rho_r;
  infile >> stringin >> S_c >> stringin >> start_D;
  infile >> stringin >> particle_spacing >> stringin >> start_depth;
  infile >> stringin >> L_H;

  // get the effective erosion (in g/cm^2/yr)
  double eff_eros_rate = start_erosion*rho_r; // YOU NEED TO UPDATE THIS WITH THE CORRECT UNIT CONVERSION

  // now set up a column 
  int start_type = 1;
  double startxloc = L_H;

  // set up a placeholder for zeta at the moment
  double zeta = 10;

  // make a LSDCRNParameters object. This holds parameters for 
  // production of cosmogenic nuclides
  LSDCRNParameters CRN_param;

  list<LSDCRNParticle> CRN_plist = initiate_SS_cosmo_column(start_type, startxloc,
		      start_depth, particle_spacing, 
		      zeta, rho_r, eff_eros_rate,
		      CRN_param);


}

