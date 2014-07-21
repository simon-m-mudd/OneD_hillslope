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

  // set the rho ratio to 1: we are not going to consider density changes
  double rho_ratio = 1;

  // get the effective erosion (in g/cm^2/yr)
  // for this to work properly rho_r needs to be in kg/m^3
  double eff_eros_rate = start_erosion*rho_r/10;
  
  cout << "Erosion rate is: " << start_erosion << " and effective erosion is: " 
       << eff_eros_rate << endl; 

  // convert the erosion rate to dimensionless units
  double U_hat_start = (start_erosion*rho_ratio*2*L_H)/(start_D*S_c);

  // now set up a column 
  int start_type = 1;
  double startxloc = L_H;

  // make a LSDCRNParameters object. This holds parameters for 
  // production of cosmogenic nuclides
  LSDCRNParameters CRN_param;

  // set the parameters to be nucleonic only
  CRN_param.set_Neutron_only_parameters();

  // now get a steady hillslope
  // create a hillslope
  OneDImplicitHillslope Hillslope;

  // set the dimensional parameters in the hillslope
  Hillslope.set_D(start_D);
  Hillslope.set_S_c(S_c);
  Hillslope.set_L_H(L_H);
  Hillslope.set_rho_ratio(rho_ratio);

  // get the dimensionless U_hat this way
  double U_hat_test = Hillslope.U_hat_from_dimensional_U(start_erosion);

  cout << "Test of dimensionalising U_hat, directly: " << U_hat_start 
    << " and from object: " << U_hat_test << endl;

  // set up a steady nondimensional hillslope
  Hillslope.set_analytical_steady(U_hat_start);

  // get the ridgetop node
  int ridgetop_node = Hillslope.get_ridgetop_nod();

  // get the dimensional zeta
  vector<double> zeta_dimen = Hillslope.dimensional_zeta_from_zeta_hat();

  // print the elevation at the ridgetop
  double ridgetop_zeta = zeta_dimen[ridgetop_node];
  cout << "The elevation at the ridgetop is: " << ridgetop_zeta << endl;

  // start with a column of particles at steady state
  list<LSDCRNParticle> CRN_plist = initiate_SS_cosmo_column(start_type, startxloc,
		      start_depth, particle_spacing, 
		      ridgetop_zeta, rho_r, eff_eros_rate,
		      CRN_param);

  // now print these particles to file to see if they are at steady state
  string CRN_suffix = "CRNpdata";
  string CRN_of_name = pathname+param_name+dot+CRN_suffix;

  ofstream CRNdataout;
  CRNdataout.open(CRN_of_name.c_str());

  double dt = 0.1;
  double t_ime = 0;
  print_particles_and_apparent_erosion_3CRN(CRN_plist,rho_r,dt,t_ime,
					    CRNdataout, CRN_param);
					    
	// now start a hillslope model: we start with a step change in D 	
  // set time parameters. Time is scaled by 
  double dt_hat =0.01;
  double endTime = 4;
  double t_ime_hat = 0;
  double dimensional_dt = 0.1;
  double tolerance = 0.000001;
  double t_ime_spacing = 10000;
  double next_print_time = t_ime_spacing;

  double D_ratio = 0.5;
  double current_D;
  double current_Uhat;
  double elev_uplift;
  double zeta_rt;
  double zeta_rt_old;
  double dzeta;
  
  list<LSDCRNParticle> eroded_list;
  //double app_10Be_eros;
  //double app_14C_eros;
  //double app_21Ne_eros;
  
  double app_D;  
  double D_err;
  vector<double> app_eros;
  
  // set the new D
  current_D = start_D*D_ratio;
  Hillslope.set_D(current_D);
  current_Uhat = Hillslope.U_hat_from_dimensional_U(start_erosion); 
  // run the loop
  while (t_ime_hat < endTime)
  {
    
    // do a timestep
    dimensional_dt = Hillslope.run_dimensional_hillslope_timestep(dt_hat,t_ime_hat,t_ime,start_erosion,tolerance);
    
    // need to update the zeta location of all the particles. CRN_funcs uses an
    // absolute coordinate system so in the advective coordinate system of the dimensional
    // hillslope we need to update the zeta locations
    elev_uplift = start_erosion*dimensional_dt;
    
    // now update the particles based on this uplift
    update_list_z_for_advective_coord_system(CRN_plist,elev_uplift);
    
    // now get the zeta new and zeta_old
    zeta_rt = Hillslope.get_current_ridgetop_dimensional_zeta();
    dzeta = Hillslope.get_dz_ridgetop();
    zeta_rt_old = zeta_rt - dzeta;
    						
    // update the particles
    eroded_list =  update_CRN_list_eros_limit_3CRN(CRN_plist,dimensional_dt, 
                     rho_r,start_type,start_depth,startxloc,
	                   zeta_rt_old, zeta_rt, particle_spacing, CRN_param);
    
    // check to see if we print
    if(t_ime >= next_print_time)
    {
      // first reset the next print time:
      next_print_time+=t_ime_spacing;
      
      // check the apparent erosion rate and the apparent D
      app_eros = calculate_apparent_erosion_3CRN_neutron(CRN_plist,
		             rho_r, CRN_param );
		  app_D = -rho_ratio*app_eros[0]/
              Hillslope.calculate_dimensional_ridgetop_curvature();
              
      D_err = fabs((app_D-current_D)/current_D);                   
		             
      cout << "Time is: " << t_ime << " years, app E (10Be) is: " << app_eros[0] 
           << " , app D is: " << app_D << " and D_err is: " << D_err << endl;
      // now print to file
      
    }
    
    
  }      
					    

}

