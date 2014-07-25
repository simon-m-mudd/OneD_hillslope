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
#include "LSDStatsTools.hpp"
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
  double rho_r;
  double S_c;
  double particle_spacing;
  double start_depth;
  double L_H;

  double start_erosion, start_D, start_Dratio;
  int N_er,N_sD, N_Dr;
  double dstart_erosion, dstart_D, dDratio;

  // read in data from parameter file
  infile >> stringin >> start_erosion >> stringin >> N_er >> stringin >> dstart_erosion;
  infile >> stringin >> start_D >> stringin >> N_sD >> stringin >> dstart_D;
  infile >> stringin >> start_Dratio >> stringin >> N_Dr >> stringin >> dDratio;
  infile >> stringin >> rho_r >> stringin >> S_c;
  infile >> stringin >> particle_spacing >> stringin >> start_depth;
  infile >> stringin >> L_H;

  cout << "start erosion: " << start_erosion << " N erosion: " << N_er << endl;
  cout << "start D: " << start_D << " N D:" << N_sD << endl;
  cout << "start D ratop: " << start_Dratio << " N Dratio:" << N_Dr << endl;
  cout << "rho_r: " << rho_r << " S_c: " << S_c << " p_spacing: " << particle_spacing
       << "start_depth: " << start_depth << " L_H: " << L_H << endl;

  // set the rho ratio to 1: we are not going to consider density changes
  double rho_ratio = 1;

  // now set up a column
  int start_type = 1;
  double startxloc = L_H;

  // make a LSDCRNParameters object. This holds parameters for
  // production of cosmogenic nuclides
  LSDCRNParameters CRN_param;

  // set the parameters to be nucleonic only
  CRN_param.set_Neutron_only_parameters();

  // now start a hillslope model: we start with a step change in D
  // set time parameters. Time is scaled by
  double dt;                  // time spacing in years
  double t_ime;               // time in years
  double dt_hat;              // dimensionless time spacing
  double endTime = 2;         // dimensionless end time
  double t_ime_hat;           // dimensionless time
  double dimensional_dt;      // another time spacing, this for particle tracking
  double tolerance = 0.000001;
  double that_print_spacing = 0.01;  // this is time spacing for printing to file
  double next_print_time = that_print_spacing;

  double current_D;           // current diffusivity
  double current_Uhat;        // the dimensionless uplift
  double elev_uplift;         // the length of uplift during a timestep
  double zeta_rt;             // elevation at the ridgetop
  double zeta_rt_old;         // old elevation at the ridgetop
  double zeta_rt_uplift_corrected;  // this is used to calculate effective erosion
  double dzeta;               // change in surface elevation
  double app_D_10Be;               // apparent diffusivity based on 10Be
  double D_err_10Be;               // error in diffusibity based on 10Be
  double app_D_14C;               // apparent diffusivity based on 14C
  double D_err_14C;               // error in diffusibity based on 14C
  vector<double> app_eros;

  double eff_eros_rate;       // the effective erosion rate in g cm^-2 yr^-1
  double U_hat_start;         // dimensionless uplift at the start

  double current_erosion;     // the current erosion rate (m/yr)
  double current_start_D;     // current D in m^2/yr
  double current_Dratio;      // ratio  of the starting D and current d

  list<LSDCRNParticle> eroded_list;  // this is just a placeholder: nothing is done with it

  // enter the loop for erosion rate
  for(int er_i = 0; er_i< N_er; er_i++)
  {
    current_erosion = pow(2,start_erosion+double(er_i)*dstart_erosion)*0.0001;
    cout << "\n\nCurrent erosion is: " << current_erosion << endl;

    // enter the loop for starting D
    for (int sD_i = 0; sD_i< N_sD; sD_i++ )
    {
      current_start_D = pow(10,start_D+double(sD_i)*dstart_D);
      cout << "Current starting diffusivity is: " << current_start_D << endl;

      // enter the loop for D_ratio
      for (int Dr_i = 0; Dr_i< N_Dr; Dr_i++)
      {

        current_Dratio = pow(2,start_Dratio+double(Dr_i)*dDratio);
        cout << "Current D ratio is: " << current_Dratio << endl;

        // get the effective erosion (in g/cm^2/yr)
        // for this to work properly rho_r needs to be in kg/m^3
        eff_eros_rate = current_erosion*rho_r/10;

        // convert the erosion rate to dimensionless units
        U_hat_start = (current_erosion*rho_ratio*2*L_H)/(current_start_D*S_c);
        cout << "Starting dimensionless U: " << U_hat_start
             << " and eff_erosion: " << eff_eros_rate << endl;

        // deal with the filenames
        string undersc = "_";
        string CRN_suffix = "CRNpdata_v3";
        string CRN_runname = dot+itoa(er_i)+undersc+itoa(sD_i)+undersc+itoa(Dr_i);
        string CRN_of_name = pathname+param_name+CRN_runname+dot+CRN_suffix;

        ofstream CRNdataout;
        CRNdataout.open(CRN_of_name.c_str());

        CRNdataout <<  current_erosion << "\t" << current_start_D << "\t"
                   <<  current_Dratio << endl;

        // create a hillslope
        OneDImplicitHillslope Hillslope;

        // set the dimensional parameters in the hillslope
        Hillslope.set_D(current_start_D);
        Hillslope.set_S_c(S_c);
        Hillslope.set_L_H(L_H);
        Hillslope.set_rho_ratio(rho_ratio);

        // set up a steady nondimensional hillslope
        Hillslope.set_analytical_steady(U_hat_start);
        Hillslope.populate_dimensional_zeta();

        // get the erosion properties
        zeta_rt = Hillslope.get_current_ridgetop_dimensional_zeta();

        // start with a column of particles at steady state
        list<LSDCRNParticle> CRN_plist = initiate_SS_cosmo_column(start_type, startxloc,
		        start_depth, particle_spacing,
		        zeta_rt, rho_r, eff_eros_rate,
		        CRN_param);

        // set the new D
        current_D = current_start_D*current_Dratio;
        Hillslope.set_D(current_D);
        current_Uhat = Hillslope.U_hat_from_dimensional_U(current_erosion);

        // reset time
        t_ime = 0;
        t_ime_hat = 0;
        dt = 1;            // this will get replaced before getting used
        dt_hat = 0.00001;   // dimensionless time spacing
        next_print_time = 0;  // reset the next time to print
        // run the loop

        while (t_ime_hat < endTime)
        //for(int i = 0; i<200; i++)
        {
          // do a timestep
          dimensional_dt = Hillslope.run_dimensional_hillslope_timestep(dt_hat,
                            t_ime_hat,t_ime,current_erosion,tolerance);

          // need to update the zeta location of all the particles. CRN_funcs uses an
          // absolute coordinate system so in the advective coordinate system of the dimensional
          // hillslope we need to update the zeta locations
          elev_uplift = current_erosion*dimensional_dt;



          // now get the zeta new and zeta_old
          zeta_rt = Hillslope.get_current_ridgetop_dimensional_zeta();
          dzeta = Hillslope.get_dz_ridgetop();
          zeta_rt_old = zeta_rt - dzeta;
          zeta_rt_uplift_corrected =  zeta_rt - elev_uplift;

          // update the particles
          eroded_list =  update_CRN_list_eros_limit_3CRN(CRN_plist,dimensional_dt,
                           rho_r,start_type,start_depth,startxloc,
      	                   zeta_rt_old, zeta_rt_uplift_corrected,
                           particle_spacing, CRN_param);

          // now update the particles based on this uplift
          // all of the particles will move toward the surface by a distance elev_uplift
          update_list_z_for_advective_coord_system(CRN_plist,elev_uplift);

          //print_particles_and_apparent_erosion_3CRN(CRN_plist,rho_r,dt,t_ime,
      	  //				    CRNdataout, CRN_param);

          // check to see if we print
	  //cout << "t_ime_hat: " << t_ime_hat << " and end_t_ime_hat is: " << endTime << endl;
	  //cout << "t_ime is: " << t_ime << " and time to next print is: " << next_print_time << endl;
          if(t_ime_hat >= next_print_time)
          {
            // first reset the next print time:
            next_print_time+=that_print_spacing;

            // check the apparent erosion rate and the apparent D
            app_eros = calculate_apparent_erosion_3CRN_neutron(CRN_plist,
      		             rho_r, CRN_param );
      	    app_D_10Be = -rho_ratio*app_eros[0]/
                    Hillslope.calculate_dimensional_ridgetop_curvature();
            app_D_14C = -rho_ratio*app_eros[2]/
                    Hillslope.calculate_dimensional_ridgetop_curvature();

            D_err_10Be = fabs((app_D_10Be-current_D)/current_D);
	    D_err_14C =  fabs((app_D_14C-current_D)/current_D);

            //cout << "Time is: " << t_ime << " and ridgetop zeta is "
            //     << Hillslope.get_current_ridgetop_dimensional_zeta()
            //     << " years, app E (10Be) is: " << app_eros[0]
            //     << " , app D is: " << app_D << " and D_err is: " << D_err << endl;


            CRNdataout << t_ime << "\t" << t_ime_hat << "\t"
                 << Hillslope.get_current_ridgetop_dimensional_zeta()
		       << "\t" << app_eros[0] << "\t" << app_eros[2]
		       << "\t" << app_D_10Be << "\t" << app_D_14C << "\t"
	               << "\t" << D_err_10Be << "\t" << D_err_14C << endl;
            // now print to file
            //print_particles_and_apparent_erosion_3CRN(CRN_plist,rho_r,dt,t_ime,
      			//		    CRNdataout, CRN_param);
          }              // end printing if statement


        }                // end model run
        CRNdataout.close();
      }                  // end D_ratio loop
    }                    // end starting_D loop
  }                      // end loop for erosion rate

}

