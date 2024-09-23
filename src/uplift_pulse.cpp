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
#include "LSDParameterParser.hpp"
#include "LSDStatsTools.hpp"
using namespace std;


int main (int nNumberofArgs,char *argv[])
{

  string version_number = "0.2d";
  string citation = "https://doi.org/10.1002/esp.3923";

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the 1D hillslope tool!                   ||" << endl;
  cout << "|| Are you ready to run some hillslopes??              ||" << endl;
  cout << "|| You better be ready! Lets go! Whoooooo              ||" << endl;  
  cout << "=========================================================" << endl;


  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);
  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

  // Check if we are doing the version or the citation
  if(f_name == "lsdtt_citation.txt")
  {

    cout << endl << endl << endl << "==============================================" << endl;
    cout << "To cite this code, please use this citation: " << endl;
    cout << citation << endl;
    cout << "Copy this url to find the full citation." << endl;
    cout << "also see above for more detailed citation information." << endl;
    cout << "=========================================================" << endl;

    ofstream ofs;
    ofs.open("./onedhillslope-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    cout << endl << endl << endl << "==============================================" << endl;    
    cout << "This is onedhillslope version number " << version_number << endl;
    cout << "If the version contains a 'd' then you are using a development version." << endl;
    cout << "=========================================================" << endl;
    ofstream ofs;
    ofs.open("./onedhillslope-version.txt");
    ofs << version_number << endl;
    ofs.close();

    exit(0);
  }

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // this will contain the help file
  map< string, vector<string> > help_map;

  //==================================================================================
  //
  // .#####....####...#####....####...##...##..######..######..######..#####....####..
  // .##..##..##..##..##..##..##..##..###.###..##........##....##......##..##..##.....
  // .#####...######..#####...######..##.#.##..####......##....####....#####....####..
  // .##......##..##..##..##..##..##..##...##..##........##....##......##..##......##.
  // .##......##..##..##..##..##..##..##...##..######....##....######..##..##...####..
  //
  //=================================================================================
  float_default_map["start_estar"] = 0.5;
  help_map["start_estar"] = { "float","0.2","Starting E* value.","Does what is says on the tin."};

  float_default_map["end_estar"] = 2.0;
  help_map["end_estar"] = {  "float","10.0","Ending E* value.","Does what is says on the tin."};

  float_default_map["end_tstar"] = 0.5;
  help_map["end_tstar"] = {  "float","0.5","End of simulation in T* units. T* scales by response time to hillslope is mostly adjusted by 1","Complete adjustment is around 3 T* but mostly adjusted by 1 T*"};

  float_default_map["dx_hat"] = 0.1;
  help_map["dx_hat"] = {  "float","0.1","Spacing of nodes in dimensionless x","Note spacing is closer near the divide."};  

  float_default_map["dimensionless_print_profile_interval"] = 0.05;
  help_map["dimensionless_print_profile_interval"] = {  "float","0.05","Interval at which the profile prints in T* units. T* scales by response time to hillslope is mostly adjusted by 1","The program solves a variable timestep so the time of printing may not be exactly this number"};
  
  float_default_map["dimensionless_print_timeseries_interval"] = 0.005;
  help_map["dimensionless_print_timeseries_interval"] = {  "float","0.005","Interval at which the timeseries prints in T* units. T* scales by response time to hillslope is mostly adjusted by 1","The program solves a variable timestep so the time of printing may not be exactly this number"};
  
  float_default_map["L_H"] = 100.0;
  help_map["L_H"] = {  "float","10.0","Hillslope length in.","For dimensionalising results."};

  float_default_map["S_c"] = 0.75;
  help_map["S_c"] = {  "float","0.75","Critical slope. Unitless.","For dimensionalising results."};

  float_default_map["D"] = 0.001;
  help_map["D"] = {  "float","0.001","Diffusion coeffient in m^2/yr.","For dimensionalising results."};

  float_default_map["rho_ratio"] = 2.0;
  help_map["rho_ratio"] = {  "float","2","Ratio between the bedrock and soil density","For dimensionalising results."};

  //=========================================================================
  //
  //.#####....####...#####....####...##...##..######..######..######..#####..
  //.##..##..##..##..##..##..##..##..###.###..##........##....##......##..##.
  //.#####...######..#####...######..##.#.##..####......##....####....#####..
  //.##......##..##..##..##..##..##..##...##..##........##....##......##..##.
  //.##......##..##..##..##..##..##..##...##..######....##....######..##..##.
  //
  //..####...##..##..######...####...##..##...####..                         
  //.##..##..##..##..##......##..##..##.##...##.....                         
  //.##......######..####....##......####.....####..                         
  //.##..##..##..##..##......##..##..##.##.......##.                         
  //..####...##..##..######...####...##..##...####..                         
  //============================================================================
  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./oned-hillslope-README.csv" << endl;
    string help_prefix = "oned-hillslope-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }
  
  // Now print the parameters for bug checking
  cout << "PRINT THE PARAMETERS..." << endl;
  LSDPP.print_parameters();

  // create a hillslope
  double dx_hat = 0.05;
  OneDImplicitHillslope thisHillslope(dx_hat);

  // some parameter that at the moment we are not changing  
  double dt_hat = 0.0005;  // this gives a reasonable starting point for iterations
  double t_ime_hat = 0;   // we always start at time 0
  double tolerance = 0.00001; // this sets convergence of model

  // user defined parameters
  double start_estar = this_float_map["start_estar"];
  double end_estar = this_float_map["end_estar"];
  double end_time_hat = double(float_default_map["end_tstar"]);
  double timeseries_print_spacing = double(float_default_map["dimensionless_print_timeseries_interval"]);
  double profile_print_spacing =  double(float_default_map["dimensionless_print_profile_interval"]);
  double next_timeseries_print = timeseries_print_spacing;
  double next_profile_print = profile_print_spacing;


  // set the hillslope to steady state
  thisHillslope.set_analytical_steady(start_estar);

  // set the dimensional values
  thisHillslope.set_D(double(this_float_map["D"]));
  thisHillslope.set_S_c(double(this_float_map["S_c"]));
  thisHillslope.set_L_H(double(this_float_map["L_H"]));
  thisHillslope.set_rho_ratio(double(this_float_map["rho_ratio"]));


  // set up outfiles
  string timeseries_outfile = "timeseries.csv";
  string profile_outfile = "profile.csv";
  cout << "Now doing timeseries and profile, filenames are: " << timeseries_outfile << " and " << profile_outfile << endl;
    
  ofstream EsRs_transient_out;
  EsRs_transient_out.open(timeseries_outfile.c_str());

  ofstream profile_transient_out;
  profile_transient_out.open(profile_outfile.c_str());

  EsRs_transient_out << "t_hat,E*,R*,analytical_R*" << endl;

  // get the x location string
  string x_string = thisHillslope.print_comma_delimited_xhat_string();


  // printing for the profile
  profile_transient_out << "t_hat,"+x_string << endl;
  // print the initial condition
  string zeta_hat_str = thisHillslope.print_comma_delimited_zetahat_string();
  profile_transient_out << "0," << zeta_hat_str << endl;


  cout << "Starting hillslope loop" << endl;
  while (t_ime_hat < end_time_hat)
  {
    thisHillslope.hillslope_timestep(dt_hat, t_ime_hat, end_estar, tolerance);
    
    if(t_ime_hat >= next_timeseries_print)
    {
      //cout << "Printing to file, time is: " << t_ime_hat << " and end time: " << end_time_hat << endl;
      double this_E_star = thisHillslope.calculate_E_star();
      EsRs_transient_out << t_ime_hat << ","
                          << thisHillslope.calculate_E_star() << "," 
                          << thisHillslope.calculate_R_star() << ","
                          << thisHillslope.analytical_R_star(this_E_star) << endl;
      next_timeseries_print += timeseries_print_spacing;                  
    }

    if(t_ime_hat >= next_profile_print)
    {
      cout << "Printing your profile, the time is: " << t_ime_hat << " and end time: " << end_time_hat << endl;
      string zeta_hat_str = thisHillslope.print_comma_delimited_zetahat_string();
      profile_transient_out << t_ime_hat << "," << zeta_hat_str << endl;
      next_profile_print += profile_print_spacing; 
    }

  }


}

