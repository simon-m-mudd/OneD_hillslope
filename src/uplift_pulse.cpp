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
  float_default_map["start_estar"] = 0.2;
  help_map["start_estar"] = { "float","0.2","Starting E* value.","Does what is says on the tin."};

  float_default_map["end_estar"] = 10.0;
  help_map["end_estar"] = {  "float","10.0","Ending E* value.","Does what is says on the tin."};
  
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
    cout << "./lsdtt-basic-metrics-README.csv" << endl;
    string help_prefix = "lsdtt-basic-metrics-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }
  
  // Now print the parameters for bug checking
  cout << "PRINT THE PARAMETERS..." << endl;
  LSDPP.print_parameters();

  // create a hillslope
  OneDImplicitHillslope thisHillslope;

  double start_estar = this_float_map["start_estar"];
  double end_estar = this_float_map["end_estar"];


  // set to steady state
  thisHillslope.set_analytical_steady(start_estar);


  string uscore = "_";
  string data_ext = ".EsRsdata";
  string this_outfile = "yoyoma"+data_ext;
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

