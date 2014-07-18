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
  
  double dlogUstart = (logUstart_high - logUstart_low)/(double(N_Ustart)-1);
  
  cout << "logUstarts are: " << endl; 
  for(int i = 0; i<N_Ustart; i++)
  {
     cout << "\t" << logUstart_low+double(i)*dlogUstart;
  }
  cout << endl;

  cout << "Ustarts are: " << endl; 
  for(int i = 0; i<N_Ustart; i++)
  {
     cout << "\t" << pow(10,logUstart_low+double(i)*dlogUstart);
  }
  cout << endl;
  
  // now the parameters for the ratio between U start and Ufinal
  double logratio_low;
  double logratio_high;
  int n_ratio;
  
  infile >> logratio_low >> logratio_high >> n_ratio;
  infile.close();
  
  double dlogR = (logratio_high - logratio_low)/(double(n_ratio)-1);
  
  cout << "logR are: " << endl; 
  for(int i = 0; i<n_ratio; i++)
  {
     cout << "\t" << logratio_low+double(i)*dlogR;
  }  
  cout << endl;
  
  cout << "R are: " << endl;
  for(int i = 0; i<n_ratio; i++)
  {
     cout << "\t" << pow(10,logratio_low+double(i)*dlogR);
  }
  cout << endl;   

	// create a hillslope
	OneDImplicitHillslope Hillslope;

	// set the tolerance for the iterations
	//double tolerance = 0.000001;

  // parameters for each iteration
  double Uhat_start;
  double Uhat_final;
  double Uhat_ratio;
  double max_Uhat = 30;  
  vector<double> step_change_data;


  
  // open the outfile2
  string data_outname1 =  param_name+"_data1.data";
  string data_out_fname1 = pathname+data_outname1;
  ofstream data_out1;
  data_out1.open(data_out_fname1.c_str());
 
  string data_outname2 =  param_name+"_data2.data";
  string data_out_fname2 = pathname+data_outname2;
  ofstream data_out2;
  data_out2.open(data_out_fname2.c_str()); 

  string data_outname3 =  param_name+"_data3.data";
  string data_out_fname3 = pathname+data_outname3;
  ofstream data_out3;
  data_out1.open(data_out_fname3.c_str()); 

  string data_outname4 =  param_name+"_data4.data";
  string data_out_fname4 = pathname+data_outname4;
  ofstream data_out4;
  data_out4.open(data_out_fname4.c_str()); 

  // print header information
  data_out1 << "Uhat_start\tUhat_ratio\tlogUhat_ratio\tUhat_final\tMax_err\tEhat_at_maxerr\tThat_at_maxerr" << endl;
  data_out2 << "Uhat_start\tUhat_ratio\tlogUhat_ratio\tUhat_final\tMax_err\tEhat_at_maxerr\tThat_at_maxerr" << endl;
  data_out3 << "Uhat_start\tUhat_ratio\tlogUhat_ratio\tUhat_final\tMax_err\tEhat_at_maxerr\tThat_at_maxerr" << endl;
  data_out4 << "Uhat_start\tUhat_ratio\tlogUhat_ratio\tUhat_final\tMax_err\tEhat_at_maxerr\tThat_at_maxerr" << endl;
  
  // loop through the starting Uis
  for (int Ui = 0; Ui<N_Ustart; Ui++)
  {
    for (int Ufi = 0; Ufi< n_ratio; Ufi++)
    {
      cout << "Ui: " << Ui+1 << " of " << N_Ustart 
           << " and Ufi: " << Ufi+1 << " of " << n_ratio << endl;
           
      // get the start   Uhat
      Uhat_start = pow(10,logUstart_low+double(Ui)*dlogUstart);
      Uhat_ratio = pow(10,logratio_low+double(Ufi)*dlogR);
      Uhat_final = Uhat_start*Uhat_ratio;
      
      // make sure you don't get a very large Uhat
      if(Uhat_final <= max_Uhat)
      {
        // now run the model
        step_change_data = Hillslope.Max_Rstart_err_after_step_change(
                                     Uhat_start, Uhat_final);
                                     
        // print these results to file
        if (Ui == 0)
        {
          data_out1 << Uhat_start << "\t" << Uhat_ratio << "\t" 
                    << logratio_low+double(Ufi)*dlogR << "\t"
                    << Uhat_final << "\t"
                    << step_change_data[0] << "\t" << step_change_data[1] << "\t" 
                    << step_change_data[2] << endl;
        }
        else if (Ui == 1)
        {
          data_out2 << Uhat_start << "\t" << Uhat_ratio << "\t" 
                    << logratio_low+double(Ufi)*dlogR << "\t"
                    << Uhat_final << "\t"
                    << step_change_data[0] << "\t" << step_change_data[1] << "\t" 
                    << step_change_data[2] << endl;
        } 
        else if (Ui == 2) 
        {
          data_out2 << Uhat_start << "\t" << Uhat_ratio << "\t" 
                    << logratio_low+double(Ufi)*dlogR << "\t"
                    << Uhat_final << "\t"
                    << step_change_data[0] << "\t" << step_change_data[1] << "\t" 
                    << step_change_data[2] << endl;
        }               
        else if (Ui == 3) 
        {
          data_out2 << Uhat_start << "\t" << Uhat_ratio << "\t" 
                    << logratio_low+double(Ufi)*dlogR << "\t"
                    << Uhat_final << "\t"
                    << step_change_data[0] << "\t" << step_change_data[1] << "\t" 
                    << step_change_data[2] << endl;
        } 
      }
    }
  }
  data_out.close();
  


}

