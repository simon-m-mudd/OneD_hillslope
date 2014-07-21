//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// oneD_implicit_nonlinear_hillslope.cpp
//
// an object for solving nonlinear hillslope evolution implicitly
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef OneDImplicitHillslope_H
#define OneDImplicitHillslope_H

class OneDImplicitHillslope
{
	public:
		OneDImplicitHillslope()				{ create(); }
		OneDImplicitHillslope(double tp_temp,double Up_temp,double Uw_temp)
											{ create(tp_temp, Up_temp, Uw_temp); }
		OneDImplicitHillslope(double tp_temp,double Up_temp,double Uw_temp, double dx)
											{ create(tp_temp, Up_temp, Uw_temp, dx); }
		OneDImplicitHillslope(int tn_nodes, int tridgetop_node, double tt_hat_peak, double tU_hat_peak,
							double tU_hat_width, Array1D<double> tzeta_hat, Array1D<double> tzeta_last_timestep,
							Array1D<double> tzeta_intermediate, Array1D<double> tf, Array2D<double> tCoeff_matrix,
							vector<double> tx_hat,vector<double> tA_hat_denom, vector<double> tB_hat_denom,
							vector<double> tA_slope_denom2, vector<double> tB_slope_denom2, 
              double tD, double tS_c, double tL_H, vector<double> tx, 
              vector<double> tzeta_dimen,vector<double> tzeta_lts_dimen, double trho_ratio)
								{ create(tn_nodes, tridgetop_node,tt_hat_peak, tU_hat_peak,tU_hat_width,
								tzeta_hat,tzeta_last_timestep,tzeta_intermediate, tf, tCoeff_matrix,
								tx_hat, tA_hat_denom, tB_hat_denom, tA_slope_denom2,tB_slope_denom2, 
                tD, tS_c, tL_H, tx, tzeta_dimen, tzeta_lts_dimen, trho_ratio); }

    // set the dimensional parameters
    void set_D(double new_D)  { D = new_D; }
    void set_S_c(double new_S_c) {S_c = new_S_c; }
    void set_L_H(double new_L_H);
    void set_rho_ratio(double new_rho_ratio) {rho_ratio = new_rho_ratio; }

    // check if the dimensional parameters have been set
    bool check_if_dimensional_parameters_set();
    
    // some tools for moving between dimensional and nondimensional space
    double dimensional_time_from_that(double dt_hat); 
    double dimensional_uplift_from_Uhat(double U_hat);
    vector<double> dimensional_zeta_from_zeta_hat();
    double U_hat_from_dimensional_U(double U);
       
    // this runs a hillslope timestep in dimensional and wraps it in dimensional
    // space
    // returns dimensional dt
    double run_dimensional_hillslope_timestep(double& dt_hat, double& t_hat_ime,
                                   double & t_ime, double uplift, double tolerance);

    // this calculates change in zeta at the ridgetop. Only works if zeta
    // has been initialised
    double get_dz_ridgetop();
    
    // this gets the current dimensional zeta at the ridgetop
    double get_current_ridgetop_dimensional_zeta();

		// reset hillslope to a new flat surface using the same spatial discritization
		void reset_hillslope(double tp_temp, double Up_temp, double Uw_temp);

		// this takes a time vector and calculates the E* and R* along this time series
		void run_based_on_data_spacing(vector<double> t_hat_data, vector<double>& E_star_modelled,
							 vector<double>& R_star_modelled, double tolerance);
		void advance_to_next_time_interval_gaussian_uplift(double& dt_hat, double& t_ime,
		                                            double target_time, double tolerance);
													// this advnaces the model to some
													// set time
		
    void hillslope_timestep(double& dt_hat, double& t_ime, double U_hat,double tolerance);
													// this does a single timestep in the
													// evolution of the hillslope
													// returns a new timestep and the new t_ime
		
    void hillslope_timestep_gaussian_uplift(double& dt_hat, double& t_ime, double tolerance);
													// this does a single timestep but solves the
													// governing equations using U_hat at
													// the future timestep determined from a
													// gaussian function
		
    int hillslope_iterator(double dt_hat, double U_hat, double tolerance);
													// this adjusts the timestep
		
    void solve_for_zeta_intermediate(double dt_hat, double U_hat);
													// this solves a timestep, replacing the
													// data array zeta_intermediate with
													// a solved zeta

		double get_zeta_rootsquare_error();			// gets the root of the square of error
													// between zeta nad zeta intermediate

		// accessor functions (used for = operator)
		int get_n_nodes() const						{return n_nodes;}
		int get_ridgetop_nod() const				{return ridgetop_node;}
		double get_t_hat_peak() const				{return t_hat_peak;}
		double get_U_hat_peak() const				{return U_hat_peak;}
		double get_U_hat_width() const				{return U_hat_width;}
		double get_rho_ratio() const      { return rho_ratio; }

		Array1D<double> get_zeta_hat() const		{ Array1D<double> yo = zeta_hat.copy();
											  		return yo; }
		Array1D<double> get_zeta_last_timestep() const
													{ Array1D<double> yo = zeta_last_timestep.copy();
											 			 return yo; }
		Array1D<double> get_zeta_intermediate() const
													{ Array1D<double> yo = zeta_intermediate.copy();
											  			return yo; }
		Array1D<double> get_f() const				{ Array1D<double> yo = f.copy();
											  			return yo; }
		Array2D<double> get_Coeff_matrix() const	{ Array2D<double> yo = Coeff_matrix.copy();
											  			return yo; }
		vector<double> get_x_hat() const			{ return x_hat; }
		vector<double> get_A_hat_denom() const		{ return A_hat_denom; }
		vector<double> get_B_hat_denom() const		{ return B_hat_denom; }
		vector<double> get_A_slope_denom2() const	{ return A_slope_denom2; }
		vector<double> get_B_slope_denom2() const	{ return B_slope_denom2; }
		double get_D() const  { return D; }
    double get_S_c() const { return S_c; }
    double get_L_H() const { return L_H; }
    vector<double> get_x() const { return x; }
    vector<double> get_zeta_dimen() const { return zeta_dimen; }
    vector<double> get_zeta_lts_dimen() const { return zeta_lts_dimen; }

		// functions for getting the anayltical solutions to the governing equations
		// these are used to test the numerical method
		Array1D<double> get_analytical_SS_Cstar(double U_star);
		Array1D<double> get_analytical_SS_Sstar(double U_star);
		Array1D<double> get_analytical_SS_zstar(double U_star);
		void print_analytical_SS(double U_star);

		// this runs the model to steady state and then tests it against
		// the analytical solutions. Done for bug checking purposes
		double test_steady(double dt_hat, double tolerance);

		// these get the calculated E* and R* values
		double calculate_R_star();
		double analytical_R_star(double U_star);
		double calculate_E_star();

		// function for caluculating uplift based on
		// a gaussian uplift curve
		double gaussian_uplift(double t_hat);

		OneDImplicitHillslope& operator=(const OneDImplicitHillslope& ODHS);

    // sets zeta to be the analytical solution for U_hat
    void set_analytical_steady(double U_hat);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // Max_Rstar_err_after_step_change
  // This tests to see what the maximum difference between predicted and measured 
  // R_star will be given a pulse of uplift
  //
  // returns a vector for holding the data
  // the first element is the maximum error in Rstar, 
  // the second element is the Estar where this occurs
  // the third element is the time after perturbation when the max occurs
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  vector<double> Max_Rstart_err_after_step_change(double Uhat_start, double Uhat_end);


	protected:

		int n_nodes;				    // number of nodes
		
		// data members for dimensional version
		double D;                       // diffusivity m^2/yr
		double S_c;                     // critical slope
		double L_H;                     // hillslope lenght (m)
		double rho_ratio;               // ratio between rock and soil density
		vector<double> x;               // distance from divide
		vector<double> zeta_dimen;      // dimensional elevation
		vector<double> zeta_lts_dimen;  // dimensional elevation from the last timestep
		
		int ridgetop_node;			// the node number of the ridgetop
		double t_hat_peak;			// time of peak uplift
		double U_hat_peak;			// maximum uplift
		double U_hat_width;			// width of gaussian uplift signal

		Array1D<double> zeta_hat;	// the elevation of the hillslope surface
		Array1D<double> zeta_last_timestep;
										// zeta (hat, the elevation) at the last timestep
		Array1D<double> zeta_intermediate;
										// intermediate zeta. This tests the convergence
										// of the matrix inversion
		Array1D<double> f;				// vector in system Coeff_matrix*zeta@t+1 = f
		Array2D<double> Coeff_matrix;
		vector<double> x_hat;			// the location of the nodes
		vector<double> A_hat_denom;		// coefficients for solving the forward
										// time problem: saves some computational time
		vector<double> B_hat_denom;		// coefficients for solving the forward
										// time problem: saves some computational time
		vector<double> A_slope_denom2;	// coefficients for solving the forward
										// time problem: saves some computational time
		vector<double> B_slope_denom2;	// coefficients for solving the forward
										// time problem: saves some computational time



	private:
		void create();
		void create(double tp_temp,double Up_temp,double Uw_temp);
		void create(double tp_temp,double Up_temp,double Uw_temp, double dx);
		void create(int tn_nodes, int tridgetop_node, double tt_hat_peak, double tU_hat_peak,
					double tU_hat_width, Array1D<double> zeta_hat, Array1D<double> tzeta_last_timestep,
					Array1D<double> tzeta_intermediate, Array1D<double> tf, Array2D<double> tCoeff_matrix,
					vector<double> tx_hat,vector<double> tA_hat_denom, vector<double> tB_hat_denom,
					vector<double> tA_slope_denom2, vector<double> tB_slope_denom2, 
          double tD, double tS_c, double tL_H, vector<double> tx, 
          vector<double> tzeta_dimen,vector<double> tzeta_lts_dimen, double rho_ratio);
};

#endif

