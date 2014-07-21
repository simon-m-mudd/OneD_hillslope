#ifndef CRN_FUNCS_H
#define CRN_FUNCS_H


#include <list>
#include <fstream>
#include "LSDParticle.hpp"
#include "LSDCRNParameters.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

// initiate a list with steady state particles
list<LSDCRNParticle> initiate_SS_cosmo_column(int start_type, double startxLoc,
		      double start_depth, double particle_spacing, 
		      double zeta, double rho_r, double eff_eros_rate,
		      LSDCRNParameters& CRN_param);

// update_list
// this is the function that inserts a particle into a list
void insert_part_into_CRN_list(list<LSDCRNParticle>& CRN_list,
						 int start_type, double startxLoc,
						 double start_depth,double zeta,
						 double rho_r);

// this is a function that takes an old and new zeta, calculates
// an erosion rate, and updates the CRN concntrations
// of all the particles in the list
// in addition, particles are popped from the list (eroded)
// if their final position is above the new surface,
// and if the lowermost particle is above the start depth by a
// distance greater than or equal to particle_spacing
// a new particle is inserted into the column
list<LSDCRNParticle> update_CRN_list(list<LSDCRNParticle>& CRN_list,
	double dt, double rho_r,
	int start_type, double start_depth,
	double startxLoc,
double zeta_old,double zeta_new,
	     double particle_spacing, LSDCRNParameters& CRN_param);

list<LSDCRNParticle> update_CRN_list_eros_limit(list<LSDCRNParticle>& CRN_list,
		double dt, double rho_r,
		int start_type, double start_depth,
		double startxLoc,
		double zeta_old,double zeta_new,
		double particle_spacing, LSDCRNParameters& CRN_param);

list<LSDCRNParticle> update_CRN_list_eros_limit_3CRN(list<LSDCRNParticle>& CRN_list,
		double dt, double rho_r,
		int start_type, double start_depth,
		double startxLoc,
		double zeta_old,double zeta_new,
						     double particle_spacing, LSDCRNParameters& CRN_param);
list<LSDCRNParticle> update_CRN_list_eros_limit_3CRN_neutron(list<LSDCRNParticle>& CRN_list,
		double dt, double rho_r,
		int start_type, double start_depth,
		double startxLoc,
		double zeta_old,double zeta_new,
	     double particle_spacing, LSDCRNParameters& CRN_param);


// this calculates the apparent erosion rates, and returns them in a vector
// erosion rates are in m/yr
// [0] = 10Be
// [1] = 21Ne
// [2] = 14C
vector<double> calculate_apparent_erosion_3CRN_neutron(list<LSDCRNParticle>& CRN_list,
		double rho_r, LSDCRNParameters& CRN_param );





// this function updates lists that are distributed in space
// effectively a distributed version of update_CRN_list
void update_CRN_list_vec( vector< list<LSDCRNParticle> >& CRN_list_vec,
							Array2D<double>& zeta,
							Array2D<double>& zeta_old,
							vector<int>& row_list,
							vector<int>& col_list,
							double dt,
							double rho_r,
							int start_type, double start_depth,
							double startxLoc,
			  double particle_spacing, LSDCRNParameters& CRN_param);

vector< list<LSDCRNParticle> > update_CRN_list_vec_with_eroded
		 				  ( vector< list<LSDCRNParticle> >& CRN_list_vec,
							Array2D<double>& zeta,
							Array2D<double>& zeta_old,
							vector<int>& row_list,
							vector<int>& col_list,
							double dt, double rho_r,
							int start_type, double start_depth,
							double startxLoc,
						    double particle_spacing, LSDCRNParameters& CRN_param);

vector< list<LSDCRNParticle> > update_CRN_list_vec_with_eroded_3CRN
		 				  ( vector< list<LSDCRNParticle> >& CRN_list_vec,
							Array2D<double>& zeta,
							Array2D<double>& zeta_old,
							vector<int>& row_list,
							vector<int>& col_list,
							double dt, double rho_r,
							int start_type, double start_depth,
							double startxLoc,
							double particle_spacing);
vector< list<LSDCRNParticle> > update_CRN_list_vec_with_eroded_3CRN_neutron
		 				  ( vector< list<LSDCRNParticle> >& CRN_list_vec,
							Array2D<double>& zeta,
							Array2D<double>& zeta_old,
							vector<int>& row_list,
							vector<int>& col_list,
							double dt, double rho_r,
							int start_type, double start_depth,
							double startxLoc,
						    double particle_spacing, LSDCRNParameters& CRN_param);


// this function initializes a CRN list vec
void initialize_CRN_list_vec( vector< list<LSDCRNParticle> >& CRN_list_vec,
							Array2D<double>& zeta, vector<int>& row_list,
							vector<int>& col_list, int start_type, double start_depth,
							double startxLoc, double rho_r);

void collect_particles(list<LSDCRNParticle>& collected_list,
						vector< list<LSDCRNParticle> >& eroded_list_vec);

void manipulate_and_print_collected_particles(double t_ime,
							double rho_r,
							list<LSDCRNParticle>& collected_list,
							vector<double>& Ne_bin_locs,
							ofstream& app_eros_out,
					      ofstream& Ne_pdf_out, LSDCRNParameters& CRN_param);

// this takes a list and chages the zeta locations of all the particles
// so that the colum can be inserted into a landscape with
// a different surface elevation than before
// the d loc data elements and eff_d remain the same
void update_list_z_location(list<LSDCRNParticle>& CRN_list,
							double zeta_new);

// This is for use with an advective coordinate system, where
// zeta_adv = zeta-zeta_bl, whith zeta_bl is the base level
//
// You give the function the uplift displacement and the zetas are all updated
// accordingly
//
// Uplift component is a distance, e.g., U*dt
void update_list_z_for_advective_coord_system(list<LSDCRNParticle>& CRN_list,
							double uplift_component);

// this prints information form a particle list
// the current function prints in the format
// n_parts t_ime d_loc1  d_loc2 ...
// n_parts t_ime C_10Be1 C_10Be2 ...
// n_parts t_ime C_14C1  C_14C2 ...
// so you have three rows, the first two data elements in each row are
// the number of particles,
// the second part is the time and the following elements are
// a series of data points,
// one for each particle
// the first row is the depths of the particles,
// second row is the 10Be conc (in atoms/g)
// and the third row is the 14C conc (in atoms/g)
void print_particle_list(list<LSDCRNParticle>& CRN_list,
							double t_ime, ofstream& particles_out);

// prints the particles and apparent erosion rates
// each particle has its own row
void print_particles_and_apparent_erosion_3CRN
		(list<LSDCRNParticle>&  eroded_part_list,
		double rho_r, double dt, double t_ime, ofstream& eroded_part_out, 
    LSDCRNParameters& CRN_param);


void print_particle_list_3CRN(list<LSDCRNParticle>& CRN_list,
							double t_ime, ofstream& particles_out);

void print_particle_list_vec(vector< list<LSDCRNParticle> >& CRN_list_vec,
							vector<int>& row_list,
							vector<int>& col_list,
							double node_spacing,
							double t_ime, ofstream& particles_out);

void print_particle_list_vec_3CRN(vector< list<LSDCRNParticle> >& CRN_list_vec,
							vector<int>& row_list,
							vector<int>& col_list,
							double node_spacing,
							double t_ime, ofstream& particles_out);

void print_eroded_list_vec(vector< list<LSDCRNParticle> >& eroded_part_list,
							vector<int>& row_list,
							vector<int>& col_list,
							double node_spacing,
							double t_ime, ofstream& eroded_part_out);

void print_eroded_list_vec_3CRN(vector< list<LSDCRNParticle> >& eroded_part_list,
							vector<int>& row_list,
							vector<int>& col_list,
							double node_spacing,
							double t_ime, ofstream& eroded_part_out);

void print_particles_for_erosion_lists_3CRN_neutron_only
						(list<LSDCRNParticle>&  eroded_part_list,
							int i,
							vector<int>& row_list,
							vector<int>& col_list,
							Array2D<double>& zeta,
							Array2D<double>& zeta_old,
							double rho_r,
							double node_spacing,
							double dt,
						 double t_ime, ofstream& eroded_part_out, LSDCRNParameters& CRN_param);

void load_particle_list_vec(vector< list<LSDCRNParticle> >& CRN_list_vec,
							vector<int>& row_list,
							vector<int>& col_list,
							Array2D<double>& zeta,
							int start_type, double startxLoc, double rho_r,
							double t_ime, ifstream& particles_in,
							ifstream& tracking_list);

void load_particle_list_vec_3CRN(vector< list<LSDCRNParticle> >& CRN_list_vec,
							vector<int>& row_list,
							vector<int>& col_list,
							Array2D<double>& zeta,
							int start_type, double startxLoc, double rho_r,
							double t_ime, ifstream& particles_in,
							ifstream& tracking_list);

#endif
