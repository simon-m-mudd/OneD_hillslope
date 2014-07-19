#ifndef CRN_FUNCS_CPP
#define CRN_FUNCS_CPP


#include <list>
#include <fstream>
#include <vector>
#include "LSDParticle.hpp"
#include "LSDCRNParameters.hpp"
#include "TNT/tnt.h"
#include "CRN_funcs.hpp"
using namespace TNT;

#endif



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function starts a column at steady state
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
list<LSDCRNParticle> initiate_SS_cosmo_column(int start_type, double startxLoc,
		      double start_depth, double particle_spacing, 
		      double zeta, double rho_r, double eff_eros_rate,
		      LSDCRNParameters& CRN_param)
{
  // get number of particles
  int N_particles = start_depth/particle_spacing;

  // create the list
  list<LSDCRNParticle> CRN_plist;

  double this_depth;

  // now loop over the depths, inserting particles and setting them to steady state
  for (int p = 0; p< N_particles; p++)
  {
    // first the depth
    this_depth = start_depth-double(p)*particle_spacing;

    // now insert the particle into the list
    insert_part_into_CRN_list(CRN_plist, start_type, startxLoc,
			      this_depth,zeta, rho_r);
  }

  // now loop through the particles, setting to steady
  list<LSDCRNParticle>::iterator part_iter;
  part_iter = CRN_plist.begin();
  while(part_iter != CRN_plist.end())
  {
    // update the CRN_concntrations
    ( *part_iter ).update_10Be_SSfull(eff_eros_rate,CRN_param);
    ( *part_iter ).update_14C_SSfull(eff_eros_rate,CRN_param);
    ( *part_iter ).update_21Ne_SSfull(eff_eros_rate,CRN_param);
 
    part_iter++;   
  }

  return CRN_plist;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// update_list
// this is the function that inserts a particle into a list
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void insert_part_into_CRN_list(list<LSDCRNParticle>& CRN_list,
	 int start_type, double startxLoc,
	 double start_depth,double zeta,
	 double rho_r)
{
	double d = start_depth;
	double eff_d = start_depth*rho_r*0.1;
	double z_p = zeta-start_depth;
	LSDCRNParticle CRN_tp(start_type,startxLoc,d,eff_d,z_p);
	CRN_list.push_back(CRN_tp);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a function that takes an old and new zeta, calculates
// an erosion rate, and updates the CRN concntrations
// of all the particles in the list
// in addition, particles are popped from the list (eroded)
// if their final position is above the new surface,
// and if the lowermost particle is above the start depth by a
// distance greater than or equal to particle_spacing
// a new particle is inserted into the column
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
list<LSDCRNParticle> update_CRN_list(list<LSDCRNParticle>& CRN_list,
	double dt, double rho_r,
	int start_type, double start_depth,
	double startxLoc,
	double zeta_old,double zeta_new,
	double particle_spacing, 
        LSDCRNParameters& CRN_param)
{
	double d,eff_d,z_p,eff_eros_rate;
	int eroded_test;

	list<LSDCRNParticle> eroded_list;
	list<LSDCRNParticle>::iterator part_iter;
	list<LSDCRNParticle>::iterator remove_iter;

	eff_eros_rate = rho_r*0.1*(zeta_new-zeta_old)/dt;

	// first go through and update the CRN concentrations
	// in the list
	part_iter = CRN_list.begin();
	while(part_iter != CRN_list.end())
	{
		z_p = ( *part_iter ).get_zetaLoc();
		d = zeta_new-z_p;
		eff_d = d*rho_r*0.1;

		// update the CRN_concntrations
		( *part_iter ).update_10Be_conc(dt,eff_eros_rate,CRN_param);
		( *part_iter ).update_14C_conc(dt,eff_eros_rate,CRN_param);

		// update the depths (note, teh z_locations arene't updated
		// because it is a rock system, no 'fluffing' of soil occurs
		( *part_iter ).update_depths(d, eff_d);
		part_iter++;
	}

	// now go through the list and see if the partiucles
	// are either eroded or a new particle needs to be added
	// particles are added to the back of the list and eroded from the front
	// so first check for erosion
	eroded_test = 0;
	do
	{
		part_iter = CRN_list.begin();
		z_p = ( *part_iter ).get_zetaLoc();

		// if the elevation of the particle is greater than the elevation
		// of the surface, it is eroded from the particle list
		// and eroded_test remains unchanged
		// if the elevation of the particle is not greater than the surface
		// elevation, then eroded_test goes to 1 and the loop is exited
		if (z_p > zeta_new)
		{
			//cout << "LINE 82 popping out!" << endl;
			eroded_list.push_back(*part_iter);
			CRN_list.pop_front();
		}
		else
		{
			eroded_test = 1;
		}
	} while (eroded_test == 0);

	// now see if we insert a particle
	z_p = (CRN_list.back()).get_zetaLoc();
	d = zeta_new-z_p;
	if ( (start_depth - d) >= particle_spacing)
	{
		//cout << "LINE 95 inserting_particle!" << endl;
		insert_part_into_CRN_list(CRN_list, start_type, startxLoc,
				 start_depth, zeta_new, rho_r);
	}


	return eroded_list;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This is almost exactly the same as the preceding function, with 
// the key difference that if a particle is eroded it only gets partial exposure 
// for the time it is in the ground. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
list<LSDCRNParticle> update_CRN_list_eros_limit(list<LSDCRNParticle>& CRN_list,
		double dt, double rho_r,
		int start_type, double start_depth,
		double startxLoc,
		double zeta_old,double zeta_new,
		double particle_spacing, LSDCRNParameters& CRN_param)
{
	double d,eff_d,z_p,eff_eros_rate;
	int eroded_test;
	double effective_dt;
	double d_frac;

	list<LSDCRNParticle> eroded_list;
	list<LSDCRNParticle>::iterator part_iter;
	list<LSDCRNParticle>::iterator remove_iter;

	eff_eros_rate = rho_r*0.1*(zeta_new-zeta_old)/dt;

	// first go through and update the CRN concentrations
	// in the list
	part_iter = CRN_list.begin();
	while(part_iter != CRN_list.end())
	{
		z_p = ( *part_iter ).get_zetaLoc();
		d = zeta_new-z_p;
		eff_d = d*rho_r*0.1;
		//cout << endl << endl << "zeta_new: " << zeta_new << " z_p: " << z_p << endl;
		//cout << "yoyoma: " << zeta_new-z_p << endl;
		//cout << "z_p: " << z_p << "  z_new: " << zeta_new << " d: " << d << endl;

		// check to see if the particle is above the surface, if
		// it is use an abbreviated exposure
		if (z_p >= zeta_new)
		{
			d_frac = (zeta_old-z_p)/(zeta_old-zeta_new);
			effective_dt = d_frac*dt;
		}
		else
		{
			effective_dt = dt;
		}

		// update the CRN_concntrations
		( *part_iter ).update_10Be_conc(effective_dt,eff_eros_rate, CRN_param);
		( *part_iter ).update_14C_conc(effective_dt,eff_eros_rate, CRN_param);

		// update the depths (note, the z_locations aren't updated
		// because it is a rock system, no 'fluffing' of soil occurs
		( *part_iter ).update_depths(d, eff_d);
		part_iter++;
	}

	// now go through the list and see if the particles
	// are either eroded or a new particle needs to be added
	// particles are added to the back of the list and eroded from the front
	// so first check for erosion
	eroded_test = 0;
	do
	{
		part_iter = CRN_list.begin();
		z_p = ( *part_iter ).get_zetaLoc();

		// if the elevation of the particle is greater than the elevation
		// of the surface, it is eroded from the particle list
		// and eroded_test remains unchanged
		// if the elevation of the particle is not greater than the surface
		// elevation, then eroded_test goes to 1 and the loop is exited
		if (z_p > zeta_new)
		{
			//cout << "LINE 82 popping out!" << endl;
			eroded_list.push_back(*part_iter);
			CRN_list.pop_front();
		}
		else
		{
			eroded_test = 1;
		}
	} while (eroded_test == 0);

	// now see if we insert a particle
	z_p = (CRN_list.back()).get_zetaLoc();
	d = zeta_new-z_p;
	if ( (start_depth - d) >= particle_spacing)
	{
		//cout << "LINE 95 inserting_particle!" << endl;
		insert_part_into_CRN_list(CRN_list, start_type, startxLoc,
				 start_depth, zeta_new, rho_r);
	}


	return eroded_list;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function updates the cosmo concentrations for
// 10Be, 21Ne and 14C
// It also calculates 'effective' concentration for particles
// leaving the surface so they only gain nuclide concetration
// while still in the regolith. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
list<LSDCRNParticle> update_CRN_list_eros_limit_3CRN(list<LSDCRNParticle>& CRN_list,
	double dt, double rho_r,
	int start_type, double start_depth,
	double startxLoc,
	double zeta_old,double zeta_new,
        double particle_spacing, LSDCRNParameters& CRN_param)
{
	double d,eff_d,z_p,eff_eros_rate;
	int eroded_test;
	double effective_dt;
	double d_frac;

	list<LSDCRNParticle> eroded_list;
	list<LSDCRNParticle>::iterator part_iter;
	list<LSDCRNParticle>::iterator remove_iter;

	eff_eros_rate = rho_r*0.1*(zeta_new-zeta_old)/dt;

	// first go through and update the CRN concentrations
	// in the list
	part_iter = CRN_list.begin();
	while(part_iter != CRN_list.end())
	{
		z_p = ( *part_iter ).get_zetaLoc();
		d = zeta_new-z_p;
		eff_d = d*rho_r*0.1;
		//cout << endl << endl << "zeta_new: " << zeta_new << " z_p: " << z_p << endl;
		//cout << "yoyoma: " << zeta_new-z_p << endl;
		//cout << "z_p: " << z_p << "  z_new: " << zeta_new << " d: " << d << endl;

		// check to see if the particle is above the surface, if
		// it is use an abbreviated exposure
		if (z_p >= zeta_new)
		{
			d_frac = (zeta_old-z_p)/(zeta_old-zeta_new);
			effective_dt = d_frac*dt;
		}
		else
		{
			effective_dt = dt;
		}

		// update the CRN_concntrations
		( *part_iter ).update_10Be_conc(effective_dt,eff_eros_rate, CRN_param);
		( *part_iter ).update_14C_conc(effective_dt,eff_eros_rate, CRN_param);
		( *part_iter ).update_21Ne_conc(effective_dt,eff_eros_rate, CRN_param);

		// update the depths (note, teh z_locations arene't updated
		// because it is a rock system, no 'fluffing' of soil occurs
		( *part_iter ).update_depths(d, eff_d);
		part_iter++;
	}

	// now go through the list and see if the partiucles
	// are either eroded or a new particle needs to be added
	// particles are added to the back of the list and eroded from the front
	// so first check for erosion
	eroded_test = 0;
	do
	{
		part_iter = CRN_list.begin();
		z_p = ( *part_iter ).get_zetaLoc();

		// if the elevation of the particle is greater than the elevation
		// of the surface, it is eroded from the particle list
		// and eroded_test remains unchanged
		// if the elevation of the particle is not greater than the surface
		// elevation, then eroded_test goes to 1 and the loop is exited
		if (z_p > zeta_new)
		{
			//cout << "LINE 82 popping out!" << endl;
			eroded_list.push_back(*part_iter);
			CRN_list.pop_front();
		}
		else
		{
			eroded_test = 1;
		}
	} while (eroded_test == 0);

	// now see if we insert a particle
	z_p = (CRN_list.back()).get_zetaLoc();
	d = zeta_new-z_p;
	if ( (start_depth - d) >= particle_spacing)
	{
		//cout << "LINE 95 inserting_particle!" << endl;
		insert_part_into_CRN_list(CRN_list, start_type, startxLoc,
				 start_depth, zeta_new, rho_r);
	}


	return eroded_list;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Same as the above function but this only updates the CRN concentration
// based on neutron spallation
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
list<LSDCRNParticle> update_CRN_list_eros_limit_3CRN_neutron(list<LSDCRNParticle>& CRN_list,
		double dt, double rho_r,
		int start_type, double start_depth,
		double startxLoc,
		double zeta_old,double zeta_new,
      	        double particle_spacing,LSDCRNParameters& CRN_param )
{
	double d,eff_d,z_p,eff_eros_rate;
	int eroded_test;
	double effective_dt;
	double d_frac;

	list<LSDCRNParticle> eroded_list;
	list<LSDCRNParticle>::iterator part_iter;
	list<LSDCRNParticle>::iterator remove_iter;

	eff_eros_rate = rho_r*0.1*(zeta_new-zeta_old)/dt;

	// first go through and update the CRN concentrations
	// in the list
	part_iter = CRN_list.begin();
	while(part_iter != CRN_list.end())
	{
		z_p = ( *part_iter ).get_zetaLoc();
		d = zeta_new-z_p;
		eff_d = d*rho_r*0.1;
		//cout << endl << endl << "zeta_new: " << zeta_new << " z_p: " << z_p << endl;
		//cout << "yoyoma: " << zeta_new-z_p << endl;
		//cout << "z_p: " << z_p << "  z_new: " << zeta_new << " d: " << d << endl;

		// check to see if the particle is above the surface, if
		// it is use an abbreviated exposure
		if (z_p >= zeta_new)
		{
			d_frac = (zeta_old-z_p)/(zeta_old-zeta_new);
			effective_dt = d_frac*dt;
		}
		else
		{
			effective_dt = dt;
		}

		// update the CRN_concntrations
		( *part_iter ).update_10Be_conc_neutron_only(effective_dt,eff_eros_rate, CRN_param);
		( *part_iter ).update_14C_conc_neutron_only(effective_dt,eff_eros_rate, CRN_param);
		( *part_iter ).update_21Ne_conc(effective_dt,eff_eros_rate, CRN_param);

		// update the depths (note, teh z_locations arene't updated
		// because it is a rock system, no 'fluffing' of soil occurs
		( *part_iter ).update_depths(d, eff_d);
		part_iter++;
	}

	// now go through the list and see if the partiucles
	// are either eroded or a new particle needs to be added
	// particles are added to the back of the list and eroded from the front
	// so first check for erosion
	eroded_test = 0;
	do
	{
		part_iter = CRN_list.begin();
		z_p = ( *part_iter ).get_zetaLoc();

		// if the elevation of the particle is greater than the elevation
		// of the surface, it is eroded from the particle list
		// and eroded_test remains unchanged
		// if the elevation of the particle is not greater than the surface
		// elevation, then eroded_test goes to 1 and the loop is exited
		if (z_p > zeta_new)
		{
			//cout << "LINE 82 popping out!" << endl;
			eroded_list.push_back(*part_iter);
			CRN_list.pop_front();
		}
		else
		{
			eroded_test = 1;
		}
	} while (eroded_test == 0);

	// now see if we insert a particle
	z_p = (CRN_list.back()).get_zetaLoc();
	d = zeta_new-z_p;
	if ( (start_depth - d) >= particle_spacing)
	{
		//cout << "LINE 95 inserting_particle!" << endl;
		insert_part_into_CRN_list(CRN_list, start_type, startxLoc,
								 start_depth, zeta_new, rho_r);
	}


	return eroded_list;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// 
// FUNCTIONS FOR SPATIALLY HETEROGENOUS COSMO PRODUCTION
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// these function updates the entire vector of CRNS at all active nodes on the
// hillslope
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void update_CRN_list_vec( vector< list<LSDCRNParticle> >& CRN_list_vec,
	Array2D<double>& zeta,
	Array2D<double>& zeta_old,
	vector<int>& row_list,
	vector<int>& col_list,
	double dt, double rho_r,
	int start_type, double start_depth,
	double startxLoc,
   	double particle_spacing, LSDCRNParameters& CRN_param)
{
	int n_columns = CRN_list_vec.size();
	int row,col;
	list<LSDCRNParticle> eroded_list;
	double z_old,z;
	for (int i = 0; i<n_columns; i++)
	{
		row = row_list[i];
		col = col_list[i];
		z = zeta[row][col];
		z_old = zeta_old[row][col];

		eroded_list = update_CRN_list(CRN_list_vec[i],dt, rho_r,
			start_type, start_depth,
			startxLoc, z_old, z,
  		      particle_spacing, CRN_param);
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// For spatially variable CRN with a 2D model
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< list<LSDCRNParticle> > update_CRN_list_vec_with_eroded
 			  ( vector< list<LSDCRNParticle> >& CRN_list_vec,
	Array2D<double>& zeta,
	Array2D<double>& zeta_old,
	vector<int>& row_list,
	vector<int>& col_list,
	double dt, double rho_r,
	int start_type, double start_depth,
	double startxLoc,
        double particle_spacing, LSDCRNParameters& CRN_param)
{
	int n_columns = CRN_list_vec.size();
	int row,col;
	list<LSDCRNParticle> eroded_list;
	vector< list<LSDCRNParticle> > elv(n_columns);

	double z_old,z;
	for (int i = 0; i<n_columns; i++)
	{
		row = row_list[i];
		col = col_list[i];
		z = zeta[row][col];
		z_old = zeta_old[row][col];

		//cout << "zeta: " << z << " zeta_old: " << z_old << endl;
		eroded_list = update_CRN_list_eros_limit(CRN_list_vec[i],dt, rho_r,
			start_type, start_depth,
			startxLoc, z_old, z,
			 particle_spacing, CRN_param);
		elv[i] = eroded_list;
	}
	return elv;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< list<LSDCRNParticle> > update_CRN_list_vec_with_eroded_3CRN
			  ( vector< list<LSDCRNParticle> >& CRN_list_vec,
			Array2D<double>& zeta,
			Array2D<double>& zeta_old,
			vector<int>& row_list,
			vector<int>& col_list,
			double dt, double rho_r,
			int start_type, double start_depth,
			double startxLoc,
		    double particle_spacing, LSDCRNParameters& CRN_param)
{
	int n_columns = CRN_list_vec.size();
	int row,col;
	list<LSDCRNParticle> eroded_list;
	vector< list<LSDCRNParticle> > elv(n_columns);

	double z_old,z;
	for (int i = 0; i<n_columns; i++)
	{
		row = row_list[i];
		col = col_list[i];
		z = zeta[row][col];
		z_old = zeta_old[row][col];

		//cout << "zeta: " << z << " zeta_old: " << z_old << endl;
		eroded_list = update_CRN_list_eros_limit_3CRN(CRN_list_vec[i],dt, rho_r,
					start_type, start_depth,
					startxLoc, z_old, z,
				      particle_spacing, CRN_param);
		elv[i] = eroded_list;
	}
	return elv;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< list<LSDCRNParticle> > update_CRN_list_vec_with_eroded_3CRN_neutron
		  ( vector< list<LSDCRNParticle> >& CRN_list_vec,
			Array2D<double>& zeta,
			Array2D<double>& zeta_old,
			vector<int>& row_list,
			vector<int>& col_list,
			double dt, double rho_r,
			int start_type, double start_depth,
			double startxLoc,
		    double particle_spacing, LSDCRNParameters& CRN_param)
{
	int n_columns = CRN_list_vec.size();
	int row,col;
	list<LSDCRNParticle> eroded_list;
	vector< list<LSDCRNParticle> > elv(n_columns);

	double z_old,z;
	for (int i = 0; i<n_columns; i++)
	{
		row = row_list[i];
		col = col_list[i];
		z = zeta[row][col];
		z_old = zeta_old[row][col];

		//cout << "zeta: " << z << " zeta_old: " << z_old << endl;
		eroded_list = update_CRN_list_eros_limit_3CRN_neutron(CRN_list_vec[i],dt, rho_r,
						start_type, start_depth,
						startxLoc, z_old, z,
						      particle_spacing, CRN_param);
		elv[i] = eroded_list;
	}
	return elv;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

void initialize_CRN_list_vec( vector< list<LSDCRNParticle> >& CRN_list_vec,
		Array2D<double>& zeta, vector<int>& row_list,
		vector<int>& col_list, int start_type, double start_depth,
		double startxLoc, double rho_r)
{
	int n_columns = row_list.size();
	//int row,col;

	double z;
	for(int i = 0; i<n_columns; i++)
	{
		z = zeta[row_list[i]][col_list[i]];
		list<LSDCRNParticle> temp_list;
		insert_part_into_CRN_list(temp_list, start_type, startxLoc,
							  start_depth,  z, rho_r);
		CRN_list_vec.push_back(temp_list);
	}
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// collect_particles
// this function collects particles that are later used to calcualte
// and apparent erosion rate and a 21Ne pdf
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void collect_particles(list<LSDCRNParticle>& collected_list,
		vector< list<LSDCRNParticle> >& eroded_list_vec)
{
	int elv_size = eroded_list_vec.size();
	list<LSDCRNParticle>::iterator part_iter;

	for (int i = 0; i< elv_size; i++)
	{

		// if anything has been eroded
		if (eroded_list_vec[i].size() >0)
		{
			// print the details to a line of the file
			part_iter = eroded_list_vec[i].begin();
			while (part_iter != eroded_list_vec[i].end())
			{
				collected_list.push_back( *part_iter );
				part_iter++;
			}
		}
	}
}


void manipulate_and_print_collected_particles(double t_ime,
			double rho_r,
			list<LSDCRNParticle>& collected_list,
			vector<double>& Ne_bin_locs,
			ofstream& app_eros_out,
		      ofstream& Ne_pdf_out, LSDCRNParameters& CRN_param)
{

	int n_parts = collected_list.size();
	list<LSDCRNParticle>::iterator part_iter;

	double sum_Ne = 0;
	double sum_Be = 0;
	double sum_C = 0;

	double avg_Ne;
	double avg_Be;
	double avg_C;

	double bin_spacing = Ne_bin_locs[1]-Ne_bin_locs[0];
	int n_bins = Ne_bin_locs.size();
	vector<int> Ne_bins_int(n_bins+1,0);
	vector<double> Pdf_Ne(n_bins+1,0);
	int bin_number;

	part_iter = collected_list.begin();
	while (part_iter != collected_list.end())
	{
		sum_Ne+= (*part_iter).getConc_21Ne();
		sum_Be+= (*part_iter).getConc_10Be();
		sum_C+= (*part_iter).getConc_14C();

		bin_number = int(double((*part_iter).getConc_21Ne()/bin_spacing));
		//cout << "bin_spacing: " << bin_spacing << " Conc: " << (*part_iter).getConc_21Ne()
		//     << " bin_number: " << bin_number << endl;
		if (bin_number > n_bins-1)
		{
			Ne_bins_int[n_bins]++;
		}
		else
		{
			Ne_bins_int[bin_number]++;
		}
		part_iter++;

	}

	avg_Ne = sum_Ne/double(n_parts);
	avg_Be = sum_Be/double(n_parts);
	avg_C  = sum_C/double(n_parts);

	//cout << "LINE 654 averages "<< endl;
	//cout << avg_Ne << " " << avg_Be << " " << avg_C << endl;

	int start_type = 0;
	double startxLoc = 0;
	double start_z = 0;
	double start_d = 0;
	double start_eff_d = 0;
	LSDCRNParticle initial_part(start_type, startxLoc, start_z,
							   start_d, start_eff_d,
							   avg_Be, avg_C, avg_Ne);

	app_eros_out << t_ime << " "
		     << initial_part.apparent_erosion_10Be_neutron_only(rho_r, CRN_param) << " "
		     << initial_part.apparent_erosion_14C_neutron_only(rho_r, CRN_param) << " "
		     << initial_part.apparent_erosion_21Ne(rho_r, CRN_param) << endl;

	//double prob;
	Ne_pdf_out << t_ime;
	for (int i = 0; i<=n_bins; i++)
	{
		Ne_pdf_out << " " << double(Ne_bins_int[i])/double(n_parts);
	}
	Ne_pdf_out << endl;

}


// this takes a list and chages the zeta locations of all the particles
// so that the colum can be inserted into a landscape with
// a different surface elevation than before
// the d loc data elements and eff_d remain the same
void update_list_z_location(list<LSDCRNParticle>& CRN_list,
							double zeta_new)
{
	list<LSDCRNParticle>::iterator part_iter;
	part_iter = CRN_list.begin();
	while (part_iter != CRN_list.end())
	{
		(*part_iter).update_zetaLoc_with_new_surface(zeta_new);
		part_iter++;
	}
}


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
			double t_ime, ofstream& particles_out)
{
	int n_particles = CRN_list.size();

	vector<double> d_vec;
	vector<double> Conc10Be_vec;
	vector<double> Conc14C_vec;

	list<LSDCRNParticle>::iterator part_iter;
	part_iter = CRN_list.begin();
	while (part_iter != CRN_list.end())
	{
		d_vec.push_back( (*part_iter).getdLoc() );
		Conc10Be_vec.push_back( (*part_iter).getConc_10Be() );
		Conc14C_vec.push_back( (*part_iter).getConc_14C() );
		part_iter++;
	}

	// print the depths
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << d_vec[i];
	}
	particles_out << endl;

	// print the berillium concentrations
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << Conc10Be_vec[i];
	}
	particles_out << endl;

	// print the berillium concentrations
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << Conc14C_vec[i];
	}
	particles_out << endl;

}


void print_particle_list_3CRN(list<LSDCRNParticle>& CRN_list,
				double t_ime, ofstream& particles_out)
{
	int n_particles = CRN_list.size();

	vector<double> d_vec;
	vector<double> Conc10Be_vec;
	vector<double> Conc14C_vec;
	vector<double> Conc21Ne_vec;

	list<LSDCRNParticle>::iterator part_iter;
	part_iter = CRN_list.begin();
	while (part_iter != CRN_list.end())
	{
		d_vec.push_back( (*part_iter).getdLoc() );
		Conc10Be_vec.push_back( (*part_iter).getConc_10Be() );
		Conc14C_vec.push_back( (*part_iter).getConc_14C() );
		Conc21Ne_vec.push_back( (*part_iter).getConc_21Ne() );
		part_iter++;
	}

	// print the depths
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << d_vec[i];
	}
	particles_out << endl;

	// print the berillium concentrations
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << Conc10Be_vec[i];
	}
	particles_out << endl;

	// print the 14C concentrations
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << Conc14C_vec[i];
	}
	particles_out << endl;

	// print the 21Ne concentrations
	particles_out << n_particles << " " << t_ime;
	for (int i = 0; i< n_particles; i++)
	{
		particles_out << " " << Conc21Ne_vec[i];
	}
	particles_out << endl;

}

void print_particle_list_vec(vector< list<LSDCRNParticle> >& CRN_list_vec,
		vector<int>& row_list,
		vector<int>& col_list,
		double node_spacing,
		double t_ime, ofstream& particles_out)
{
	int n_bins = CRN_list_vec.size();
	//cout << "LINE 252 CRN_funcs.cpp, n_bins: " << n_bins << endl;
	for (int i = 0; i<n_bins; i++)
	{
		particles_out << double(row_list[i])*node_spacing << " "
					  << double(col_list[i])*node_spacing << endl;
		print_particle_list(CRN_list_vec[i], t_ime, particles_out);
	}
}

void print_particle_list_vec_3CRN(vector< list<LSDCRNParticle> >& CRN_list_vec,
		vector<int>& row_list,
		vector<int>& col_list,
		double node_spacing,
		double t_ime, ofstream& particles_out)
{
	int n_bins = CRN_list_vec.size();
	//cout << "LINE 252 CRN_funcs.cpp, n_bins: " << n_bins << endl;
	for (int i = 0; i<n_bins; i++)
	{
		particles_out << double(row_list[i])*node_spacing << " "
					  << double(col_list[i])*node_spacing << endl;
		print_particle_list_3CRN(CRN_list_vec[i], t_ime, particles_out);
	}
}

void print_eroded_list_vec(vector< list<LSDCRNParticle> >& eroded_part_list,
		vector<int>& row_list,
		vector<int>& col_list,
		double node_spacing,
		double t_ime, ofstream& eroded_part_out)
{
	int size_epl = eroded_part_list.size();
	list<LSDCRNParticle>::iterator part_iter;

	// loop through all the possible nodes
	for (int i  = 0; i<size_epl; i++)
	{
		// if anything has been eroded
		if (eroded_part_list[i].size() >0)
		{
			// print the details to a line of the file


			part_iter = eroded_part_list[i].begin();
			while (part_iter != eroded_part_list[i].end())
			{
				eroded_part_out << t_ime << " "
			              << double(row_list[i])*node_spacing << " "
					      << double(col_list[i])*node_spacing << " "
					      << (*part_iter).getConc_10Be() << " "
						  << (*part_iter).getConc_14C() << " " << endl;
				part_iter++;
			}
		}
	}
}

void print_eroded_list_vec_3CRN(vector< list<LSDCRNParticle> >& eroded_part_list,
		vector<int>& row_list,
		vector<int>& col_list,
		double node_spacing,
		double t_ime, ofstream& eroded_part_out)
{
	int size_epl = eroded_part_list.size();
	list<LSDCRNParticle>::iterator part_iter;

	// loop through all the possible nodes
	for (int i  = 0; i<size_epl; i++)
	{
		// if anything has been eroded
		if (eroded_part_list[i].size() >0)
		{
			// print the details to a line of the file


			part_iter = eroded_part_list[i].begin();
			while (part_iter != eroded_part_list[i].end())
			{
				eroded_part_out << t_ime << " "
			              << double(row_list[i])*node_spacing << " "
					      << double(col_list[i])*node_spacing << " "
					      << (*part_iter).getConc_10Be() << " "
						  << (*part_iter).getConc_14C() << " "
						  << (*part_iter).getConc_21Ne() << " "
						  << endl;
				part_iter++;
			}
		}
	}
}


void print_eroded_list_vec_3CRN(vector< list<LSDCRNParticle> >& eroded_part_list,
		vector<int>& row_list,
		vector<int>& col_list,
		Array2D<double>& zeta,
		Array2D<double>& zeta_old,
		double rho_r,
		double node_spacing,
		double dt,
	double t_ime, ofstream& eroded_part_out, LSDCRNParameters& CRN_param)
{
	int size_epl = eroded_part_list.size();
	list<LSDCRNParticle>::iterator part_iter;

	// loop through all the possible nodes
	for (int i  = 0; i<size_epl; i++)
	{
		// if anything has been eroded
		if (eroded_part_list[i].size() >0)
		{
			// print the details to a line of the file


			part_iter = eroded_part_list[i].begin();
			while (part_iter != eroded_part_list[i].end())
			{
				eroded_part_out << t_ime << " "
			              << double(row_list[i])*node_spacing << " "
					      << double(col_list[i])*node_spacing << " "
					      << (*part_iter).getConc_10Be() << " "
						  << (*part_iter).getConc_14C() << " "
						  << (*part_iter).getConc_21Ne() << " "
						  << (zeta[ row_list[i] ][ col_list[i] ]
						       - zeta_old[ row_list[i] ][ col_list[i] ])/dt << " "
						<< (*part_iter).apparent_erosion_10Be_neutron_only(rho_r, CRN_param) << " "
						<< (*part_iter).apparent_erosion_14C_neutron_only(rho_r, CRN_param) << " "
						<< (*part_iter).apparent_erosion_21Ne(rho_r, CRN_param) << " "
						  << endl;
				part_iter++;
			}
		}
	}
}

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
	 double t_ime, ofstream& eroded_part_out, LSDCRNParameters& CRN_param)
{
	list<LSDCRNParticle>::iterator part_iter;
	part_iter = eroded_part_list.begin();
	while (part_iter != eroded_part_list.end())
	{
		eroded_part_out << t_ime << " "
		  << double(row_list[i])*node_spacing << " "
		  << double(col_list[i])*node_spacing << " "
		<< (*part_iter).getConc_10Be() << " "
		  << (*part_iter).getConc_14C() << " "
		  << (*part_iter).getConc_21Ne() << " "
		  << (zeta[ row_list[i] ][ col_list[i] ]
		   - zeta_old[ row_list[i] ][ col_list[i] ])/dt << " "
		<< (*part_iter).apparent_erosion_10Be_neutron_only(rho_r, CRN_param) << " "
		<< (*part_iter).apparent_erosion_14C_neutron_only(rho_r, CRN_param) << " "
		<< (*part_iter).apparent_erosion_21Ne(rho_r, CRN_param) << " "
		  << endl;
		part_iter++;
	}
}


void load_particle_list_vec(vector< list<LSDCRNParticle> >& CRN_list_vec,
			vector<int>& row_list,
		vector<int>& col_list,
		Array2D<double>& zeta,
		int start_type, double startxLoc, double rho_r,
		double t_ime, ifstream& particles_in,
		ifstream& tracking_list)
{
	double start_z;
	double start_eff_d;
	//double row_dist, col_dist;
	double temp_double;
	int temp_int;
	double node_spacing;
	int n_p_columns;
	double z;
	int n_parts;
	vector<double> d_locs;
	vector<double> Conc_10Be;
	vector<double> Conc_14C;


	// first load the row and column list
	tracking_list >> n_p_columns >> node_spacing;
	cout << "LINE 285 CRN_funcs, n_list: " << n_p_columns << endl;
	for (int i = 0; i<n_p_columns; i++)
	{
		tracking_list >> temp_int;
		row_list.push_back(temp_int);
		tracking_list >> temp_int;
		col_list.push_back(temp_int);
	}

	// initialize a list vec of the appropriate number of lists
	vector< list<LSDCRNParticle> > temp_CRN_list_vec(n_p_columns);

	// now load the particle data
	for (int i = 0; i<n_p_columns; i++)
	{
		//cout <<  "LINE 300 looping through column  " << i << endl;

		// generate a temporary list
		list<LSDCRNParticle> temp_CRN_list;

		// vectors for storing the data
		vector<double> d_locs;
		vector<double> Conc_10Be;
		vector<double> Conc_14C;

		particles_in >> temp_double;
		//cout << "LINE 311 CRN_funcs.cpp yloc: " << temp_double << endl;
		particles_in >> temp_double;
		//cout << "LINE 313 CRN_funcs.cpp xloc: " << temp_double << endl;


		z = zeta[row_list[i]][col_list[i]];
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			d_locs.push_back(temp_double);
		}
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			Conc_10Be.push_back(temp_double);
		}
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			Conc_14C.push_back(temp_double);
		}

		// now loop through the nodes
		// as you read in teh data from teh file,
		// the first node was nearest to the surface
		for (int j = 0; j<n_parts; j++)
		{
			start_z = z-d_locs[j];
			start_eff_d = rho_r*d_locs[j]*0.1;	// 0.1 is to convert to g/cm^2

			// create the particle
			LSDCRNParticle initial_part(start_type, startxLoc, start_z,
								  d_locs[j], start_eff_d,
									Conc_10Be[j],Conc_14C[j]);

			// add it to the temporary list
			temp_CRN_list.push_back(initial_part);
		}
		// add the temp list to the list vec
		temp_CRN_list_vec[i] = temp_CRN_list;

	}

	// now reset the main list_vec
	CRN_list_vec = temp_CRN_list_vec;
}

void load_particle_list_vec_3CRN(vector< list<LSDCRNParticle> >& CRN_list_vec,
			vector<int>& row_list,
			vector<int>& col_list,
			Array2D<double>& zeta,
			int start_type, double startxLoc, double rho_r,
			double t_ime, ifstream& particles_in,
			ifstream& tracking_list)
{
	double start_z;
	double start_eff_d;
	//double row_dist, col_dist;
	double temp_double;
	int temp_int;
	double node_spacing;
	int n_p_columns;
	double z;
	int n_parts;
	vector<double> d_locs;
	vector<double> Conc_10Be;
	vector<double> Conc_14C;
	vector<double> Conc_21Ne;


	// first load the row and column list
	tracking_list >> n_p_columns >> node_spacing;
	cout << "LINE 285 CRN_funcs, n_list: " << n_p_columns << endl;
	for (int i = 0; i<n_p_columns; i++)
	{
		tracking_list >> temp_int;
		row_list.push_back(temp_int);
		tracking_list >> temp_int;
		col_list.push_back(temp_int);
	}

	// initialize a list vec of the appropriate number of lists
	vector< list<LSDCRNParticle> > temp_CRN_list_vec(n_p_columns);

	// now load the particle data
	for (int i = 0; i<n_p_columns; i++)
	{
		//cout <<  "LINE 300 looping through column  " << i << endl;

		// generate a temporary list
		list<LSDCRNParticle> temp_CRN_list;

		// vectors for storing the data
		vector<double> d_locs;
		vector<double> Conc_10Be;
		vector<double> Conc_14C;
		vector<double> Conc_21Ne;

		particles_in >> temp_double;
		//cout << "LINE 311 CRN_funcs.cpp yloc: " << temp_double << endl;
		particles_in >> temp_double;
		//cout << "LINE 313 CRN_funcs.cpp xloc: " << temp_double << endl;


		z = zeta[row_list[i]][col_list[i]];
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			d_locs.push_back(temp_double);
		}
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			Conc_10Be.push_back(temp_double);
		}
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			Conc_14C.push_back(temp_double);
		}
		particles_in >> n_parts >> temp_double;
		for (int j = 0; j<n_parts; j++)
		{
			particles_in >> temp_double;
			Conc_21Ne.push_back(temp_double);
		}


		// now loop through the nodes
		// as you read in teh data from teh file,
		// the first node was nearest to the surface
		for (int j = 0; j<n_parts; j++)
		{
			start_z = z-d_locs[j];
			start_eff_d = rho_r*d_locs[j]*0.1;	// 0.1 is to convert to g/cm^2

			// create the particle
			LSDCRNParticle initial_part(start_type, startxLoc, start_z,
								  d_locs[j], start_eff_d,
									Conc_10Be[j],Conc_14C[j], Conc_21Ne[j]);
			//LSDCRNParticle initial_part(start_type, startxLoc, start_z,
			//				  d_locs[j], start_eff_d,
			//					Conc_10Be[j],Conc_14C[j]);


			// add it to the temporary list
			temp_CRN_list.push_back(initial_part);
		}
		// add the temp list to the list vec
		temp_CRN_list_vec[i] = temp_CRN_list;

	}

	// now reset the main list_vec
	CRN_list_vec = temp_CRN_list_vec;
}



