#ifndef DATA_OPS_H
#define DATA_OPS_H

#include "snapshot_io.h"

/////////////////////////////////////////////////// DEFINITIONS

/**
*	order type enumeration to specify how to sort the particles with the reordering function
*/
enum OrderType
{
     SORT_ID,		///<sort the particles by their ID
     SORT_TYPE,		///<sort the particles by their type (0..gas, 1..halo...)
     SORT_TYPE_ID	///<sort the particles first by their type, then by their id
};

/**
*	struct to store ratios
*/
typedef struct 
{
     float time;		///<time of the snapshot
     float ratio_ns_wd;		///<ratio of new formed stars wake / disk (w/o IC new stars)
     float ratio_ns_wt;		///<ratio of new formed stars wake / total (w/o IC new stars)
     float ratio_ns_wd_ss;	///<ratio of new formed stars wake / disk from a singel snapshot
     float ratio_gas;		///<ratio of gas wake / total
     float ratio_gastemp;	///<ratio of mean gas temperature wake / disk
} InvestData;

typedef struct 
{
     double temperature;
     double u;
     double density;
     double rho;
     int counter;
} GridData;

#define HOT_GAS_TEMP	5E6

////////////////////////////////////////////// TYPE DEFINITIONS

//////////////////////////////////////////// FUNCTION PROTOTYPES
/**
*	reorders the particles in the array
*
*	\param pd		pointer to the set of particle data
*	\param size		number of particles
*	\param ot		sorting criterion, use enum OrderType
*	\return			zero
*/
int reordering(particle_data *pd, int size, enum OrderType ot);
/**
*	filters the particles from the ICM by sorting the particles by id and cutting off
*
*	\param P		pointer to the set of particle data
*	\param size		number of particles
*	\param nP		pointer to a pointer of particle date where the filtered particles are stored,
*					the memory is allocated in the function
*	\param numICM		maximum particle id till which the particles belong to the ICM
*	\param masses		array of particle masses to search for
*	\param numMasses	the number of masses
*	\return			number of particles which don't belong to the ICM
*/
int filterICM(particle_data *const P, const int size, particle_data **nP, const int numICM, const float *const masses, const int numMasses);
/** 	this template shows how one may convert from Gadget's units
* 	to cgs units.
* 	In this example, the temperature of the gas is computed.
* 	(assuming that the electron density in units of the hydrogen density
*	was computed by the code. This is done if cooling is enabled.)
*
*	\param P	pointer to the set of particle data
*	\param size	number of particles
*	\return
*/
int unitConversion(particle_data *P, const int size);
/**
*	calculates the radial density profile,
*	the output is in 1E10 sun masses/kpc²(y) and kpc(x)
*
*	\param P	pointer to the set of particle data
*	\param size	number of particles
*	\param binSize	the size of the bins to calculate the histogram
*	\param numBins	number of bins
*	\param centroid a centroid can be defined here, if NULL, the centroid of the particles P is used
*	\param data	the results are returned here
*	\return		0
*/
int calcRadialDensityProfile(const particle_data *const P, const int size, const float binSize, const int numBins, float *centroid, float *data);
/**
*	calculates some star formation and gas ratios
*
*	\param P		pointer to the set of particle data
*	\param header		pointer to an io_header structure, corresponding to the set of particles
*	\param wakeDist		distance from the centroid of the old star disk, from which of the particles are considered to be in the wake [kpc]
*	\param idata		pointer to the structure where the calculated ratios are stored
*	\param initNsMass	initial mass of new stars in the galaxy
*	\param phi		inclination angle of the galaxy relative to the wind tunnel
*	\param h		height (thickness of the galaxy, respectively the distance where particles are considered to belong to the wake)
*	\return			0
*/
int calcWakeDiskRatios(const particle_data *const P, const io_header *const header, const float wakeDist, InvestData *const idata, const float initNsMass, const float phi, const float h, const int first_run);
#ifdef HSML
/**
*	calculates the mean smoothing length of the gas particles
*
*	\param P	pointer to the set of particle data
*	\param size	number of particles
*	\return		the mean smoothing length
*/
float calcMeanSL(const particle_data *const P, const int size);
#endif
/**
*
*
*/
GridData* CICmethod(const particle_data *const P, const int size, int type, GridData **gd_p);

int inDisk(const double centroid[3], const float r, const float h, float phi, const float pos[3]);


#endif
