#ifndef SNAPSHOT_IO_H
#define SNAPSHOT_IO_H

/////////////////////////////////////////////////// DEFINITIONS

#define HEADER_SIZE	256	///size of the GADGET header in bytes

#define P_GAS			1
#define P_HALO		2
#define P_DISK		4
#define P_BULGE		8
#define P_STARS		16
#define P_BNDRY		32

////////////////////////////////////////////// TYPE DEFINITIONS

#define IO_NBLOCKS 17   /*!< total number of defined information blocks for snapshot files.
                             Must be equal to the number of entries in "enum iofields" */

enum iofields /*!< this enumeration lists the defined output blocks in snapshot files. Not all of them need to be present. */
{
// 	IO_HEADER,
	IO_POS, IO_VEL, IO_ID, IO_MASS, IO_U, IO_RHO,
	IO_NE, //COOLING
	IO_NH, //COOLING
	IO_HSML,
	IO_SFR, //SFR
	IO_Z, //METALS
	IO_STELLARAGE, //STELLARAGE
	IO_POT,
	IO_ACCEL,
	IO_DTENTR,
	IO_TSTP,
	IO_COOR, //OUTPUTCOOLRATE (cooling rate)
};

typedef struct  
{
     //all particles
     float	Pos[3];
     float 	Vel[3];
     float	Mass;
     float 	Timestep;
     int 	Id;     	//particle id
     float 	Pot;	//gravitational potential 
     float 	Accel;
     
     //only gas
     float Entropy;  //change of entropy
     float Rho; 	//density
     float U; 	//internal energy
     float Ne;	//electron abundance
     float Nnh;	//neutral hydrogen abundance
     float CoolRate; //current cooling rate
     float Hsml;	//smoothing length
     float sfr;	//star formation rate (only gas 0), according to G2 code, unit is M_s/yr
     
     //only stellar particles (type 4)
     float age;	//mean stellar age (in units from timeBegin-timeMax in .tex file) (only type 4)
     
     //0 & 4 gas and stars
     float metals;	//Z metallicity (only for gas 0 and newly formed stars 4 particles)
     
     //additional, not in snapshot file
     int   	Type;
     float 	Temp;	//temperature
} particle_data;



typedef struct 
{
     unsigned int npart[6];
     double mass[6];
     double time;
     double redshift;
     int flag_sfr;
     int flag_feedback;
     int npartTotal[6];
     int flag_cooling;
     int num_files;
     double BoxSize;
     double Omega0;
     double OmegaLambda;
     double HubbleParam;
     int flagAge;	//unused
     int flagMetals;	//unused
     int nAllHW[6];
     int flagEntrICs;
     //char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} io_header;


//////////////////////////////////////////// FUNCTION PROTOTYPES


/**
*	reads a GADGET snapshot file
*	
*	\param filename	name of the file
*	\param header	header of the GADGET file to read
*	\param P	pointer to the set of particle data, the memory is allocated in the function
*	\return		number of loaded particles on success, negative number else
*
*	\todo check block sizes instead of SKIP
*/
int load_snapshot(const char *const filename, io_header *const header, particle_data **const P);
/**
*	reads a ICFormat 2 GADGET snapshot file 
*	
*	\param filename	name of the file
*	\param header	header of the GADGET file to read
*	\param P	pointer to the set of particle data, the memory is allocated in the function
*	\return		number of loaded particles on success, negative number else
*
*	\todo check block sizes instead of SKIP
*/
int load_snapshotF2(const char *const filename, io_header *const header, particle_data **const P);
int load_snapshot_format2(const char *const filename, io_header *const header, particle_data **const P);
/**
*	writes a GADGET snapshot file
*	
*	\param filename	name of the file
*	\param header	header of the GADGET file to write
*	\param P	pointer to the set of particle data
*	\return		0 on success, negative number else
*/
int write_snapshot(const char *const filename, const io_header *const header, particle_data *const P);
/**
*	counts the individual particles and sets the appropriate fields in the header
*
*	\param header	pointer to an io_header structure, only the fields npart, mass and npartTotal are set
*	\param P	pointer to the set of particle data
*	\param size	number of particles
*	\param useMass 	if 1 use definitely the mass of each particle 
*	\return 	number of particles (size)
*/

int write_snapshot2(const char *const filename, const io_header *const header, particle_data *const P);

int generateHeader(io_header *const header, const particle_data *const P, const int size, const int useMass);
/**
*	write an Ifrit particle file
*	use {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0} as bounding box to keep coordinates as written in the file
*
*	\param filename	name of the file
*	\param P			pointer to the set of particle data
*	\param size		number of particles
*	\param box		boundary box of the particles (xl,yl,zl,xh,yh,zh)
*	\return			0
*/
int writeIfritParticleFile(const char *const filename, const particle_data *const P, const int size, const float box[6]);
/**
*	calculates the bounding box of a data set, as edge length, the maximum distance from
*	two particles is used to avoid stretching in ifrit particle files
*
*	\todo not done yet
*
*	\param P	pointer to the set of particle data
*	\param size	number of particles
*	\param box	boundary box of the particles (xl,yl,zl,xh,yh,zh)
*	\return		0
*/
int getBoundingBox(const particle_data *const P, const int size, float box[6]);



/**
*	load stellar age file from vacuum stellar age particles
*	loads particle id's and age
*	
*	\param filename	filename of the stellar particle filename
*	\param sp	stellar particles 
*	\param size 	number of stellar particles loaded
*	\return		0 on success
*/
int loadStellarAge(const char *const filename, particle_data * sp, int *const size);

/**
*	write stellar age file from vacuum stellar age particles
*	write particle id's and age
*	
*	\param filename	filename of the stellar particle filename
*	\param sp	stellar particles 
*	\param size 	number of stellar particles to be written
*	\return		0 on success
*/
int writeStellarAge(const char *const filename, particle_data *sp, int size);

#endif
