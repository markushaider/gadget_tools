#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "data_ops.h"
#include "snapshot_io.h"

////////////////////////////////////////// LOCAL TYPE DEFINITIONS

//////////////////////////////////////////////// LOCAL PROTOTYPES
static int compareId(const void *a, const void *b);
static int compareType(const void *a, const void *b);
static int compareTypeId(const void *a, const void *b);

////////////////////////////////////////////////// IMPLEMENTATION

static int compareId(const void *a, const void *b)
{
     int idA = ((particle_data*)a)->Id;
     int idB = ((particle_data*)b)->Id;
     
     return ((idA < idB) ? -1 : ((idA == idB) ? 0 : 1));
}

static int compareType(const void *a, const void *b)
{
     int typeA = ((particle_data*)a)->Type;
     int typeB = ((particle_data*)b)->Type;
     
     return ((typeA < typeB) ? -1 : ((typeA == typeB) ? 0 : 1));
}

static int compareTypeId(const void *a, const void *b)
{
     int idA = ((particle_data*)a)->Id;
     int idB = ((particle_data*)b)->Id;
     int typeA = ((particle_data*)a)->Type;
     int typeB = ((particle_data*)b)->Type;
     
     return ((typeA < typeB) ? -1 : ((typeA == typeB) ? ((idA < idB) ? -1 : ((idA == idB) ? 0 : 1)) : 1));
}

int reordering(particle_data *pd, int size, enum OrderType ot)
{
     switch(ot)
     {
	  default:
	  case SORT_ID:
	  {
	       qsort(pd, size, sizeof(particle_data), compareId);
	       break;
	  }
	  case SORT_TYPE:
	  {
	       qsort(pd, size, sizeof(particle_data), compareType);
	       break;
	  }
	  case SORT_TYPE_ID:
	  {
	       qsort(pd, size, sizeof(particle_data), compareTypeId);
	       break;
	  }
     }
     
     return 0;
}


int filterICM(particle_data *const P, const int size, particle_data **nP, const int numICM, const float *const masses, const int numMasses)
{
     int i, j;
     int nsz;//, pc;
     particle_data *nP_p;
     
     if((numICM && numMasses) || !(numICM || numMasses))
     {
	  fprintf(stderr, "filterICM: use either masses or the number of particles");
	  return -1;
     }
     
     if(numICM)
     {
	  if((nsz = size - numICM) <= 0)
	  {
	       fprintf(stderr, "error in number of ICM particles\n");
	       return -1;
	  }
     }
     else
     {
	  nsz = 0;
	  for(j = 0; j < size; j++)
	  {
	       if(P[j].Type > 0)
		    nsz++;
	       else
	       {
		    for(i = 0; i < numMasses; i++)
			 if(P[j].Mass == masses[i])
			      break;
		    if(i == numMasses)
			 nsz++;
	       }
	       /*
	       for(i = 0; i < numMasses; i++)
		    if(P[j].Type == 0 && P[j].Mass == masses[i])
		    {
			 nsz++;
			 break;
		    }
	       */
	  }
     }
     
     //reordering(P, size, SORT_ID); //sort by Id
     
     if(!(nP_p = *nP = (particle_data *)malloc(sizeof(particle_data) * nsz)))
     {
	  fprintf(stderr, "error: can't allocate memory! \n");
	  exit(-2);
     }
     
     /*
     pc = 0;
     if(numICM)
     {
	  for(j = 0; j < size; j++)
	       if(numICM && P[j].Id > numICM)
		    (*nP)[pc++] = P[j];
     }
     else
     {
	  for(j = 0; j < size; j++)
	       for(i = 0; i < numMasses; i++)
		    if(P[j].Mass == masses[i])
		    {
			 (*nP)[pc++] = P[j];
			 break;
		    }
     }
     */
    
     if(numICM)
     {
	  for(j = 0; j < size; j++)
	       if(P[j].Id > numICM)
		    *nP_p++ = P[j];
     }
     else
     {
	  for(j = 0; j < size; j++)
	  {
	       if(P[j].Type == 0)	//check only for gas
	       {
		    for(i = 0; i < numMasses; i++)
		    {
			 if(P[j].Mass == masses[i])
			 {
			      //*nP_p++ = P[j];
			      break;
			 }
		    }
		    if(i == numMasses)
			 *nP_p++ = P[j];
	       }
	       else
		    *nP_p++ = P[j];
	  }
     }
    
     //for(j = numICM; j < size; j++)
	  //*nP_p++ = P[j];
     
     
     return nsz;
}

int unitConversion(particle_data *P, int size)
{
     double GRAVITY, BOLTZMANN, PROTONMASS;
     double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
     double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
     double G, Xh, HubbleParam;

     int i;
     double MeanWeight, u, gamma;
     
     //printf("unit conversion:\n");

     /* physical constants in cgs units */
     GRAVITY   = 6.672e-8;
     BOLTZMANN = 1.3806e-16;
     PROTONMASS = 1.6726e-24;

     /* internal unit system of the code */
     UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
     UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
     UnitVelocity_in_cm_per_s= 1.0e5;

     UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
     UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
     UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
     UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

     G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);


     Xh= 0.76;  /* mass fraction of hydrogen */
     HubbleParam= 0.65;

     gamma= 5.0/3;

     for(i=0; i < size; i++)
     {
	  if(P[i].Type==0)  /* gas particle */
	  {
	       MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;

	       /* convert internal energy to cgs units */

	       u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	       /* get temperature in Kelvin */

	       P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  }
     }
     
     return 0;
}

//2D & 3D
#define DISK_THICK	5.0f	//in kpc

int calcRadialDensityProfile(const particle_data *const P, const int size, const float binSize, const int numBins, float *centroid, float *data)
{
     int i, j;
     int ctr = 0;
     //float data[numBins];
     //float centroid[3] = {0.0, 0.0, 0.0};
     float totalMass = 0.0;
     float r = 0.0;
     
     
     if(centroid == NULL)
     {
	  ctr = 1;
	  centroid = (float *)malloc(sizeof(float) * 3);
	  //calc centroid -> gaseous centroid (or better old disk centroid???)
	  for(j = 0; j < size; j++)
	  {
	       //if(P[j].Type == 0)
	       {
		    centroid[0] += P[j].Pos[0] * P[j].Mass;
		    centroid[1] += P[j].Pos[1] * P[j].Mass;
		    centroid[2] += P[j].Pos[2] * P[j].Mass;
		    totalMass += P[j].Mass;
	       }
	  }
	  centroid[0] /= totalMass;
	  centroid[1] /= totalMass;
	  centroid[2] /= totalMass;
	  //centroid[0] /= size;
	  //centroid[1] /= size;
	  //centroid[2] /= size;
	  
	  //centroid[0] = 0.2;
	  //centroid[1] = 3.9;
	  
	  printf("gaseous centroid: %.3f %.3f %.3f tot mass: %f\n", centroid[0], centroid[1], centroid[2], totalMass);
	  //return 0;
     }

     /*
     for(j = 0; j < numBins; j++)
     {
	  float mass = 0.0;
	  float dl = j*binSize, dh = (j+1)*binSize;
	  float area = M_PI * (pow(dh, 2) -pow(dl, 2));
	  
	  dl = pow(dl, 2);
	  dh = pow(dh, 2);
	  
	  for(i = 0; i < size; i++)
	  {
	       float dist = pow(centroid[0]-P[i].Pos[0], 2) + pow(centroid[1]-P[i].Pos[1], 2);
	       if(dist >= dl && dist < dh)
		    mass += P[i].Mass;// *1E10;
	  }
	  data[j] = mass / area;
	  //printf("bin: %d mass: %f area: %f\n", j, mass, area);
     }
     */
     
     totalMass = 0.0;
     
     r = 0.0;
     for(j = 0; j < numBins; j++)
     {
	  float mass = 0.0;
	  int dN = 0;
	  float dl ,dh;
	  
	  dh = r + binSize/2.0;
	  dl = (j) ? r - binSize/2.0 : 0;
	  	  
	  
	  for(i = 0; i < size; i++)
	  {
	       float dist;
	       /*float dist = sqrt(pow(centroid[0]-P[i].Pos[0], 2) + pow(centroid[1]-P[i].Pos[1], 2) + pow(centroid[2]-P[i].Pos[2], 2));
	       if(dist > dl && dist < dh)
	       {
		    dN++;
		    mass += P[i].Mass;
	       }
	       */
	       
	       //only use the gas content in the disk
	       if(fabs(P[i].Pos[2]-centroid[2]) > DISK_THICK/2.0f) //not for inclined disk
		    continue;
	       
	       dist = sqrt(pow(centroid[0]-P[i].Pos[0], 2) + pow(centroid[1]-P[i].Pos[1], 2));
	       
	       if(dist >= dl && dist < dh)
	       {
		    dN++;
		    mass += P[i].Mass;
		    totalMass += P[i].Mass;
	       }
	  }
	  
	  //data[j] = dN/(4*M_PI*r*r);
	  
	  data[j] = mass*1e10/((dh*dh - dl*dl) * M_PI); //surface density
	  //data[j] = mass*1e10/((dh*dh*dh - dl*dl*dl)*4/3*M_PI); //volumetric density in M_sun/kpc^3
	  
	  r += binSize;
	  
	  	  
     }
     
     //write to file
     /*FILE *fd = fopen("radDensProf.txt", "w");
     //fprintf(fd, "# radial density profile\n");
     for(j = 0; j < numBins; j++)
     {
	  fprintf(fd, "%f %e\n", j*binSize, data[j]);
     }
     fclose(fd);
     */
     
     printf("overall total mass in the disk: %f\n", totalMass);
     
     if(ctr)
	  free(centroid);
     
     return 0;
}

int inDisk(const double centroid[3], const float r, const float h, float phi, const float pos[3])
{
     float d[3];
     
     phi *= M_PI/180.0 * (-1.0);
     d[0] = pos[0] - centroid[0]; //x
     d[1] = (pos[1] - centroid[1])*cos(phi) - (pos[2] - centroid[2])*sin(phi); //y
     d[2] = (pos[1] - centroid[1])*sin(phi) + (pos[2] - centroid[2])*cos(phi); //z
     
     
     //common case
     /*tmp[0] = point[0]*cos(alpha)*cos(beta) +
		    point[1]*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) +
		    point[2]*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
     tmp[1] = point[0]*sin(alpha)*cos(beta) +
		    point[1]*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) +
		    point[2]*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
     tmp[2] = point[0]*(-sin(beta)) + point[1]*cos(beta)*sin(gamma) + point[2]*cos(beta)*cos(gamma);
     */
     
     //printf("d: %f %f %f\n", d[0], d[1], d[2]);
     return ((d[0]*d[0]+d[1]*d[1] <= r*r) && (fabs(d[2]) <= h/2.0));
}

/**	used for the binary search algorithm from the c standard library
*	\param key_p	a pointer to the key (element to search for)
*	\param el_p	a pointer to an array element
*	\return		a value less, equal or greater than 0 if key is less, equal or greater than el
*/
static int bs_compare(const void *key_p, const void *el_p)
{
     int key = *((const int *)key_p);
     int el = *((const int *)el_p);
     
     return ((key > el) ? 1 : ((key < el) ? -1 : 0));
}


//assume ICM already filtered
#define CENTROID_TYPE		2	///0 gas, 1 halo, 2 stars
int calcWakeDiskRatios(const particle_data *const P, const io_header *const header, const float wakeDist, InvestData *const idata, const float initNsMass, const float phi, const float h, const int first_run)
{
     int i, j, num_part = 0;
     
     static int *id = 0;
     static int size = 0;
     
     int *id_p;
     
     //get total number of particles
     for(j = 0; j < 6; j++)
	  num_part += header->npartTotal[j];
     //printf("number of particles: %d\n", num_part);
     
     if(first_run)
     {
	  /*if((id_p = id = malloc(sizeof(int)*header->npartTotal[4])) <= 0)
	  {
	       fprintf(stderr, "calcWakeDiskRatio: unable to malloc!\n");
	       exit(-1);
	  }
	  */
	       
	  idata->time = header->time;
	  idata->ratio_ns_wd = 0.0;
	  idata->ratio_ns_wt = 0.0;
	  idata->ratio_gas = 0.0;
	  idata->ratio_gastemp = 0.0;
     }
     else
     {
	  double total_mass[5] = {0.0};
	  //double total_ns_mass = 0.0;
	  //double wake_ns_mass = 0.0;
	  double wake_mass[5] = {0.0};
	  double centroid[3] = {0.0};
	  //int num_newpart = 0;
	  
	  //get centroid (gas or old stars? -> use def)
	  //get total masses per type
	  for(j = 0; j < num_part; j++)
	  {
	       if(P[j].Type == CENTROID_TYPE)
	       {
		    for(i = 0; i < 3; i++)
			 centroid[i] += P[j].Pos[i]*P[j].Mass;
	       }
	       
	       total_mass[P[j].Type] += P[j].Mass;
	  }
	  for(i = 0; i < 3; i++)
	       centroid[i] /= total_mass[CENTROID_TYPE];
	  
	  //get total mass of newly formed stars in the last snapshot
	  //get total - " - in the wake
	  for(j = 0; j < num_part; j++)
	  {
	       /*if((P[j].Type == 4) && (bsearch(&(P[j].Id), id, size, sizeof(int), bs_compare) == (void *)0))	//a new particle?
	       {
		    total_ns_mass += P[j].Mass;
		    if(fabs(P[j].Pos[2]-centroid[2]) > 10.0)
			 wake_ns_mass += P[j].Mass;
		    num_newpart++;
	       }*/
	       //if(fabs(P[j].Pos[2]-centroid[2]) > 10.0) //////// nur für face on!!!
	       //	wake_mass[P[j].Type] += P[j].Mass;
	       if(!inDisk(centroid, 150, h, phi, P[j].Pos))
		    wake_mass[P[j].Type] += P[j].Mass;
	       
	  }
	  
	  idata->time = header->time;
	  
	  //idata->ratio_ns_wd = wake_ns_mass/(total_ns_mass - wake_ns_mass);
	  //idata->ratio_ns_wd = wake_ns_mass/total_ns_mass;
	  //printf("centroid x,y,z [kpc]: %f, %f, %f\n", centroid[0], centroid[1], centroid[2]);
	  //printf("nnpart: %d total: %d\n", num_newpart, header->npartTotal[4]);
	  //printf("tot: %f wake: %f\n\n", total_ns_mass, wake_ns_mass);
	  
	  idata->ratio_ns_wd = wake_mass[4]/(total_mass[4]-wake_mass[4]);
	  idata->ratio_ns_wt = wake_mass[4]/total_mass[4];
	  idata->ratio_gas = wake_mass[0]/total_mass[0];
	  idata->ratio_gastemp = 0; 
	  
	  printf("time: %f wake_mass: %e\n", header->time, wake_mass[4]);
	  
	  /*
	  for(j = 0; j < size; j++)
	  {
	       //if(!inDisk(centroid, 100, h, phi, P[j].Pos)) //TODO automatisch mit Dichteprofil?
	       if(fabs(P[j].Pos[2]-centroid[2]) > 10.0)
	       wake_mass[P[j].Type] += P[j].Mass;
	  }
	  
	  idata->time = header->time;
	  idata->ratio_ns_wd = wake_mass[4]/(total_mass[4]-wake_mass[4]);
	  idata->ratio_ns_wt = wake_mass[4]/total_mass[4];
	  idata->ratio_gas = wake_mass[0]/total_mass[0];
	  idata->ratio_gastemp = 0; 
	  */
     }
     
     //save new particles
     size = header->npartTotal[4];
     if((id_p = id = realloc(id, sizeof(int)*size)) <= 0)
     {
	  fprintf(stderr, "calcWakeDiskRatio: unable to realloc memory!\n");
	  exit(-1);
     }
     for(j = 0; j < num_part; j++)
     {
	  if(P[j].Type == 4)
	       *id_p++ = P[j].Id;
     }
	       
     qsort(id, size, sizeof(int), bs_compare);
     
     //printf("wd: %f :: %f %f\n", idata->ratio_ns_wd, wake_mass[4], total_mass[4]);
     //printf("centroid: %f %f %f\n", centroid[0], centroid[1], centroid[2]);
     
     return 0;
     
     ////////////////////////////////////////////////
/*     
     float total_gas_mass = 0.0;	///the total gas mass
     float total_star_mass = 0.0;	///the total mass of stars
     float total_newstar_mass = 0.0;	///the total mass of new formed stars
     
     float centroid[3] = {0.0, 0.0, 0.0};	///old disk population centroid
     
     float mass_star_disk = 0.0;
     float mass_star_wake = 0.0;
     float mass_gas_wake = 0.0;
     float mass_gas_disk = 0.0;
     
     float mass_hotgas = 0.0;
     
     float mean_gastemp_wake = 0.0;
     float mean_gastemp_disk = 0.0;
     int num_gas_disk = 0, num_gas_wake = 0;
     
     static float nsmass_disk_prevSnap = 0.0;
     static float nsmass_wake_prevSnap = 0.0;
     
     ////////////////////////////////////////////////

     for(j = 0; j < 6; j++)
	  size += header->npartTotal[j];
     
     //get total masses and old stellar centroid
     for(j = 0; j < size; j++)
     {
	  //get mass
	  switch(P[j].Type)
	  {
	       case 0:
	       {
		    total_gas_mass += P[j].Mass;
		    break;
	       }
	       case 2:
	       {
		    int i;
		    total_star_mass += P[j].Mass;
		    for(i = 0; i < 3; i++)
			 centroid[i] += P[j].Pos[i]*P[j].Mass;
		    break;
	       }
	       case 4:
	       {
		    total_newstar_mass += P[j].Mass;
		    break;
	       }
	       default:
		    break;
	  }
     }
     

     if(total_star_mass == 0.0)
	  return -1;
     for(j = 0; j < 3; j++)
	  centroid[j] /= total_star_mass;
     
     //total_star_mass += total_newstar_mass;
     
     //use a cylinder to distinguish wake and disk
     
     
     //printf("total star mass [1E10 sun masses]: %f\n", total_star_mass);
     //printf("centroid x,y,z [kpc]: %f, %f, %f\n", centroid[0], centroid[1], centroid[2]);
     
     
     for(j = 0; j < size; j++)
     {
	  if(P[j].Type == 4)
	  {
	       if(fabs(P[j].Pos[2] - centroid[2]) < wakeDist)
		    mass_star_disk += P[j].Mass;
	       else
		    mass_star_wake += P[j].Mass;
	  }
	  else if(P[j].Type == 0)
	  {
	       if(fabs(P[j].Pos[2] - centroid[2]) < wakeDist)
	       {
		    mass_gas_disk += P[j].Mass;
		    mean_gastemp_disk += P[j].Temp;
		    num_gas_disk++;
	       }
	       else
	       {
		    mass_gas_wake += P[j].Mass;
		    mean_gastemp_wake += P[j].Temp;
		    num_gas_wake++;
	       }
	       
	       if(P[j].Temp > HOT_GAS_TEMP)
		    mass_hotgas += P[j].Mass;
	  }
     }
     
     if(num_gas_disk)
	  mean_gastemp_disk /= (float)num_gas_disk;
     else 
	  mean_gastemp_disk = 0.0;
     
     if(num_gas_wake)
	  mean_gastemp_wake /= (float)num_gas_wake;
     else
	  mean_gastemp_wake = 0.0;
     
     //printf("t:%f msd:%f msw:%f totgas:%f mgas:%f hotgas:%f\n", header->time, mass_star_disk, mass_star_wake, total_gas_mass, mass_gas_wake, mass_hotgas);
     
     //set calculated data
     idata->time = header->time;     
     idata->ratio_ns_wd = (mass_star_disk - initNsMass != 0.0) ? mass_star_wake / (mass_star_disk - initNsMass) : 0.0;
     idata->ratio_ns_wt = (total_newstar_mass - initNsMass != 0.0) ? (mass_star_wake) / (total_newstar_mass - initNsMass) : 0.0;
     //idata->ratio_ns_wd_ss = (total_star_mass - nsmass_disk_prevSnap != 0.0) ? (mass_star_wake - nsmass_wake_prevSnap) / (total_star_mass - nsmass_disk_prevSnap) : 0.0;
     idata->ratio_ns_wd_ss = total_star_mass - nsmass_disk_prevSnap - nsmass_wake_prevSnap;
     idata->ratio_gas = (total_gas_mass != 0.0) ? mass_gas_wake / total_gas_mass : 0.0;
     idata->ratio_gastemp = (mean_gastemp_disk != 0.0) ?mean_gastemp_wake / mean_gastemp_disk : 0.0;
     
     //
     nsmass_disk_prevSnap = mass_star_disk;
     nsmass_wake_prevSnap = mass_star_wake;
     //

     return 0;
     
*/
}

#undef CENTROID_TYPE


//#ifdef HSML
float calcMeanSL(const particle_data *const P, const int size)
{
     int j, count = 0;
     float sk = 0.0;
     
     for(j = 0; j < size; ++j)
	  if(P[j].Type == 0)
	  {
	       sk += P[j].Hsml;
	       ++count;
	  }
	  
     return (count) ? sk/count : 0.0;
}
//#endif

static inline void getWeights(const particle_data *const P, const float *const pos, const float *const s, float *const weights)
{
     int i;
     
     for(i = 0; i < 3; i++)
	  weights[i] = 1.0 - fabs(P->Pos[i] - pos[i])/s[i];
     
     if((weights[3] = weights[0]*weights[1]*weights[2]) < 0.0)
	  weights[3] = 0.0;
}

GridData* CICmethod(const particle_data *const P, const int size, int type, GridData **gd_p)
{
     float centroid[3] = {0.0, 0.0, 0.0};
     float totStarMass = 0.0;
     int i, j, k;//, num = 0;
     char fn_gastemp[200], fn_nsdens[200];
     FILE *fd_gastemp, *fd_nsdens;
     //particle_data *tmp;
     
     //
     //int resolution = 400;
     float sizeCube[3] = {300.0, 300.0, 300.0};	//size of the cube that contains the data
     int sz = 1;
     float center[3] = {sizeCube[0]*0.5, sizeCube[1]*0.5, sizeCube[2]*0.9};
     int n[3] = {300, 300, 300};	//number of elements in the cube
     float s[3];			//size of a single cube element
     //float weight[4];
     
     //float min = FLT_MAX, max = FLT_MIN;
     
     GridData *gridData;
     //
     
     //init
     for(i = 0; i < 3; i++)
     {
	  //n[i] = resolution;
	  s[i] = sizeCube[i]/(float)n[i];
	  //sz *= sizeCube[i];
	  sz *= n[i];
     }
     *gd_p = gridData = (GridData *)calloc(sz, sizeof(GridData));
     if(!gridData)
     {
	  fprintf(stderr, "unable to allocate memory!\n");
	  return NULL;
     }
     //memset(gridData, 0, sizeof(struct GridData)*sz);
     //printf("gd : %f %f %d\n", gridData[20].temperature, gridData[20].rho, gridData[20].counter);
     
     //get centroid
     for(j = 0; j < size; j++)
	  if(P[j].Type == 2)
	  {
	       totStarMass += P[j].Mass;
	       for(i = 0; i < 3; i++)
		    centroid[i] += P[j].Pos[i] * P[j].Mass;
	  }
     for(i = 0; i < 3; i++)
	  centroid[i] /= totStarMass;
     
     //calc
     //tmp = (particle_data *)malloc(sizeof(particle_data)*size);
     for(j = 0; j < size; j++)
     {
	  particle_data tmp = P[j];
	  //set particles
	  for(i = 0; i < 3; i++)
	       tmp.Pos[i] = P[j].Pos[i] - centroid[i] + center[i];
	  
	  if ((tmp.Mass == 0.0) && (tmp.Pos[0] < 0.0) && 
	       (tmp.Pos[1] < 0.0) && (tmp.Pos[2] < 0.0)) 
	       continue;
	  
	  if(tmp.Type == type)
	  {
	       int pp[3];
	       float pos[3];
	       float weights[4];
	       int index;
	       
	       for(i = 0; i < 3; i++)
	       {
		    pp[i] = (int)(floor(((tmp.Pos[i])/s[i])));
		    pos[i] = s[i]*((float)pp[i]);
	       }
	       if ((pp[0] < 0) || (pp[0] >= (n[0]-1)) 
		    || (pp[1] < 0) || (pp[1] >= (n[1]-1))
                    || (pp[2] < 0) || (pp[2] >= (n[2]-1))) 
		    continue; //not in the grid
	       //
	       getWeights(&tmp, pos, s, weights);
	       index = pp[0] + pp[1]*n[0] + pp[2]*n[0]*n[1];
	       gridData[index].temperature += (double)tmp.Temp; //if(tmp.Temp < 0.0) printf("neg temp.!\n");
	       gridData[index].u += (double)tmp.U; 
	       //gridData[index].newStarDensity += (double)((tmp.Mass * weights[3])/(s[0]*s[1]*s[2]));
	       gridData[index].density += (double)(tmp.Mass);// /(s[0]*s[1]*s[2]));
	       gridData[index].counter++;
	       //gridData[index].rho += (double)(tmp.Mass*weights[3]/(s[0]*s[1]*s[2]));
	  }
     }
     
     //
     for(j = 0; j < sz; j++)
     {
	  gridData[j].temperature = ((gridData[j].counter == 0) ? 0.0 : (gridData[j].temperature/gridData[j].counter));
	  gridData[j].u = ((gridData[j].counter == 0) ? 0.0 : (gridData[j].u/gridData[j].counter));
	  gridData[j].density /= (s[0]*s[1]*s[2]);
     }
	  
//////////////////////////////
//#if defined(CIC_2D) || defined(CIC_3D)
     sprintf(fn_gastemp, "gastemp_part.raw");
     sprintf(fn_nsdens, "ns_density.raw");
     if(!(fd_gastemp = fopen(fn_gastemp, "w")))
     {
	  fprintf(stderr, "Unable to open file!\n");
	  free(gridData);
	  return NULL;
     }
     if(!(fd_nsdens = fopen(fn_nsdens, "w")))
     {
	  fprintf(stderr, "Unable to open file!\n");
	  fclose(fd_gastemp);
	  free(gridData);
	  return NULL;
     }
//#endif
     
//#ifdef CIC_2D
     //2D projection on y-z
     printf("calc CIC_2D!\n");
     int cmax = 0, cother = 0;
     for(i = 0; i < n[1]; i++)
     {
	  for(j = 0; j < n[2]; j++)
	  {
	       double tmpTemp = 0.0;
	       double tmpNsd = 0.0;
	       //tempField[i][j] = 0.0;
	       int index = n[1]*i + n[2]*n[1]*j;
	       for(k = 0; k < sizeCube[0]; k++)
	       {
		    //tempField[i][j] += gridData[index++].temperature;
		    tmpTemp += gridData[index].temperature;
		    tmpNsd += gridData[index++].density;
	       }
	       
	       //
	       //if(tmpNsd > 0) printf("nsd: %lf\n", tmpNsd);
	       
	       //TODO nur für die neu geformten Sterne des letzten snapshots (oder mehrerer snapshots)
		    
	       if((tmpNsd *= 1000) > 0.04)
		    ;//tmpNsd = 0.04, cmax++;
	       else
		    if(tmpNsd > 0.0)
			 cother++;
	       
	       /*
	       if(tmpNsd > max)
		    max = tmpNsd;
	       if(tmpNsd < min && tmpNsd > 0.0)
		    min = tmpNsd;
	       */
	       
	       //
		    
	       fwrite(&tmpTemp, sizeof(double), 1, fd_gastemp);
	       fwrite(&tmpNsd, sizeof(double), 1, fd_nsdens);
	       //printf("write: %d\n", j);
	  }
     }
     printf("cmax: %d coth: %d\n", cmax, cother);
     //printf("min: %.10f max: %.10f\n", min, max);
     
     
//#elif CIC_3D 
/*
     //3D
     printf("\t\tsz: %d\n", sz);
     //double tempField[n[0]][n[1]][n[2]]; //TODO malloc
     printf("\t\t--->\n");
     for(i = 0; i < n[0]; i++)
	  for(j = 0; j < n[1]; j++)
	       for(k = 0; k < n[2]; k++)
	       {
		    int index = i + j*n[1] + k*n[1]*n[2];
		    if(index >= sz)
			 printf("xxx\n"), exit(-1);
		    //tempField[i][j][k] = gridData[index].temperature;
		    //fwrite(&(gridData[index].temperature), sizeof(double), 1, fd);
	       }
	       
     printf("\t\tend CIC_3d\n");
#endif
*/

//#if defined(CIC_2D) || defined(CIC_3D)
     fclose(fd_nsdens);
     fclose(fd_gastemp);
//#endif

     //write to file
/*     sprintf(fn, "gastemp.raw");
     if(!(fd = fopen(fn, "w")))
	  fprintf(stderr, "Unable to open file!\n");
     else
     {
#ifdef CIC_2D	  
	  printf("\twrite file <%s> of size %ldB\n", "gastemp.raw", sizeof(float)*n[1]*n[2]);
	  for(i = 0; i < n[1]; i++)
	       for(j = 0; j < n[2]; j++)
		    fwrite(&(tempField[i][j]), sizeof(double), 1, fd);
	       
	  fclose(fd);
#elif CIC_3D
	  //3D
	  //printf("\twrite file <%s> of size %ldB\n", "gastemp3D.raw", sizeof(float)*n[0]*n[1]*n[2]);
	  
#endif
     }

*/


     //clean up -> done in calling function
     //free(gridData);
     
     return gridData;
}
