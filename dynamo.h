//////////////////////////////////////////////////////////////////////////////////
//EXTERNAL STANDARD LIBRARIES
//////////////////////////////////////////////////////////////////////////////////
#ifndef __DYNAMO__
#define __DYNAMO__

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstddef>
#include<cmath>
#include<cstring>
#include<unistd.h>

//////////////////////////////////////////////////////////////////////////////////
//USEFUL LIBRARIES
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//GSL MODULES
//////////////////////////////////////////////////////////////////////////////////
#include<gsl/gsl_math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_statistics_int.h>
#include<gsl/gsl_statistics_int.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>

//////////////////////////////////////////////////////////////////////////////////
//GLOBAL MACROS
//////////////////////////////////////////////////////////////////////////////////

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//NUMERICAL
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define POW gsl_pow_int

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//PHYSICAL
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define CLIGHT 3E10 //cm/s
#define YEAR 365.0*24*60*60 //s
#define MSUN 1.989E33 //g
#define KPC 3.2616E3*CLIGHT*YEAR //cm
#define GYR 1E9*YEAR //s
#define KPS 1E5 //cm/s

#define UM 1E10*MSUN 
#define UL 1*KPC 
#define UV 1*KPS
#define UT 0.9785*GYR
#define GGAD 43019.1 //UL^3/(UM UT^2)
#define GCONST GGAD*UL*UL*UL/(UM*UT*UT)

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//FILES
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define FSIZE 100 //MAXIMUM SIZE OF THE NAME OF A FILE
#define LSIZE 500 //MAXIMUM SIZE OF THE LINE OF A DATAFILE
#define NSUB 50 //MAXIMUM NUMBER OF SUBSTRUCTURES
#define DUMMY(stream) fread(&dummy,sizeof(dummy),1,stream);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//VERBOSITY
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define BAR_SIZE 80
#define BAR(stream,character) for(int i=0;i<BAR_SIZE;i++){fprintf(stream,"%c",character);}fprintf(stream,"\n");
#define EXIT exit(1);
#define VERBOSE(n) (VERBOSITY>=n?1:0)
#define TITLE(stream,character,title) BAR(stream,character);fprintf(stream,title);fprintf(stream,"\n");BAR(stream,character);
#define VPRINTF(n) if(VERBOSE(n)) printf
#define STPRINTF fprintf(stdout,"--> ");printf

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//OPTIONS MACROS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define SET_OPTIONS(opt) strcpy(OPTIONS,opt)
#define SET_USAGE(opt) strcpy(USAGE,opt)
#define ITEROPTIONS ((OPTION=getopt(argc,argv,OPTIONS))!=-1)
#define OPTION_ERRORS \
case 'h':\
PRINT_USAGE;\
exit(0);\
break;\
case ':':\
fprintf(stderr,"Option -%c requires and argument\n",optopt);\
PRINT_USAGE;\
exit(1);\
break;\
default:\
fprintf(stderr,"Unrecognized option: -%c\n",optopt);\
PRINT_USAGE;\
exit(1);\
break;
#define PRINT_USAGE fprintf(stderr,"%s",USAGE);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//MISCELLANEOUS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define DEFVAL(x,val,def) ((x)==(val)?(def):(x))
#define SET_MIN(x,min) (x=x<min?x:min)
#define SET_MAX(x,max) (x=x>max?x:max)

#define MINREAL -1E-38
#define MAXREAL +1E+38
#define EPSREAL 1E-7

#define MAXLAYERS 10
#define MAXGRID 101

#define EQUAL(x,y) (gsl_fcmp(x,y,EPSREAL)==0)

//COMPUTE THE OFFSET OF A MEMBER IN STRUCTURE PARTICLE
#define PART_MEMBER(member) offsetof(particle,member)
//GET A MEMBER OF THE STRUCTURE object OF TYPE type FROM ITS OFFSET offset
#define GET_MEMBER(object,type,offset) (*(type*)((char*)&object+offset))

//////////////////////////////////////////////////////////////////////////////////
//CUSTOM DATA TYPES AND ENUMERATORS
//////////////////////////////////////////////////////////////////////////////////
typedef float real;
typedef double real2;
typedef double* darray;
typedef char* cadena;
typedef FILE* file;
typedef FILE* File;
typedef void* pars;

typedef struct{
  real x,y,z;
} vector;

typedef struct
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npart_total[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam; 
  char fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];/*fills to 256 Bytes*/
} GHeader;

typedef struct
{
  real mtot;
  int ntot;
  int npartg[6];
  int nsubs,stype[NSUB];
  char subtype[NSUB][100];
  int nparts[NSUB];
  real mtotg[6];
  real mtots[NSUB];
  vector jtot;
  vector rcm,vcm;
  real z,omegaL,Ho,boxsize; 
  bool qupdate;//TRUE IF IT REQUIRES AN UPDATE
  int rf;//REFERENCE FRAME: 0-Cosmological simulation, 1-Rotated and translated
} SHeader;

typedef struct {
  int sid;//SORT ID (WHEN PARTICLES ARE SORTED)
  int pid;//GADGET ID
  int tid;//TYPE ID ACCORDING TO SUBSTRUCTURE
  int gtype;//GADGET TYPE:0,1,2,3,4,5
  real m,h;//MASS AND SOFTENING FACTOR
  vector pos,vel;//POSITION AND VELOCITY IN
  vector posc,velc;//POSITION AND VELOCITIES IN CYLINDRICAL COORDINATES
  vector j;//ANGULAR MOMENTUM RESPECT CENTER OF BOX
  real r,vcirc,jcirc,epsilon;//KINEMATIC INFORMATION
} particle;

typedef particle* particles;

//////////////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES FROM OTHER MODULES
//////////////////////////////////////////////////////////////////////////////////
extern char *optarg;

//////////////////////////////////////////////////////////////////////////////////
//MISCELLANEA
//////////////////////////////////////////////////////////////////////////////////
using namespace std;

//////////////////////////////////////////////////////////////////////////////////
//CUSTOM FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////
#include <functions.h>

//////////////////////////////////////////////////////////////////////////////////
//ROUTINE PROTOTYPES AND GLOBAL VARIABLES
//////////////////////////////////////////////////////////////////////////////////

extern char OPTIONS[100];
extern char OPTION;
extern char USAGE[1000];
extern int VERBOSITY;
extern real2 *XFit,*YFit,*EFit;
extern int NFit;
extern real2 (*FFit)(real2,pars);


particles readGadget(char simulation[],char snapshot[],SHeader* sheader,bool qsub=false);

int updateParticles(particles parts,SHeader sheader) ;

int updateHeader(particles parts,SHeader *sheader) ;

int checkSimulation(SHeader sheader);

int writeSimulation(particles parts,SHeader sheader,char simfile[],int frac=1) ;

particles readSimulation(char simfile[],SHeader* sheader) ;

cadena str2particle(char str[],particle *part) ;

cadena particle2str(particle part) ;

int readGadgetParam(char paramfile[],const char variable[],char type,void *result);

int transformSimulation(particles parts,SHeader *sheader,vector ro,vector vo,vector jo);

real* sliceSimulation(particles parts,int ini,int end,size_t member);

int constVector(vector* v,real k);

int scaleVector(vector* v,real k,vector u);

int sumVector(vector* v,real k,vector u1,vector u2);

int crossVector(vector* v,vector u1,vector u2);

real distVector(vector v1,vector v2);

int copyVector(vector* v,vector u);

real magVector(vector v);

int cylVector(vector *v,vector u);

cadena vector2str(vector r);

int arrayPrintf(file stream,const char frm[],real v[],int n);

int arrayVector(vector *v,real a[]);

real arrayMin(real* array,int n);

real arrayMax(real* array,int n);

int fileExists(char *filename);

int isBlank(char *string);

int isSpace(char *string);

int pause(void);

file fileOpen(const char filename[],const char mode[]);

int splitString(char str[],const char sep[],real list[]);

int splitString(char str[],const char sep[],real2 list[]);

int splitString(char str[],const char sep[],char list[][LSIZE]);

int grepString(char str[],const char search[]);

int grepArray(real value,real array[],int n,real tol/*=1E-5*/);

int readLine(char row[],real array[],int nfields);

int readLine(char row[],real2 array[],int nfields);

int systemOutput(const char cmd[],char out[]);

int countLines2(const char file[]);

int countLines(const char file[]);

real2* readColumn(const char fname[],int nlines,int ncols,int col);

real* readColumns(const char fname[],int nlines,int ncols,int col);

real2 chiSquare(const gsl_vector* param,pars tmp);

int fitData(real2* X,real2* Y,real2* E,int Ndata,real2 (*func)(real2,pars),int Npars,real2 Params[],real2& cs,real2& pval);

real gravSoft(real x,real h);

#endif