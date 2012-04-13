/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Purpose: Transform simulation data to another reference system

  Changes log:

  - Apr/7-JZ: Checked

  * I have performed checks with a mock disk generated with
    gendata.py, i.e. TEST_000.sim.  We recommend to check when changes
    be introduced in code

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#include <dynamo.h>

int main(int argc,char *argv[])
{
  //////////////////////////////////////////////////////////////
  //PROGRAM VARIBALES
  //////////////////////////////////////////////////////////////
  char simfile[FSIZE]="",transfile[FSIZE]="";
  char origin[130]="",velocity[130]="",orientation[130]="";
  real list[3];
  vector ro,vo,jo;
  int rflabel;
  //real list[3];

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"TRANSFORM TO ANOTHER COORDINATE SYSTEM");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVs:t:o:e:j:l:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -s <simulation_file> [-t <transformation_file>]\n"
"\t                               -o <xo>,<yo>,<zo> -e <vox>,<voy>,<voz>\n"
"\t                               -j <jx>,<jy>,<jz> [-l <rf_label>]\n"
"\n"
"Read the data from <simulation_file>_and transform it to the reference frame centered at\n"
"<xo>,<yo>,<zo> with velocity <vox>,<voy>,<voz> and orientation according to the set of\n"
"axis <jx>,<jy>,<jz>.\n"
"A reference frame label <rf_label> can be set.  Default 1.\n"
"=======================================================================================\n"
);

  //////////////////////////////////////////////////////////////
  //READ OPTIONS
  //////////////////////////////////////////////////////////////
  while(ITEROPTIONS){
    switch(OPTION){
    case 's':
      strcpy(simfile,optarg);
      break;
    case 't':
      strcpy(transfile,optarg);
      break;
    case 'o':
      strcpy(origin,optarg);
      break;
    case 'e':
      strcpy(velocity,optarg);
      break;
    case 'j':
      strcpy(orientation,optarg);
      break;
    case 'l':
      rflabel=atoi(optarg);
      break;
    //========================================
    //COMMON
    //========================================
    case 'v':
      VERBOSITY=1;
      break;
    case 'V':
      VERBOSITY=2;
      break;
    //DETECT ERRORS
    OPTION_ERRORS;
    }
  }

  //////////////////////////////////////////////////////////////
  //VALIDATE OPTIONS
  //////////////////////////////////////////////////////////////
  if(isBlank(simfile)){
    fprintf(stderr,"Error: No simulation file name was provided\n");
    PRINT_USAGE;
    EXIT;
  }
  if(!fileExists(simfile)){
    fprintf(stderr,"Error: Datafile '%s' does not exist\n",simfile);
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(transfile)){
    strcpy(transfile,simfile);
  }
  if(isBlank(origin))
    strcpy(origin,"0.0,0.0,0.0");
  if(isBlank(velocity))
    strcpy(velocity,"0.0,0.0,0.0");
  if(isBlank(orientation))
    strcpy(orientation,"0.0,0.0,1.0");

  //ORIGIN AND ORIENTATION VECTOR
  splitString(origin,",",list);
  arrayVector(&ro,list);
  splitString(velocity,",",list);
  arrayVector(&vo,list);
  splitString(orientation,",",list);
  arrayVector(&jo,list);
  
  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Simulation file: %s\n",simfile);
    fprintf(stdout,"Origin: %s\n",vector2str(ro));
    fprintf(stdout,"Velocity: %s\n",vector2str(vo));
    fprintf(stdout,"Orientation: %s\n",vector2str(jo));
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////
  
  //==================================================
  //READ THE SIMULATION DATA
  //==================================================
  particles parts;
  SHeader sheader;
  STPRINTF("Reading simulation data...\n");
  parts=readSimulation(simfile,&sheader);
  if(VERBOSE(2)) checkSimulation(sheader);
  STPRINTF("\t%d particles read...\n",sheader.ntot);

  //==================================================
  //TRANSFORM
  //==================================================
  STPRINTF("Transforming data to new reference frame...\n");

  VPRINTF(2)("Particle 0 before transformation: %s\n",particle2str(parts[0]));
  transformSimulation(parts,&sheader,ro,vo,jo);
  sheader.rf=rflabel;
  VPRINTF(2)("Particle 0 after transformation: %s\n",particle2str(parts[0]));
  if(VERBOSE(2)) checkSimulation(sheader);

  //==================================================
  //WRITE TRANSFORMATION
  //==================================================
  STPRINTF("Writing %d particles in file '%s'...",sheader.ntot,transfile);
  writeSimulation(parts,sheader,transfile);
  fprintf(stdout,"Done.\n");

  return 0;
}
