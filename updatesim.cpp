#include <dynamo.h>

int main(int argc,char *argv[])
{
  //////////////////////////////////////////////////////////////
  //PROGRAM VARIBALES
  //////////////////////////////////////////////////////////////
  char simfile[FSIZE]="",nsimfile[FSIZE]="";
  int i;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"UPDATE SIMFILE");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:s:n:c:b:t:N:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <simfile> [-s <new_simfile>]\n"
"\n"
"Update particle properties in simulation file <simfile>.  If no new simulation file is\n"
"provided the updated data will be stored in <simfile>.\n"
"=======================================================================================\n"
);

  //////////////////////////////////////////////////////////////
  //READ OPTIONS
  //////////////////////////////////////////////////////////////
  while(ITEROPTIONS){
    switch(OPTION){
    case 'f':
      strcpy(simfile,optarg);
      break;
    case 's':
      strcpy(nsimfile,optarg);
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
    fprintf(stderr,"Error: No simfile was provided\n");
    PRINT_USAGE;
    EXIT;
  }
  if(!fileExists(simfile)){
    fprintf(stderr,"Error: Simfile '%s' does not exist\n",simfile);
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(nsimfile)){
    strcpy(nsimfile,simfile);
  }

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Simulation file: %s\n",simfile);
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////
  particles parts;
  SHeader sheader;
  int ntot;

  //Read simulation data
  ntot=countLines(simfile);
  STPRINTF("Reading %d particles in simulation data...\n",ntot);
  parts=readSimulation(simfile,&sheader);
  sheader.ntot=ntot;

  //Update particles and header
  STPRINTF("Updating particles and header...\n");
  updateParticles(parts,sheader);
  updateHeader(parts,&sheader);
  STPRINTF("\t%d particles updated...\n",sheader.ntot);

  //Write updated simulation
  STPRINTF("Saving simulation file '%s'...\n",nsimfile);
  writeSimulation(parts,sheader,nsimfile);

  STPRINTF("Done.\n");
  return 0;
}
