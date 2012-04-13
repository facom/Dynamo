/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Purpose: Read the data from a Gadget file and write it in a
  simulation file in ASCII format.

  Changes log:

  - Apr/7-JZ: Checked
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#include <dynamo.h>

int main(int argc,char *argv[])
{
  //////////////////////////////////////////////////////////////
  //PROGRAM VARIBALES
  //////////////////////////////////////////////////////////////
  char simulation[FSIZE]="",snapshot[FSIZE]="";
  char binfile[FSIZE]="",simfile[FSIZE]="";
  bool qsub=false;
  
  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"CONVERSION FROM BINARY TO ASCII SIMULATION");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVSs:n:d:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -s <simulation_name> -n <number_snapshot> [-d <simulation_datafile>]\n"
"\t          [-S]\n"
"\n"
"Read the data from <simulation_name>_<number_snapshot> and convert to\n"
"ascii in <simulation_datafile>\n"
"\n"
"Optional:\n"
"\n"
"\t-S : include index from substructure description file\n"
"\t     <simulation_name>.sub\n"
"=======================================================================================\n"
);
  
  //////////////////////////////////////////////////////////////
  //READ OPTIONS
  //////////////////////////////////////////////////////////////
  while(ITEROPTIONS){
    switch(OPTION){
    case 's':
      strcpy(simulation,optarg);
      break;
    case 'n':
      strcpy(snapshot,optarg);
      break;
    case 'd':
      strcpy(simfile,optarg);
      break;
    case 'S':
      qsub=true;
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
      OPTION_ERRORS;
    }
  }

  //////////////////////////////////////////////////////////////
  //VALIDATE OPTIONS
  //////////////////////////////////////////////////////////////
  if(isBlank(simulation)){
    fprintf(stderr,"Error: No simulation name was provided\n");
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(snapshot)){
    fprintf(stderr,"Error: No snapshot number provided\n");
    PRINT_USAGE;
    EXIT;
  }
  sprintf(binfile,"%s_%s",simulation,snapshot);
  if(!fileExists(binfile)){
    fprintf(stderr,"Error: Binary file '%s' does not exist\n",
	    binfile);
    PRINT_USAGE;
    EXIT;
  }
  if(qsub){
    char subfile[FSIZE];
    sprintf(subfile,"%s.sub",simulation);
    if(!fileExists(subfile)){
      fprintf(stderr,"Error: Substructure file '%s' does not exist\n",
	      subfile);
      PRINT_USAGE;
      EXIT;
    }
  }
  if(isBlank(simfile)){
    sprintf(simfile,"%s.sim",binfile);
  }
  /*
  if(fileExists(simfile)){
    fprintf(stderr,"Simulation file '%s' already exist.  Do you want to continue? ",simfile);
    pause();
  }
  */

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Simulation: %s\n",simulation);
    fprintf(stdout,"Snapshot: %s\n",snapshot);
    fprintf(stdout,"Simulation file: %s\n",simfile);
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////

  //========================================
  //READ GADGET
  //========================================
  STPRINTF("Reading data from gadget file...\n");
  SHeader sheader;
  particles parts;
  parts=readGadget(simulation,snapshot,&sheader,qsub);
  STPRINTF("\t%d particles read...\n",sheader.ntot);
  if(VERBOSE(2)) checkSimulation(sheader);

  //========================================
  //UPDATE HEADER
  //========================================
  STPRINTF("Computing properties from data read from file... ");
  updateParticles(parts,sheader);
  updateHeader(parts,&sheader);
  fprintf(stdout,"Done.\n");
  if(VERBOSE(2)) checkSimulation(sheader);
  
  //========================================
  //WRITE SIMULATION FILE
  //========================================
  STPRINTF("Writing simulation file...");
  writeSimulation(parts,sheader,simfile);
  fprintf(stdout,"Done.\n");

  return 0;
}
