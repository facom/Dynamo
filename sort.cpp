/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Sort datafile according to a given column

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
  char datafile[FSIZE]="",sortfile[FSIZE]="";
  char reverse[5]="",options[LSIZE]="";
  int col=0;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"SORT DATA FILE");
  
  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:s:c:r");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <datafile> [-s <sortedfile>] -c <col>\n"
"\t                        [-r] [-o '<options>']\n"
"\n"
"Sort <datafile> according to column <col>. Use option -r if want a reverse order\n"
"Optative <options> are additional options for the bash sort command\n"
"\n"
"=======================================================================================\n"
);

  //////////////////////////////////////////////////////////////
  //READ OPTIONS
  //////////////////////////////////////////////////////////////
  while(ITEROPTIONS){
    switch(OPTION){
    case 'f':
      strcpy(datafile,optarg);
      break;
    case 's':
      strcpy(sortfile,optarg);
      break;
    case 'c':
      col=atoi(optarg);
      break;
    case 'r':
      strcpy(reverse,"-r");
      break;
    case 'o':
      strcpy(options,optarg);
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
  if(isBlank(datafile)){
    fprintf(stderr,"Error: No datafile was provided\n");
    PRINT_USAGE;
    EXIT;
  }
  if(!fileExists(datafile)){
    fprintf(stderr,"Error: Datafile '%s' does not exist\n",datafile);
    PRINT_USAGE;
    EXIT;
  }
  if(col==0){
    col=1;
  }
  if(isBlank(sortfile)){
    sprintf(sortfile,"%s.srt",datafile);
  }

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Datafile: %s\n",datafile);
    fprintf(stdout,"Sorted file: %s\n",sortfile);
    fprintf(stdout,"Selected column: %d\n",col);
    fprintf(stdout,"Reverse: %s\n",reverse);
    fprintf(stdout,"Other options: %s\n",options);
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////
  char cmd[LSIZE];
  
  fprintf(stdout,
	  "Sorting file '%s' according to column %d...\n",
	  datafile,col);
  sprintf(cmd,"grep '#' %s > %s",datafile,sortfile);
  system(cmd);
  sprintf(cmd,
	  "grep -v '#' %s | sort -g -k %d %s %s >> %s",
	  datafile,col,reverse,options,sortfile);
  system(cmd);
  fprintf(stdout,
	  "Sorted file '%s' saved...\n",
	  sortfile);
  return 0;
}
