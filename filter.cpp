/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Filter data file according to numeric criteria

  Changes log:

  - Apr/7-JZ: 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#include <dynamo.h>

int main(int argc,char *argv[])
{
  //////////////////////////////////////////////////////////////
  //PROGRAM VARIBALES
  //////////////////////////////////////////////////////////////
  char datafile[FSIZE]="",filterfile[FSIZE]="",eq[100]="",cmd[LSIZE]="";
  real min=-MAXREAL,max=MAXREAL,eqval;
  int numcols=0,col=0;
  bool eqcmp=true;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"FILTER DATA FILE");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:s:n:c:m:M:e:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <datafile> [-s <filterfile>] -n <numcols>\n"
"\t                         -c <col> [-m <min>] [-M <Max>] | [-e <equal>]\n"
"\n"
"Filter <datafile> according to column <col> with the limits defined by <min>\n"
"and <max> or with the exact value <equal>.\n"
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
      strcpy(filterfile,optarg);
      break;
    case 'n':
      numcols=atoi(optarg);
      break;
    case 'c':
      col=atoi(optarg);
      break;
    case 'm':
      min=atof(optarg);
      break;
    case 'M':
      max=atof(optarg);
      break;
    case 'e':
      strcpy(eq,optarg);
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
  if(numcols<1){
    fprintf(stderr,"Error: The number of columns should be different from 0\n");
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(filterfile)){
    sprintf(filterfile,"%s.fil",datafile);
  }
  if(col==0){
    fprintf(stderr,"Error: You must provide a column number\nColumns available:");
    sprintf(cmd,"grep '#1:' %s",datafile);
    system(cmd);
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(eq)){
    eqcmp=false;
  }else{
    eqval=atof(eq);
    if(eqval<min || eqval>max){
      fprintf(stderr,"Error: Value to compare '%g' is outside the input range [%g,%g]\n",
	      eqval,min,max);
      PRINT_USAGE;
      EXIT;
    }
  }

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Datafile: %s\n",datafile);
    fprintf(stdout,"Filterfile: %s\n",filterfile);
    fprintf(stdout,"Number of columns: %d\n",numcols);
    fprintf(stdout,"Selected column: %d\n",col);
    fprintf(stdout,"Minimum: %+14.7e\n",min);
    fprintf(stdout,"Maximum: %+14.7e\n",max);
    if(eqcmp)
      fprintf(stdout,"Equal: %+14.7e\n",eqval);
  }
  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////
  STPRINTF("Filtering datafile '%s'...\n",datafile);
  file fd=fileOpen(datafile,"r");
  file fs=fileOpen(filterfile,"w");

  real value,*line;
  char linea[LSIZE];

  //SELECT DATA FROM DATAFILE
  line=(real*)malloc(numcols*sizeof(real));

  //Read file content
  int r=0;
  int i=0;
  int c=col-1;
  while(1){
    fgets(linea,sizeof linea,fd);
    if(feof(fd)) break;

    //Check line with comment
    if(linea[0]=='#'){
      fprintf(fs,"%s",linea);
      continue;
    }
    r++;
    //fprintf(stdout,"LINEA %d: %s",r,linea);

    //Read values in line
    readLine(linea,line,numcols);
    value=line[c];
    //fprintf(stdout,"\tCOLUMN %d VALUE: %+14.7e\n",col,value);

    //Check if filter is for an equal value
    if(!eqcmp) eqval=value;
   
    //Filter
    if(EQUAL(value,eqval) && value>=min && value<=max){
      fprintf(fs,"%s",linea);
      //fprintf(stdout,"ROW SELECTED\n");
      i++;
    }else{
      //fprintf(stdout,"VALUE NOT MATCHING (eq:%e,min:%e,max:%e) = %e\n",eqval,min,max,value);
    }
  }
  STPRINTF("\t%d lines read...\n",r);
  STPRINTF("\t%d lines selected...\n",i);
  fclose(fs);
  fclose(fd);
  STPRINTF("Filtered file '%s' saved...\n",filterfile);
  return 0;
}
