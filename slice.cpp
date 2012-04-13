/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Slice a data file choosing a given set of rows and columns

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
  char datafile[FSIZE]="",dataslice[FSIZE]="";
  int numcols=0;
  char cols[FSIZE]="";
  char rows[FSIZE]="";
  char colspec[FSIZE]="";
  real columns[100];
  real list[100];
  real selrows[3];
  int i,c;
  bool allcols=false,allrows=false;
  
  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"SLICE DATA FILE");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:s:n:c:r:C:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <datafile> [-s <dataslice>] -n <numcols>\n"
"\t                         -c <col1>,<col2>|<col1>-<col2>\n"
"\t                        [-C '<column_specification>']\n"
"\t                         -r <row1>:<row2>:<freq>\n\n"
"Take a slice of a data file with <numcols> columns. If <row1> and/or <row2>\n"
"are left empty they will be assummed equal to the minimum and maximum value\n"
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
    case 'n':
      numcols=atoi(optarg);
      break;
    case 'c':
      strcpy(cols,optarg);
      break;
    case 'r':
      strcpy(rows,optarg);
      break;
    case 'C':
      strcpy(colspec,optarg);
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
  if(isBlank(dataslice)){
    sprintf(dataslice,"%s.slc",datafile);
  }
  if(isBlank(cols)){
    sprintf(cols,"%d-%d",1,numcols);
    allcols=true;
  }
  if(isBlank(rows)){
    sprintf(rows,"-:-:1");
    allrows=true;
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COLUMNS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int ncolsp;
  if(grepString(cols,",")){
    ncolsp=splitString(cols,",",columns);
  }
  else if(grepString(cols,"-")){
    ncolsp=splitString(cols,"-",list);
    ncolsp=list[1]-list[0]+1;
    for(i=0;i<ncolsp;i++)
      columns[i]=list[0]+i;
  }
  else{
    ncolsp=1;
    columns[0]=atoi(cols);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //ROWS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int rini,rend,rfreq;
  splitString(rows,":",selrows);
  rini=DEFVAL((int)selrows[0],0,1);
  rend=DEFVAL((int)selrows[1],0,2E9);
  rfreq=DEFVAL((int)selrows[2],0,1);

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Datafile: %s\n",datafile);
    fprintf(stdout,"Dataslice: %s\n",dataslice);
    fprintf(stdout,"Number of columns: %d\n",numcols);
    fprintf(stdout,"Selected columns: %s\n",cols);
    fprintf(stdout,"\tColumns (%d): ",ncolsp);
    for(i=0;i<ncolsp;i++)
      fprintf(stdout,"%.0f ",columns[i]);
    fprintf(stdout,"\n");
    fprintf(stdout,"Rows specification: %s\n",rows);
    fprintf(stdout,"\tinit.row, final row, step row: %d, %d, %d\n",
	    rini,rend,rfreq);
    if(!isBlank(colspec))
      fprintf(stdout,"Columns specification: %s\n",colspec);
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////

  fprintf(stdout,"Slicing datafile '%s'...\n",datafile);

  file fd=fileOpen(datafile,"r");
  file fs=fileOpen(dataslice,"w");
  real value,*line;
  char linea[LSIZE],tmp[LSIZE],*pos;
  int r,rr;

  //SELECT DATA FROM DATAFILE
  line=(real*)malloc(numcols*sizeof(real));
  r=0;
  rr=0;//READING ROW

  //Header information
  fprintf(fs,"#Source file:%s\n",datafile);
  fprintf(fs,"#Columns:%s\n",cols);
  fprintf(fs,"#Rows:%s\n",rows);
  if(isBlank(colspec)){
    fprintf(fs,"#");
    for(c=0;c<ncolsp;c++){
      sprintf(colspec,"%d:C%d",(int)columns[c],(int)columns[c]);
      fprintf(fs,"%-14s ",colspec);
    }
    fprintf(fs,"\n");
  }else{
    char colnames[100][LSIZE];
    splitString(colspec," ",colnames);
    fprintf(fs,"#");
    for(c=0;c<ncolsp;c++){
      sprintf(colspec,"%d:%s",(int)columns[c],colnames[c]);
      fprintf(fs,"%-14s ",colspec);
    }
    fprintf(fs,"\n");
  }

  //Read file content
  while(1){
    fgets(linea,sizeof linea,fd);
    if(feof(fd)) break;

    //Check line with comment
    if(linea[0]=='#') continue;
    r++;
    //fprintf(stdout,"LINEA %d: %s\n",r,linea);

    //Check row
    if((r>=rini && r<=rend) && (((r-rini)%rfreq)==0)){
      rr++;
      //fprintf(stdout,"READING ROW %d:\n\t%s\n",r,linea);

      //Read all columns from line
      readLine(linea,line,numcols);

      //Save columns selected by user
      //fprintf(stdout,"\tCOLUMNS: ");
      for(i=0;i<ncolsp;i++){
	c=(int)columns[i];
	fprintf(fs,"%+14.7e ",line[c-1]);
	//fprintf(stdout,"%d:%.2e ",c,line[c-1]);
      }
      fprintf(fs,"\n");
      //fprintf(stdout,"\n");
    }else{
      //fprintf(stdout,"SKIPPING ROW %d\n",r);
    }
    //End if row is larger than rend
    if(r>rend) break;
  }
  fprintf(stdout,"\t%d lines read...\n",rr);
  fclose(fs);
  fclose(fd);
  fprintf(stdout,"Slice file '%s' saved...\n",dataslice);
  return 0;
}
