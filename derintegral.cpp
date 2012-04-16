/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Given a set of data compute the first and second derivative
  and the integral.

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
  char datafile[FSIZE]="",derintfile[FSIZE]="";
  int numcols=0,colx=0,coly=0;
  int i,ndata;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"COMPUTE DERIVATIVES AND INTEGRAL OF DATA");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:d:n:x:y:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <datafile> [-d <derint_file>] -n <numcols> -x <colx> -y <coly>\n"
"\n"
"Take columns <colx> and <coly> from <datafile> and computes the first\n"
"and second derivative but also the integral in the interval.  The\n"
"results are stored in file <derint_file>\n"
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
    case 'd':
      strcpy(derintfile,optarg);
      break;
    case 'n':
      numcols=atoi(optarg);
      break;
    case 'x':
      colx=atoi(optarg);
      break;
    case 'y':
      coly=atoi(optarg);
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
  if(isBlank(derintfile)){
    sprintf(derintfile,"%s.der",datafile);
  }
  if((ndata=countLines(datafile))==0){
    fprintf(stderr,"Error: Datafile '%s' seems empty\n",datafile);
    PRINT_USAGE;
    EXIT;
  }
  if(numcols<1){
    fprintf(stderr,"Error: The number of columns should be different from 0\n");
    PRINT_USAGE;
    EXIT;
  }
  if(colx==0){
    colx=0;
  }
  if(coly==0){
    coly=1;
  }

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    BAR(stdout,'O');
    fprintf(stdout,"Datafile: %s\n",datafile);
    fprintf(stdout,"Der. Int. file: %s\n",derintfile);
    fprintf(stdout,"Number of columns: %d\n",numcols);
    fprintf(stdout,"Columns x,y: %d,%d\n",colx,coly);
    fprintf(stdout,"Number of points in datafile: %d\n",ndata);
    BAR(stdout,'O');
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING DATA
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  STPRINTF("Reading %d data points from datafile '%s'...\n",ndata,datafile);
  real *X=readColumns(datafile,ndata,numcols,colx);
  real *Y=readColumns(datafile,ndata,numcols,coly);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COMPUTING DERIVATIVES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  STPRINTF("Computing derivatives and integral...\n");
  file fd=fileOpen(derintfile,"w");
  fprintf(fd,"#Datafile: %s\n",datafile);
  fprintf(fd,"#NumberColumns: %d\n",numcols);
  fprintf(fd,"#ColX,ColY: %d %d\n",colx,coly);
  fprintf(fd,"%-14s %-14s %-14s %-14s %-14s\n",
	  "#1:X","2:Y","3:dydx","4:d2ydx2","5:integ");

  real dx,dydx,d2ydx2,integ=0;
  for(i=0;i<ndata;i++){
    if(i==0){
      dx=X[i+1]-X[i];
      dydx=(Y[i+1]-Y[i])/dx;
      d2ydx2=(Y[2]-2*Y[1]+Y[0])/((X[2]-X[0])/2*(X[2]-X[0])/2);
      integ+=(Y[i]+Y[i+1])/2*dx;
    }else if(i==ndata-1){
      dx=X[i]-X[i-1];
      dydx=(Y[i]-Y[i-1])/dx;
      d2ydx2=(Y[ndata-1]-2*Y[ndata-2]+Y[ndata-3])/
	((X[ndata-1]-X[ndata-3])/2*(X[ndata-1]-X[ndata-3])/2);
      integ+=(Y[i]+Y[i-1])/2*dx;
    }else{
      dx=X[i+1]-X[i-1];
      dydx=(Y[i+1]-Y[i-1])/dx;
      d2ydx2=(Y[i+1]-2*Y[i]+Y[i-1])/(dx/2*dx/2);
      integ+=(Y[i+1]+Y[i-1])/2*(dx/2);
    }
    fprintf(fd,"%+14.7e %+14.7e %+14.7e %+14.7e %+14.7e\n",
	    X[i],Y[i],dydx,d2ydx2,integ);
  }
  fclose(fd);
  STPRINTF("Done.\n");
  
  return 0;
}
