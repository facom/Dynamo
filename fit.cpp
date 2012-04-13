/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Fit a given set of data with a function defined in the
  functions.cpp file.

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
  char datafile[FSIZE]="",fitfile[FSIZE]="";
  char fitfunc[FSIZE]="",inipars[FSIZE]="";
  real2 params[MAXPARAMS];
  int numcols=0,colx=0,coly=0,cole=0;
  int nfac=2;
  int i,ndata;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"FIT DATA TO A GIVEN FUNCTION");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:F:u:p:n:x:y:e:N:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <datafile> [-F <fitfile>] -u <fit_function> [-p <initial_params>]\n"
"\t          [-N <num_sampling_fitting_func>]\n"
"\t          -n <numcols> -x <colx> -y <coly> [-e <cole>]\n"
"\n"
"Fit the 2D data described by <colx> and <coly> with errors <cole> using fitting function\n"
"<fit_function> as given by function file functions.hpp.  The initial set of parameters\n"
"<initial_params> should be provided as a list of comma separated real values.\n"
"The result of the fit is stored in file <fitfile> where the fitted parameters, the \n"
"chisquare and the p-value are stored in the header.  <fitfile> is a two column file with\n"
"the value of the fitted function at intermediate points (<num_sampling_fitting_func> x \n"
"number of fitted points).\n"
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
    case 'F':
      strcpy(fitfile,optarg);
      break;
    case 'u':
      strcpy(fitfunc,optarg);
      break;
    case 'p':
      strcpy(inipars,optarg);
      break;
    case 'n':
      numcols=atoi(optarg);
      break;
    case 'N':
      nfac=atoi(optarg);
      break;
    case 'x':
      colx=atoi(optarg);
      break;
    case 'y':
      coly=atoi(optarg);
      break;
    case 'e':
      cole=atoi(optarg);
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
  if(isBlank(fitfile)){
    sprintf(fitfile,"%s.fit",datafile);
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
    colx=1;
  }
  if(coly==0){
    coly=2;
  }
  if(isBlank(fitfunc)){
    fprintf(stderr,"Error: No fit function was provided\n");
    PRINT_USAGE;
    EXIT;
  }else{
    getFunction(fitfunc);
  }
  if(isBlank(inipars)){
    for(i=0;i<FNpars;i++){
      if(i==FNpars-1)
	strcat(inipars,"0.0");
      else
	strcat(inipars,"0.0,");
    }
  }
  splitString(inipars,",",params);

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    fprintf(stdout,"Datafile: %s\n",datafile);
    fprintf(stdout,"Fitfile: %s\n",fitfile);
    fprintf(stdout,"Fit function: %s\n",fitfunc);
    fprintf(stdout,"Fit parameters: %s\n",inipars);
    fprintf(stdout,"Test call: %s(%+14.7e;%s) = %+14.7e\n",
	    fitfunc,FXtest,inipars,FFunc(FXtest,params));
    fprintf(stdout,"Number of columns: %d\n",numcols);
    fprintf(stdout,"Columns x,y: %d,%d\n",colx,coly);
    if(cole!=0)
      fprintf(stdout,"Column error: %d\n",cole);
    else
      fprintf(stdout,"Errors assumed 1\n");
    fprintf(stdout,"Factor of over sampling: %d\n",nfac);
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING DATA
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Reading %d data points from datafile '%s'...\n",ndata,datafile);
  real2 *X=readColumn(datafile,ndata,numcols,colx);
  real2 *Y=readColumn(datafile,ndata,numcols,coly);
  real2 *E;
  if(cole)
    E=readColumn(datafile,ndata,numcols,cole);
  else{
    E=(real2*)calloc(ndata,sizeof(real2));
    for(int i=0;i<ndata;i++) E[i]=1.0;
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //FITTING DATA
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real2 chis,pval;
  fprintf(stdout,"Fitting %d data points with function '%s'...\n",
	  ndata,fitfunc);
  fitData(X,Y,E,ndata,FFunc,FNpars,params,chis,pval);

  fprintf(stdout,"Fit succesful:\n");
  fprintf(stdout,"\tBest fit parameters: ");
  for(i=0;i<FNpars;i++)
    fprintf(stdout,"%+14.7e ",params[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"\tChisquare: %+14.7e\n",chis);
  fprintf(stdout,"\tP-val (nu = %d): %+14.7e\n",ndata,pval);
  if(pval>0.05)
    fprintf(stdout,"\tData is compatible with fitting function\n");
  else
    fprintf(stdout,"\tData cannot be fitted with function\n");

  file fs=fileOpen(fitfile,"w");
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //STORING 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Storing fit results...\n",fitfile);
  
  fprintf(fs,"#Datafile: %s\n",datafile);
  fprintf(fs,"#TotCols,ColX,ColY,ColE: %d %d %d %d\n",
	  numcols,colx,coly,cole);
  fprintf(fs,"#FitFunction: %s\n",fitfunc);
  fprintf(fs,"#Initial parameters: %s\n",inipars);
  fprintf(fs,"#NumberDataPoints: %d\n",ndata);
  fprintf(fs,"#FitParameters: ");
  for(i=0;i<FNpars;i++)
    fprintf(fs,"%+14.7e ",params[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#Chisquare: %+14.7e\n",chis);
  fprintf(fs,"#P-val: %+14.7e\n",pval);
  fprintf(fs,"%-14s %-14s\n","#1:X","2:Y");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //SAMPLING FUNCTION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Sampling fitting function...\n",fitfile);
  
  real2 x,y;
  real2 xini=X[0];
  real2 dx=(X[ndata-1]-X[0])/(ndata*nfac);
  for(i=0;i<=ndata*nfac;i++){
    x=xini+i*dx;
    y=FFunc(x,params);
    fprintf(fs,"%+14.7e %+14.7e\n",x,y);
  }

  fprintf(stdout,"Fitting result stored in %s...\n",fitfile);
  fclose(fs);
  return 0;
}
