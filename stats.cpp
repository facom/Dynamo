/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Compute statistics on a column of a given datafile

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
  char datafile[FSIZE]="",statsfile[FSIZE]="",typebin='l';
  int numcols=0,col=0,colw=0;
  int numbin=10,npb=0;
  int numlines;
  int i;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"STATISTICS FOR DATA FILE");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:s:n:c:b:t:N:");
  SET_USAGE(
"=======================================================================================\n"
"Usage:\n\n"
"\t./program -f <datafile> [-s <statsfile>] -n <numcols> -c <col> [-w <col_weight>]\n"
"\t                         -b <numbins> [-t <type_binning> -N <min_number_per_bin>]\n"
"\n"
"Compute basic statistics on column <col> in <datafile> and generate an histogram.\n"
"Optionally a column <col_weight> with weights used to make weight statistics, could\n"
"be probided (example mass of particles)\n"
"\n"
"Basic statistics are stored in the header of the <statsfile>.\n"
"The type of bining <type_binning> is selected between (l)inear, lo(g)arithmic,\n"
"(a)daptative (with <min_number_per_bin> threshold\n"
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
      strcpy(statsfile,optarg);
      break;
    case 'n':
      numcols=atoi(optarg);
      break;
    case 'c':
      col=atoi(optarg);
      break;
    case 'w':
      colw=atoi(optarg);
      break;
    //-------------------
    //OPTIONAL 
    //-------------------
    case 'b':
      numbin=atoi(optarg);
      break;
    case 't':
      typebin=optarg[0];
      break;
    case 'N':
      npb=atoi(optarg);
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
  if((numlines=countLines(datafile))==0){
    fprintf(stderr,"Error: Datafile '%s' seems empty\n",datafile);
    PRINT_USAGE;
    EXIT;
  }
  if(numcols<1){
    fprintf(stderr,"Error: The number of columns should be different from 0\n");
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(statsfile)){
    sprintf(statsfile,"%s.his",datafile);
  }
  if(col==0){
    col=1;
  }else if(col>numcols){
    fprintf(stderr,"Error: Column %d out of range (%d)\n",col,numcols);
    PRINT_USAGE;
    EXIT;
  }
  if(colw>numcols){
    fprintf(stderr,"Error: Column %d out of range (%d)\n",col,numcols);
    PRINT_USAGE;
    EXIT;
  }

  if(npb==0){
    npb=numlines/numbin;
  }

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    BAR(stdout,'O');
    fprintf(stdout,"Datafile: %s\n",datafile);
    fprintf(stdout,"Statistics file: %s\n",statsfile);
    fprintf(stdout,"Number of columns: %d\n",numcols);
    fprintf(stdout,"Selected column: %d\n",col);
    if(colw)
      fprintf(stdout,"Weight column: %d\n",colw);
    else
      fprintf(stdout,"Constant weight");
    fprintf(stdout,"Number of bins: %d\n",numbin);
    fprintf(stdout,"Type of binning: %c\n",typebin);
    fprintf(stdout,"Minimum per bin: %d\n",npb);
    BAR(stdout,'O');
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING DATA
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Reading data from datafile '%s'...\n",datafile);
  real2 *values=readColumn(datafile,numlines,numcols,col);
  if(colw){
    real2 *weights=readColumn(datafile,numlines,numcols,colw);
    
  }


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //STATISTICS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int ic;
  real2 min=MAXREAL,max=MINREAL,mean=0,median=0,rms=0,sd=0;
  real2 q25,q50,q75,mode=MINREAL,nmode=MINREAL;

  fprintf(stdout,"Computing basic statistics...\n",numlines);
  //Sort
  gsl_sort(values,1,numlines);

  //Basic statistics
  mean=gsl_stats_mean(values,1,numlines);
  sd=gsl_stats_sd(values,1,numlines);
  rms=sqrt((numlines-1)/numlines*sd*sd+mean*mean);
  gsl_stats_minmax(&min,&max,values,1,numlines);
  q25=gsl_stats_quantile_from_sorted_data(values,1,numlines,0.25);
  median=q50=gsl_stats_quantile_from_sorted_data(values,1,numlines,0.50);
  q75=gsl_stats_quantile_from_sorted_data(values,1,numlines,0.75);

  //Binning
  real2* bins=(real2*)calloc(numbin+1,sizeof(real2));
  real2* nhis=(real2*)calloc(numbin,sizeof(real2));
  real2 xini,x,dx,xmed,xend;

  switch(typebin){
  case 'l':
    xini=min;
    dx=(max-min)/numbin;
    break;
  case 'g':
    if(gsl_isinf(log10(min))!=0 || gsl_isinf(log10(min))!=0){
      fprintf(stderr,"Error: A logarithmic binning cannot be used for this data (min:%g,max:%g).\n\tTry with linear or adaptative binning\n",min,max);
      PRINT_USAGE;
      EXIT;
    }
    xini=log10(min);
    dx=(log10(max)-log10(min))/numbin;
    break;
  default:
    bins[0]=xini=min;
    dx=0;
    break;
  }

  if(typebin!='a'){
    for(ic=0;ic<=numbin;ic++){
      bins[ic]=xini+ic*dx;
    }
    xend=bins[1];
  }else{
    bins[0]=min;
  }

  //Compute the histogram
  fprintf(stdout,"Computing histogram...\n",numlines);
  ic=0;
  for(i=0;i<numlines;i++){
    //Value of x
    if(typebin=='g') x=log10(values[i]);
    else x=values[i];
    //Test change of bin
    if(typebin=='a'){
      if(nhis[ic]>npb){
	ic++;
	bins[ic]=x;
      }
      if(ic>numbin){
	fprintf(stderr,"Error: The number of bins is not enough for the amount of data\n");
	PRINT_USAGE;
	EXIT;
      }
    }else if(gsl_fcmp(x,xend,EPSREAL)>0){
      ic++;
      xend=bins[ic+1];
    }
    //Increment
    nhis[ic]++;
  }
  if(typebin=='a'){
    numbin=ic+1;
    bins[numbin]=max;
    fprintf(stdout,"Number of detected bins: %d\n",numbin);
  }

  //Compute derivative quantities
  for(ic=0;ic<numbin;ic++)
    if(nhis[ic]>nmode){
      mode=(bins[ic]+bins[ic+1])/2;
      nmode=nhis[ic];
    }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //STORING HISTOGRAM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file fs=fileOpen(statsfile,"w");

  //Header
  fprintf(fs,"#Data file: %s, Column: %d\n",datafile,col);
  fprintf(fs,"#Typebin: %c\n",typebin);
  fprintf(fs,"#Numbins: %d\n",numbin);
  fprintf(fs,"#NPB: %d\n",npb);
  fprintf(fs,"#Tot.stats(Ntot,Min,Max): %d,%+14.7e %+14.7e\n",numlines,min,max);
  fprintf(fs,"#Stats(mean,median,mode,rms,disp): %+14.7e %+14.7e %+14.7e %+14.7e %+14.7e\n",
	  mean,median,mode,rms,sd);
  fprintf(fs,"#Quartiles(25,50,75):%+14.7e,%+14.7e,%+14.7e\n",q25,q50,q75);
  fprintf(fs,"%-14s %-14s %-14s %-7s %-14s %-14s ",
	  "#1:xini","2:xmed","3:xend","4:n","5:h","6:f");
  fprintf(fs,"%-14s %-14s %-14s ","7:dn","8:dh","9:df");
  fprintf(fs,"%-7s %-14s ","10:F","11:P");
  fprintf(fs,"\n");

  //Histogram
  int nh,F;
  real2 dn,hhis,dh,fhis,df,P;
  F=0;
  P=0;
  for(ic=0;ic<numbin;ic++){
    if(typebin=='g'){
      xini=pow10(bins[ic]);
      xend=pow10(bins[ic+1]);
    }else{
      xini=bins[ic];
      xend=bins[ic+1];
    }
    xmed=(xini+xend)/2;
    dx=(xend-xini);

    nh=(int)nhis[ic];
    hhis=nhis[ic]/numlines;
    fhis=hhis/dx;

    dn=sqrt(nh);
    dh=dn/numlines;
    df=dh/dx;

    F+=nh;
    P+=hhis;

    //xini xend nhis hhis fhis
    fprintf(fs,"%+14.7e %+14.7e %+14.7e %-7d %+14.7e %+14.7e ",
	    xini,xmed,xend,nh,hhis,fhis);

    //Errors
    fprintf(fs,"%+14.7e %+14.7e %+14.7e ",dn,dh,df);

    //Cummulative
    fprintf(fs,"%-7d %+14.7e ",F,P);

    fprintf(fs,"\n");

  }
  
  fclose(fs);
  fprintf(stdout,"Histogram file '%s' saved...\n",statsfile);
  return 0;
}
