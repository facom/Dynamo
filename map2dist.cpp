/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Generate a scalar map of a given set of data

  Changes log:

  - Apr/7-JZ: Checked

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#include <dynamo.h>

int positionDist(int ix,int iy,int iz,int i,int j,int k,int nx,int ny,int nz)
{
  int p;
  int ireal,jreal,kreal;
  //fprintf(stdout,"Retrieved: %d,%d,%d - ",i,j,k);
  switch(ix){
  case 0:ireal=i;break;
  case 1:jreal=i;break;
  case 2:kreal=i;break;
  }
  switch(iy){
  case 0:ireal=j;break;
  case 1:jreal=j;break;
  case 2:kreal=j;break;
  }
  switch(iz){
  case 0:ireal=k;break;
  case 1:jreal=k;break;
  case 2:kreal=k;break;
  }
  //fprintf(stdout,"Real: %d,%d,%d - ",ireal,jreal,kreal);
  p=ireal+jreal*nx+kreal*nx*ny;
  //fprintf(stdout,"p = %d\n",p);
  return p;
}

int main(int argc,char *argv[])
{
  //////////////////////////////////////////////////////////////
  //PROGRAM VARIBALES
  //////////////////////////////////////////////////////////////
  char mapfile[FSIZE]="",distfile[FSIZE]="",type[2]="";
  char ctype,coper='s';
  int cx=0,cy=1;
  
  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"COMPUTE A SCALAR MAP");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVm:d:t:o:x:y:");
  SET_USAGE(
"============================================================================\n"
"Usage:\n"
"\n"
"\t./program -m <mapfile> [-d <distfile>] -t <type_dist> [-o <type_oper>]\n"
"\t                        -x <direction_x> [-y <direction_y>]\n"
/*

Collapse the information in the <mapfile>.  Collapsing a map means to
produce a 2-D or a 1-D (<type_dist>: 1D or 2D) from a full 3-D map.
What collapse do is to take a given direction in the grid and apply
over all the cells in that direction an operation (<type_oper>: sum,
average).  

You can indicate the directions that should be preserved
(<direction_x> and/or <direction_y>).  

 */
"\n"
"=========================================================================\n"
);

  //////////////////////////////////////////////////////////////
  //READ OPTIONS
  //////////////////////////////////////////////////////////////
  while(ITEROPTIONS){
    switch(OPTION){
    case 'm':
      strcpy(mapfile,optarg);
      break;
    case 'd':
      strcpy(distfile,optarg);
      break;
    case 't':
      strcpy(type,optarg);
      break;
    case 'o':
      coper=optarg[0];
      break;
    case 'x':
      cx=atoi(optarg);
      break;
    case 'y':
      cy=atoi(optarg);
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
  if(isBlank(mapfile)){
    fprintf(stderr,"Error: No mapfile was provided\n");
    PRINT_USAGE;
    EXIT;
  }
  if(!fileExists(mapfile)){
    fprintf(stderr,"Error: Mapfile '%s' does not exist\n",mapfile);
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(distfile)){
    sprintf(distfile,"%s.dst",mapfile);
  }
  if(isBlank(type)){
    strcpy(type,"1D");
  }
  ctype=type[0];
  if(ctype=='2'){
    if(cx==cy){
      fprintf(stderr,"Error: x(%d) and y(%d) cannot be equal\n",cx,cy,mapfile);
      PRINT_USAGE;
      EXIT;
    }
  }else{
    cy=(cx+1)%3;
  }
  
  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    BAR(stdout,'O');
    fprintf(stdout,"Map file: %s\n",mapfile);
    fprintf(stdout,"Distributuion file: %s\n",distfile);
    fprintf(stdout,"Type of distribution: %s\n",type);
    fprintf(stdout,"Type of operation: %c\n",coper);
    fprintf(stdout,"X column: %d\n",cx);
    if(ctype!='2')
      fprintf(stdout,"Y column: %d\n",cy);
    BAR(stdout,'O');
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////
  int i,j,k,n,p;
  char linea[LSIZE],ctmp[LSIZE];
  int nx,ny,nz;
  real *fieldmap;
  real x,y,z;
  real C[3][MAXGRID];
  int ix,iy,iz;
  int ndx,ndy,ndz;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READ MAP FILE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  STPRINTF("Getting data from mapfile '%s'...\n",mapfile);
  file fm=fopen(mapfile,"r");
  
  //Read header
  for(i=1;i<=6;i++)
    fgets(linea,sizeof linea,fm);

  //Get number of data points in the grid
  fgets(linea,sizeof linea,fm);
  sscanf(linea,"%s %d %d %d",ctmp,&nx,&ny,&nz);
  fieldmap=(real*)calloc(nx*ny*nz,sizeof(real));
  
  //Read data points
  fgets(linea,sizeof linea,fm);
  n=0;
  for(k=0;k<nz;k++){
    for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
	p=i+j*nx+k*nx*ny;
	fgets(linea,sizeof linea,fm);n++;
	//fprintf(stdout,"DATA LINEA %d: %s",n,linea);
	sscanf(linea,"%f %f %f %f",
	       &x,&y,&z,&fieldmap[p]);
	//fprintf(stdout,"%d %d %d : %d = %g\n",i,j,k,p,fieldmap[p]);
	C[1][j]=y;
      }
      C[0][i]=x;
      //Blank line
      fgets(linea,sizeof linea,fm);n++;
      //fprintf(stdout,"LINEA %d: %s",n,linea);
    }
    C[2][k]=z;
    //Blank line
    fgets(linea,sizeof linea,fm);n++;
    //fprintf(stdout,"LINEA %d: %s",n,linea);
  }
  fclose(fm);
  
  /*
  fprintf(stdout,"X:");
  arrayPrintf(stdout,"%.2e ",C[0],nx);
  fprintf(stdout,"\n");

  fprintf(stdout,"Y:");
  arrayPrintf(stdout,"%.2e ",C[1],ny);
  fprintf(stdout,"\n");

  fprintf(stdout,"Z:");
  arrayPrintf(stdout,"%.2e ",C[2],nz);
  fprintf(stdout,"\n");
  */

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COLLAPSE DATA
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  STPRINTF("Collapsing data...\n");
  real XDIST[MAXGRID],YDIST[MAXGRID];
  real distmap[MAXGRID][MAXGRID];
  
  switch(cx){
  case 0:
    ndx=nx;ix=0;
    switch(cy){
    case 1:
      ndy=ny;iy=1;
      ndz=nz;iz=2;
      break;
    case 2:
      ndy=nz;iy=2;
      ndz=ny;iz=1;
      break;
    }
    break;
  case 1:
    ndx=ny;ix=1;
    switch(cy){
    case 0:
      ndy=nx;iy=0;
      ndz=nz;iz=2;
      break;
    case 2:
      ndy=nz;iy=2;
      ndz=nx;iz=0;
      break;
    }
    break;
  case 2:
    ndx=nz;ix=2;
    switch(cy){
    case 0:
      ndy=nx;iy=0;
      ndz=ny;iz=1;
      break;
    case 1:
      ndy=ny;iy=1;
      ndz=nx;iz=0;
      break;
    }
    break;
  }
  /*
  fprintf(stdout,"ix,iy,iz: %d,%d,%d\n",ix,iy,iz);
  fprintf(stdout,"nx,ny,nz: %d,%d,%d\n",ndx,ndy,ndz);
  */

  /*
  i=1;j=3;k=2;
  p=positionDist(ix,iy,iz,i,j,k,nx,ny,nz);
  fprintf(stdout,"ENTRY (%d,%d,%d:%d): %e\n",i,j,k,p,fieldmap[p]);
  exit(0);
  //*/

  switch(ctype){
  case '1':
    for(i=0;i<ndx;i++){
      XDIST[i]=C[ix][i];
      //fprintf(stdout,"CONST %d: %g\n",i,XDIST[i]);
      distmap[i][0]=0;
      for(j=0;j<ndy;j++){
	for(k=0;k<ndz;k++){
	  p=positionDist(ix,iy,iz,i,j,k,nx,ny,nz);
	  /*
	  fprintf(stdout,"\tSumming up %d:%d,%d:%d,%d:%d = %d: %g\n",
		  ix,i,iy,j,iz,k,p,fieldmap[p]);
	  */
	  distmap[i][0]+=fieldmap[p];
	}
      }
      if(coper=='a') distmap[i][0]/=(ndy*ndz);
    }
    break;
  case '2':
    for(i=0;i<ndx;i++){
      XDIST[i]=C[ix][i];
      for(j=0;j<ndy;j++){
	YDIST[j]=C[iy][j];
	distmap[i][j]=0;
	for(k=0;k<ndz;k++){
	  p=positionDist(ix,iy,iz,i,j,k,nx,ny,nz);
	  distmap[i][j]+=fieldmap[p];
	}
	if(coper=='a') distmap[i][j]/=ndz;
      }
    }
    break;
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //STORE DISTRIBUTION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  STPRINTF("Storing distribution in file '%s'...\n",distfile);
  file fd=fileOpen(distfile,"w");
  
  fprintf(fd,"#MapFile: %s\n",mapfile);
  fprintf(fd,"#TypeDist: %s\n",type);

  switch(ctype){
  case '1':
    fprintf(fd,"#ColapseDirections: %d\n",cx);
    fprintf(fd,"%-14s %-14s\n","#1:X","2:DIST");
    for(i=0;i<ndx;i++){
      //fprintf(fd,"%+14.7e %+14.7e\n",XDIST[i],distmap[i][0]);
      fprintf(fd,"%+14.7e %g\n",XDIST[i],distmap[i][0]);
    }
    break;
  case '2':
    fprintf(fd,"#ColapseDirections: %d %d\n",cx,cy);
    fprintf(fd,"%-14s %-14s\t%-14s\n","#1:X","2:Y","3:DIST");
    for(i=0;i<ndx;i++){
      for(j=0;j<ndy;j++){
	//fprintf(fd,"%+14.7e %+14.7e\t%+14.7e\n",XDIST[i],YDIST[j],distmap[i][j]);
	fprintf(fd,"%+14.7e %+14.7e\t%g\n",XDIST[i],YDIST[j],distmap[i][j]);
      }
      fprintf(fd,"\n");
    }
    break;
  }  

  return 0;
}
