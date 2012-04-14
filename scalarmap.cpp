/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Purpose: Generate a scalar map of a given set of data

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
  char simfile[FSIZE]="",mapfile[FSIZE]="";
  char coordinates[FSIZE]="",types[FSIZE]="",type[3][LSIZE],field[FSIZE]="";
  char cini[FSIZE]="",cend[FSIZE]="",ngrid[FSIZE]="";
  int i,j,k,p,ndata;
  particles parts;
  SHeader sheader;
  char ctype[3],ccoord,cfield;

  //////////////////////////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////////////////////////
  TITLE(stdout,'*',"COMPUTE A SCALAR MAP");

  //////////////////////////////////////////////////////////////
  //SET OPTIONS AND USAGE
  //////////////////////////////////////////////////////////////
  SET_OPTIONS(":hvVf:m:s:c:d:i:e:n:");
  SET_USAGE(
"============================================================================\n"
"Usage:\n"
"\n"
"\n"
"\t./program -f <simfile> [-m <mapfile>] -s <scalar_field>\n"
"\t          [-c <coordinates>] [-d <type_of_grid>]\n"
"\t          -i <cini> -e <cend> -n <ngrid>\n"
"\n"
"Create a map of the <scalar_field> associated fo the simulation\n"
"<simfile> and stores the\n results in <mapfile>.\n"
"\n"
"<scalar_field>: number_density, mass_density, gravitational_potential\n"
"\n"
"The grid for the scalar map should be described using the following\n"
"options:\n"
"\n"
"<coordinates>: (car)tesian, (cyl)indrical, (sph)rical (default: car)\n"
"<type_of_grid>: (lin)ear, (log)arithmic (default: lin)\n"
"<cini>: xini,yini,zini\n"
"<cend>: xend,yend,zend\n"
"<ngrid>: nx,ny,nz\n"
"\n"
"Where x, y and z stand also for (rho,phi,z) or (r,phi,theta)\n"
"accrodingly to the coordinate system.  If you use '-' instead of\n"
"values for the cini and cend specification the minimum and maximum\n"
"values of the coordinates will be used.\n"
"=========================================================================\n"
);

  //////////////////////////////////////////////////////////////
  //READ OPTIONS
  //////////////////////////////////////////////////////////////
  while(ITEROPTIONS){
    switch(OPTION){
    case 'f':
      strcpy(simfile,optarg);
      break;
    case 'm':
      strcpy(mapfile,optarg);
      break;
    case 's':
      strcpy(field,optarg);
      break;
    case 'c':
      strcpy(coordinates,optarg);
      break;
    case 't':
      strcpy(types,optarg);
      break;
    case 'i':
      strcpy(cini,optarg);
      break;
    case 'e':
      strcpy(cend,optarg);
      break;
    case 'n':
      strcpy(ngrid,optarg);
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
  if(isBlank(field)){
    fprintf(stderr,"Error: You must specify a scalar field\n");
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(mapfile)){
    sprintf(mapfile,"%s.map",simfile);
  }
  if((ndata=countLines(simfile))==0){
    fprintf(stderr,"Error: Simfile '%s' seems empty\n",simfile);
    PRINT_USAGE;
    EXIT;
  }
  if(isBlank(coordinates)){
    strcpy(coordinates,"car");
  }
  if(isBlank(types)){
    strcpy(types,"lin,lin,lin");
  }
  if(isBlank(cini)){
    strcpy(cini,"-,-,-");
  }
  if(isBlank(cend)){
    strcpy(cend,"-,-,-");
  }
  if(isBlank(ngrid)){
    fprintf(stderr,"Error: you must specify number of points for the grid\n");
    PRINT_USAGE;
    EXIT;
  }

  //Assign symbols to coordinates, type of grid and scalar field
  if(strcmp(coordinates,"car")==0) ccoord='c';
  else if(strcmp(coordinates,"cyl")==0) ccoord='y';
  else if(strcmp(coordinates,"sph")==0) ccoord='s';
  else{
    fprintf(stderr,"Error: Coordinate system '%s' not recognized\n",coordinates);
    PRINT_USAGE;
    EXIT;
  }

  splitString(types,",",type);
  for(i=0;i<3;i++){
    if(strcmp(type[i],"lin")==0) ctype[i]='l';
    else if(strcmp(type[i],"log")==0) ctype[i]='g';
    else{
      fprintf(stderr,"Error: Type of grid '%s' not recognized\n",type[i]);
      PRINT_USAGE;
      EXIT;
    }
  }

  if(strcmp(field,"number_density")==0) cfield='n';
  else if(strcmp(field,"mass_density")==0) cfield='m';
  else if(strcmp(field,"gravitational_potential")==0) cfield='g';
  else{
    fprintf(stderr,"Error: Field '%s' not recognized\n",field);
    PRINT_USAGE;
    EXIT;
  }

  //CHECK DEPTH OF THE GRID
  int nx,ny,nz;
  real list_n[3];
  splitString(ngrid,",",list_n);
  nx=DEFVAL(list_n[0],0,10);
  ny=DEFVAL(list_n[1],0,10);
  nz=DEFVAL(list_n[2],0,10);
  /*
  if(nx>MAXGRID){
    fprintf(stdout,"Error the depth of the grid in the x direction '%d' is larger than maximum %d\n",nx,MAXGRID);
    PRINT_USAGE;
    EXIT;
  }
  if(ny>MAXGRID){
    fprintf(stdout,"Error the depth of the grid in the y direction '%d' is larger than maximum %d\n",ny,MAXGRID);
    PRINT_USAGE;
    EXIT;
  }
  if(nz>MAXGRID){
    fprintf(stdout,"Error the depth of the grid in the z direction '%d' is larger than maximum %d\n",nz,MAXGRID);
    PRINT_USAGE;
    EXIT;
  }
  */

  //////////////////////////////////////////////////////////////
  //REPORT INPUT INFORMATION
  //////////////////////////////////////////////////////////////
  if(VERBOSE(1)){
    BAR(stdout,'O');
    fprintf(stdout,"Simulation file: %s\n",simfile);
    fprintf(stdout,"Map file: %s\n",mapfile);
    fprintf(stdout,"Scalar field: %s\n",field);
    fprintf(stdout,"Coordinates: %s\n",coordinates);
    fprintf(stdout,"Type of grid: %s,%s,%s\n",type[0],type[1],type[2]);
    fprintf(stdout,"Initial points: %s\n",cini);
    fprintf(stdout,"Final points: %s\n",cend);
    fprintf(stdout,"Number of points: %s\n",ngrid);
    BAR(stdout,'O');
  }

  //////////////////////////////////////////////////////////////
  //PROGRAM
  //////////////////////////////////////////////////////////////

  //========================================
  //READ DATA AND GET SOME BASIC PROPERTIES
  //========================================
  STPRINTF("Reading simulation data...\n");
  parts=readSimulation(simfile,&sheader);
  //checkSimulation(sheader);
  
  real *X,*Y,*Z,*R,*r;
  real xmin,xmax,ymin,ymax,Rmin,Rmax,rmin,rmax;

  X=sliceSimulation(parts,0,ndata-1,PART_MEMBER(pos.x));
  Y=sliceSimulation(parts,0,ndata-1,PART_MEMBER(pos.y));
  Z=sliceSimulation(parts,0,ndata-1,PART_MEMBER(pos.z));
  R=sliceSimulation(parts,0,ndata-1,PART_MEMBER(posc.x));
  r=sliceSimulation(parts,0,ndata-1,PART_MEMBER(r));

  //========================================
  //COMPUTE PROPERTIES OF THE GRID
  //========================================
  STPRINTF("Preparing the grid (%d,%d,%d)...\n",nx,ny,nz);
  real list_ini[3],list_end[3];
  real xini,xend,yini,yend,zini,zend;
  real dx,dy,dz;
  size_t xoff,yoff,zoff;
  real x,y,z;
  splitString(cini,",",list_ini);
  splitString(cend,",",list_end);

  switch(ccoord){
  case 'c':
    xini=DEFVAL(list_ini[0],0,arrayMin(X,ndata));
    yini=DEFVAL(list_ini[1],0,arrayMin(Y,ndata));
    zini=DEFVAL(list_ini[2],0,arrayMin(Z,ndata));

    xend=DEFVAL(list_end[0],0,arrayMax(X,ndata));
    yend=DEFVAL(list_end[1],0,arrayMax(Y,ndata));
    zend=DEFVAL(list_end[2],0,arrayMax(Z,ndata));

    xoff=PART_MEMBER(pos.x);
    yoff=PART_MEMBER(pos.y);
    zoff=PART_MEMBER(pos.z);
    break;
  case 'y':
    xini=DEFVAL(list_ini[0],0,arrayMin(R,ndata));
    yini=DEFVAL(list_ini[1],0,0);
    zini=DEFVAL(list_ini[2],0,arrayMin(Z,ndata));

    xend=DEFVAL(list_end[0],0,arrayMax(R,ndata));
    yend=DEFVAL(list_end[1],0,2*M_PI);
    zend=DEFVAL(list_end[2],0,arrayMax(Z,ndata));

    xoff=PART_MEMBER(posc.x);
    yoff=PART_MEMBER(posc.y);
    zoff=PART_MEMBER(posc.z);
    break;
  case 's':
    xini=DEFVAL(list_ini[0],0,arrayMin(r,ndata));
    yini=DEFVAL(list_ini[1],0,0);
    zini=DEFVAL(list_ini[2],0,0);

    xend=DEFVAL(list_end[0],0,arrayMax(r,ndata));
    yend=DEFVAL(list_end[1],0,2*M_PI);
    zend=DEFVAL(list_end[2],0,M_PI);

    xoff=PART_MEMBER(r);
    yoff=PART_MEMBER(posc.y);
    zoff=PART_MEMBER(posc.x);//THETA IS NOT COMPUTED SO R IS RETURNED
    break;
  }
  dx=(xend-xini)/nx;
  dy=(yend-yini)/ny;
  dz=(zend-zini)/nz;

  if(EQUAL(dx,0)) dx=1;
  if(EQUAL(dy,0)) dy=1;
  if(EQUAL(dz,0)) dz=1;


  VPRINTF(2)("xini,yini,zini = %e,%e,%e\n",xini,yini,zini);
  VPRINTF(2)("xend,yend,zend = %e,%e,%e\n",xend,yend,zend);
  VPRINTF(2)("dx,dy,dz = %e,%e,%e\n",dx,dy,dz);

  //========================================
  //COMPUTE THE SCALAR MAP
  //========================================
  STPRINTF("Computing the scalar map...\n");
  int ibox,jbox,kbox;
  real *fieldmap;
  fieldmap=(real*)calloc(nx*ny*nz,sizeof(real));

  switch(cfield){
  case 'n':case 'm':
    //==================================================
    //NUMBER AND MASS DENSITY
    //==================================================
    for(i=0;i<ndata;i++){
      x=GET_MEMBER(parts[i],real,xoff);
      y=GET_MEMBER(parts[i],real,yoff);
      z=GET_MEMBER(parts[i],real,zoff);
      if(x<xini || x>xend) continue;
      if(y<yini || y>yend) continue;
      if(z<zini || z>zend) continue;

      if(ccoord=='s') z=asin(z/x);//COMPUTE THETA
      ibox=(int)(x-xini)/dx;
      jbox=(int)(y-yini)/dy;
      kbox=(int)(z-zini)/dz;
      p=ibox+jbox*nx+kbox*nx*ny;
      /*
      if(cfield=='n') fieldmap[ibox][jbox][kbox]+=1;
      else fieldmap[ibox][jbox][kbox]+=parts[i].m;
      */
      if(cfield=='n') fieldmap[p]+=1;
      else fieldmap[p]+=parts[i].m;
    }
  case 'g':
    for(k=0;k<nz;k++){
      z=zini+k*dz+dz/2;
      for(i=0;i<nx;i++){
	x=xini+i*dx+dx/2;
	for(j=0;j<ny;j++){
	  y=yini+j*dy+dy/2;
	  
	  p=i+j*nx+k*nx*ny;
	  fieldmap[p]=phi;
	}
      }
    }
    for(i=0;i<ndata;i++){
      x=GET_MEMBER(parts[i],real,xoff);
      y=GET_MEMBER(parts[i],real,yoff);
      z=GET_MEMBER(parts[i],real,zoff);
      if(x<xini || x>xend) continue;
      if(y<yini || y>yend) continue;
      if(z<zini || z>zend) continue;

      if(ccoord=='s') z=asin(z/x);//COMPUTE THETA
      ibox=(int)(x-xini)/dx;
      jbox=(int)(y-yini)/dy;
      kbox=(int)(z-zini)/dz;
      p=ibox+jbox*nx+kbox*nx*ny;
      /*
      if(cfield=='n') fieldmap[ibox][jbox][kbox]+=1;
      else fieldmap[ibox][jbox][kbox]+=parts[i].m;
      */
      if(cfield=='n') fieldmap[p]+=1;
      else fieldmap[p]+=parts[i].m;
    }
    break;
  }

  //========================================
  //STORE THE SCALAR MAP
  //========================================
  STPRINTF("Storing the scalar map in '%s'...\n",mapfile);
  file fs=fileOpen(mapfile,"w");
  
  //Header of the map file
  fprintf(fs,"#Simfile: %s\n",simfile);
  fprintf(fs,"#Field: %s\n",field);
  fprintf(fs,"#Coordinates: %s\n",coordinates);
  fprintf(fs,"#TypeGrid: %s,%s,%s\n",type[0],type[1],type[2]);
  fprintf(fs,"#Initial: %+14.7e %+14.7e %+14.7e\n",
	  xini,yini,zini);
  fprintf(fs,"#Final: %+14.7e %+14.7e %+14.7e\n",
	  xend,yend,zend);
  fprintf(fs,"#NumberPoints: %d %d %d\n",nx,ny,nz);
  fprintf(fs,"%-14s %-14s %-14s\t%-14s\n","#1:X","2:Y","3:Z","4:field");

  for(k=0;k<nz;k++){
    z=zini+k*dz;
    for(i=0;i<nx;i++){
      x=xini+i*dx;
      for(j=0;j<ny;j++){
	y=yini+j*dy;
	/*
	fprintf(fs,"%+14.7e %14.7e %+14.7e\t%+14.7e\n",
		x,y,z,fieldmap[i][j][k]);
	*/
	p=i+j*nx+k*nx*ny;
	/*
	fprintf(fs,"%+14.7e %14.7e %+14.7e\t%g\n",
		x,y,z,fieldmap[i][j][k]);
	*/
	fprintf(fs,"%+14.7e %14.7e %+14.7e\t%g\n",
		x,y,z,fieldmap[p]);
      }
      fprintf(fs,"\n");
    }
    fprintf(fs,"\n");
  }
      
  fclose(fs);

  return 0;
}
