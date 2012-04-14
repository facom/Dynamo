#include <dynamo.h>

//////////////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////////////////////////////////////////////

/*D*/

char OPTIONS[100];
char OPTION;
char USAGE[1000];
int VERBOSITY;

real2 *XFit,*YFit,*EFit;
int NFit;
real2 (*FFit)(real2,pars);

/*D*/

//////////////////////////////////////////////////////////////////////////////////
//ROUTINES
//////////////////////////////////////////////////////////////////////////////////

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//I/O
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
/*P*/
particles readGadget(char simulation[],char snapshot[],SHeader* sheader,
		     bool qsub/*=false=*/)
{
  int i,j;
  file fd,fp;
  int dummy;
  GHeader header;
  char fsnap[FSIZE],fparam[FSIZE];

  sprintf(fsnap,"%s_%s",simulation,snapshot);

  ////////////////////////////////////////////////////////////////////////
  //READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  fd=fileOpen(fsnap,"r");

  ////////////////////////////////////////////////////////////////////////
  //READ HEADER
  ////////////////////////////////////////////////////////////////////////
  DUMMY(fd);fread(&header,sizeof(GHeader),1,fd);DUMMY(fd);
  sheader->z=header.redshift;
  sheader->omegaL=header.OmegaLambda;
  sheader->Ho=header.HubbleParam;
  sheader->boxsize=header.BoxSize;

  ////////////////////////////////////////////////////////////////////////
  //READ DATA FOR PARTICLES
  ////////////////////////////////////////////////////////////////////////
  int ntot=0;
  for(i=0;i<6;i++){
    ntot+=header.npart_total[i];
    sheader->npartg[i]=header.npart_total[i];
  }
  sheader->ntot=ntot;
  particles parts;
  parts=(particles)malloc(ntot*sizeof(particle));

  //==================================================
  //POSITION AND VELOCITIES
  //==================================================
  DUMMY(fd);
  for(i=0;i<ntot;i++)
    fread(&parts[i].pos,sizeof(real),3,fd);
  DUMMY(fd);DUMMY(fd);
  for(i=0;i<ntot;i++)
    fread(&parts[i].vel,sizeof(real),3,fd);
  DUMMY(fd);DUMMY(fd);

  //==================================================
  //IDS
  //==================================================
  for(i=0;i<ntot;i++)
    fread(&parts[i].pid,sizeof(int),1,fd);
  DUMMY(fd);

  //==================================================
  //MASSES
  //==================================================
  real mass,mtot,m;
  i=0;
  mtot=0;
  for(int gtid=0;gtid<6;gtid++){
    mass=0;
    for(j=0;j<header.npart_total[gtid];j++){
      if(header.mass[gtid]==0)
	fread(&parts[i].m,sizeof(float),1,fd);
      else
	parts[i].m=header.mass[gtid];
      parts[i].gtype=gtid;
      parts[i].tid=gtid;
      parts[i].sid=i;
      m=parts[i].m;
      mass+=m;
      i++;
    }
    sheader->mtotg[gtid]=mass;
    printf("Type %d, mass = %d x %e = %e\n",
	   gtid,header.npart_total[gtid],m,mass);
    mtot+=mass;
  }
  printf("M.tot. = %e\n",
	 mtot);
  sheader->mtot=mtot;
  DUMMY(fd);

  //BY DEFAULT SUBSTRUCTURES ARE THE SAME AS GADGET TYPES
  sheader->nsubs=6;
  for(i=0;i<6;i++){
    sheader->stype[i]=i;
    sprintf(sheader->subtype[i],"Type%d",i);
    sheader->nparts[i]=sheader->npartg[i];
    sheader->mtots[i]=sheader->mtotg[i];
  }

  //==================================================
  //TYPES
  //==================================================
  if(qsub){
    //****************************************
    //READ VECTOR OF SUBTYPES
    //****************************************
    int nsub[NSUB],stype[NSUB],nsubs;
    char line[LSIZE],subtype[NSUB][LSIZE];
    sprintf(fsnap,"%s.sub",simulation);
    file fs=fileOpen(fsnap,"r");
    i=0;
    while(fgets(line,sizeof line,fs)!=NULL){
      sscanf(line,"%d %d %s",&nsub[i],&stype[i],subtype[i]);
      sheader->nparts[i]=0;
      sheader->mtots[i]=0;
      sheader->stype[i]=stype[i];
      strcpy(sheader->subtype[i],subtype[i]);
      i++;
    }
    nsubs=i;
    sheader->nsubs=nsubs;
    fclose(fs);
    
    //****************************************
    //ASSIGN SUBTYPES TO PARTICLES
    //****************************************
    for(i=0;i<ntot;i++){
      int nacum=0;
      for(j=0;j<nsubs;j++){
	if(parts[i].sid>=nacum &&
	   parts[i].sid<nacum+nsub[j]){
	  parts[i].tid=stype[j];
	  sheader->nparts[j]++;
	  sheader->mtots[j]+=parts[i].m;
	  break;
	}
	nacum+=nsub[j];
      }
      //fprintf(stdout,"Particle %d of substructure %d\n",i,parts[i].tid);
    }
    //i=4;fprintf(stdout,"Total particles in substructure %d: %d\n",stype[i],sheader->nparts[i]);
  }

  //==================================================
  //SOFTENING LENGTHS
  //==================================================
  real sg,sh,sd,sb,ss,sy;
  sprintf(fparam,"%s.param",simulation);
  readGadgetParam(fparam,"SofteningGas",'f',&sh);
  readGadgetParam(fparam,"SofteningHalo",'f',&sh);
  readGadgetParam(fparam,"SofteningDisk",'f',&sd);
  readGadgetParam(fparam,"SofteningBulge",'f',&sb);
  readGadgetParam(fparam,"SofteningStars",'f',&ss);
  readGadgetParam(fparam,"SofteningBndry",'f',&sy);
  for(i=0;i<ntot;i++){
    switch(parts[i].gtype){
    case 0:parts[i].h=sg;break;
    case 1:parts[i].h=sh;break;
    case 2:parts[i].h=sd;break;
    case 3:parts[i].h=sb;break;
    case 4:parts[i].h=ss;break;
    case 5:parts[i].h=sy;break;
    }
  }

  //HEADER REQUIRES AN UPDATE: JTOT, RCM, VCM
  sheader->qupdate=true;
  sheader->rf=0;
  constVector(&sheader->jtot,0.0);
  constVector(&sheader->rcm,0.0);
  constVector(&sheader->vcm,0.0);

  return parts;
}

/*
  This routine update derivative properties of particles: position in
  cylindrical coordinates, angular momentum, radial distance, circular
  velocity, circular angular momentum, epsilon parameter
 */
/*P*/
int updateParticles(particles parts,SHeader sheader) 
{
  int i,j;
  size_t* sids;
  double* rs;
  sids=(size_t*)malloc(sheader.ntot*sizeof(size_t));
  rs=(double*)malloc(sheader.ntot*sizeof(double));

  for(i=0;i<sheader.ntot;i++){
    //ANGULAR MOMENTUM
    crossVector(&parts[i].j,parts[i].pos,parts[i].vel);
    scaleVector(&parts[i].j,parts[i].m,parts[i].j);

    //COMPUTE THE RADIAL POSITION
    parts[i].r=magVector(parts[i].pos);

    //CONVERT TO CYLINDRICAL COORDINATES
    cylVector(&parts[i].posc,parts[i].pos);
    cylVector(&parts[i].velc,parts[i].vel);
    /*
    printf("Part.%d: Cartesian %s -> Cylindrical %s\n",
	   i,vector2str(parts[i].pos),vector2str(parts[i].posc));
    //*/
    /*
    printf("Part.%d: Cartesian %s -> Cylindrical %s\n",
	   i,vector2str(parts[i].vel),vector2str(parts[i].velc));
    //*/

    //STORE VECTORS OF RADIAL DISTANCES
    rs[i]=parts[i].r;
  }  

  //SORT BY r
  gsl_sort_index(sids,rs,1,sheader.ntot);

  //COMPUTE THE CIRCULAR VELOCITY, MOMENTUM AND EPSILON PARAMETER

  //First particle
  parts[sids[0]].vcirc=0;
  parts[sids[0]].jcirc=0;
  parts[sids[0]].epsilon=0;

  size_t sid;
  real mr;

  //Content mass
  mr=parts[sids[0]].m;
  for(i=1;i<sheader.ntot;i++){
    sid=sids[i];
    parts[sid].vcirc=sqrt(GGAD*mr/parts[sid].r);
    parts[sid].jcirc=parts[sid].m*parts[sid].vcirc*parts[sid].r;
    parts[sid].epsilon=parts[sid].j.z/parts[sid].jcirc;
    /*
    printf("Particle %lu: mr=%e, vcirc=%e, jcirc=%e, jz=%e, epsilon=%e\n",
	   sid,mr,parts[sid].vcirc,parts[sid].jcirc,parts[sid].j.z,parts[sid].epsilon);
    //*/
    mr+=parts[sid].m;
  }

  return 0;
}

/*
  This routine recounts the number of particles and compute the total
  angular momentum, the position and velocity of the center of mass
  and update other important information about the simulation
 */
/*P*/
int updateHeader(particles parts,SHeader *sheader) 
{
  int st,stype,gtype;

  for(st=0;st<sheader->nsubs;st++){
    stype=sheader->stype[st];
    sheader->nparts[st]=0;
    sheader->mtots[st]=0;
  }
  for(gtype=0;gtype<6;gtype++){
    sheader->npartg[gtype]=0;
    sheader->mtotg[gtype]=0;
  }
  constVector(&sheader->jtot,0.0);
  constVector(&sheader->rcm,0.0);
  constVector(&sheader->vcm,0.0);

  /*
  printf("Total number of particles %d\n",sheader->ntot);
  printf("Total number of subtypes %d\n",sheader->nsubs);
  */
  //TOTAL MASS
  sheader->mtot=0;for(int i=0;i<sheader->ntot;i++) sheader->mtot+=parts[i].m;
  for(int i=0;i<sheader->ntot;i++){
    //============================================================
    //DETERMINE THE TYPE OF PARTICLE AND ACCUMULATE NUM. AND MASS
    //============================================================
    for(st=0;st<sheader->nsubs;st++){
      stype=sheader->stype[st];
      if(parts[i].tid==stype){
	sheader->nparts[st]++;
	sheader->mtots[st]+=parts[i].m;
	break;
      }
    }
    for(gtype=0;gtype<6;gtype++){
      if(parts[i].gtype==gtype){
	sheader->npartg[gtype]++;
	sheader->mtotg[gtype]+=parts[i].m;
	break;
      }
    }
    
    //============================================================
    //COMPUTE THE TOTAL ANGULAR MOMENTUM,CENTER OF MASS PROPERTIES
    //============================================================
    sumVector(&sheader->jtot,1,parts[i].j,sheader->jtot);
    sumVector(&sheader->rcm,parts[i].m/sheader->mtot,parts[i].pos,sheader->rcm);
    sumVector(&sheader->vcm,parts[i].m/sheader->mtot,parts[i].vel,sheader->vcm);
  }

  /*
  printf("\nTotal mass: %e\n",sheader->mtot);
  printf("\nTotal angular momentum: %s\n",vector2str(sheader->jtot));
  printf("rcm: %s\n",vector2str(sheader->rcm));
  printf("vcm: %s\n",vector2str(sheader->vcm));
  
  for(st=0;st<sheader->nsubs;st++){
    stype=sheader->stype[st];
    printf("Simulation type %d:\n",stype);
    printf("\tTotal number\t%d\n",sheader->nparts[st]);
    printf("\tTotal mass of sim.type\t%e\n",sheader->mtots[st]);
  }
  printf("\n\n");
  for(gtype=0;gtype<6;gtype++){
    printf("Gadget type %d:\n",gtype);
    printf("\tTotal number\t%d\n",sheader->npartg[gtype]);
    printf("\tTotal mass of sim.type\t%e\n",sheader->mtotg[gtype]);
  }
  printf("\n");
  //*/

}

/*P*/
int checkSimulation(SHeader sheader)
{
  int i;
  TITLE(stdout,'=',"SIMULATION CHECK");
  fprintf(stdout,"#NTOT,MTOT\t%-7d %+14.7e\n",sheader.ntot,sheader.mtot);
  fprintf(stdout,"#NPART.G\t");
  for(i=0;i<=5;i++)
    fprintf(stdout,"%-7d\t",sheader.npartg[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#MTOT.G.\t");
  for(i=0;i<=5;i++)
    fprintf(stdout,"%+14.7e\t",sheader.mtotg[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#NSUBS\t\t%d\n",sheader.nsubs);
  fprintf(stdout,"#SUBS.ID.\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(stdout,"%-2d ",sheader.stype[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#SUBS.STR.\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(stdout,"%-s ",sheader.subtype[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#NPART.SIM.\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(stdout,"%-7d ",sheader.nparts[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#MTOT.SIM.\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(stdout,"%+14.7e\t",sheader.mtots[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#ANG.MOM.TOT\t%44s\n",vector2str(sheader.jtot));
  fprintf(stdout,"#POSIT.CM\t%44s\n",vector2str(sheader.rcm));
  fprintf(stdout,"#VELOC.CM\t%44s\n",vector2str(sheader.vcm));
  fprintf(stdout,"#Z,OL,HO,BOX\t%+14.7e %+14.7e %+14.7e %+14.7e\n",
	  sheader.z,sheader.omegaL,sheader.Ho,sheader.boxsize);
  fprintf(stdout,"#REF.FRAME\t%-2d\n",sheader.rf);
  BAR(stdout,'=');
  return 0;
}

/*P*/
int writeSimulation(particles parts,SHeader sheader,char simfile[],int frac/*=1=*/) 
{
  file fs=fileOpen(simfile,"w");
  int i,j;

  //======================================================================
  //WRITE HEADER
  //======================================================================
  fprintf(fs,"#NTOT,MTOT:\t%-7d %+14.7e\n",sheader.ntot,sheader.mtot);
  fprintf(fs,"#NPART.G:\t");
  for(i=0;i<=5;i++)
    fprintf(fs,"%-7d\t",sheader.npartg[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#MTOT.G.:\t");
  for(i=0;i<=5;i++)
    fprintf(fs,"%+14.7e\t",sheader.mtotg[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#NSUBS:\t\t%d\n",sheader.nsubs);
  fprintf(fs,"#SUBS.ID.:\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(fs,"%-2d ",sheader.stype[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#SUBS.STR.:\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(fs,"%-s ",sheader.subtype[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#NPART.SIM.:\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(fs,"%-7d ",sheader.nparts[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#MTOT.SIM.:\t");
  for(i=0;i<sheader.nsubs;i++)
    fprintf(fs,"%+14.7e ",sheader.mtots[i]);
  fprintf(fs,"\n");
  fprintf(fs,"#ANG.MOM.TOT:\t%44s\n",vector2str(sheader.jtot));
  fprintf(fs,"#POSIT.CM:\t%44s\n",vector2str(sheader.rcm));
  fprintf(fs,"#VELOC.CM:\t%44s\n",vector2str(sheader.vcm));
  fprintf(fs,"#Z,OL,HO,BOX:\t%+14.7e %+14.7e %+14.7e %+14.7e\n",
	  sheader.z,sheader.omegaL,sheader.Ho,sheader.boxsize);
  fprintf(fs,"#REF.FRAME:\t%-2d\n",sheader.rf);

  //======================================================================
  //DESCRIBE FIELDS
  //======================================================================
  fprintf(fs,	  
	  "%-7s %-7s %-7s %-7s "
	  "%-14s %-14s "
	  "%-44s %-44s %-44s %-44s "
	  "%-44s "
	  "%-14s %-14s %-14s %-14s\n",
	  "#1:sid","2:pid","3:tid","4:gtype",
	  "5:m","6:h",
	  "7:pos","10:vel","13:posc","16:velc",
	  "19:j",
	  "22:r","23:vc","24:jc","25:eps");

  //======================================================================
  //WRITE DATA
  //======================================================================
  for(i=0;i<sheader.ntot;i++){
    if((i%frac)==0)
      fprintf(fs,"%s\n",particle2str(parts[i]));
  }

  fclose(fs);
  return 0;
}

/*P*/
particles readSimulation(char simfile[],SHeader* sheader) 
{
  file fs=fileOpen(simfile,"r");
  char ctmp[500];
  int i,j;
  particles parts;

  //======================================================================
  //READ HEADER
  //======================================================================
  fscanf(fs,"%s %d %f",ctmp,&sheader->ntot,&sheader->mtot);
  fscanf(fs,"%s",ctmp);
  for(i=0;i<=5;i++)
    fscanf(fs,"%d",&sheader->npartg[i]);
  fscanf(fs,"%s",ctmp);
  for(i=0;i<=5;i++)
    fscanf(fs,"%f",&sheader->mtotg[i]);
  fscanf(fs,"%s %d",ctmp,&sheader->nsubs);
  fscanf(fs,"%s",ctmp);
  for(i=0;i<sheader->nsubs;i++)
    fscanf(fs,"%d",&sheader->stype[i]);
  fscanf(fs,"%s",ctmp);
  for(i=0;i<sheader->nsubs;i++)
    fscanf(fs,"%s",sheader->subtype[i]);
  fscanf(fs,"%s",ctmp);
  for(i=0;i<sheader->nsubs;i++)
    fscanf(fs,"%d",&sheader->nparts[i]);
  fscanf(fs,"%s",ctmp);
  for(i=0;i<sheader->nsubs;i++)
    fscanf(fs,"%f",&sheader->mtots[i]);
  fscanf(fs,"%s %f %f %f",ctmp,&sheader->jtot.x,&sheader->jtot.y,&sheader->jtot.z);
  fscanf(fs,"%s %f %f %f",ctmp,&sheader->rcm.x,&sheader->rcm.y,&sheader->rcm.z);
  fscanf(fs,"%s %f %f %f",ctmp,&sheader->vcm.x,&sheader->vcm.y,&sheader->vcm.z);
  fscanf(fs,"%s %f %f %f %f",ctmp,
	 &sheader->z,&sheader->omegaL,&sheader->Ho,&sheader->boxsize);
  fscanf(fs,"%s %d",ctmp,&sheader->rf);
  fgets(ctmp,sizeof ctmp,fs);fgets(ctmp,sizeof ctmp,fs);

  //======================================================================
  //INFORMATION IN HEADER FILE COULD BE DIFFERENT THAT ACTUAL DATA
  //======================================================================
  sheader->ntot=countLines(simfile);

  //======================================================================
  //READ PARTICLE DATA
  //======================================================================
  parts=(particles)malloc(sheader->ntot*sizeof(particle));
  for(i=0;i<sheader->ntot;i++){
    fgets(ctmp,sizeof ctmp,fs);
    str2particle(ctmp,&parts[i]);
    //fprintf(stdout,"Particle %d: %s\n",i,particle2str(parts[i]));
    //break;
  }
  fclose(fs);

  //Update particle and header information
  updateParticles(parts,*sheader);
  updateHeader(parts,sheader);

  return parts;
}

/*P*/
cadena str2particle(char str[],particle *part) 
{
  sscanf(str,
	 "%d %d %d %d "
	 "%f %f "
	 "%f %f %f "
	 "%f %f %f "
	 "%f %f %f "
	 "%f %f %f "
	 "%f %f %f "
	 "%f %f %f %f ",
	 &part->sid,&part->pid,&part->tid,&part->gtype,
	 &part->m,&part->h,
	 &part->pos.x,&part->pos.y,&part->pos.z,
	 &part->vel.x,&part->vel.y,&part->vel.z,
	 &part->posc.x,&part->posc.y,&part->posc.z,
	 &part->velc.x,&part->velc.y,&part->velc.z,
	 &part->j.x,&part->j.y,&part->j.z,
	 &part->r,&part->vcirc,&part->jcirc,&part->epsilon
	 );
  return 0;
}

/*P*/
cadena particle2str(particle part) 
{
  cadena str;
  str=(cadena)malloc(500*sizeof(char));
  sprintf(str,
	  "%-7d %-7d %-7d %-7d "
	  "%+14.7e %+14.7e "
	  "%44s %44s %44s %44s "
	  "%44s "
	  "%+14.7e %+14.7e %+14.7e %+14.7e ",
	  part.sid,part.pid,part.tid,part.gtype,
	  part.m,part.h,
	  vector2str(part.pos),vector2str(part.vel),
	  vector2str(part.posc),vector2str(part.velc),
	  vector2str(part.j),
	  part.r,part.vcirc,part.jcirc,part.epsilon
	  );
  return str;
}

/*
  Usage:
      readGadgetParam("SNAP.param","SofteningDisk",'f',&T);
      readGadgetParam("SNAP.param","ICFormat",'d',&d);
      readGadgetParam("SNAP.param","InitCondFile",'s',cadena);
 */
/*P*/
int readGadgetParam(char paramfile[],const char variable[],char type,void *result)
{
  char tmpfile[]="/tmp/readgad";
  char cmd[FSIZE];
  char value[100];
  float tmp;

  sprintf(cmd,"grep '^%s ' %s | awk '{print $2}' > %s",variable,paramfile,tmpfile);
  system(cmd);
  FILE *fp=fopen(tmpfile,"r");
  fscanf(fp,"%s",value);
  fclose(fp);
  if(VERBOSE(2)) 
    fprintf(stdout,"Value read from param. file : %s\n",value);
  switch(type){
  case 'f':{
    float *valread=(float*)result;
    *valread=atof(value);
    break;
  }
  case 'd':{
    int *valread=(int*)result;
    *valread=atoi(value);
    break;
  }
  case 's':{
    char *valread=(char*)result;
    strcpy(valread,value);
    break;
  }
  }
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//ANALYSIS
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
/*P*/
int transformSimulation(particles parts,SHeader *sheader,
			vector ro,vector vo,vector jo)
{
  int i;
  real M1[3][3],M2[3][3],M3[3][3],M4[3][3];
  vector vec1,vec2;
  real alpha,beta,phi,theta,chi,delta;

  //////////////////////////////////////////
  //TRANSLATE TO RO
  //////////////////////////////////////////
  for(i=0;i<sheader->ntot;i++){
    sumVector(&parts[i].pos,-1,ro,parts[i].pos);
    sumVector(&parts[i].vel,-1,vo,parts[i].vel);
  }

  //////////////////////////////////////////
  //ROTATE TO JO
  //////////////////////////////////////////
  if(jo.x!=0 || jo.y!=0){
    //GET THE ROTATION ANGLES
    phi = atan(fabs(jo.y)/fabs(jo.x));
    theta = atan(fabs(jo.z)/fabs(jo.x));

    phi = phi*180.0/M_PI;
    theta = theta*180.0/M_PI;
    //printf("phi=%lf theta=%lf\n",phi,theta);

    //Z-ANGLE
    if((jo.x>=0.0) && (jo.y>=0.0))
      alpha = fabs(phi);
    if((jo.x<0.0) && (jo.y>=0.0))
      alpha = 180.0-fabs(phi);
    if((jo.x<0.0) && (jo.y<0.0))
      alpha = 180.0+fabs(phi);
    if((jo.x>=0.0) && (jo.y<0.0))
      alpha = 360.0-fabs(phi);

    //Y-ANGLE
    if((jo.x>=0.0) && (jo.z>=0.0))
      beta = 270.0+fabs(theta);
    if((jo.x<0.0) && (jo.z>=0.0))
      beta = 90.0-fabs(theta);
    if((jo.x<0.0) && (jo.z<0.0))
      beta = 90.0+fabs(theta);
    if((jo.x>=0.0) && (jo.z<0.0))
      beta = 270.0-fabs(theta);
  
    alpha = alpha*M_PI/180.0;
    beta = beta*M_PI/180.0;
    chi = chi*M_PI/180.0;
    delta = delta*M_PI/180.0;

    //rotate about Z 
    M1[0][0]=cos(alpha);
    M1[0][1]=sin(alpha);
    M1[0][2]=0.0;
    M1[1][0]=-sin(alpha);
    M1[1][1]=cos(alpha);
    M1[1][2]=0.0;
    M1[2][0]=0.0;
    M1[2][1]=0.0;
    M1[2][2]=1.0;
  
    //rotate about Y
    M2[0][0]=cos(beta);
    M2[0][1]=0.0;
    M2[0][2]=sin(beta);
    M2[1][0]=0.0;
    M2[1][1]=1.0;
    M2[1][2]=0.0;
    M2[2][0]=-sin(beta);
    M2[2][1]=0.0;
    M2[2][2]=cos(beta);
  
    //rotate about Z 
    M3[0][0]=cos(chi);
    M3[0][1]=sin(chi);
    M3[0][2]=0.0;
    M3[1][0]=-sin(chi);
    M3[1][1]=cos(chi);
    M3[1][2]=0.0;
    M3[2][0]=0.0;
    M3[2][1]=0.0;
    M3[2][2]=1.0;
  
    //rotate about X 
    M4[0][0]=1.0;
    M4[0][1]=0.0;
    M4[0][2]=0.0;
    M4[1][0]=0.0;
    M4[1][1]=cos(delta);
    M4[1][2]=sin(delta);
    M4[2][0]=0.0;
    M4[2][1]=-sin(delta);
    M4[2][2]=cos(delta);
  
    for(i=0;i<sheader->ntot;i++){
      vec1.x=parts[i].pos.x;
      vec1.y=parts[i].pos.y;
      vec1.z=parts[i].pos.z;      
      vec2.x=parts[i].vel.x;
      vec2.y=parts[i].vel.y;
      vec2.z=parts[i].vel.z;
      
      parts[i].pos.x = M1[0][0]*vec1.x + M1[0][1]*vec1.y + M1[0][2]*vec1.z;
      parts[i].pos.y = M1[1][0]*vec1.x + M1[1][1]*vec1.y + M1[1][2]*vec1.z;
      parts[i].pos.z = M1[2][0]*vec1.x + M1[2][1]*vec1.y + M1[2][2]*vec1.z;
      parts[i].vel.x = M1[0][0]*vec2.x + M1[0][1]*vec2.y + M1[0][2]*vec2.z;
      parts[i].vel.y = M1[1][0]*vec2.x + M1[1][1]*vec2.y + M1[1][2]*vec2.z;
      parts[i].vel.z = M1[2][0]*vec2.x + M1[2][1]*vec2.y + M1[2][2]*vec2.z; 
    
      vec1.x=parts[i].pos.x;
      vec1.y=parts[i].pos.y;
      vec1.z=parts[i].pos.z;      
      vec2.x=parts[i].vel.x;
      vec2.y=parts[i].vel.y;
      vec2.z=parts[i].vel.z; 
   
      parts[i].pos.x = M2[0][0]*vec1.x + M2[0][1]*vec1.y + M2[0][2]*vec1.z;
      parts[i].pos.y = M2[1][0]*vec1.x + M2[1][1]*vec1.y + M2[1][2]*vec1.z;
      parts[i].pos.z = M2[2][0]*vec1.x + M2[2][1]*vec1.y + M2[2][2]*vec1.z;
      parts[i].vel.x = M2[0][0]*vec2.x + M2[0][1]*vec2.y + M2[0][2]*vec2.z;
      parts[i].vel.y = M2[1][0]*vec2.x + M2[1][1]*vec2.y + M2[1][2]*vec2.z;
      parts[i].vel.z = M2[2][0]*vec2.x + M2[2][1]*vec2.y + M2[2][2]*vec2.z;
    
      vec1.x=parts[i].pos.x;
      vec1.y=parts[i].pos.y;
      vec1.z=parts[i].pos.z;      
      vec2.x=parts[i].vel.x;
      vec2.y=parts[i].vel.y;
      vec2.z=parts[i].vel.z; 
    
      parts[i].pos.x = M3[0][0]*vec1.x + M3[0][1]*vec1.y + M3[0][2]*vec1.z;
      parts[i].pos.y = M3[1][0]*vec1.x + M3[1][1]*vec1.y + M3[1][2]*vec1.z;
      parts[i].pos.z = M3[2][0]*vec1.x + M3[2][1]*vec1.y + M3[2][2]*vec1.z;
      parts[i].vel.x = M3[0][0]*vec2.x + M3[0][1]*vec2.y + M3[0][2]*vec2.z;
      parts[i].vel.y = M3[1][0]*vec2.x + M3[1][1]*vec2.y + M3[1][2]*vec2.z;
      parts[i].vel.z = M3[2][0]*vec2.x + M3[2][1]*vec2.y + M3[2][2]*vec2.z;
    
      vec1.x=parts[i].pos.x;
      vec1.y=parts[i].pos.y;
      vec1.z=parts[i].pos.z;      
      vec2.x=parts[i].vel.x;
      vec2.y=parts[i].vel.y;
      vec2.z=parts[i].vel.z; 
    
      parts[i].pos.x = M4[0][0]*vec1.x + M4[0][1]*vec1.y + M4[0][2]*vec1.z;
      parts[i].pos.y = M4[1][0]*vec1.x + M4[1][1]*vec1.y + M4[1][2]*vec1.z;
      parts[i].pos.z = M4[2][0]*vec1.x + M4[2][1]*vec1.y + M4[2][2]*vec1.z;
      parts[i].vel.x = M4[0][0]*vec2.x + M4[0][1]*vec2.y + M4[0][2]*vec2.z;
      parts[i].vel.y = M4[1][0]*vec2.x + M4[1][1]*vec2.y + M4[1][2]*vec2.z;
      parts[i].vel.z = M4[2][0]*vec2.x + M4[2][1]*vec2.y + M4[2][2]*vec2.z;
    } 
  }//END ROTATION

  //UPDATE PARTICLE AND HEADER
  updateParticles(parts,*sheader);
  updateHeader(parts,sheader);

  return 0;
}

/*P*/
real* sliceSimulation(particles parts,int ini,int end,size_t member)
{
  int i;
  real* vals=(real*)calloc((end-ini+1),sizeof(real));
  
  for(i=ini;i<=end;i++)
    vals[i-ini]=GET_MEMBER(parts[i],real,member);
  
  return vals;
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//VECTOR
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
/*P*/
int constVector(vector* v,real k)
{
  v->x=k;
  v->y=k;
  v->z=k;
}

/*P*/
int scaleVector(vector* v,real k,vector u)
{
  v->x=k*u.x;
  v->y=k*u.y;
  v->z=k*u.z;
  
}

/*P*/
int sumVector(vector* v,real k,vector u1,vector u2)
{
  v->x=k*u1.x+u2.x;
  v->y=k*u1.y+u2.y;
  v->z=k*u1.z+u2.z;
}

/*P*/
int crossVector(vector* v,vector u1,vector u2)
{
  v->x=+(u1.y*u2.z-u1.z*u2.y);
  v->y=-(u1.x*u2.z-u1.z*u2.x);
  v->z=+(u1.x*u2.y-u1.y*u2.x);
}

/*P*/
int copyVector(vector* v,vector u)
{
  v->x=u.x;
  v->y=u.y;
  v->z=u.z;
}

/*P*/
real magVector(vector v)
{
  return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

/*
  Transform the components of a vector to cyclindrical coordinates
 */
/*P*/
int cylVector(vector *v,vector u)
{
  real r;
  //v->x stands for R
  v->x=sqrt(u.x*u.x+u.y*u.y);
  r=sqrt(v->x*v->x+u.z*u.z);

  //v->y stands for phi
  v->y=atan2(u.y,u.x);
  v->y=v->y<0?2*M_PI+v->y:v->y;

  //v->z stands for theta
  //v->z=u.z;
  v->z=acos(u.z/r);
}

/*P*/
cadena vector2str(vector r)
{
  cadena str=(cadena)malloc(100*sizeof(char));
  sprintf(str,"%+14.7e %+14.7e %+14.7e",r.x,r.y,r.z);
  return str;
}

/*P*/
int arrayPrintf(file stream,const char frm[],real v[],int n)
{
  for(int i=0;i<n;i++)
    fprintf(stream,frm,v[i]);
  return 0;
}

/*P*/
int arrayVector(vector *v,real a[])
{
  v->x=a[0];
  v->y=a[1];
  v->z=a[2];
}

/*P*/
real arrayMin(real* array,int n)
{
  int i;
  real min=MAXREAL;
  for(i=0;i<n;i++)
    SET_MIN(min,array[i]);
  return min;
}

/*P*/
real arrayMax(real* array,int n)
{
  int i;
  real max=MINREAL;
  for(i=0;i<n;i++)
    SET_MAX(max,array[i]);
  return max;
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//UTILITARY
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
/*P*/
int fileExists(char *filename)
{
  FILE *f;
  if((f=fopen(filename,"r"))==NULL) return 0;
  else{
    fclose(f);
    return 1;
  }
}

/*P*/
int isBlank(char *string)
{
  if(strlen(string)) return 0;
  else return 1;
}

/*P*/
int isSpace(char *string)
{
  char *pos,tmp[LSIZE];
  strcpy(tmp,string);
  pos=strtok(tmp,"1234567890");
  if(pos==NULL)
    return 1;
  else if(strcmp(pos,"\n")==0 || strcmp(pos," ")==0)
    return 1;
  else
    return 0;
}

/*P*/
int pause(void)
{
  fprintf(stderr,"<ENTER TO CONTINUE/CTRL+C TO CANCEL>");
  getc(stdin);
}

/*P*/
file fileOpen(const char filename[],const char mode[])
{
  file f;

  if((f=fopen(filename,mode))==NULL){
    fprintf(stderr,"Error DYNAMO: File '%s' cannot be openned as '%s'\n",
	    filename,mode);
    exit(1);
  }
  return f;
}

/*P*/
int splitString(char str[],const char sep[],real list[])
{
  char wstr[LSIZE];
  char tmp[50];
  char *pos;
  strcpy(wstr,str);
  pos=strtok(wstr,sep);

  int i=0;
  while(pos!=NULL){
    strcpy(tmp,pos);
    list[i]=atof(tmp);
    pos=strtok(NULL,sep);
    i++;
  }
  return i;
}

/*P*/
int splitString(char str[],const char sep[],real2 list[])
{
  char wstr[LSIZE];
  char tmp[50];
  char *pos;
  strcpy(wstr,str);
  pos=strtok(wstr,sep);

  int i=0;
  while(pos!=NULL){
    strcpy(tmp,pos);
    list[i]=atof(tmp);
    pos=strtok(NULL,sep);
    i++;
  }
  return i;
}

/*P*/
int splitString(char str[],const char sep[],char list[][LSIZE])
{
  char wstr[LSIZE];
  char tmp[50];
  char *pos;
  strcpy(wstr,str);
  pos=strtok(wstr,sep);

  int i=0;
  while(pos!=NULL){
    strcpy(tmp,pos);
    strcpy(list[i],tmp);
    pos=strtok(NULL,sep);
    i++;
  }
  return i;
}

/*P*/
int grepString(char str[],const char search[])
{
  if(strstr(str,search)==NULL) return 0;
  else return 1;
}

/*P*/
int grepArray(real value,real array[],int n,real tol/*=1E-5*/)
{
  for(int i=0;i<n;i++){
    if(fabs(value-array[i])<tol) return 1;
  }
  return 0;
}

/*P*/
int readLine(char row[],real array[],int nfields)
{
  char *pos,str[LSIZE],tmp[LSIZE];
  strcpy(str,row);
  pos=strtok(str," \t");
  for(int i=0;i<nfields;i++){
    strcpy(tmp,pos);
    array[i]=strtof(tmp,NULL);
    pos=strtok(NULL," \t");
  }
  return 0;
}

/*P*/
int readLine(char row[],real2 array[],int nfields)
{
  char *pos,str[LSIZE],tmp[LSIZE];
  strcpy(str,row);
  pos=strtok(str," \t");
  for(int i=0;i<nfields;i++){
    strcpy(tmp,pos);
    array[i]=strtod(tmp,NULL);
    pos=strtok(NULL," \t");
  }
  return 0;
}

/*P*/
int systemOutput(const char cmd[],char out[])
{
  char comm[FSIZE]="",tmpf[FSIZE]="",line[LSIZE]="";
  sprintf(tmpf,"/tmp/sysout");
  sprintf(comm,"%s > %s",cmd,tmpf);
  system(comm);
  file f=fileOpen(tmpf,"r");
  while(1){
    fgets(line,sizeof line,f);
    if(feof(f)) break;
    strcat(out,line);
  }
  fclose(f);
  return 0;
}

/*P*/
int countLines(const char file[])
{
  int numlines=0;
  char cmd[LSIZE]="",out[LSIZE]="";
  
  sprintf(cmd,"grep -v '#' %s | grep -v '^$' | wc -l",file);
  systemOutput(cmd,out);
  numlines=atoi(out);
  return numlines;
}

/*P*/
real2* readColumn(const char fname[],int nlines,int ncols,int col)
{
  real2 *line,*values;
  char linea[LSIZE]="";
  file fd=fileOpen(fname,"r");
 
  line=(real2*)malloc(ncols*sizeof(real2));
  values=(real2*)malloc(nlines*sizeof(real2));

  //Read values
  int i=-1;
  while(1){
    fgets(linea,sizeof linea,fd);
    //fprintf(stdout,"LINE: %s",linea);
    if(feof(fd)) break;
    //Check for comments
    if(linea[0]=='#' ||
       isSpace(linea)){
      //fprintf(stdout,"EXCLUDED LINE\n");
      continue;
    }
    //Check for blank lines
    //if() continue;
    i++;
    readLine(linea,line,ncols);
    values[i]=line[col-1];
  }

  fclose(fd);
  return values;
}

/*P*/
real* readColumns(const char fname[],int nlines,int ncols,int col)
{
  real *line,*values;
  char linea[LSIZE]="";
  file fd=fileOpen(fname,"r");
 
  line=(real*)malloc(ncols*sizeof(real2));
  values=(real*)malloc(nlines*sizeof(real2));

  //Read values
  int i=-1;
  while(1){
    fgets(linea,sizeof linea,fd);
    if(feof(fd)) break;
    if(linea[0]=='#' ||
       isSpace(linea)) continue;
    i++;
    readLine(linea,line,ncols);
    values[i]=line[col-1];
  }

  fclose(fd);
  return values;
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//ADVANCED NUMERICAL ROUTINES
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
/*P*/
real2 chiSquare(const gsl_vector* param,pars tmp)
{
  int i;
  real2 chis;
  real2 xo,yo,s,yt;
  real2* ps;
  ps=param->data;
  
  chis=0;
  for(i=0;i<NFit;i++){
    //Calculo valores teoricos
    xo=XFit[i];
    yo=YFit[i];
    s=EFit[i];
    yt=FFit(xo,ps);
    if(EQUAL(s,0)) s=1;
    chis+=(yt-yo)*(yt-yo)/(s*s);
  }
  return chis;
}

/*P*/
int fitData(real2* X,real2* Y,real2* E,int Ndata,
	    real2 (*func)(real2,pars),int Npars,real2 Params[],
	    real2& cs,real2& pval)
{
  #ifndef MAXITER
  #define MAXITER 1000
  #endif

  //==================================================
  //Initialize global fitting variables
  //==================================================
  XFit=X;
  YFit=Y;
  EFit=E;
  NFit=Ndata;
  FFit=func;

  //==================================================
  //Initialize parameters
  //==================================================
  int i;
  gsl_vector *params=gsl_vector_alloc(Npars),*rparams;
  gsl_vector *dparams=gsl_vector_alloc(Npars);
  for(i=0;i<Npars;i++)
    gsl_vector_set(params,i,Params[i]);
  gsl_vector_set_all(dparams,1E-1);

  //==================================================
  //Declare and initialize minimizer
  //==================================================
  gsl_multimin_fminimizer*
    mins=gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex,FNpars);
  gsl_multimin_function function;
  function.f=&chiSquare;
  function.n=Npars;
  function.params=NULL;
  gsl_multimin_fminimizer_set(mins,&function,params,dparams);
  
  //==================================================
  //Minimize chisquare
  //==================================================
  int status;
  real2 ftom;
  i=0;
  do{
    gsl_multimin_fminimizer_iterate(mins);
    ftom=gsl_multimin_fminimizer_size(mins);
    status=gsl_multimin_test_size(ftom,1E-5);
    if((i%(MAXITER/10))==0)
      VPRINTF(2)("\tIteration %d: %e\n",i,ftom);
    i++;
  }while(status!=GSL_SUCCESS && i<MAXITER);
  VPRINTF(2)("Fitting results:\n");
  VPRINTF(2)("\tNumber of iterations: %d\n",i);
  VPRINTF(2)("\tDistance to minimum: %g\n",ftom);
  if(i>=MAXITER)
    VPRINTF(2)("\tMaximum number of iterations reached\n",ftom);

  //==================================================
  //Return the best fit parameters
  //==================================================
  rparams=gsl_multimin_fminimizer_x(mins);
  for(i=0;i<Npars;i++)
    Params[i]=gsl_vector_get(rparams,i);

  //==================================================
  //Return the chisquare and the p-value
  //==================================================
  cs=gsl_multimin_fminimizer_minimum(mins);
  
  //Degrees of freedom are assumed equal to the number of points 
  int nu=Ndata;
  pval=gsl_cdf_chisq_Q(cs,nu);

  gsl_multimin_fminimizer_free(mins);
  gsl_vector_free(params);
  gsl_vector_free(dparams);

  return 0;
}

