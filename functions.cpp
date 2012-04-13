/*
  Set of custom functions
*/
#include <dynamo.h>
/*D*/
real2 (*FFunc)(real2,pars);
int FNpars;
real2 FXtest;
/*D*/

/*
  ////////////////////////////////////////////////////////////////////////////////
  For each function you add, add the respective conditional here
  ////////////////////////////////////////////////////////////////////////////////
  */
/*P*/
int getFunction(const char funcname[])
{
  /*
    FFunc: Function
    FNpars: Number of free parameters
    FXtest: X value to test function
   */

  if(strcmp(funcname,"lineFunction")==0){
    FFunc=lineFunction;
    FNpars=2;
    FXtest=1.0;
  }

  if(strcmp(funcname,"asymBell")==0){
    FFunc=asymBell;
    FNpars=3;
    FXtest=1.0;
  }

  if(strcmp(funcname,"gaussianBell")==0){
    FFunc=gaussianBell;
    FNpars=3;
    FXtest=1.0;
  }

  if(strcmp(funcname,"exponentialDistribution")==0){
    FFunc=exponentialDistribution;
    FNpars=3;
    FXtest=1.0;
  }
}

//////////////////////////////////////////////////////////////////////////////////
//CUSTOM FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////
/*P*/
real2 lineFunction(real2 x,pars params)
{
  /*
    ps[0]: steep
    ps[1]: intercept
   */
  real2* ps=(real2*)params;
  real2 y;

  y=ps[0]*x+ps[1];

  return y;
}

/*P*/
real2 asymBell(real2 x,pars params)
{
  /*
    ps[0] x * exp[ - (x-ps[1])/ps[2] ]
  */
  real2* ps=(real2*)params;
  real2 y;
  
  y=ps[0]*x*exp(-(x-ps[1])/ps[2]);

  return y;
}

/*P*/
real2 gaussianBell(real2 x,pars params)
{
  /*
    ps[0] exp [ - (x-ps[1])^2 / ps[2]^2 ]
  */
  real2* ps=(real2*)params;
  real2 y;
  
  y=ps[0]*exp(-(x-ps[1])*(x-ps[1])/(ps[2]*ps[2]));

  return y;
}

/*P*/
real2 exponentialDistribution(real2 x,pars params)
{
  /*
    ps[0] * exp ( - |x - ps[1]| / ps[2] )
   */
  real2* ps=(real2*)params;
  real2 y;

  y=ps[0]*exp(-fabs(x-ps[1])/ps[2]);

  return y;
}
