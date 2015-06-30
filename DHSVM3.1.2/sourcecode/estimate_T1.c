/*
 * SUMMARY:      estimate_T1.c - Calculate soil temperature
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Keith Cherkauer
 * ORG:          University of Washington, Department of Civil Engineering
 * ORIG-DATE:    Jul-1998   
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
/**********************************************************************
  estimate_T1                

  uses Xu Liangs 3-layer energy balance formulation to estimate the 
  temperature between the first and second layers.  Formerly calculated
  independently in each of the surface energy balance equation routines.

  Modifications:
  01-20-00 removed from end of func_surf_energy_bal.c and put into a
           separate file                                           KAC

**********************************************************************/
float estimate_T1(float Ts, 
		   float T1_old,
		   float T2,
		   float D1, 
		   float D2, 
		   float kappa1, 
		   float kappa2, 
		   float Cs1, 
		   float Cs2, 
		   float dp,
		   int delta_t) {

  float C1;
  float C2;
  float C3;
  float T1;

  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);

  T1 = (kappa1/2./D1/D2*(Ts) + C1/delta_t*T1_old
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);

  return(T1);

}
