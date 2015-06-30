/*
 * SUMMARY:      SurfaceEnergyBalance.c - Calculate surface energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculate surface energy balance.  This group of functions
 *               is used by the iterative Brent method to determine the
 *               surface temperature 
 * DESCRIP-END.
 * FUNCTIONS:    SurfaceEnergyBalance()
 * COMMENTS:
 * $Id: SurfaceEnergyBalance.c,v 1.4 2003/07/01 21:26:26 olivier Exp $     
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  Function name: SurfaceEnergyBalance()

  Purpose      : Calculate the surface energy balance in the absence of snow

  Required     :
    float TSurf           - new estimate of effective surface temperature
    va_list ap            - Argument list initialized by va_start().  For
                            elements of list and order, see beginning of
                            routine

  Returns      :
    float RestTerm        - Rest term in the energy balance

  Modifies     : none

  Comments     :
*****************************************************************************/
float SurfaceEnergyBalance(float TSurf, va_list ap)
{
  /* start of list of arguments in variable argument list */

  int Dt;				/* Model time step (seconds) */
  float F;              /* Fractional coverage (%) */
  float Ra;				/* Aerodynamic resistance (s/m) */
  float Z;				/* Reference height (m) */
  float Displacement;	/* Displacement height (m) */
  float Z0;				/* Surface roughness (m) */
  float Wind;			/* Wind speed (m/s) */
  float ShortRad;		/* Net incident shortwave radiation (W/m2) */
  float LongRadIn;		/* Incoming longwave radiation (W/m2) */
  float AirDens;		/* Density of air (kg/m3) */
  float Lv;			    /* Latent heat of vaporization (J/kg3) */
  float ETot;			/* Total evapotranspiration (m) */


  float kappa1;			/* Top layer effective soil thermal conductivity (W/(m*K)) */
  float kappa2;         /* Second layer effective soil thermal conductivity (W/(m*K)) */
  float Cs1;		    /* Top layer soil thermal capacity (J/(m3*C)) */
  float Cs2;		    /* Second layer soil thermal capacity (J/(m3*C)) */
  float tau;            /* Transmittance for overstory vegetation layer */
  float dp;			    /* Damping depth (m) */
  float D1;				/* Depth of 1st soil heat profile (m) */
  float D2;				/* Depth of 2nd soil heat profile (m) */
  float Tair;			/* Air temperature (C) */
  float *T1;            /* Soil temperature of 1st soil layer (C) */
  float T2;       		/* Average soil temperature at the deep soil layer (C) */
  float T1_old;			/* Soil temperature at pervious time step (C) */
  float OldTSurf;		/* Surface temperature during previous time step */
  float MeltEnergy;		/* Energy used to melt/refreeze snow pack (W/m2) */
  float *deltaH;
  float *grnd_flux;
  /* end of list of arguments in variable argument list */

  float LatentHeat;		    /* latent heat exchange at surface (W/m2) */
  float LongRadOut;			/* long wave radiation emitted by surface (W/m2) */
  float NetRad;				/* net radiation exchange at surface (W/m2) */
  float RestTerm;			/* rest term in surface energy balance (W/m2) */
  float SensibleHeat;		/* sensible heat exchange at surface (W/m2) */
  float TMean = 0.;			/* Mean temperature during interval (C) */
  double Tmp;				/* temporary variable */

  /* Assign the elements of the array to the appropriate variables.  The list
     is traversed as if the elements are doubles, because:
     In the variable-length part of variable-length argument lists, the old
     ``default argument promotions'' apply: arguments of type float are
     always promoted (widened) to type double, and types char and short int
     are promoted to int. Therefore, it is never correct to invoke
     va_arg(argp, float); instead you should always use va_arg(argp,
     double). 
     (quoted from the comp.lang.c FAQ list)
   */

  Dt = va_arg(ap, int);
  F = (float)va_arg(ap, double);
  Ra = (float) va_arg(ap, double);
  Z = (float) va_arg(ap, double);
  Displacement = (float) va_arg(ap, double);
  Z0 = (float) va_arg(ap, double);
  Wind = (float) va_arg(ap, double);
  ShortRad = (float) va_arg(ap, double);
  LongRadIn = (float) va_arg(ap, double);
  AirDens = (float) va_arg(ap, double);
  Lv = (float) va_arg(ap, double);
  ETot = (float) va_arg(ap, double);
  kappa1 = (float) va_arg(ap, double);
  kappa2 = (float) va_arg(ap, double);
  Cs1 = (float) va_arg(ap, double);
  Cs2 = (float) va_arg(ap, double);
  tau = (float) va_arg(ap, double);
  dp = (float) va_arg(ap, double);
  D1 = (float) va_arg(ap, double);
  D2 = (float) va_arg(ap, double);
  Tair = (float) va_arg(ap, double);
  T1  = (float *) va_arg(ap, double*);
  T2  = (float) va_arg(ap, double);
  T1_old  = (float) va_arg(ap, double);
  OldTSurf = (float) va_arg(ap, double);
  MeltEnergy = (float) va_arg(ap, double);
  grnd_flux  = (float *) va_arg(ap, double*);
  deltaH  = (float *) va_arg(ap, double*);

  /* In this routine transport of energy to the surface is considered positive */
  TMean = TSurf;
  //TMean = 0.5*(TSurf+OldTSurf);
  Tmp = TMean + 273.15;
 
  /* Apply the stability correction to the aerodynamic resistance */
  if (Wind > 0.0)
    Ra /= StabilityCorrection(Z, Displacement, TMean, Tair, Wind, Z0);
  else
    Ra = DHSVM_HUGE;

  /* Use Liang et al. 1999 Equations to Calculate Ground Heat Flux.
     NOTE: T2 is not the temperature of the 2nd soil layer, nor at depth dp;
     T2 is the constant temperature at depths (much) greater than dp.*/
  *T1 = (float) estimate_T1(TMean, T1_old, T2, D1, D2, kappa1, kappa2, Cs1, Cs2, dp, Dt);

  /* Calculate the ground heat flux at the depth D1.
     The equation is modified from the original Liang's equation. */
  *grnd_flux = ((1 - F) + F * tau) * (kappa1 / D1 * ((*T1) - TMean));
  /**grnd_flux = ((1 - F) + F * tau) * (kappa1 / D1 * ((*T1) - TMean)
	  + (kappa2 / D2 * ( 1. - exp( -D2 / dp )) * (T2 - (*T1)))) / 2.;*/

  /* Calculate the heat storage change in the ground heat storage in the upper soil */
  *deltaH = (Cs1 * ((OldTSurf + T1_old) - (TMean + *T1)) * D1 / Dt / 2.);

  /* Compute net surface radiation for evaporation estimates */
  LongRadOut = STEFAN * (Tmp * Tmp * Tmp * Tmp);
  NetRad = ShortRad + LongRadIn - LongRadOut + *grnd_flux + *deltaH;

  /* Calculate the sensible heat flux */
  SensibleHeat = AirDens * CP * (Tair - TMean) / Ra;

  /* Calculate the latent heat flux */
  LatentHeat = -(Lv * ETot) / Dt * WATER_DENSITY;

  /* Calculate the net energy exchange at the surface.  The left hand side of 
     the equation should go to zero for the balance to close, so we want to 
     minimize the absolute value of the left hand side */
  RestTerm =
    MeltEnergy + NetRad + SensibleHeat + LatentHeat;

  return RestTerm;
}

