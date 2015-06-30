/*
 * SUMMARY:      SensibleHeatFlux.c - Calculate sensible heat flux
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate sensible heat flux
 * DESCRIP-END.
 * FUNCTIONS:    SensibleHeatFlux()
 *               NoSensibleHeatFlux()
 * COMMENTS:
 * $Id: SensibleHeatFlux.c,v 1.4 2003/07/01 21:26:24 olivier Exp $     
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "brent.h"
#include "functions.h"

/*****************************************************************************
  SensibleHeatFlux()

  Source: Liang, X., 1995, Modeling ground heat flux in land surface 
          parameterization schemes

  This function calculates outgoing longwave, sensible heat flux, ground heat flux, 
  and storage of heat in the thin upper layer, based on 3 soil layer temperature.
  The first two layers are defined by the soil layer depth of D1 and D2. The 3rd 
  layer is defined by the damping depth (set to a constant 4.0 m) at which soil 
  temperature variation is negligible.  
*****************************************************************************/
void SensibleHeatFlux(int y, int x, int Dt, float Displacement, 
			  MET_MAP_PIX *LocalMet, EVAPPIX *LocalEvap, SOILTABLE *SoilType,
		      SOILPIX *LocalSoil, SNOWPIX *LocalSnow, VEGTABLE *VType, 
			  PRECIPPIX *LocalPrecip)
{
  int NSoilLayers;
  int Tsurf_fbflag;     /* Flag to indicate when any temperature iterations fail to converge */ 
  int Tsurf_fbcount;    /* count of numbers of Tsurf_fbflag */
  float F;              /* fractional coverage of overstory (%) */
  float dp;		        /* the depth at which the interannual change in soil temperature 
						   is negligible (m) */
  float LowerRa;		/* Aerodynamic resistance for lower layer (s/m) */
  float LowerWind;	    /* Wind for lower layer (m/s) */
  float MaxTSurf;		/* Upper bracket for effective surface temperature (C) */
  float MinTSurf;		/* Lower bracket for effective surface temperature (C) */
  float OldTSurf;		/* Effective surface temperature at the end of the last timestep (C) */
  float TMean;			/* Average surface temperature (C) */
  float D1;             /* Depth of the first thermal soil layer (m) */
  float D2;             /* Depth of the second thermal soil layer (m) */
  float D3;             /* Depth of the third thermal soil layer (m) */
  float kappa1;         /* Soil thermal conductivity of the top thin layer (W/m/K) */
  float kappa2;         /* Soil thermal conductivity of the layer of depth [D1, D1+D2] (W/m/K) */
  float Cs1;            /* Soil heat volumetric capacity oof the layer of depth [D1, D1+D2] (W/m^3/K) */
  float Cs2;            /* Soil heat volumetric capacity of the top thin layer (W/m^3/K) */
  float T1_old;         /* Soil temperature at the boundary between 1st and 2nd layers (C) */
  float T1;             /* Soil temperature at the boundary between 1st and 2nd layers (C) */
  float T2;             /* Soil temperature at the damping depth where the soil temperature variation
						   is negligible (C). (>> dp; *NOT* at depth D2) */
  float deltaH;         /* Change in heat capacity */
  float grnd_flux;      /* ground heat flux at the depth D1	*/
  float UpperRa;		/* Aerodynamic resistance for upper layer (s/m) */
  float UpperWind;	    /* Wind for upper layer (m/s) */
  float ZRef;		    /* Reference height for sensible heat calculation (m) */
  float Z0;			    /* Roughness length (m) */
  double Tmp;			/* Temporary value */

  /* Initialize T_fbflag */
  Tsurf_fbflag = 0;
  Tsurf_fbcount = 0;

  /* Layers of the soil */
  NSoilLayers = SoilType->NLayers;
  
  if (LocalSnow->HasSnow == TRUE) {
    ZRef = 2. + Z0_SNOW;
    Z0 = Z0_SNOW;
  }
  else {
    ZRef = 2. + Z0_GROUND;
    Z0 = Z0_GROUND;
  }
  
  /* calculate the actual aerodynamic resistances and wind speeds */
  UpperWind = VType->U[0] * LocalMet->wind_speed;
  UpperRa = VType->Ra[0] / LocalMet->wind_speed;
  if (VType->OverStory == TRUE) {
    LowerWind = VType->U[1] * LocalMet->wind_speed;
    LowerRa = VType->Ra[1] / LocalMet->wind_speed;
  }
  else {
    LowerWind = UpperWind;
    LowerRa = UpperRa;
  }
  
  /**************************************************
    Set All Variables For Use
  **************************************************/
  OldTSurf = LocalSoil->TSurf;
  T1_old = LocalSoil->Temp[1];
  T1 = 0.;
  T2 = avg_temp; // a constant lower boundary tempeature
  D1 = VType->RootDepth[0];
  D2 = VType->RootDepth[1];
  D3 = VType->RootDepth[2];

  if (VType->OverStory == TRUE) 
	F = VType->Fract[0];
  else 
	F = 0.;
  
  /* Set the upper and lower boundary */
  MaxTSurf = 0.5*(LocalSoil->TSurf+LocalMet->air_temp)+DELTAT;
  MinTSurf = 0.5*(LocalSoil->TSurf+LocalMet->air_temp)-DELTAT;

  /* The damping depth is defined as the depth at which the interannual soil 
     temperature variation is negligible. According to VIC-based studies, by using 
	 a deeper dp, Liang's approach produces comparable results to the analytical solution */
  dp = 4.0;

  /* Calculate the effective thermal conductivity of the top and second layer */
  CalcEffectiveKh(NSoilLayers, dp, VType->RootDepth,
	  SoilType->KhDry, SoilType->KhSol, LocalSoil->Moist,
	  SoilType->Porosity, LocalSoil->Temp, LocalSoil->LayerKh);
  /* Comments: the variable LocalSoil->Temp is valid to be used as the criteria here????
     unless the new Ts and Ts values are transferred back?? */
  kappa1 = LocalSoil->LayerKh[0];
  kappa2 = LocalSoil->LayerKh[1];

  /* Calculate volumetric heat capacities in J/m^3/K */
  CalcHeatCapacity(NSoilLayers, dp, VType->RootDepth, SoilType->Porosity, 
	  LocalSoil->Moist,LocalSoil->Temp, LocalSoil->Ch);
  Cs1 = LocalSoil->Ch[0];
  Cs2 = LocalSoil->Ch[1];

  /* Find surface temperature using root brent method */
  grnd_flux = 0.;
  deltaH = 0.;
  LocalSoil->TSurf =
    RootBrent(y, x, MinTSurf, MaxTSurf, SurfaceEnergyBalance, Dt, F, LowerRa, 
	      ZRef, Displacement, Z0, LocalMet->wind_speed, LocalSoil->NetShortIn, 
		  LocalSoil->NetLongIn, LocalMet->air_dens, LocalMet->Lv, 
		  LocalEvap->ETot, kappa1, kappa2, Cs1, Cs2, LocalSoil->tau,
		  dp, D1, D2, LocalMet->air_temp, &T1, T2, T1_old, OldTSurf, 
		  LocalSoil->MeltEnergy, &grnd_flux, &deltaH);

  if (LocalSoil->TSurf <= -998) {
	LocalSoil->TSurf = OldTSurf;
	Tsurf_fbflag = 1;
    Tsurf_fbcount++;
  }
  
  LocalSoil->Temp[0] = LocalSoil->TSurf;
  LocalSoil->Temp[1] = T1;
  LocalSoil->Temp[2] = T2+(T1-T2)*exp(-(D2+D3)/dp);
  
  /* Recalculate the terms of the energy balance. This is similar to the
     code in SurfaceEnergyBalance.c */
  TMean = LocalSoil->TSurf;

  if (LocalMet->wind_speed > 0.0)
	 LowerRa /= StabilityCorrection(ZRef, Displacement, TMean, LocalMet->air_temp,
			    LocalMet->wind_speed, Z0);
  else
    LowerRa = DHSVM_HUGE;
  LocalSoil->Ra = LowerRa;

  /**/
  Tmp = TMean + 273.15;
  LocalSoil->Qnet = LocalSoil->NetShortIn + LocalSoil->NetLongIn 
	  - STEFAN * (Tmp * Tmp * Tmp * Tmp);
  
  /* Calculate the effective thermal conductivity of the top and second layer */
  CalcEffectiveKh(NSoilLayers, dp, VType->RootDepth,
	  SoilType->KhDry, SoilType->KhSol, LocalSoil->Moist,
	  SoilType->Porosity, LocalSoil->Temp, LocalSoil->LayerKh);
  kappa1 = LocalSoil->LayerKh[0];
  kappa2 = LocalSoil->LayerKh[1];

  /* Calculate volumetric heat capacities in J/m^3/K */
  CalcHeatCapacity(NSoilLayers, dp, VType->RootDepth, SoilType->Porosity, 
	  LocalSoil->Moist, LocalSoil->Temp, LocalSoil->Ch);
  Cs1 = LocalSoil->Ch[0];
  Cs2 = LocalSoil->Ch[1];

  /* Compute the sensible heat flux from the surface */
  LocalSoil->Qs = LocalMet->air_dens * CP * (LocalMet->air_temp - TMean) / LowerRa;

  /* Compute the latent heat from the Surface and covering Vegetation */
  LocalSoil->Qe = -(LocalMet->Lv * LocalEvap->ETot) / Dt * WATER_DENSITY;
  
  /* Compute the ground heat flux from the top soil layer */
  LocalSoil->Qg = ((1 - F) + F * LocalSoil->tau) * (kappa1 / D1 * (T1 - TMean));
  /*LocalSoil->Qg = ((1 - F) + F * LocalSoil->tau) * (kappa1 / D1 * (T1 - TMean)
	  + (kappa2 / D2 * ( 1. - exp( -D2 / dp )) * (T2 - T1))) / 2.;*/
  
  LocalSoil->Qst = (Cs1 * ((OldTSurf + T1_old) - (TMean + T1)) * D1 / Dt / 2.);

  /* Calculate the energy balance error term */
  LocalSoil->Qrest = LocalSoil->Qnet + LocalSoil->Qs + LocalSoil->Qe +
    LocalSoil->Qg + LocalSoil->Qst + LocalSoil->MeltEnergy;

  /* debug */
  if (y == 395 && x > 109 && x < 180)
    printf("%2.0f%% F=%.2f t=%4.1f %4.1f %4.1f %4.1f Qnet=%3.0f Qs=%4.0f Qg=%4.0f Qe=%4.0f Qst=%3.0f\n",
	VType->ImpervFrac*100, F, LocalSoil->TSurf, OldTSurf, T1, LocalMet->air_temp, LocalSoil->Qnet, LocalSoil->Qs, 
	  LocalSoil->Qg, LocalSoil->Qe, LocalSoil->Qst);
}

/*****************************************************************************
  Function name: NoSensibleHeatFlux()

  Purpose      : Calculate latent heat flux in W/m2

  Required     : 
    int Dt             - Model timestep (seconds)
    PIXMET LocalMet    - Met data for the current pixel 
    float ETot         - Total vapor flux (mm/timestep) 
    SOILPIX *LocalSoil - Structure with soil moisture data and energy data
                         for the current pixel

  Returns      :
    void

  Modifies     :
    members of LocalSoil

  Comments     : This function sets all the energy fluxes at the pixel level 
                 to 0.0, but calculates the evapotranspiration in W/m2
*****************************************************************************/
void NoSensibleHeatFlux(int Dt, PIXMET * LocalMet, float ETot,
			SOILPIX * LocalSoil)
{
  LocalSoil->TSurf = 0.0;
  LocalSoil->Ra = 0.0;
  LocalSoil->Qnet = 0.0;
  LocalSoil->Qs = 0.0;
  LocalSoil->Qe = -(LocalMet->Lv * ETot) / Dt * WATER_DENSITY;
  LocalSoil->Qg = 0.0;
  LocalSoil->Qst = 0.0;
  LocalSoil->Qrest = 0.0;
}
