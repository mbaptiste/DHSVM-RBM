/*
 * SUMMARY:      CalcEffectiveKh.c - Calculate effective thermal conducitivity
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the effective thermal conductivity of a soil under
 *               dry conditions
 * DESCRIP-END.
 * FUNCTIONS:    CalcEffectiveKh()
 * COMMENTS:
 * $Id: CalcEffectiveKh.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "DHSVMerror.h"
#include "functions.h"

/*****************************************************************************
  CalcEffectiveKh()

  Source: Farouki, O. T., 1986, Thermal properties of soils, 
                          Trans Tech Publications

  This function calculates the effective thermal conductivity of each soil 
  layer based on the thermal conductivity under dry conditions, KhDry, and
  the thermal conductivity under saturated conditions, KhSat.  The latter
  differs for frozen and unfrozen soils, and is function of the effective
  solids thermal conductivity, KhSol, and the saturates soil moisture or 
  ice content.  

  The method followed here is Johansen's method, section 7.11 [Farouk, 1986]
 
  First the effective conductivity is calculated for each layer, after which
  the total effective thermal conductivity is calculated for the specified 
  depth.
*****************************************************************************/
void CalcEffectiveKh(int NSoilLayers, float Bottom, float *SoilDepth, 
			float *KhDry, float *KhSol, float *Moisture, float *Porosity, 
			float *TSoil, float *LayerKh)
{
  int i;			    /* counter */
  int NLayers;			/* Number of soil layers in depth interval */
  int StartLayer = 0;	/* First layer below top */
  char NoEndLayer;		/* flag to indicate whether an end layer has been determined */
  char NoStartLayer;	/* flag to indicate whether a start layer has been determined */
  float Dz;			    /* Depth from soil surface (m) */
  float Ke;			    /* Kersten number (see reference) */
  float KhSat;			/* Thermal conductivity for saturated soils (W/(m*K)) */
  float Sr;			    /* degree of saturation */
  float TotalDepth;		/* Depth of soil column for which to  calculate the effective thermal 
				           conductivity (m) */

  TotalDepth = Bottom;

  Dz = 0.0;
  NLayers = 0;
  NoStartLayer = TRUE;
  NoEndLayer = TRUE;

  for (i = 0; i < NSoilLayers && NoEndLayer; i++) {
    Dz += SoilDepth[i];
	NLayers++;
	if (NoStartLayer) {
	  StartLayer = i;
	  NoStartLayer = FALSE;
	}
    else 
	  NoEndLayer = FALSE;
  }

  for (i = 0; i < NLayers; i++) {
    Sr = Moisture[i + StartLayer] / Porosity[i + StartLayer];

    /* Assume for now that either all the water is either frozen or unfrozen */

    /* frozen soil */
    if (TSoil[i + StartLayer] < 0) {
      Ke = Sr;
      KhSat = pow((double)KhSol[i + StartLayer], (double)(1-Porosity[i + StartLayer])) *
		  pow((double)KhICE, (double)Porosity[i + StartLayer]);
	}
    /* unfrozen soil */
    else {
      if (Sr > 0.1)
		Ke = log10((double) Sr) + 1.0;
      else
		Ke = 0.0;
	  KhSat = pow((double) KhSol[i + StartLayer], (double)(1-Porosity[i + StartLayer])) *
		  pow((double)KhH2O, (double)Porosity[i + StartLayer]);
    }
    LayerKh[i] = (KhSat-KhDry[i + StartLayer])*Ke + KhDry[i + StartLayer];
  }
}

/**********************************************************************
  This function calculates the soil volumetric heat capacity based 
  on the fractional volume of its component parts.

  Constant values are volumetric heat capacities in J/m^3/K
**********************************************************************/
void CalcHeatCapacity(int NSoilLayers, float Bottom, float *SoilDepth,
					  float *Porosity, float *Moisture, float *TSoil,
					  float *Ch)
{
  int i;			    /* counter */
  int NLayers;			/* Number of soil layers in depth interval */
  int StartLayer = 0;	/* First layer below top */
  char NoEndLayer;		/* flag to indicate whether an end layer has been determined */
  char NoStartLayer;	/* flag to indicate whether a start layer has been determined */
  float Dz;			    /* Depth from soil surface (m) */
  float TotalDepth;		/* Depth of soil column for which to  calculate the effective thermal 
				           conductivity (m) */

  TotalDepth = Bottom;

  Dz = 0.0;
  NLayers = 0;
  NoStartLayer = TRUE;
  NoEndLayer = TRUE;

  for (i = 0; i < NSoilLayers && NoEndLayer; i++) {
    Dz += SoilDepth[i];
	NLayers++;
	if (NoStartLayer) {
	  StartLayer = i;
	  NoStartLayer = FALSE;
	}
    else 
	  NoEndLayer = FALSE;
  }

  for (i = 0; i < NLayers; i++) {
    Ch[i] = 2.0e6 * Porosity[i + StartLayer];
    if (TSoil[i + StartLayer] >= 0) 
      Ch[i] += Moisture[i + StartLayer] * CH_WATER;
    else 
      Ch[i] += Moisture[i + StartLayer] * CH_WATER;
    }  
}
