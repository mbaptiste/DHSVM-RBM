/*
 * SUMMARY:      globals.c - global constants for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    29-May-97 at 20:27:40
 * $Id: globals.c,v 1.4 2003/07/01 21:26:30 olivier Exp $
 */

float LAI_SNOW_MULTIPLIER;	/* multiplier to calculate the amount of 
				   available snow interception as a function of
				   LAI */
float LAI_WATER_MULTIPLIER;	/* multiplier to determine maximum interception 
				   storage as a function of LAI  */
float LIQUID_WATER_CAPACITY;	/* water holding capacity of snow as a fraction
				   of snow-water-equivalent */
float MAX_SNOW_TEMP;		/* maximum temperature at which snow can 
				   occur (C) */
float MIN_INTERCEPTION_STORAGE;	/* the amount of snow on the canopy that can 
				   only be melted off. (m) */
float MIN_RAIN_TEMP;		/* minimum temperature at which rain can 
				   occur (C) */
int NWINDMAPS;			/* Number of wind maps in case the wind source
				   is model */
unsigned char OUTSIDEBASIN;	/* Mask value indicating outside the basin */
float PRECIPLAPSE;		/* Precipitation lapse rate in m/timestep / m */
float PRECIPMULTIPLIER;
float MINELEV;          /* Smallest elevation of all grid cells (m) */
float TEMPLAPSE;		/* Temperature lapse rate in C/m */
float Z0_GROUND;		/* Roughness length for bare soil (m) */
float Z0_SNOW;			/* Roughness length for snow (m) */
float Zref;			    /* Reference height (m) */
float MASSITER;         /* Maximum number of iterations for mass wasting*/
float DEBRISd50;
float DEBRISd90;
float CHANNELd50;
float CHANNELd90;