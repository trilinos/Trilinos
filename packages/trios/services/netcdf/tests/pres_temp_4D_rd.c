/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is an example which reads some 4D pressure and
   temperatures. The data file read by this program is produced by the
   companion program pres_temp_4D_wr.c. It is intended to illustrate
   the use of the netCDF C API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: pres_temp_4D_rd.c,v 1.5 2006/06/26 20:37:31 russ Exp $
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "pres_temp_4D.nc"

/* We are reading 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
   timesteps of data. */
#define NDIMS 4
#define NLAT 6
#define NLON 12
#define LAT_NAME "latitude"
#define LON_NAME "longitude"
#define NREC 2
#define REC_NAME "time"
#define LVL_NAME "level"
#define NLVL 2

/* Names of things. */
#define PRES_NAME "pressure"
#define TEMP_NAME "temperature"
#define UNITS "units"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"

/* These are used to calculate the values we expect to find. */
#define SAMPLE_PRESSURE 900
#define SAMPLE_TEMP 9.0
#define START_LAT 25.0
#define START_LON -125.0

/* For the units attributes. */
#define UNITS "units"
#define PRES_UNITS "hPa"
#define TEMP_UNITS "celsius"
#define LAT_UNITS "degrees_north"
#define LON_UNITS "degrees_east"
#define MAX_ATT_LEN 80

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

int
main()
{
   int ncid, pres_varid, temp_varid;
   int lat_varid, lon_varid;

   /* The start and count arrays will tell the netCDF library where to
      read our data. */
   size_t start[NDIMS], count[NDIMS];

   /* Program variables to hold the data we will read. We will only
      need enough space to hold one timestep of data; one record. */
   float pres_in[NLVL][NLAT][NLON];
   float temp_in[NLVL][NLAT][NLON];

   /* These program variables hold the latitudes and longitudes. */
   float lats[NLAT], lons[NLON];

   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   /* Open the file. */
   if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
      ERR(retval);

   /* Get the varids of the latitude and longitude coordinate
    * variables. */
   if ((retval = nc_inq_varid(ncid, LAT_NAME, &lat_varid)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, LON_NAME, &lon_varid)))
      ERR(retval);

   /* Read the coordinate variable data. */
   if ((retval = nc_get_var_float(ncid, lat_varid, &lats[0])))
      ERR(retval);
   if ((retval = nc_get_var_float(ncid, lon_varid, &lons[0])))
      ERR(retval);

   /* Check the coordinate variable data. */
   for (lat = 0; lat < NLAT; lat++)
      if (lats[lat] != START_LAT + 5.*lat)
	 return 2;
   for (lon = 0; lon < NLON; lon++)
      if (lons[lon] != START_LON + 5.*lon)
	 return 2;

   /* Get the varids of the pressure and temperature netCDF
    * variables. */
   if ((retval = nc_inq_varid(ncid, PRES_NAME, &pres_varid)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, TEMP_NAME, &temp_varid)))
      ERR(retval);

   /* Read the data. Since we know the contents of the file we know
    * that the data arrays in this program are the correct size to
    * hold one timestep. */
   count[0] = 1;
   count[1] = NLVL;
   count[2] = NLAT;
   count[3] = NLON;
   start[1] = 0;
   start[2] = 0;
   start[3] = 0;

   /* Read and check one record at a time. */
   for (rec = 0; rec < NREC; rec++)
   {
      start[0] = rec;
      if ((retval = nc_get_vara_float(ncid, pres_varid, start, 
				      count, &pres_in[0][0][0])))
	 ERR(retval);
      if ((retval = nc_get_vara_float(ncid, temp_varid, start,
				      count, &temp_in[0][0][0])))
	 ERR(retval);

      /* Check the data. */
      i = 0;
      for (lvl = 0; lvl < NLVL; lvl++)
	 for (lat = 0; lat < NLAT; lat++)
	    for (lon = 0; lon < NLON; lon++)
	    {
	       if (pres_in[lvl][lat][lon] != SAMPLE_PRESSURE + i) 
		  return 2;
	       if (temp_in[lvl][lat][lon] != SAMPLE_TEMP + i) 
		  return 2;
	       i++;
	    }

   } /* next record */

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

   printf("*** SUCCESS reading example file pres_temp_4D.nc!\n");
   return 0;
}
