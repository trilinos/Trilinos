/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is an example program which writes some 4D pressure and
   temperatures. It is intended to illustrate the use of the netCDF
   C API. The companion program pres_temp_4D_rd.c shows how
   to read the netCDF data file created by this program.
   
   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
   
   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: pres_temp_4D_wr.c,v 1.4 2006/06/29 20:45:58 russ Exp $
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

/* This is the name of the data file we will create. */
#define FILE_NAME "pres_temp_4D.nc"

/* We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
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

/* These are used to construct some example data. */
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
   /* IDs for the netCDF file, dimensions, and variables. */
   int ncid, lon_dimid, lat_dimid, lvl_dimid, rec_dimid;
   int lat_varid, lon_varid, pres_varid, temp_varid;
   int dimids[NDIMS];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[NDIMS], count[NDIMS];

   /* Program variables to hold the data we will write out. We will only
      need enough space to hold one timestep of data; one record. */
   float pres_out[NLVL][NLAT][NLON];
   float temp_out[NLVL][NLAT][NLON];

   /* These program variables hold the latitudes and longitudes. */
   float lats[NLAT], lons[NLON];

   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   /* Create some pretend data. If this wasn't an example program, we
    * would have some real data to write, for example, model
    * output. */
   for (lat = 0; lat < NLAT; lat++)
      lats[lat] = START_LAT + 5.*lat;
   for (lon = 0; lon < NLON; lon++)
      lons[lon] = START_LON + 5.*lon;
   
   for (lvl = 0; lvl < NLVL; lvl++)
      for (lat = 0; lat < NLAT; lat++)
	 for (lon = 0; lon < NLON; lon++)
	 {
	    pres_out[lvl][lat][lon] = SAMPLE_PRESSURE + i;
	    temp_out[lvl][lat][lon] = SAMPLE_TEMP + i++;
	 }

   /* Create the file. */
   if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
   if ((retval = nc_def_dim(ncid, LVL_NAME, NLVL, &lvl_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, LAT_NAME, NLAT, &lat_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, LON_NAME, NLON, &lon_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
      ERR(retval);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&lat_dimid) and
      similarly for (&lon_dimid). */
   if ((retval = nc_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid, 
			    &lat_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid, 
			    &lon_varid)))
      ERR(retval);

   /* Assign units attributes to coordinate variables. */
   if ((retval = nc_put_att_text(ncid, lat_varid, UNITS, 
				 strlen(DEGREES_NORTH), DEGREES_NORTH)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, lon_varid, UNITS, 
				 strlen(DEGREES_EAST), DEGREES_EAST)))
      ERR(retval);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = rec_dimid;
   dimids[1] = lvl_dimid;
   dimids[2] = lat_dimid;
   dimids[3] = lon_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, PRES_NAME, NC_FLOAT, NDIMS, 
			    dimids, &pres_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, TEMP_NAME, NC_FLOAT, NDIMS, 
			    dimids, &temp_varid)))
      ERR(retval);

   /* Assign units attributes to the netCDF variables. */
   if ((retval = nc_put_att_text(ncid, pres_varid, UNITS, 
				 strlen(PRES_UNITS), PRES_UNITS)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, temp_varid, UNITS, 
				 strlen(TEMP_UNITS), TEMP_UNITS)))
      ERR(retval);

   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_float(ncid, lat_varid, &lats[0])))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid, lon_varid, &lons[0])))
      ERR(retval);

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = 1;
   count[1] = NLVL;
   count[2] = NLAT;
   count[3] = NLON;
   start[1] = 0;
   start[2] = 0;
   start[3] = 0;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
   for (rec = 0; rec < NREC; rec++)
   {
      start[0] = rec;
      if ((retval = nc_put_vara_float(ncid, pres_varid, start, count, 
				      &pres_out[0][0][0])))
	 ERR(retval);
      if ((retval = nc_put_vara_float(ncid, temp_varid, start, count, 
				      &temp_out[0][0][0])))
	 ERR(retval);
   }

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** SUCCESS writing example file %s!\n", FILE_NAME);
   return 0;
}
