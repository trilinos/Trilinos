/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is an example which reads some surface pressure and
   temperatures. The data file read by this program is produced
   companion program sfc_pres_temp_wr.c. It is intended to illustrate
   the use of the netCDF C API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: sfc_pres_temp_rd.c,v 1.3 2006/06/25 11:40:00 ed Exp $
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "sfc_pres_temp.nc"

/* We are reading 2D data, a 6 x 12 lat-lon grid. */
#define NDIMS 2
#define NLAT 6
#define NLON 12

#define LAT_NAME "latitude"
#define LON_NAME "longitude"
#define PRES_NAME "pressure"
#define TEMP_NAME "temperature"

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

   /* We will read surface temperature and pressure fields. */
   float pres_in[NLAT][NLON];
   float temp_in[NLAT][NLON];

   /* For the lat lon coordinate variables. */
   float lats_in[NLAT], lons_in[NLON];

   /* To check the units attributes. */
   char pres_units_in[MAX_ATT_LEN], temp_units_in[MAX_ATT_LEN];
   char lat_units_in[MAX_ATT_LEN], lon_units_in[MAX_ATT_LEN];

   /* We will learn about the data file and store results in these
      program variables. */
   int ndims_in, nvars_in, ngatts_in, unlimdimid_in;

   /* Loop indexes. */
   int lat, lon;

   /* Error handling. */
   int retval;

   /* Open the file. */
   if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
      ERR(retval);

   /* There are a number of inquiry functions in netCDF which can be
      used to learn about an unknown netCDF file. NC_INQ tells how
      many netCDF variables, dimensions, and global attributes are in
      the file; also the dimension id of the unlimited dimension, if
      there is one. */
   if ((retval = nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in, 
			&unlimdimid_in)))
      ERR(retval);

   /* In this case we know that there are 2 netCDF dimensions, 4
      netCDF variables, no global attributes, and no unlimited
      dimension. */
   if (ndims_in != 2 || nvars_in != 4 || ngatts_in != 0 || 
       unlimdimid_in != -1) return 2;

   /* Get the varids of the latitude and longitude coordinate
    * variables. */
   if ((retval = nc_inq_varid(ncid, LAT_NAME, &lat_varid)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, LON_NAME, &lon_varid)))
      ERR(retval);

   /* Read the coordinate variable data. */
   if ((retval = nc_get_var_float(ncid, lat_varid, &lats_in[0])))
      ERR(retval);
   if ((retval = nc_get_var_float(ncid, lon_varid, &lons_in[0])))
      ERR(retval);

   /* Check the coordinate variable data. */
   for (lat = 0; lat < NLAT; lat++)
      if (lats_in[lat] != START_LAT + 5.*lat)
	 return 2;
   for (lon = 0; lon < NLON; lon++)
      if (lons_in[lon] != START_LON + 5.*lon)
	 return 2;

   /* Get the varids of the pressure and temperature netCDF
    * variables. */
   if ((retval = nc_inq_varid(ncid, PRES_NAME, &pres_varid)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, TEMP_NAME, &temp_varid)))
      ERR(retval);

   /* Read the data. Since we know the contents of the file we know
    * that the data arrays in this program are the correct size to
    * hold all the data. */
   if ((retval = nc_get_var_float(ncid, pres_varid, &pres_in[0][0])))
      ERR(retval);
   if ((retval = nc_get_var_float(ncid, temp_varid, &temp_in[0][0])))
      ERR(retval);

   /* Check the data. */
   for (lat = 0; lat < NLAT; lat++)
      for (lon = 0; lon < NLON; lon++)
	 if (pres_in[lat][lon] != SAMPLE_PRESSURE + (lon * NLAT + lat) ||
	     temp_in[lat][lon] != SAMPLE_TEMP + .25 * (lon * NLAT + lat))
	    return 2;

   /* Each of the netCDF variables has a "units" attribute. Let's read
      them and check them. */
   if ((retval = nc_get_att_text(ncid, lat_varid, UNITS, lat_units_in)))
      ERR(retval);
   if (strncmp(lat_units_in, LAT_UNITS, strlen(LAT_UNITS))) 
      return 2;

   if ((retval = nc_get_att_text(ncid, lon_varid, UNITS, lon_units_in)))
      ERR(retval);
   if (strncmp(lon_units_in, LON_UNITS, strlen(LON_UNITS))) 
      return 2;

   if ((retval = nc_get_att_text(ncid, pres_varid, UNITS, pres_units_in)))
      ERR(retval);
   if (strncmp(pres_units_in, PRES_UNITS, strlen(PRES_UNITS)))
      return 2;

   if ((retval = nc_get_att_text(ncid, temp_varid, UNITS, temp_units_in)))
      ERR(retval);
   if (strncmp(temp_units_in, TEMP_UNITS, strlen(TEMP_UNITS))) return 2;

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

   printf("*** SUCCESS reading example file sfc_pres_temp.nc!\n");
   return 0;
}
