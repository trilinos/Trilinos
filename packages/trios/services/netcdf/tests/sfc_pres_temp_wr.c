/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This example writes some surface pressure and temperatures. It is
   intended to illustrate the use of the netCDF C API. The companion
   program sfc_pres_temp_rd.c shows how to read the netCDF data file
   created by this program.
   
   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
   
   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: sfc_pres_temp_wr.c,v 1.3 2006/06/13 20:46:16 ed Exp $
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

/* This is the name of the data file we will create. */
#define FILE_NAME "sfc_pres_temp.nc"

/* We are writing 2D data, a 6 x 12 lat-lon grid. We will need two
 * netCDF dimensions. */
#define NDIMS 2
#define NLAT 6
#define NLON 12
#define LAT_NAME "latitude"
#define LON_NAME "longitude"

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

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

int
main()
{
   int ncid, lon_dimid, lat_dimid, pres_varid, temp_varid;

/* In addition to the latitude and longitude dimensions, we will also
   create latitude and longitude netCDF variables which will hold the
   actual latitudes and longitudes. Since they hold data about the
   coordinate system, the netCDF term for these is: "coordinate
   variables." */
   int lat_varid, lon_varid;

   int dimids[NDIMS];

   /* We will write surface temperature and pressure fields. */
   float pres_out[NLAT][NLON];
   float temp_out[NLAT][NLON];
   float lats[NLAT], lons[NLON];

   /* It's good practice for each netCDF variable to carry a "units"
    * attribute. */
   char pres_units[] = "hPa";
   char temp_units[] = "celsius";

   /* Loop indexes. */
   int lat, lon;
   
   /* Error handling. */
   int retval;

   /* Create some pretend data. If this wasn't an example program, we
    * would have some real data to write, for example, model
    * output. */
   for (lat = 0; lat < NLAT; lat++)
      lats[lat] = START_LAT + 5.*lat;
   for (lon = 0; lon < NLON; lon++)
      lons[lon] = START_LON + 5.*lon;
   
   for (lat = 0; lat < NLAT; lat++)
      for (lon = 0; lon < NLON; lon++)
      {
	 pres_out[lat][lon] = SAMPLE_PRESSURE + (lon * NLAT + lat);
	 temp_out[lat][lon] = SAMPLE_TEMP + .25 * (lon * NLAT + lat);
      }

   /* Create the file. */
   if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. */
   if ((retval = nc_def_dim(ncid, LAT_NAME, NLAT, &lat_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, LON_NAME, NLON, &lon_dimid)))
      ERR(retval);

   /* Define coordinate netCDF variables. They will hold the
      coordinate information, that is, the latitudes and longitudes. A
      varid is returned for each.*/
   if ((retval = nc_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid, 
			    &lat_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid, 
			    &lon_varid)))
      ERR(retval);

   /* Define units attributes for coordinate vars. This attaches a
      text attribute to each of the coordinate variables, containing
      the units. Note that we are not writing a trailing NULL, just
      "units", because the reading program may be fortran which does
      not use null-terminated strings. In general it is up to the
      reading C program to ensure that it puts null-terminators on
      strings where necessary.*/
   if ((retval = nc_put_att_text(ncid, lat_varid, UNITS, 
				 strlen(DEGREES_NORTH), DEGREES_NORTH)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, lon_varid, UNITS, 
				 strlen(DEGREES_EAST), DEGREES_EAST)))
      ERR(retval);

   /* Define the netCDF variables. The dimids array is used to pass
      the dimids of the dimensions of the variables.*/
   dimids[0] = lat_dimid;
   dimids[1] = lon_dimid;
   if ((retval = nc_def_var(ncid, PRES_NAME, NC_FLOAT, NDIMS, 
			    dimids, &pres_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, TEMP_NAME, NC_FLOAT, NDIMS, 
			    dimids, &temp_varid)))
      ERR(retval);

   /* Define units attributes for vars. */
   if ((retval = nc_put_att_text(ncid, pres_varid, UNITS, 
				 strlen(pres_units), pres_units)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, temp_varid, UNITS, 
				 strlen(temp_units), temp_units)))
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

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays of data are the same size
      as the netCDF variables we have defined. */
   if ((retval = nc_put_var_float(ncid, pres_varid, &pres_out[0][0])))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid, temp_varid, &temp_out[0][0])))
      ERR(retval);

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** SUCCESS writing example file sfc_pres_temp.nc!\n");
   return 0;
}
