#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <mpi.h>

void
check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
           (void) fprintf(stderr, "line %d of %s: %s\n", line, file, ncmpi_strerror(stat));
        exit(1);
    }
}

int
main(int argc, char **argv) {			/* create foo.nc */

   int  stat;			/* return status */
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int lat_dim;
   int lon_dim;
   int time_dim;

   /* dimension lengths */
   size_t lat_len = 10;
   size_t lon_len = 5;
   size_t time_len = NC_UNLIMITED;

   /* variable ids */
   int lat_id;
   int lon_id;
   int time_id;
   int z_id;
   int t_id;
   int p_id;
   int rh_id;

   /* rank (number of dimensions) for each variable */
#  define RANK_lat 1
#  define RANK_lon 1
#  define RANK_time 1
#  define RANK_z 3
#  define RANK_t 3
#  define RANK_p 3
#  define RANK_rh 3

   /* variable shapes */
   int lat_dims[RANK_lat];
   int lon_dims[RANK_lon];
   int time_dims[RANK_time];
   int z_dims[RANK_z];
   int t_dims[RANK_t];
   int p_dims[RANK_p];
   int rh_dims[RANK_rh];

   /* attribute vectors */
   double z_valid_range[2];
   double p__FillValue[1];
   int rh__FillValue[1];

  int stat=0;
   MPI_Init(&argc, &argv);
   /* enter define mode */
   stat = ncmpi_create(MPI_COMM_WORLD, "foo.nc", NC_CLOBBER, MPI_INFO_NULL, &ncid);
   check_err(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = ncmpi_def_dim(ncid, "lat", lat_len, &lat_dim);
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_def_dim(ncid, "lon", lon_len, &lon_dim);
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_def_dim(ncid, "time", time_len, &time_dim);
   check_err(stat,__LINE__,__FILE__);

   /* define variables */

   lat_dims[0] = lat_dim;
   stat = ncmpi_def_var(ncid, "lat", NC_INT, RANK_lat, lat_dims, &lat_id);
   check_err(stat,__LINE__,__FILE__);

   lon_dims[0] = lon_dim;
   stat = ncmpi_def_var(ncid, "lon", NC_INT, RANK_lon, lon_dims, &lon_id);
   check_err(stat,__LINE__,__FILE__);

   time_dims[0] = time_dim;
   stat = ncmpi_def_var(ncid, "time", NC_INT, RANK_time, time_dims, &time_id);
   check_err(stat,__LINE__,__FILE__);

   z_dims[0] = time_dim;
   z_dims[1] = lat_dim;
   z_dims[2] = lon_dim;
   stat = ncmpi_def_var(ncid, "z", NC_FLOAT, RANK_z, z_dims, &z_id);
   check_err(stat,__LINE__,__FILE__);

   t_dims[0] = time_dim;
   t_dims[1] = lat_dim;
   t_dims[2] = lon_dim;
   stat = ncmpi_def_var(ncid, "t", NC_FLOAT, RANK_t, t_dims, &t_id);
   check_err(stat,__LINE__,__FILE__);

   p_dims[0] = time_dim;
   p_dims[1] = lat_dim;
   p_dims[2] = lon_dim;
   stat = ncmpi_def_var(ncid, "p", NC_DOUBLE, RANK_p, p_dims, &p_id);
   check_err(stat,__LINE__,__FILE__);

   rh_dims[0] = time_dim;
   rh_dims[1] = lat_dim;
   rh_dims[2] = lon_dim;
   stat = ncmpi_def_var(ncid, "rh", NC_INT, RANK_rh, rh_dims, &rh_id);
   check_err(stat,__LINE__,__FILE__);

   /* assign attributes */
   stat = ncmpi_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_put_att_text(ncid, time_id, "units", 7, "seconds");
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_put_att_text(ncid, z_id, "units", 6, "meters");
   check_err(stat,__LINE__,__FILE__);
   z_valid_range[0] = 0;
   z_valid_range[1] = 5000;
   stat = ncmpi_put_att_double(ncid, z_id, "valid_range", NC_DOUBLE, 2, z_valid_range);
   check_err(stat,__LINE__,__FILE__);
   p__FillValue[0] = -9999;
   stat = ncmpi_put_att_double(ncid, p_id, "_FillValue", NC_DOUBLE, 1, p__FillValue);
   check_err(stat,__LINE__,__FILE__);
   rh__FillValue[0] = -1;
   stat = ncmpi_put_att_int(ncid, rh_id, "_FillValue", NC_INT, 1, rh__FillValue);
   check_err(stat,__LINE__,__FILE__);

   /* leave define mode */
   stat = ncmpi_enddef (ncid);
   check_err(stat,__LINE__,__FILE__);

   {			/* store lat */
    static int lat[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    ncmpi_begin_indep_data(ncid);
    stat = ncmpi_put_var_int(ncid, lat_id, lat);
    ncmpi_end_indep_data(ncid);
    check_err(stat,__LINE__,__FILE__);
   }

   {			/* store lon */
    static int lon[] = {-140, -118, -96, -84, -52};
    ncmpi_begin_indep_data(ncid);
    stat = ncmpi_put_var_int(ncid, lon_id, lon);
    ncmpi_end_indep_data(ncid);
    check_err(stat,__LINE__,__FILE__);
   }
   stat = ncmpi_close(ncid);
   check_err(stat,__LINE__,__FILE__);
   MPI_Finalize();
   return 0;
}
