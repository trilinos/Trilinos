/*
 * netcdf_client_private.h
 *
 *  Created on: Mar 16, 2009
 *      Author: thkorde
 */

#ifndef NETCDF_CLIENT_PRIVATE_H_
#define NETCDF_CLIENT_PRIVATE_H_


#ifdef __cplusplus
extern "C" {
#endif

    extern int nc_begin_indep_data(int ncid);
    extern int nc_end_indep_data(int ncid);
    extern int nc_set_file_state(int ncid);
    extern int nc_sync_wait(int ncid);
    extern int nc_close_wait(int ncid);

    extern int netcdf_client_fini(void);

#ifdef __cplusplus
}
#endif


#endif /* NETCDF_CLIENT_PRIVATE_H_ */
