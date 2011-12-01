/**
 *   @file netcdf_config.h
 *
 *   @brief  Parse the nssi configuration XML file.
 *
 *   @author Todd Kordenbrock (thkorde\@sandia.gov)
 *
 */

#ifndef _NETCDF_CONFIG_H_
#define _NETCDF_CONFIG_H_

#include "Trios_nssi_types.h"
#include "netcdf_args.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * @brief A structure to represent the configuration of
     * NSSI core services.
     */
    struct netcdf_config {

        /** @brief Number of available storage servers */
        int num_servers;

        /** @brief storage service IDs */
        char **netcdf_server_urls;

        /** @brief number of clients the server can expect */
        unsigned int num_participants;

        /** @brief The type of write operation the client wishes the server to perform */
        enum write_type write_type;

        /** @brief The number of bytes the server can cache for aggregation */
        size_t bytes_per_server;

        /** @brief The number of bytes the server can cache for aggregation */
        size_t use_subchunking;

    };

    /**
     * Options and arguments passed to the client driver.
     */
    struct netcdf_args {
            std::string server_url;
            std::string url_file;
            log_level verbose;
            std::string logfile;
            bool daemon_flag;
    };


#ifdef __cplusplus
}
#endif

#endif
