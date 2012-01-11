/**
 *   @file config_parser.h
 *
 *   @brief  Parse the nssi configuration XML file.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 1.23 $
 *   $Date: 2005/11/09 20:15:51 $
 *
 */

#ifndef _NETCDF_CONFIG_PARSER_H_
#define _NETCDF_CONFIG_PARSER_H_

#include "Trios_nssi_types.h"
#include "netcdf_args.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern log_level config_debug_level;

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



#if defined(__STDC__) || defined(__cplusplus)

    extern int parse_netcdf_config_file(const char *fname,
            struct netcdf_config *config);

    extern void netcdf_config_free(
            struct netcdf_config *config);

#endif

#ifdef __cplusplus
}
#endif

#endif
