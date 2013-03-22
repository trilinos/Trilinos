/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */

#include "Trios_config.h"
#ifdef HAVE_TRIOS_PNETCDF

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

#include "ezxml.h"
#include "Trios_logger.h"
#include "Trios_nssi_types.h"
#include "Trios_nssi_fprint_types.h"

#include "netcdf_config_parser.h"
#include "netcdf_debug.h"


/* debug level for the configuration parser */
log_level netcdf_config_debug_level = LOG_UNDEFINED;


/* The config parser parses an NSSI config file (xml)
 * and returns the nssi_config data structure.
 */
static int
parse_service(ezxml_t node, char *url)
{
    int rc = NSSI_OK;
    const char *attr;

    url[0]='\0';
    /* get the URL attribute */
    attr = ezxml_attr(node, "url");
    if (attr != NULL) {
        strncpy(url, attr, NNTI_URL_LEN);
    }

    log_debug(netcdf_config_debug_level, "url=%s", url);

    return rc;
}


static int
parse_write_type(ezxml_t node, write_type *write_type)
{
    int rc = NSSI_OK;
    const char *attr;

    attr = ezxml_attr(node, "default");
    if (!strcmp(attr, "WRITE_DIRECT")) {
        *write_type = WRITE_DIRECT;
        log_debug(netcdf_config_debug_level, "using %s", attr);
    } else if (!strcmp(attr, "WRITE_AGGREGATE_INDEPENDENT")) {
        *write_type = WRITE_AGGREGATE_INDEPENDENT;
        log_debug(netcdf_config_debug_level, "using %s", attr);
    } else if (!strcmp(attr, "WRITE_AGGREGATE_COLLECTIVE")) {
        *write_type = WRITE_AGGREGATE_COLLECTIVE;
        log_debug(netcdf_config_debug_level, "using %s", attr);
    } else if (!strcmp(attr, "WRITE_CACHING_INDEPENDENT")) {
        *write_type = WRITE_CACHING_INDEPENDENT;
        log_debug(netcdf_config_debug_level, "using %s", attr);
    } else if (!strcmp(attr, "WRITE_CACHING_COLLECTIVE")) {
        *write_type = WRITE_CACHING_COLLECTIVE;
        log_debug(netcdf_config_debug_level, "using %s", attr);
    }

    return rc;
}

static int
parse_num_participants(ezxml_t node, unsigned int *num_participants)
{
    int rc = NSSI_OK;
    const char *attr;

    attr = ezxml_attr(node, "default");
    sscanf(attr, "%u", num_participants);

    return rc;
}

static int
parse_bytes_per_server(ezxml_t node, size_t *bytes_per_server)
{
    int rc = NSSI_OK;
    const char *attr;

    attr = ezxml_attr(node, "default");
    sscanf(attr, "%ld", bytes_per_server);
    log_debug(netcdf_debug_level, "bytes_per_server %ld", *bytes_per_server);

    return rc;
}

static int
parse_use_subchunking(ezxml_t node, size_t *use_subchunking)
{
    int rc = NSSI_OK;
    const char *attr;

    attr = ezxml_attr(node, "default");
    sscanf(attr, "%ld", use_subchunking);
    log_debug(netcdf_debug_level, "use_subchunking %ld", *use_subchunking);

    return rc;
}

static int
parse_netcdf(
    ezxml_t node,
    const int index,
    struct netcdf_config *config)
{
    int rc = NSSI_OK;
    ezxml_t server_id, write_type, num_participants, bytes_per_server, use_subchunking;

    /* server list */
    server_id = ezxml_child(node, "server-id");
    if (server_id) {
        rc = parse_service(server_id, (config->netcdf_server_urls[index]));
        if (rc != NSSI_OK) {
            log_error(netcdf_config_debug_level,
                    "error parsing serverlist: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    write_type = ezxml_child(node, "write-type");
    if (write_type) {
        rc = parse_write_type(write_type, &config->write_type);
        if (rc != NSSI_OK) {
            log_error(netcdf_config_debug_level,
                    "error parsing write_type: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    num_participants = ezxml_child(node, "num-participants");
    if (num_participants) {
        rc = parse_num_participants(num_participants, &config->num_participants);
        if (rc != NSSI_OK) {
            log_error(netcdf_config_debug_level,
                    "error parsing num_participants: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    bytes_per_server = ezxml_child(node, "bytes-per-server");
    if (bytes_per_server) {
        rc = parse_bytes_per_server(bytes_per_server, &config->bytes_per_server);
        if (rc != NSSI_OK) {
            log_error(netcdf_config_debug_level,
                    "error parsing bytes_per_server: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    use_subchunking = ezxml_child(node, "use-subchunking");
    if (use_subchunking) {
        rc = parse_use_subchunking(use_subchunking, &config->use_subchunking);
        if (rc != NSSI_OK) {
            log_error(netcdf_config_debug_level,
                    "error parsing use_subchunking: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    return rc;
}



static int
parse_config(
    ezxml_t node,
    struct netcdf_config *config)
{
    int rc;
    int i;
    int count=0;
    ezxml_t service;

    /* process all the storage servers */

    /* first count the children */
    for (service = ezxml_child(node, "netcdf"); service; service = service->next)
    {
        count++;
    }

    config->num_servers = count;

    config->netcdf_server_urls = (char **)calloc(count, sizeof(char *));
    if (!config->netcdf_server_urls) {
        log_error(netcdf_config_debug_level, "could not allocate netcdf services");
        return NSSI_ENOMEM;
    }

    for (i=0;i<count;i++) {
        config->netcdf_server_urls[i] = (char *)calloc(NNTI_URL_LEN, sizeof(char));
        if (!config->netcdf_server_urls[i]) {
            log_error(netcdf_config_debug_level, "could not allocate netcdf services");
            return NSSI_ENOMEM;
        }
    }


    /* initialize count again */
    count = 0;

    /* parse the service list */
    for (service = ezxml_child(node, "netcdf"); service; service = service->next)
    {
        if (service) {
            rc = parse_netcdf(service, count++, config);
            if (rc != NSSI_OK) {
                log_error(netcdf_config_debug_level, "error parsing netcdf service: %s",
                        nssi_err_str(rc));
                return rc;
            }
        }
    }

    return rc;
}


/*------------------ EXTERNAL APIs ------------------ */

int
parse_netcdf_config_file(
    const char *docname,
    struct netcdf_config *netcdf_cfg)
{
    int rc = NSSI_OK;
    int fd = 0;
    ezxml_t doc, config;

    log_debug(netcdf_config_debug_level, "entered parse_config_file");

    if (!docname) {
        rc = NSSI_EINVAL;
        return rc;
    }

    memset(netcdf_cfg, 0, sizeof(struct netcdf_config));

    fd = open(docname, O_RDONLY, 0);
    if (fd < 0) {
        log_error(netcdf_config_debug_level, "failed to open config file (%s): %s", docname, strerror(errno));
        return NSSI_EIO;
    }
    doc = ezxml_parse_fd(fd);
    close(fd);
    if (doc == NULL) {
        log_error(netcdf_config_debug_level, "failed to parse config file (%s)", docname);
        return NSSI_EINVAL;
    }

    config = ezxml_child(doc, "config");
    if (config) {
        rc = parse_config(config, netcdf_cfg);
        if (rc != NSSI_OK) {
            log_error(netcdf_config_debug_level,
                    "could not parse config file: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    ezxml_free(doc);

    log_debug(netcdf_config_debug_level, "finished parse_config_file");

    return rc;
}


void netcdf_config_free(struct netcdf_config *netcdf_cfg)
{
    int i;

    /* release the space allocated for the ss_server_ids */
    for (i=0;i<netcdf_cfg->num_servers;i++) {
        free(netcdf_cfg->netcdf_server_urls[i]);
    }
    free(netcdf_cfg->netcdf_server_urls);
}

#endif // HAVE_TRIOS_PNETCDF
