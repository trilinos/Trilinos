/*
 * nnti_utils.c
 *
 *  Created on: Feb 16, 2011
 *      Author: thkorde
 */


#include <string.h>
#include <errno.h>

#include "nnti_utils.h"

#include "nnti_internal.h"
#include "Trios_logger.h"

NNTI_result_t nnti_url_get_transport(const char *url, char *outstr, const int maxlen)
{
    int translen=0;
    char *sep=NULL;

    sep=strstr(url, "://");
    if (sep == NULL) {
        /* invalid URL */
        return(NNTI_EINVAL);
    } else {
        translen = sep-url;
        if (translen >= maxlen) {
            /*  */
            memcpy(outstr, url, maxlen-1);
            outstr[maxlen-1]='\0';
        } else {
            memcpy(outstr, url, translen);
            outstr[translen]='\0';
        }
    }

    return(NNTI_OK);
}

NNTI_result_t nnti_url_get_address(const char *url, char *outstr, const int maxlen)
{
    int addrlen=0;
    char *sep=NULL;
    char *address=NULL;

    sep=strstr(url, "://");
    if (sep == NULL) {
        /* invalid URL */
        return(NNTI_EINVAL);
    } else {
        address=sep+3;
        sep=strchr(address, '/');
        if (sep == NULL) {
            /* no trailing slash */
            addrlen=strlen(address);
        } else {
            addrlen=sep-address;
        }
    }
    if (addrlen == 0) {
        /* invalid URL */
        return(NNTI_EINVAL);
    }

    if (addrlen >= maxlen) {
        /*  */
        memcpy(outstr, address, maxlen-1);
        outstr[maxlen-1]='\0';
    } else {
        memcpy(outstr, address, addrlen);
        outstr[addrlen]='\0';
    }

    return(NNTI_OK);
}

NNTI_result_t nnti_url_get_memdesc(const char *url, char *outstr, const int maxlen)
{
    int mdlen=0;
    char *sep=NULL;
    char *memdesc=NULL;

    sep=strstr(url, "://");
    if (sep == NULL) {
        /* invalid URL */
        return(NNTI_EINVAL);
    } else {
        char *address=sep+3;
        sep=strchr(address, '/');
        if (sep == NULL) {
            /* no trailing slash, so no memory descriptor */
            outstr[0]='\0';
            return(NNTI_OK);
        } else {
            memdesc=sep+1;
            mdlen=strlen(memdesc);
        }
    }
    if (mdlen == 0) {
        /* nothing after trailing slash, so no memory descriptor */
        outstr[0]='\0';
        return(NNTI_OK);
    }

    if (mdlen >= maxlen) {
        /*  */
        memcpy(outstr, memdesc, maxlen-1);
        outstr[maxlen-1]='\0';
    } else {
        memcpy(outstr, memdesc, mdlen);
        outstr[mdlen]='\0';
    }

    return(NNTI_OK);
}

NNTI_result_t nnti_url_get_params(const char *url, char *outstr, const int maxlen)
{
    int plen=0;
    char *sep=NULL;
    char *params=NULL;

    sep=strstr(url, "://");
    if (sep == NULL) {
        /* invalid URL */
        return(NNTI_EINVAL);
    } else {
        char *address=sep+3;
        sep=strchr(address, '?');
        if (sep == NULL) {
            /* no question mark, so no parameters */
            outstr[0]='\0';
            return(NNTI_OK);
        } else {
            params=sep+1;
            plen=strlen(params);
        }
    }
    if (plen == 0) {
        /* nothing after question mark, so no params */
        outstr[0]='\0';
        return(NNTI_OK);
    }

    if (plen >= maxlen) {
        /*  */
        memcpy(outstr, params, maxlen-1);
        outstr[maxlen-1]='\0';
    } else {
        memcpy(outstr, params, plen);
        outstr[plen]='\0';
    }

    return(NNTI_OK);
}

NNTI_result_t nnti_sleep(const uint64_t msec)
{
    int rc=0;
    struct timespec ts, rmtp;

    ts.tv_sec=0;
    if (msec < 1000) {
        ts.tv_nsec=msec*1000*1000; /* 1sec == 1ns*1000*1000*1000; */
    } else {
        uint64_t sec=msec/1000;
        ts.tv_sec=sec;
        ts.tv_nsec=(msec-(sec*1000))*1000*1000;
    }

    rc=nanosleep(&ts, &rmtp);
    if (rc!=0) log_error(nnti_debug_level, "nanosleep failed: %s\n", strerror(errno));

    return(rc);
}
