/*
 * nnti_utils.h
 *
 *  Created on: Feb 16, 2011
 *      Author: thkorde
 */

#ifndef NNTI_UTILS_H_
#define NNTI_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <Trios_nnti_xdr.h>

NNTI_result_t nnti_url_get_transport(const char *url, char *outstr, const int maxlen);
NNTI_result_t nnti_url_get_address(const char *url, char *outstr, const int maxlen);
NNTI_result_t nnti_url_get_memdesc(const char *url, char *outstr, const int maxlen);
NNTI_result_t nnti_url_get_params(const char *url, char *outstr, const int maxlen);


NNTI_result_t nnti_sleep(const uint64_t msec);

#ifdef __cplusplus
}
#endif

#endif /* NNTI_UTILS_H_ */
