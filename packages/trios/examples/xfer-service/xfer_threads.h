/*
 * xfer_threads.h
 *
 *  Created on: Aug 20, 2012
 *      Author: raoldfi
 */

#ifndef XFER_THREADS_H_
#define XFER_THREADS_H_

#include "Trios_nssi_server.h"

int xfer_start_server_threads(const int num_threads, const int max_reqs);
int xfer_enqueue_rpc_request(nssi_svc_rpc_request *req);
int xfer_cancel_server_threads();


#endif /* XFER_THREADS_H_ */
