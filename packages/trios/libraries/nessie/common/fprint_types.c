/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/**
 *   @file fprint_types.c
 *
 *   @brief Pretty print the different NSSI types.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 1640 $.
 *   $Date: 2007-11-28 11:59:53 -0700 (Wed, 28 Nov 2007) $.
 *
 */

#include <stdio.h>
#include <string.h> /* find memcpy() */

#include "Trios_logger.h"
#include "Trios_nssi_types.h"
#include "Trios_nssi_fprint_types.h"

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include <unistd.h>



/* ----------- Utility functions to output types ---- */


/**
 * @brief Print out a return code.
 *
 * @param fp  The output file.
 * @param name The name of the variable.
 * @param prefix Text that precedes the variable.
 * @param The return code.
 */
void fprint_nssi_return_code(
                FILE *fp,
                const char *name,
                const char *prefix,
                const int rc)
{
        logger_mutex_lock();
        fprintf(fp, "%s %s = %s\n", prefix, name, nssi_err_str(rc));
        logger_mutex_unlock();
}

static const char *myitoa(int val) {
        static char buf[32];

        snprintf(buf, sizeof(buf), "%d", val);
        return buf;
}

/**
 * @brief Output a string associated with a return code.
 *
 * @param rc @input the return code.
 *
 * @returns A string associated with the return code.
 */
const char *nssi_err_str(int rc)
{
    switch (rc) {

    case NSSI_OK:
        return "NSSI_OK";

    /** @brief Unspecified I/O error. */
    case NSSI_EIO:
        return "NSSI_EIO";

    /** @brief The size of the message is larger than the supported maximum size. */
    case NSSI_EMSGSIZE:
        return "NSSI_EMSGSIZE";

    /** @brief The operation or process has been canceled. */
    case NSSI_ECANCELED:
        return "NSSI_ECANCELED";

    /** @brief An operation timed out. */
    case NSSI_ETIMEDOUT:
        return "NSSI_ETIMEDOUT";

    /** @brief The value of the variable or parameter is invalid. */
    case NSSI_EINVAL:
        return "NSSI_EINVAL";

    /** @brief  No memory is available. */
    case NSSI_ENOMEM:
        return "NSSI_ENOMEM";

    /** @brief No such entry. */
    case NSSI_ENOENT:
        return "NSSI_ENOENT";

    /** @brief Unsupported operation. */
    case NSSI_ENOTSUP:
        return "NSSI_ENOTSUP";

    /** @brief The item already exists. */
    case NSSI_EEXIST:
        return "NSSI_EEXIST";

    /** @brief Unsuccessful RPC operation. */
    case NSSI_EBADRPC:
        return "NSSI_EBADRPC";

    /** @brief Not initialized. */
    case NSSI_ENOTINIT:
        return "NSSI_ENOTINIT";

    /** @brief Insufficient priveleges to perform operation. */
    case NSSI_EPERM:
        return "NSSI_EPERM";

    /** @brief Error decoding an RPC request. */
    case NSSI_EDECODE:
        return "NSSI_EDECODE";

    /** @brief Error encoding an RPC request. */
    case NSSI_EENCODE:
        return "NSSI_EENCODE";

    default:
        return myitoa(rc);
    }
}


/**
 * @brief Output a string associated with a return code.
 *
 * @param rc @input the return code.
 *
 * @returns A string associated with the return code.
 */
const char *nnti_err_str(int rc)
{
    switch (rc) {

    case NNTI_OK:
        return "NNTI_OK";

    /** @brief Unspecified I/O error. */
    case NNTI_EIO:
        return "NNTI_EIO";

    /** @brief The size of the message is larger than the supported maximum size. */
    case NNTI_EMSGSIZE:
        return "NNTI_EMSGSIZE";

    /** @brief The operation or process has been canceled. */
    case NNTI_ECANCELED:
        return "NNTI_ECANCELED";

    /** @brief An operation timed out. */
    case NNTI_ETIMEDOUT:
        return "NNTI_ETIMEDOUT";

    /** @brief The value of the variable or parameter is invalid. */
    case NNTI_EINVAL:
        return "NNTI_EINVAL";

    /** @brief  No memory is available. */
    case NNTI_ENOMEM:
        return "NNTI_ENOMEM";

    /** @brief No such entry. */
    case NNTI_ENOENT:
        return "NNTI_ENOENT";

    /** @brief Unsupported operation. */
    case NNTI_ENOTSUP:
        return "NNTI_ENOTSUP";

    /** @brief The item already exists. */
    case NNTI_EEXIST:
        return "NNTI_EEXIST";

    /** @brief Unsuccessful RPC operation. */
    case NNTI_EBADRPC:
        return "NNTI_EBADRPC";

    /** @brief Not initialized. */
    case NNTI_ENOTINIT:
        return "NNTI_ENOTINIT";

    /** @brief Insufficient priveleges to perform operation. */
    case NNTI_EPERM:
        return "NNTI_EPERM";

    default:
        return myitoa(rc);
    }
}


/**
 * @brief Output the contents of a Portals remote memory descriptor.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param addr    The remote memory address.
 */
void fprint_NNTI_remote_addr(
                FILE *fp,
                const char *name,
                const char *prefix,
                const NNTI_remote_addr_t *addr)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s   ", prefix);

        logger_mutex_lock();
        if (addr == NULL) {
                fprintf(fp, "%s %s = NULL\n", prefix, name);
                return;
        }

        /* header */
        fprintf(fp, "%s %s = {\n", prefix, name);

        /* contents */
        fprintf(fp, "%s    transport_id = %lld,\n", subprefix, (long)addr->transport_id);
        switch (addr->transport_id) {
            case NSSI_RPC_PTL:
                fprintf(fp, "%s    buffer_id  = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.portals.buffer_id);
                fprintf(fp, "%s    match_bits = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.portals.match_bits);
                fprintf(fp, "%s    size       = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.portals.size);
                break;
            case NSSI_RPC_IB:
                fprintf(fp, "%s    buf      = %p,\n",   subprefix, addr->NNTI_remote_addr_t_u.ib.buf);
                fprintf(fp, "%s    key      = %x,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.ib.key);
                fprintf(fp, "%s    size     = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.ib.size);
                fprintf(fp, "%s    ack_buf  = %p,\n",   subprefix, addr->NNTI_remote_addr_t_u.ib.ack_buf);
                fprintf(fp, "%s    ack_key  = %x,\n",   subprefix, addr->NNTI_remote_addr_t_u.ib.ack_key);
                fprintf(fp, "%s    ack_size = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.ib.ack_size);
                break;
            case NSSI_RPC_LUC:
                fprintf(fp, "%s    buf  = %p,\n", subprefix, addr->NNTI_remote_addr_t_u.luc.buf);
                fprintf(fp, "%s    size = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.luc.size);
                break;
            case NSSI_RPC_GEMINI:
                fprintf(fp, "%s    type       = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.gni.type);
                fprintf(fp, "%s    buf        = %p,\n", subprefix, addr->NNTI_remote_addr_t_u.gni.buf);
                fprintf(fp, "%s    mem_hdl    = (%llu, %llu),\n",
                        subprefix, addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword1, addr->NNTI_remote_addr_t_u.gni.mem_hdl.qword2);
                fprintf(fp, "%s    size       = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.gni.size);
                fprintf(fp, "%s    wc_addr    = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.gni.wc_addr);
                fprintf(fp, "%s    wc_mem_hdl = (%llu, %llu),\n",
                        subprefix, addr->NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1, addr->NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2);
                break;
            case NSSI_RPC_MPI:
                fprintf(fp, "%s    rtr_tag    = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.mpi.rtr_tag);
                fprintf(fp, "%s    rts_tag    = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.mpi.rts_tag);
                fprintf(fp, "%s    data_tag   = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.mpi.data_tag);
                fprintf(fp, "%s    size       = %llu,\n", subprefix, (unsigned long long)addr->NNTI_remote_addr_t_u.mpi.size);
                break;
        }


        /* footer */
        fprintf(fp, "%s }\n", subprefix);
        logger_mutex_unlock();
}

/**
 * @brief Output the contents of a Portals remote memory descriptor.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param addr    The remote memory address.
 */
void fprint_NNTI_peer(
                FILE *fp,
                const char *name,
                const char *prefix,
                const NNTI_peer_t *addr)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s   ", prefix);

        logger_mutex_lock();
        if (addr == NULL) {
                fprintf(fp, "%s %s = NULL\n", prefix, name);
                return;
        }

        /* header */
        fprintf(fp, "%s %s = {\n", prefix, name);

        /* contents */
        fprintf(fp, "%s    url = %s,\n", subprefix, addr->url);
        fprintf(fp, "%s    transport_id = %lld,\n", subprefix, (long)addr->peer.transport_id);
        switch (addr->peer.transport_id) {
            case NSSI_RPC_PTL:
                fprintf(fp, "%s    nid = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.portals.nid);
                fprintf(fp, "%s    pid = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.portals.pid);
                break;
            case NSSI_RPC_IB:
                fprintf(fp, "%s    addr = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.ib.addr);
                fprintf(fp, "%s    port = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.ib.port);
                fprintf(fp, "%s    qpn  = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.ib.qpn);
                break;
            case NSSI_RPC_LUC:
                fprintf(fp, "%s    nid = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.portals.nid);
                fprintf(fp, "%s    pid = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.portals.pid);
                break;
            case NSSI_RPC_GEMINI:
                fprintf(fp, "%s    addr    = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.gni.addr);
                fprintf(fp, "%s    port    = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.gni.port);
                fprintf(fp, "%s    inst_id = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.gni.inst_id);
                break;
            case NSSI_RPC_MPI:
                fprintf(fp, "%s    rank    = %llu,\n", subprefix, (unsigned long long)addr->peer.NNTI_remote_process_t_u.mpi.rank);
                break;
        }

        /* footer */
        fprintf(fp, "%s }\n", subprefix);
        logger_mutex_unlock();
}

/**
 * @brief Output the contents of a Portals remote memory descriptor.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param addr    The remote memory address.
 */
void fprint_NNTI_buffer(
                FILE *fp,
                const char *name,
                const char *prefix,
                const NNTI_buffer_t *addr)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s   ", prefix);

        logger_mutex_lock();
        if (addr == NULL) {
                fprintf(fp, "%s %s = NULL\n", prefix, name);
                return;
        }

        /* header */
        fprintf(fp, "%s %s (%p) = {\n", prefix, name, addr);

        /* contents */
        fprintf(fp, "%s    transport_id = %ld,\n", subprefix, (long)addr->transport_id);
        fprint_NNTI_peer(fp, "buffer_owner", subprefix, &addr->buffer_owner);
        fprint_NNTI_remote_addr(fp, "buffer_addr", subprefix, &addr->buffer_addr);
        fprintf(fp, "%s    ops = %ld,\n", subprefix, (long)addr->ops);

        fprintf(fp, "%s    payload_size = %llu,\n", subprefix, (unsigned long long)addr->payload_size);
        fprintf(fp, "%s    payload = %p,\n", subprefix, addr->payload);
        fprintf(fp, "%s    transport_private = %p,\n", subprefix, addr->transport_private);

        /* footer */
        fprintf(fp, "%s }\n", subprefix);
        logger_mutex_unlock();
}

/**
 * @brief Output the contents of a Portals remote memory descriptor.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param addr    The remote memory address.
 */
void fprint_NNTI_status(
                FILE *fp,
                const char *name,
                const char *prefix,
                const NNTI_status_t *status)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s   ", prefix);

        logger_mutex_lock();
        if (status == NULL) {
                fprintf(fp, "%s %s = NULL\n", prefix, name);
                return;
        }

        /* header */
        fprintf(fp, "%s %s = {\n", prefix, name);

        /* contents */
        fprintf(fp, "%s    op     = %llu,\n", subprefix, (unsigned long long)status->op);
        fprintf(fp, "%s    result = %llu,\n", subprefix, (unsigned long long)status->result);
        fprintf(fp, "%s    start  = %p,\n",   subprefix, (void *)status->start);
        fprintf(fp, "%s    offset = %llu,\n", subprefix, (unsigned long long)status->offset);
        fprintf(fp, "%s    length = %llu,\n", subprefix, (unsigned long long)status->length);

        fprint_NNTI_peer(fp, "src", subprefix, &status->src);
        fprint_NNTI_peer(fp, "dest", subprefix, &status->dest);

        /* footer */
        fprintf(fp, "%s }\n", subprefix);
        logger_mutex_unlock();
}

/**
 * @brief Output the contents of a request header.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param hdr     The request header.
 */
void fprint_nssi_request_header(
                FILE *fp,
                const char *name,
                const char *prefix,
                const nssi_request_header *hdr)
{
        char subprefix[100];
        snprintf(subprefix, 100, "%s   ", prefix);

        logger_mutex_lock();
        /* header */
        fprintf(fp, "%s %s = {\n", prefix, name);

        /* contents */
        fprintf(fp, "%s    id = %lu,\n", subprefix, (long unsigned int)hdr->id);
        fprintf(fp, "%s    opcode = %u,\n", subprefix, hdr->opcode);
        fprintf(fp, "%s    fetch_args = %d,\n", subprefix, hdr->fetch_args);
        fprint_NNTI_buffer(fp, "args_addr", subprefix, &(hdr->args_addr));
        fprint_NNTI_buffer(fp, "data_addr", subprefix, &hdr->data_addr);
        fprint_NNTI_buffer(fp, "res_addr", subprefix, &(hdr->res_addr));

        /* footer */
        fprintf(fp, "%s }\n", subprefix);
        logger_mutex_unlock();
}

/**
 * @brief Output the contents of a result header.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param hdr     The result header.
 */
void fprint_nssi_result_header(
                FILE *fp,
                const char *name,
                const char *prefix,
                const nssi_result_header *hdr)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s   ", prefix);

        logger_mutex_lock();
        /* header */
        fprintf(fp, "%s %s = {\n", prefix, name);

        /* contents */
        fprintf(fp, "%s    id = %lu,\n", subprefix, (long unsigned int)hdr->id);
        fprintf(fp, "%s    opcode = %lu,\n", subprefix, (long unsigned int)hdr->opcode);
        fprintf(fp, "%s    fetch_result = %d,\n", subprefix, hdr->fetch_result);
        fprintf(fp, "%s    result_size = %d,\n", subprefix, hdr->result_size);
        fprint_NNTI_buffer(fp, "res_addr", subprefix, &hdr->result_addr);
        fprintf(fp, "%s    rc = %d,\n", subprefix, hdr->rc);

        /* footer */
        fprintf(fp, "%s }\n", subprefix);
        logger_mutex_unlock();
}



/**
 * @brief Print out nssi_rpc_encode.
 *
 * @param fp  The output file.
 * @param name The name of the variable.
 * @param prefix Text that precedes the variable.
 * @param The encoding.
 */
void fprint_nssi_rpc_encode(
                FILE *fp,
                const char *name,
                const char *prefix,
                const nssi_rpc_encode *rpc_encode)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s", prefix);

        /* contents */
        switch (*rpc_encode) {
                case NSSI_RPC_XDR:
                        fprintf(fp, "%s    %s = NSSI_RPC_XDR,\n", subprefix, name);
                        break;

                default:
                        fprintf(fp, "%s    %s = UNDEFINED,\n", subprefix, name);
                        break;
        }
}




/**
 * @brief Print out an nssi service descriptor.
 *
 * @param fp  The output file.
 * @param name The name of the variable.
 * @param prefix Text that precedes the variable.
 * @param The service.
 */
void fprint_nssi_service(
                FILE *fp,
                const char *name,
                const char *prefix,
                const nssi_service *svc)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s\t", prefix);

        logger_mutex_lock();
        /* header */
        fprintf(fp, "%s %s = {\n", prefix, name);

        /* contents */
        fprint_nssi_rpc_encode(fp, "rpc_encode", subprefix, &(svc->rpc_encode));
        fprint_NNTI_peer(fp, "svc_host", subprefix, &(svc->svc_host));
        fprint_NNTI_buffer(fp, "req_addr", subprefix, &(svc->req_addr));
        fprintf(fp, "%s    max_reqs = %u,\n", subprefix, svc->max_reqs);

        /* footer */
        fprintf(fp, "%s }\n", prefix);
        logger_mutex_unlock();
}


/**
 * @brief Print out an ssize.
 *
 * @param fp  The output file.
 * @param name The name of the variable.
 * @param prefix Text that precedes the variable.
 * @param The oid.
 */
void fprint_nssi_ssize(
                FILE *fp,
                const char *name,
                const char *prefix,
                const nssi_ssize *ssize)
{
        char subprefix[100];

        snprintf(subprefix, 100, "%s\t", prefix);

        logger_mutex_lock();
        /* contents */
        fprintf(fp, "%s %s = %lu\n", prefix, name, (unsigned long)*ssize);
        logger_mutex_unlock();
}
