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
#include <sstream>
#include <ostream>

#include "Trios_logger.h"
#include "Trios_nssi_types.h"
#include "Trios_nnti_fprint_types.h"
#include "Trios_nssi_fprint_types.h"

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include <unistd.h>



/* ----------- Utility functions to output types ---- */



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
 * @brief Print out a return code.
 *
 * @param fp  The output file.
 * @param name The name of the variable.
 * @param prefix Text that precedes the variable.
 * @param The return code.
 */
void fprint_nssi_return_code(
        std::ostream &out,
        const char *name,
        const char *prefix,
        const int rc)
{
    out << prefix << " " << name << " = " << nssi_err_str(rc) << std::endl;
}

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
    std::stringstream out(std::stringstream::out);

    fprint_nssi_return_code(out, name, prefix, rc);

    fprintf(fp, "%s", out.str().c_str());
}


/**
 * @brief Output the contents of a request header.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param hdr     The request header.
 */
void fprint_nssi_request_header(
        std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_request_header *hdr)
{
    std::string subprefix(prefix);
    subprefix.append("  ");  // add two spaces

    if (hdr == NULL) {
        out << prefix << " " << name << " = NULL " << std::endl;
        return;
    }

    /* header */
    out << prefix << " " << name << " = {" << std::endl;

    /* contents */
    out << subprefix << " id = " << hdr->id << std::endl;
    out << subprefix << " opcode = " << hdr->opcode << std::endl;
    out << subprefix << " fetch_args = " << hdr->fetch_args << std::endl;
    fprint_NNTI_buffer(out, "args_addr", subprefix.c_str(), &hdr->args_addr);
    fprint_NNTI_buffer(out, "data_addr", subprefix.c_str(), &hdr->data_addr);
    fprint_NNTI_buffer(out, "res_addr", subprefix.c_str(), &hdr->res_addr);

    /* footer */
    out << subprefix << " }" << std::endl;
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
    std::stringstream out(std::stringstream::out);

    fprint_nssi_request_header(out, name, prefix, hdr);

    fprintf(fp, "%s", out.str().c_str());
}


/**
 * @brief Output the contents of a result header.
 *
 * @param fp      File pointer (where to send output)
 * @param prefix  Text to put on every line before the output.
 * @param hdr     The result header.
 */
void fprint_nssi_result_header(
        std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_result_header *hdr)
{
    std::string subprefix(prefix);
    subprefix.append("  ");  // add two spaces

    if (hdr == NULL) {
        out << prefix << " " << name << " = NULL " << std::endl;
        return;
    }

    /* header */
    out << prefix << " " << name << " = {" << std::endl;

    /* contents */
    out << subprefix << " id = " << hdr->id << std::endl;
    out << subprefix << " opcode = " << hdr->opcode << std::endl;
    out << subprefix << " fetch_result = " << hdr->fetch_result << std::endl;
    out << subprefix << " result_size = " << hdr->result_size << std::endl;
    fprint_NNTI_buffer(out, "res_addr", subprefix.c_str(), &hdr->result_addr);
    out << subprefix << " rc = " << hdr->rc << std::endl;

    /* footer */
    out << subprefix << " }" << std::endl;
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
    std::stringstream out(std::stringstream::out);

    fprint_nssi_result_header(out, name, prefix, hdr);

    fprintf(fp, "%s", out.str().c_str());
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
        std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_rpc_encode *rpc_encode)
{
    std::string subprefix(prefix);
    subprefix.append("  ");  // add two spaces

    /* contents */
    switch (*rpc_encode) {
    case NSSI_RPC_XDR:
        out << subprefix << " " << name << " = NSSI_RPC_XDR" << std::endl;
        break;

    default:
        out << subprefix << " " << name << " = UNDEFINED" << std::endl;
        break;
    }
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
    std::stringstream out(std::stringstream::out);

    fprint_nssi_rpc_encode(out, name, prefix, rpc_encode);

    fprintf(fp, "%s", out.str().c_str());
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
        std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_service *svc)
{
    std::string subprefix(prefix);
    subprefix.append("  ");  // add two spaces

    if (svc == NULL) {
        out << prefix << " " << name << " = NULL " << std::endl;
        return;
    }

    /* header */
    out << prefix << " " << name << " = {" << std::endl;

    /* contents */
    fprint_nssi_rpc_encode(out, "rpc_encode", subprefix.c_str(), &(svc->rpc_encode));
    fprint_NNTI_peer(out, "svc_host", subprefix.c_str(), &(svc->svc_host));
    fprint_NNTI_buffer(out, "req_addr", subprefix.c_str(), &(svc->req_addr));
    out << subprefix << " max_reqs = " << svc->max_reqs << std::endl;

    /* footer */
    out << subprefix << " }" << std::endl;
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
    std::stringstream out(std::stringstream::out);

    fprint_nssi_service(out, name, prefix, svc);

    fprintf(fp, "%s", out.str().c_str());
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
        std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_ssize *ssize)
{
    std::string subprefix(prefix);
    subprefix.append("  ");  // add two spaces

    out << subprefix << " " << name << " = " << *ssize << std::endl;
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
    std::stringstream out(std::stringstream::out);

    fprint_nssi_ssize(out, name, prefix, ssize);

    fprintf(fp, "%s", out.str().c_str());
}
