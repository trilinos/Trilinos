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
/**
 *   @file Trios_nssi_fprint_types.h
 *
 *   @brief Pretty print the different NSSI types.
 *
 */

#ifndef _TRIOS_NSSI_FPRINT_TYPES_H_
#define _TRIOS_NSSI_FPRINT_TYPES_H_

#include <stdio.h>
#include <ostream>
#include "Trios_nssi_types.h"


#if defined(__STDC__) || defined(__cplusplus)

extern const char* nssi_err_str(const int rc);

    /**
     * @brief Print out a return code.
     *
     * @param fp @input_type The output file.
     * @param name @input_type The name of the variable.
     * @param prefix @input_type Text that precedes the variable.
     * @param rc @input_type The return code.
     */
    void fprint_nssi_return_code(
        FILE *fp,
        const char *name,
        const char *prefix,
        const int rc);
    void fprint_nssi_return_code(
        std::ostream &out,
        const char *name,
        const char *prefix,
        const int rc);


    /**
     * @brief Output the contents of a request header.
     *
     * @param fp      @input_type File pointer (where to send output)
     * @param name    @input_type The name of the variable.
     * @param prefix  @input_type Text to put on every line before the output.
     * @param hdr     @input_type The request header.
     */
    extern void fprint_nssi_request_header(
        FILE *fp,
        const char *name,
        const char *prefix,
        const nssi_request_header *hdr);
    extern void fprint_nssi_request_header(
            std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_request_header *hdr);

    /**
     * @brief Output the contents of a result header.
     *
     * @param fp      @input_type File pointer (where to send output)
     * @param name    @input_type The name of the variable.
     * @param prefix  @input_type Text to put on every line before the output.
     * @param hdr     @input_type The result header.
     */
    extern void fprint_nssi_result_header(
        FILE *fp,
        const char *name,
        const char *prefix,
        const nssi_result_header *hdr);
    extern void fprint_nssi_result_header(
            std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_result_header *hdr);


    extern void fprint_nssi_rpc_encode(
        FILE *fp,
        const char *name,
        const char *prefix,
        const nssi_rpc_encode *rpc_encode);
    extern void fprint_nssi_rpc_encode(
            std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_rpc_encode *rpc_encode);

    extern void fprint_nssi_service(
        FILE *fp,
        const char *name,
        const char *prefix,
        const nssi_service *svc);
    extern void fprint_nssi_service(
            std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_service *svc);

    extern void fprint_nssi_ssize(
        FILE *fp,
        const char *name,
        const char *prefix,
        const nssi_ssize *ssize);
    extern void fprint_nssi_ssize(
            std::ostream &out,
        const char *name,
        const char *prefix,
        const nssi_ssize *ssize);

#else /* K&R C */


#endif /* K&R C */

#endif
