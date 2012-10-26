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
/*
 * ncAttr.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#include <string.h>
#include <Trios_nssi_client.h>
#include <pnetcdf.h>
#include "NcAttInfo.h"
#include "netcdf_args.h"
#include "netcdf_debug.h"


#if USE_NC_TYPE
NcAttInfo::NcAttInfo(const char *nm, const nc_type xt, const size_t l) :
            name(nm), xtype(xt), len(l)
#else
NcAttInfo::NcAttInfo(const char *nm, const int xt, const size_t l) :
            name(nm), xtype(xt), len(l)
#endif
{
}

NcAttInfo::NcAttInfo(const struct nc_att &att) :
    name(att.name), xtype(att.xtype), len(att.len)
{
    log_debug(netcdf_debug_level, "Created attribute (%s, type=%d, len=%d)",
            name.c_str(), (int)xtype, (int)len);
}

NcAttInfo::~NcAttInfo()
{
}

int NcAttInfo::copyTo(struct nc_att &att)
{
    att.xtype = this->xtype;
    att.name = strdup(this->name.c_str());
    att.len = this->len;

    return NC_NOERR;
}

#if USE_NC_TYPE
int NcAttInfo::inq_att(char *name, nc_type *xtypep, size_t *lenp)
#else
int NcAttInfo::inq_att(char *name, int *xtypep, size_t *lenp)
#endif
{
    if (name != NULL) {
        strcpy(name, this->name.c_str());
    }

    *xtypep = static_cast<nc_type>(this->xtype);
    *lenp = this->len;

    return NC_NOERR;
}


int NcAttInfo::inq_attname(char *name)
{
    if (name != NULL) {
        strcpy(name, this->name.c_str());
    }
    return NC_NOERR;
}


#if USE_NC_TYPE
int NcAttInfo::inq_atttype(nc_type *xtypep)
#else
int NcAttInfo::inq_atttype(int *xtypep)
#endif
{
    *xtypep = static_cast<nc_type>(this->xtype);
    return NC_NOERR;
}


int NcAttInfo::inq_attlen(size_t *attlenp)
{
    *attlenp = this->len;
    return NC_NOERR;
}





