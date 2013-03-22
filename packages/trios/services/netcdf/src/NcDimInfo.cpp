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
 * ncDim.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#include "Trios_config.h"
#ifdef HAVE_TRIOS_PNETCDF

#include <string>
using namespace std;

#include "Trios_nssi_types.h"

#include <memory.h>
#include <pnetcdf.h>
#include "netcdf_args.h"
#include "NcDimInfo.h"

NcDimInfo::NcDimInfo(
        const int dimid,
        const char *name,
        const size_t len) :
            _dimid(dimid), _name(name), _len(len)
{
}

NcDimInfo::NcDimInfo(const nc_dim &dim) :
    _dimid(dim.dimid), _name(dim.name), _len(dim.len)
{
}

NcDimInfo::~NcDimInfo() {

}

/**
 * Convert NcDimInfo into struct nc_dim.
 */
int NcDimInfo::copyTo(struct nc_dim &dim)
{
    dim.dimid = this->_dimid;
    dim.name = strdup(this->_name.c_str());
    dim.len = this->_len;

    return NC_NOERR;
}


/** Get information about the dimension. */
int NcDimInfo::inq_dim(char *name, size_t *lengthp)
{
    int rc = NC_NOERR;
    inq_dimname(name);
    inq_dimlen(lengthp);
    return rc;
}

/** Get the name of the dimension. */
int NcDimInfo::inq_dimname(char *name)
{
    int rc = NC_NOERR;
    strcpy(name, this->_name.c_str());
    return rc;
}

/** Get the length of the dimension. */
int NcDimInfo::inq_dimlen(size_t *lengthp)
{
    int rc = NC_NOERR;
    *lengthp = this->_len;
    return rc;
}

/** Get the ID of the dimension. */
int NcDimInfo::inq_dimid(int *dimid)
{
    int rc = NC_NOERR;
    *dimid = this->_dimid;
    return rc;
}

/** Rename a dimension. */
int NcDimInfo::rename_dim(const char *newname)
{
    int rc = NC_NOERR;
    this->_name = string(newname);
    return rc;
}

#endif // HAVE_TRIOS_PNETCDF
