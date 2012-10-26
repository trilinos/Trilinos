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
 * ncDataset.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */
#include <pnetcdf.h>
#include "NcFileInfo.h"
#include "netcdf_args.h"

NcFileInfo::NcFileInfo(
        const char *path,
        const int mode,
        const size_t initialsz,
        const size_t chunksize)
: path(path), mode(mode), initialsz(initialsz), chunksize(chunksize)
{
}

NcFileInfo::NcFileInfo(const NcFileInfo &copy)
: path(copy.path), mode(copy.mode), initialsz(copy.initialsz), chunksize(copy.chunksize), format(copy.format)
{
}


NcFileInfo::~NcFileInfo() {
}


/** Set default creation format. */
int NcFileInfo::set_format(const int format)
{
    int rc = NC_NOERR;
    this->format = format;
    return rc;
}

/** Return the format of the dataset. */
int NcFileInfo::inq_format(int *formatp)
{
    int rc = NC_NOERR;
    *formatp = this->format;
    return rc;
}
