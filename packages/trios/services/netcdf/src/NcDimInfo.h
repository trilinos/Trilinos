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
 * ncDim.h
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#ifndef NCDIM_H_
#define NCDIM_H_

#include <string>

class NcGroupInfo;

/**
 *  Dimensions are defined while the dataset is in define mode.
 *  A netCDF dimension has a name and length.
 */
class NcDimInfo {

public:

    /** Get information about the dimension. */
    int inq_dim(char *name, size_t *lengthp);

    /** Get the name of the dimension. */
    int inq_dimname(char *name);

    /** Get the length of the dimension. */
    int inq_dimlen(size_t *lengthp);

    /** Get the ID of the dimension. */
    int inq_dimid(int *dimid);

    /** Rename a dimension. */
    int rename_dim(const char *name);

    /** Convert to struct nc_dim */
    int copyTo(struct nc_dim &);


public:
    /** Create a new dimension. */
    NcDimInfo(const int dimid, const char *name, const size_t len);

    /** Creat a new dimension from struct nc_dim */
    NcDimInfo(const struct nc_dim &dim);

    virtual ~NcDimInfo();

private:
    const int dimid;
    std::string name;
    size_t len;

    friend class NcGroupInfo;

};

#endif /* NCDIM_H_ */
