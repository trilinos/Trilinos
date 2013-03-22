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
 * ncVar.h
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#ifndef NCVAR_H_
#define NCVAR_H_

#include <string>
#include <vector>
#include <map>

using namespace std;

#include "netcdf_args.h"


class NcFileInfo;
class NcAttInfo;

/* whether or not to use the nc_type since netcdf and pnetcdf conflict */
#undef USE_NC_TYPE
#define USE_NC_TYPE 1

class NcVarInfo {

public:

    /** Convert to struct nc_var */
    int copyTo(struct nc_var &);

    /** Get the variable ID. */
    int inq_varid(int *varidp);

    /** Get information about a variable. */
#if USE_NC_TYPE
    int inq_var(char *name, nc_type *xtypep, int *ndmisp,
            int dimids[], int *nattsp);
#else
    int inq_var(char *name, int *xtypep, int *ndmisp,
            int dimids[], int *nattsp);
#endif

    /** Get name of variable. */
    int inq_varname(char *name);

    /** Get type of variable. */
#if USE_NC_TYPE
    int inq_vartype(nc_type *xtypep);
#else
    int inq_vartype(int *xtypep);
#endif

    /** Get the number of dimensions used for this variable. */
    int inq_varndims(int *ndimsp);

    /** Get the dimension ids for this variable. */
    int inq_vardimid(int dimids[]);

    /** Get the number of attributes for this variable. */
    int inq_varnatts(int *nattsp);


    /** Set the chunking parameters for netCDF-4 files. */
    int def_var_chunking(const int contiguous, int *chunksizep);

    /** Inquire about chunking paramenters for this variable. */
    int inq_var_chunking(int *contiguousp, int *chunksizep);

    /** Define fill parameters for a variable. */
    int def_var_fill(const int no_fill, void *fill_value);

    /** Inquire about fill parameters. */
    int inq_var_fill(int *no_fill, void *fill_value);

    /** Define compression parameters. */
    int def_var_deflate(const int shuffle, const int deflate,
            const int deflate_level);

    /** Inquire about compression parameters. */
    int inq_var_deflate(int *shufflep, int *deflatep, int *deflate_levelp);

    /** Define fletcher32 parameters for netcdf-4 variable. */
    int def_var_fletcher32(const int fletcher32);

    /** Inquire about fletcher32 parameter. */
    int inq_var_fletcher32(int *fletcher32p);

    /** Define endianness of variable. NC_ENDIAN_NATIVE, NC_ENDIAN_LITTLE, NC_ENDIAN_BIG */
    int def_var_endian(const int endian);

    /** Inquire about endiannes of variable. */
    int inq_var_endian(int *endianp);



    /** Variable attributes
     *
     *  This family of functions returns information about a netCDF attribute.
     *  All but one of these functions require the variable ID and attribute
     *  name; the exception is nc_inq_attname. Information about an attribute
     *  includes its type, length, name, and number. See the nc_get_att family
     *  for getting attribute values.
     *
     *  The function inq_attname gets the name of an attribute, given its
     *  variable ID and number. This function is useful in generic applications
     *  that need to get the names of all the attributes associated with a
     *  variable, since attributes are accessed by name rather than number
     *  in all other attribute functions. The number of an attribute is more
     *  volatile than the name, since it can change when other attributes of
     *  the same variable are deleted. This is why an attribute number is not
     *  called an attribute ID.
     *
     *  The function nc_inq_att returns the attribute's type and length.
     *  The other functions each return just one item of information about
     *  an attribute.
     */
#if USE_NC_TYPE
    int inq_att    (const char *name,
            nc_type *xtypep, size_t *lenp);
    int inq_atttype(const char *name,
            nc_type *xtypep);
#else
    int inq_att    (const char *name,
            int *xtypep, size_t *lenp);
    int inq_atttype(const char *name,
            int *xtypep);
#endif
    int inq_attlen  (const char *name, size_t *lenp);
    int inq_attname(int attnum, char *name);
    int inq_attid   (const char *name, int *attnump);

    /** Define an attribute for this variable. */
#if USE_NC_TYPE
    int def_att(const char *name, const nc_type xtype, const size_t len);
#else
    int def_att(const char *name, const int xtype, const size_t len);
#endif

    /** Delete an attribute. */
    int del_att (const char* name);



public:

    /** Create a new variable for a netcdf dataset. */
#if USE_NC_TYPE
    NcVarInfo(const int varid, const char *name, const nc_type xtype,
            const int ndims, const int dimids[]);
#else
    NcVarInfo(const int varid, const char *name, const int xtype,
            const int ndims, const int dimids[]);
#endif

    NcVarInfo(const struct nc_var &var);

    virtual ~NcVarInfo();


private:
    const int _varid;

    string _name;
#if USE_NC_TYPE
    const nc_type _xtype;
#else
    const int _xtype;
#endif


    /* chunking parameters */
    int _contiguous;
    int _chunksize;

    /* fill parameters */
    int   _no_fill;
    void *_fill_value;

    /* compression parameters */
    int _shuffle;
    int _deflate;
    int _deflate_level;

    /* fletcher32 parameters */
    int _fletcher32;

public:

    /* Map of attributes, indexed by name. */
    map<string, NcAttInfo *> _atts;

    /* Vector of dimension IDs used by this variable.  The actual
     * dimensions are stored in the NcDataset structure for this ncid.
     */
    vector<int> _dimids;   /* dimension IDS for this variable */
};

#endif /* NCVAR_H_ */
