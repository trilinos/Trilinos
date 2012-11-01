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
 * ncVar.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#include "Trios_config.h"
#ifdef HAVE_TRIOS_PNETCDF

#include <Trios_nssi_client.h>
#include <vector>
using namespace std;

#include <memory.h>
#include <pnetcdf.h>
#include "NcVarInfo.h"
#include "NcAttInfo.h"
#include "netcdf_debug.h"



#if USE_NC_TYPE
NcVarInfo::NcVarInfo(
        const int varid,
        const char *nm,
        const nc_type xt,
        const int nd,
        const int dids[]) :
            varid(varid), name(nm), xtype(xt), dimids(dids, dids+nd)
#else
NcVarInfo::NcVarInfo(
        const int varid,
        const char *nm,
        const int xt,
        const int nd,
        const int dids[]) :
            varid(varid), name(nm), xtype(xt), dimids(dids, dids+nd)
#endif
{ }

NcVarInfo::NcVarInfo(const nc_var &var) :
    varid(var.varid), name(var.name), xtype((nc_type)var.xtype),
    dimids(var.dimids.dimids_val, var.dimids.dimids_val + var.dimids.dimids_len)
{
    for (int i=0;i<var.atts.atts_len;i++) {
        this->atts[var.atts.atts_val[i].name] = new NcAttInfo(var.atts.atts_val[i]);
    }
}



NcVarInfo::~NcVarInfo()
{
}

/**
 * Convert to struct nc_var.
 */
int NcVarInfo::copyTo(struct nc_var &var)
{
    int natts, ndims;
    log_level debug_level = netcdf_debug_level;

    memset(&var, 0, sizeof(struct nc_var));

    var.name = strdup(this->name.c_str());
    var.varid = this->varid;
    var.xtype = this->xtype;

    /* copy attributes */
    natts = this->atts.size();
    var.atts.atts_len = natts;
    log_debug(debug_level, "copy %d atts", natts);
    if (natts) {
        map<string, NcAttInfo *>::iterator att_iter;
        int i=0;
        var.atts.atts_val = (struct nc_att *)calloc(natts, sizeof(struct nc_att));
        for (att_iter = this->atts.begin(); att_iter != this->atts.end(); att_iter++) {
            att_iter->second->copyTo(var.atts.atts_val[i++]);
        }
    }

    /* copy dimids */
    ndims = this->dimids.size();
    var.dimids.dimids_len = ndims;
    log_debug(debug_level, "copy %d dimids", ndims);
    if (ndims) {
        int i=0;
        vector<int>::iterator dim_iter;
        var.dimids.dimids_val = (int *)calloc(ndims, sizeof(int));
        for (dim_iter = this->dimids.begin(); dim_iter != this->dimids.end(); dim_iter++) {
            var.dimids.dimids_val[i++] = *dim_iter;
        }
    }

    return NC_NOERR;
}


/** Get the variable ID. */
int NcVarInfo::inq_varid(int *varidp)
{
    int rc = NC_NOERR;
    *varidp = this->varid;
    return rc;
}


/** Get information about a variable. */
#if USE_NC_TYPE
int NcVarInfo::inq_var(char *name, nc_type *xtypep, int *ndimsp,
        int dimids[], int *nattsp)
#else
int NcVarInfo::inq_var(char *name, int *xtypep, int *ndimsp,
        int dimids[], int *nattsp)
#endif
{
    int rc = NC_NOERR;
    if (name != NULL) {
        strcpy(name, this->name.c_str());
    }

    *xtypep = this->xtype;
    *ndimsp = this->dimids.size();

    if (dimids != NULL) {
        std::copy(this->dimids.begin(), this->dimids.end(), dimids);
    }

    *nattsp = this->atts.size();

    return rc;
}


/** Get name of variable. */
int NcVarInfo::inq_varname(char *name)
{
    int rc = NC_NOERR;

    if (name != NULL) {
        strcpy(name, this->name.c_str());
    }

    return rc;
}

/** Get type of variable. */
#if USE_NC_TYPE
int NcVarInfo::inq_vartype(nc_type *xtypep)
#else
int NcVarInfo::inq_vartype(int *xtypep)
#endif
{
    int rc = NC_NOERR;
    *xtypep = this->xtype;
    return rc;
}

/** Get the number of dimensions used for this variable. */
int NcVarInfo::inq_varndims(int *ndimsp)
{
    int rc = NC_NOERR;
    *ndimsp = this->dimids.size();
    return rc;
}

/** Get the dimension ids for this variable. */
int NcVarInfo::inq_vardimid(int dimids[])
{
    int rc = NC_NOERR;
    if (dimids != NULL) {
        std::copy(this->dimids.begin(), this->dimids.end(), dimids);
    }
    return rc;
}


/** Get the number of attributes for this variable. */
int NcVarInfo::inq_varnatts(int *nattsp)
{
    int rc = NC_NOERR;
    *nattsp = this->atts.size();
    return rc;
}


/** Set the chunking parameters for netCDF-4 files. */
int NcVarInfo::def_var_chunking(const int contiguous, int *chunksizep)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Inquire about chunking paramenters for this variable. */
int NcVarInfo::inq_var_chunking(int *contiguousp, int *chunksizep)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Define fill parameters for a variable. */
int NcVarInfo::def_var_fill(const int no_fill, void *fill_value)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Inquire about fill parameters. */
int NcVarInfo::inq_var_fill(int *no_fill, void *fill_value)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Define compression parameters. */
int NcVarInfo::def_var_deflate(const int shuffle, const int deflate,
        const int deflate_level)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Inquire about compression parameters. */
int NcVarInfo::inq_var_deflate(int *shufflep, int *deflatep, int *deflate_levelp)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Define fletcher32 parameters for netcdf-4 variable. */
int NcVarInfo::def_var_fletcher32(const int fletcher32)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Inquire about fletcher32 parameter. */
int NcVarInfo::inq_var_fletcher32(int *fletcher32p)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Define endianness of variable. NC_ENDIAN_NATIVE, NC_ENDIAN_LITTLE, NC_ENDIAN_BIG */
int NcVarInfo::def_var_endian(const int endian)
{
    int rc = NC_ENOTSUPP;
    return rc;
}

/** Inquire about endiannes of variable. */
int NcVarInfo::inq_var_endian(int *endianp)
{
    int rc = NC_ENOTSUPP;
    return rc;
}



/** Get attributes for a variable.
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
/** Define an attribute for this variable. */
#if USE_NC_TYPE
int NcVarInfo::def_att(
        const char *name,
        const nc_type xtype,
        const size_t len)
#else
int NcVarInfo::def_att(
        const char *name,
        const int xtype,
        const size_t len)
#endif
{
    int rc = NC_NOERR;

    if (atts.find(name) == atts.end()) {
        atts[name] = new NcAttInfo(name, xtype, len);
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}

#if USE_NC_TYPE
int NcVarInfo::inq_att    (const char *name,
        nc_type *xtypep, size_t *lenp)
#else
int NcVarInfo::inq_att    (const char *name,
        int *xtypep, size_t *lenp)
#endif
{
    int rc = NC_NOERR;

    if (atts.find(name) != atts.end()) {
        atts[name]->inq_atttype(xtypep);
        atts[name]->inq_attlen(lenp);
    }
    else {
        rc = NC_ENOTATT;
    }

    return rc;
}

#if USE_NC_TYPE
int NcVarInfo::inq_atttype(const char *name,
        nc_type *xtypep)
#else
int NcVarInfo::inq_atttype(const char *name,
        int *xtypep)
#endif
{
    int rc = NC_NOERR;
    if (atts.find(name) != atts.end()) {
        atts[name]->inq_atttype(xtypep);
    }
    else {
        rc = NC_ENOTATT;
    }
    return rc;
}


int NcVarInfo::inq_attlen  (const char *name, size_t *lenp)
{
    int rc = NC_NOERR;
    if (atts.find(name) != atts.end()) {
        atts[name]->inq_attlen(lenp);
    }
    else {
        rc = NC_ENOTATT;
    }
    return rc;
}


int NcVarInfo::inq_attname(int attnum, char *name)
{
    int rc = NC_NOERR;
    std::map<std::string, NcAttInfo *>::iterator iter;

    iter=atts.begin();
    for (int i=0;i<attnum && iter!=atts.end();i++) iter++;

    (*iter).second->inq_attname(name);

    return rc;
}


int NcVarInfo::inq_attid   (const char *name, int *attnump)
{
    int rc = NC_ENOTSUPP;
    return rc;
}


#endif // HAVE_TRIOS_PNETCDF

