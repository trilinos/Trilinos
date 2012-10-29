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
 * ncGroup.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#include "Trios_config.h"
#ifdef HAVE_TRIOS_PNETCDF

#include <Trios_nssi_client.h>
#include <assert.h>
#include <string.h>
#include <iostream>

#include <pnetcdf.h>
#include <map>
#include <string>
#include <algorithm>
using namespace std;


#include "netcdf_args.h"
#include "netcdf_debug.h"
#include "NcGroupInfo.h"
#include "NcVarInfo.h"
#include "NcDimInfo.h"
#include "NcFileInfo.h"
#include "NcAttInfo.h"

#ifndef NC_ENOGRP
#define NC_ENOGRP  (-125)  /* no group found */
#endif

/** Global map of open groups (used by ncmpi_client and netcdf_client) */
map<int, NcGroupInfo *> group_map;


/** Create a new root group */
NcGroupInfo::NcGroupInfo(const int ncid, const NcFileInfo &finfo) :
    ncid(ncid), name("/"), parent(NULL), fileInfo(finfo), unlimdimid(-1)
{ }


/** Create a new child group. */
NcGroupInfo::NcGroupInfo(const int ncid, const char *name, NcGroupInfo &parent) :
    ncid(ncid), name(name), parent(&parent), fileInfo(parent.fileInfo), unlimdimid(-1)
{
}

NcGroupInfo::NcGroupInfo(const struct nc_group &group, NcGroupInfo &parent) :
    ncid(group.ncid), parent(&parent), fileInfo(parent.fileInfo)
{
    copyFrom(group);
}


/** Create a group class from the struct nc_group data structure. */
NcGroupInfo::NcGroupInfo(const struct nc_group &group, const NcFileInfo &finfo) :
    ncid(group.ncid), parent(NULL), fileInfo(finfo)
{
    /* fill contents of this class with info from nc_group */
    copyFrom(group);
}


int NcGroupInfo::copyFrom(const struct nc_group &group)
{
    log_level debug_level = netcdf_debug_level;

    int ndims, nvars, natts, ngrps;
    int i;

    this->unlimdimid = group.unlimdimid;

    /* copy dimensions */
    log_debug(debug_level, "copy %d dims", group.dims.dims_len);
    for (i=0; i<group.dims.dims_len; i++) {
        nc_dim *dim = &group.dims.dims_val[i];
        this->dims[dim->dimid] = new NcDimInfo(*dim);
    }

    /* copy vars */
    log_debug(debug_level, "copy %d vars", group.vars.vars_len);
    for (i=0; i<group.vars.vars_len; i++) {
        nc_var *var = &group.vars.vars_val[i];
        this->vars[var->varid] = this->varsByName[var->name] = new NcVarInfo(*var);
    }

    /* copy global attributes */
    log_debug(debug_level, "copy %d atts", group.atts.atts_len);
    for (i=0; i<group.atts.atts_len; i++) {
        nc_att *att = &group.atts.atts_val[i];
        this->atts[att->name] = new NcAttInfo(*att);
    }


    /* copy subgroups */
    log_debug(debug_level, "copy %d groups", group.groups.groups_len);
    for (i=0; i<group.groups.groups_len; i++) {
        nc_group *subgrp = &group.groups.groups_val[i];
        this->children[subgrp->ncid] = new NcGroupInfo(*subgrp, *this);
    }

    return 0;
}




struct delete_object
{
    template <typename keyT, typename dataT>
    void operator()(pair<keyT, dataT *> &p) {
        dataT *data = p.second;
        if (data) {
            delete data;
            data = 0;
        }
    }
};

NcGroupInfo::~NcGroupInfo()
{
    log_level debug_level = netcdf_debug_level;

    /* Delete variables */
    log_debug(debug_level, "Deleting vars (%d)", vars.size());
    for_each(vars.begin(), vars.end(), delete_object());

    /* Delete dimensions */
    log_debug(debug_level, "Deleting dims (%d)", dims.size());
    for_each(dims.begin(), dims.end(), delete_object());

    /* Delete attributes */
    log_debug(debug_level, "Deleting atts (%d)", atts.size());
    for_each(atts.begin(), atts.end(), delete_object());

    /* Delete subgroups */
    log_debug(debug_level, "Deleting children (%d)", children.size());
    for_each(children.begin(), children.end(), delete_object());
}

/**
 * Convert to a serializable struct nc_group.
 */
int NcGroupInfo::copyTo(struct nc_group &group)
{
    int ndims, nvars, natts, ngrps;
    log_level debug_level = netcdf_debug_level;

    memset(&group, 0, sizeof(struct nc_group));

    group.ncid = this->ncid;
    group.parent_ncid = (this->parent)? this->parent->ncid: -1;
    group.name = strdup(this->name.c_str());


    /* copy dimensions */
    ndims = this->dims.size();
    group.dims.dims_len = ndims;
    log_debug(debug_level, "copy %d dims", ndims);
    if (ndims) {
        std::map<int, NcDimInfo *>::iterator dim_iter;
        int i=0;
        group.dims.dims_val = (struct nc_dim *)calloc(ndims, sizeof(struct nc_dim));
        for (dim_iter = this->dims.begin(); dim_iter != this->dims.end(); dim_iter++) {
            dim_iter->second->copyTo(group.dims.dims_val[i++]);
        }
    }

    /* copy vars */
    nvars = this->vars.size();
    group.vars.vars_len = nvars;
    log_debug(debug_level, "copy %d vars", nvars);
    if (nvars) {
        std::map<int, NcVarInfo *>::iterator var_iter;
        int i=0;
        group.vars.vars_val = (struct nc_var *)calloc(nvars, sizeof(struct nc_var));
        for (var_iter = this->vars.begin(); var_iter != this->vars.end(); var_iter++) {
            var_iter->second->copyTo(group.vars.vars_val[i++]);
        }
    }

    /* copy global attributes */
    natts = this->atts.size();
    group.atts.atts_len = natts;
    log_debug(debug_level, "copy %d atts", natts);
    if (natts) {
        std::map<std::string, NcAttInfo *>::iterator att_iter;
        int i=0;
        group.atts.atts_val = (struct nc_att *)calloc(natts, sizeof(struct nc_att));
        for (att_iter = this->atts.begin(); att_iter != this->atts.end(); att_iter++) {
            att_iter->second->copyTo(group.atts.atts_val[i++]);
        }
    }

    /* copy subgroups */
    ngrps = this->children.size();
    group.groups.groups_len = ngrps;
    log_debug(debug_level, "copy %d subgroups", ngrps);
    if (ngrps) {
        std::map<int, NcGroupInfo *>::iterator grp_iter;
        int i=0;
        group.groups.groups_val = (struct nc_group *)calloc(ngrps, sizeof(struct nc_group));
        for (grp_iter = this->children.begin(); grp_iter != this->children.end(); grp_iter++) {
            grp_iter->second->copyTo(group.groups.groups_val[i++]);
        }
    }

    return NC_NOERR;
}

/**
 * Return ndims, nvars, natts, and unlimdimp.
 */
int NcGroupInfo::inq(
        int *ndimsp,
        int *nvarsp,
        int *nattsp,
        int *unlimdimidp)
{
    int rc = NC_NOERR;

    *ndimsp = this->dims.size();
    *nvarsp = this->vars.size();
    *nattsp = this->atts.size();
    *unlimdimidp = this->unlimdimid;

    return rc;
}



/** Return the unlimdimid of this group. */
int NcGroupInfo::inq_unlimdimid(int *unlimdimidp)
{
    int rc = NC_NOERR;
    *unlimdimidp = this->unlimdimid;
    return rc;
}

/** Return the ncid of this group. */
int NcGroupInfo::inq_ncid(int *grp_ncid)
{
    int rc = NC_NOERR;
    *grp_ncid = this->ncid;
    return rc;
}

/** Return the ncid of a named group */
int NcGroupInfo::inq_ncid(const char *name, int *grp_ncid)
{
    int rc = NC_NOERR;
    return rc;
}

/** Return the number of groups and copy the ncids of
 *  each group in a previously allocated array of ints.
 *  If array is NULL, return the number of groups.
 */
int NcGroupInfo::inq_grps(int *numgrps, int *ncids)
{
    int rc = NC_NOERR;

    *numgrps = children.size();

    if (ncids != NULL) {
        int count = 0;
        map<int, NcGroupInfo *>::iterator iter;

        for (iter = children.begin(); iter != children.end(); iter++) {
            ncids[count++] = iter->first;
        }
    }

    return rc;
}

/** Return the length of the dimension (dimid) */
int NcGroupInfo::inq_dimlen (const int dimid, size_t *lenp)
{
    int rc = NC_NOERR;

    if (dims.find(dimid) != dims.end()) {
        *lenp = dims[dimid]->len;
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}

/** Return the name and length of the dimension (dimid) */
int NcGroupInfo::inq_dim (const int dimid, char *name, size_t *lenp)
{
    int rc = NC_NOERR;

    if (dims.find(dimid) != dims.end()) {
        *lenp = dims[dimid]->len;
        if (name != NULL) {
            strcpy(name, dims[dimid]->name.c_str());
        }
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}

/** Return the dimid of the dimension named 'name' */
int NcGroupInfo::inq_dimid (const char *name, int *dimid)
{
    int rc = NC_EBADDIM;
    std::map<int, NcDimInfo *>::iterator iter;

    if (dims.empty()) {
        rc=NC_EBADDIM;
    } else {
        iter = dims.begin();
        for (;iter != dims.end(); iter++) {
            if (!strcmp((*iter).second->name.c_str(), name)) {
                *dimid=(*iter).first;
                rc = NC_NOERR;
                break;
            }
        }
    }

    return rc;
}

/** Return the dimension IDs */
int NcGroupInfo::inq_vardimid (const int varid, int dimids[])
{
    int rc = NC_NOERR;

    if (vars.find(varid) != vars.end()) {
        NcVarInfo *varInfo = vars[varid];
        copy(varInfo->dimids.begin(), varInfo->dimids.end(), dimids);
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}

/** Return the dimension IDs */
int NcGroupInfo::inq_varndims (const int varid, int *ndimsp)
{
    int rc = NC_NOERR;
    *ndimsp = 0;

    if (vars.find(varid) != vars.end()) {
        *ndimsp = vars[varid]->dimids.size();
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}


/** Add a new variable to this group. */
#if USE_NC_TYPE
int NcGroupInfo::def_var(
        const int varid,
        const char *name,
        const nc_type xtype,
        const int ndims,
        const int dimids[])
#else
int NcGroupInfo::def_var(
        const int varid,
        const char *name,
        const int xtype,
        const int ndims,
        const int dimids[])
#endif
{
    int rc = NC_NOERR;

    if (vars.find(varid) == vars.end()) {
        NcVarInfo *varInfo = new NcVarInfo(varid, name, xtype, ndims, dimids);
        vars[varid] = varInfo;
        varsByName[name] = varInfo;
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}

/** Return the number of variables and copy the varids
 *  of each variable into the previously allocated array
 *  of ints.  If the array is NULL, just return the number
 *  of variables.
 */
int NcGroupInfo::inq_vars(int *numvars, int *varids)
{
    int rc = NC_NOERR;
    map<int, NcVarInfo *>::iterator iter;

    *numvars = this->vars.size();

    if (varids != NULL) {
        int count = 0;

        for (iter = vars.begin(); iter != vars.end(); iter++) {
            varids[count++] = iter->first;
        }

        assert(*numvars == count);
    }

    return rc;
}

/** Get info about a variable, referenced by varid. */
#if USE_NC_TYPE
int NcGroupInfo::inq_var (
        const int varid,
        char *name,
        nc_type *xtypep,
        int *ndimsp,
        int dimids[],
        int *nattsp)
#else
int NcGroupInfo::inq_var (
        const int varid,
        char *name,
        int *xtypep,
        int *ndimsp,
        int dimids[],
        int *nattsp)
#endif
{
    int rc = NC_NOERR;

    if (vars.find(varid) != vars.end()) {
        rc = vars[varid]->inq_var(name, xtypep, ndimsp, dimids, nattsp);
    }
    else {
        rc = NC_ENOTVAR;
    }

    return rc;
}

int NcGroupInfo::inq_varid(
        const char *name,
        int *varidp)
{
    int rc = NC_NOERR;

    if (varsByName.find(name) != varsByName.end()) {
        rc = varsByName[name]->inq_varid(varidp);
    }
    else {
        rc = NC_ENOTVAR;
    }
    return rc;
}

/**
 * Returns the number of values in a variable by summing
 * up the dimension lengths of each dimension.
 *
 * TODO: cache total dimension sum.
 *
 * We need to traverse the data structure every time because one of the
 * dimensions may not be fixed.  There is room for optimization here.
 */
int NcGroupInfo::inq_varsize(const int varid, size_t *countp)
{
    int rc = NC_NOERR;

    *countp = 1;

    if (vars.find(varid) != vars.end()) {
        int i;
        NcVarInfo *varinfo = vars[varid];

        for (i=0; i<varinfo->dimids.size(); i++) {
            int dimid = varinfo->dimids[i];

            /* find the dim */
            if (dims.find(dimid) != dims.end()) {
                NcDimInfo *dim = dims[dimid];
                *countp *= dim->len;
            }
            else {
                log_error(netcdf_debug_level, "Could not find dimid=%d", dimid);
                rc = NC_ENOTVAR;
                goto cleanup;
            }
        }
    }
    else {
        log_error(netcdf_debug_level, "Could not find varid=%d", varid);
        rc = NC_ENOTVAR;
        goto cleanup;
    }

cleanup:
    return rc;
}

/**
 * Returns the type of a specified variable.
 *
 */
#if USE_NC_TYPE
int NcGroupInfo::inq_vartype(const int varid, nc_type *xtypep)
#else
int NcGroupInfo::inq_vartype(const int varid, int *xtypep)
#endif
{
        int rc = NC_NOERR;

    if (vars.find(varid) != vars.end()) {
            rc = vars[varid]->inq_vartype(xtypep);
    }

    return rc;
}



/** Define a new demension */
int NcGroupInfo::def_dim(
        const int dimid,
        const char *name,
        const size_t len)
{
    int rc = NC_NOERR;

    if (dims.find(dimid) == dims.end()) {
        dims[dimid] = new NcDimInfo(dimid, name, len);
        if (len == NC_UNLIMITED) {
            this->unlimdimid = dimid;
        }
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}

/** Return the dimension ids for all dimensions in a group, or
 *  any of its parents.  If the array is NULL, return the number
 *  of dimensions defined for all associated groups.
 */
int NcGroupInfo::inq_dims(int *numdims, int *dimids)
{
    int rc = NC_NOERR;
    map<int, NcDimInfo *>::iterator iter;

    *numdims = this->dims.size();

    if (dimids != NULL) {
        int count = 0;

        for (iter = dims.begin(); iter != dims.end(); iter++) {
            dimids[count++] = iter->first;
        }

        assert(*numdims == count);
    }

    return rc;
}

/**
 * As records are written in the unlimited dimension, we have to increase
 * the length of the unlimited dimension.
 */
int NcGroupInfo::update_unlimdim_dimlen(size_t new_len)
{
    int rc = NC_NOERR;

    if (this->unlimdimid == -1) {
        rc = NC_EBADDIM;
    } else {
        if (dims.find(this->unlimdimid) != dims.end()) {
            if (dims[this->unlimdimid]->len < new_len) {
                dims[this->unlimdimid]->len=new_len;
            }
        }
        else {
            rc = NC_EBADDIM;
        }
    }

    return rc;
}


/** Copy the name of this group into a previously allocated
 *  char array.  If NULL is passed as the array, return only
 *  the length of the array.
 */
int NcGroupInfo::inq_name(int *namelen, char *name)
{
    int rc = NC_NOERR;

    *namelen = this->name.size();

    if (name != NULL) {
        strcpy(name, this->name.c_str());
    }

    return rc;
}

/** Return the id of the groups parent.
 *
 *  @returns NC_ENOGRP if this is the root group.
 */
int NcGroupInfo::inq_grp_parent(int *parent_ncid)
{
    int rc = NC_NOERR;
    if (this->parent) {
        rc = this->parent->inq_ncid(parent_ncid);
    }
    else {
        rc = NC_ENOGRP;
    }
    return rc;
}

/** Create a new child of this group.
 */
int NcGroupInfo::def_grp(const int new_ncid, const char *name)
{
    int rc = NC_NOERR;

    if (children.find(new_ncid) == children.end()) {
        children[new_ncid] = new NcGroupInfo(new_ncid, name, *this);
    }
    else {
        rc = NC_EEXIST;
    }

    return rc;
}



/** Global attributes
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
int NcGroupInfo::def_att(
        const char *name,
        const nc_type xtype,
        const size_t len)
#else
int NcGroupInfo::def_att(
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
int NcGroupInfo::inq_att(
        const int varid,
        const char *name,
        nc_type *xtypep,
        size_t *lenp)
#else
int NcGroupInfo::inq_att(
        const int varid,
        const char *name,
        int *xtypep,
        size_t *lenp)
#endif
{
    int rc = NC_NOERR;

    if (varid == NC_GLOBAL) {
        if (atts.find(name) != atts.end()) {
            atts[name]->inq_atttype(xtypep);
            atts[name]->inq_attlen(lenp);
        }
        else {
            rc = NC_ENOTATT;
        }
    }

    else {
        if (vars.find(varid) != vars.end()) {
            vars[varid]->inq_att(name, xtypep, lenp);
        }
        else {
            rc = NC_ENOTVAR;
        }
    }
    return rc;
}


#if USE_NC_TYPE
int NcGroupInfo::inq_atttype(
        const int varid,
        const char *name,
        nc_type *xtypep)
#else
int NcGroupInfo::inq_atttype(
        const int varid,
        const char *name,
        int *xtypep)
#endif
{
    int rc = NC_NOERR;

    if (varid == NC_GLOBAL) {
        if (atts.find(name) != atts.end()) {
            atts[name]->inq_atttype(xtypep);
        }
        else {
            rc = NC_ENOTATT;
        }
    }

    else {
        if (vars.find(varid) != vars.end()) {
            vars[varid]->inq_atttype(name, xtypep);
        }
        else {
            rc = NC_ENOTVAR;
        }
    }
    return rc;
}


int NcGroupInfo::inq_attlen(
        const int varid,
        const char *name,
        size_t *lenp)
{
    int rc = NC_NOERR;

    if (varid == NC_GLOBAL) {
        if (atts.find(name) != atts.end()) {
            atts[name]->inq_attlen(lenp);
        }
        else {
            rc = NC_ENOTATT;
        }
    }

    else {
        if (vars.find(varid) != vars.end()) {
            vars[varid]->inq_attlen(name, lenp);
        }
        else {
            rc = NC_ENOTVAR;
        }
    }
    return rc;
}

int NcGroupInfo::inq_attname(
        const int varid,
        const int attnum,
        char *name)
{
    int rc = NC_NOERR;

    if (varid == NC_GLOBAL) {
        std::map<std::string, NcAttInfo *>::iterator iter;

        iter=atts.begin();
        for (int i=0;i<attnum && iter!=atts.end();i++) iter++;

        (*iter).second->inq_attname(name);
    }

    else {
        if (vars.find(varid) != vars.end()) {
            vars[varid]->inq_attname(attnum, name);
        }
        else {
            rc = NC_ENOTVAR;
        }
    }
    return rc;
}

int NcGroupInfo::inq_attid(
        const int varid,
        const char *name,
        int *attnump)
{
    int rc = NC_ENOTSUPP;
    return rc;
}



/** Delete an attribute. */
int NcGroupInfo::del_att (const char* name)
{
    int rc = NC_NOERR;
    return rc;
}

#endif // HAVE_TRIOS_PNETCDF
