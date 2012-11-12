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
 * ncGroup.h
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#ifndef NCGROUP_H_
#define NCGROUP_H_

#include <string>
#include <map>

#include "NcFileInfo.h"

class NcDimInfo;
class NcAttInfo;
class NcVarInfo;

#undef USE_NC_TYPE
#define USE_NC_TYPE 1

/**
 * The ncGroup represent a netCDF group using the netCDF-4 definition.
 *
 * Excerpt from netcdf page...
 *
 * Groups are identified with a ncid, which identifies both the open file, and the
 * group within that file. When a file is opened with nc_open or nc_create,
 * the ncid for the root group of that file is provided. Using that as a
 * starting point, users can add new groups, or list and navigate existing groups.
 *
 * All netCDF calls take a ncid which determines where the call will take its
 * action. For example, the nc_def_var function takes a ncid as its first parameter.
 * It will create a variable in whichever group its ncid refers to. Use the root
 * ncid provided by nc_create or nc_open to create a variable in the root group.
 * Or use nc_def_grp to create a group and use its ncid to define a variable in
 * the new group.
 *
 * Variable are only visible in the group in which they are defined. The
 * same applies to attributes. �Global� attributes are associated with
 * the group whose ncid is used.
 *
 * Dimensions are visible in their groups, and all child groups.
 *
 * Group operations are only permitted on netCDF-4 files - that is,
 * files created with the HDF5 flag in nc_create. (see nc_create).
 * Groups are not compatible with the netCDF classic data model, so
 * files created with the NC_CLASSIC_MODEL file cannot contain groups
 * (except the root group).
 */

class NcGroupInfo {

public:

    /** Create a new root group */
    NcGroupInfo(const int ncid, const NcFileInfo &finfo);


    /** Create a new child group. */
    NcGroupInfo(const int ncid, const char *name, NcGroupInfo &parent);

    /** Create a new child group from the info in struct nc_group. */
    NcGroupInfo(const struct nc_group &group,  NcGroupInfo &parent);


    /** Create a root group from the struct nc_group. */
    NcGroupInfo(const struct nc_group &group, const NcFileInfo &finfo);


    virtual ~NcGroupInfo();

    /** Construct from a struct nc_group */
    NcGroupInfo(const struct nc_group &group, NcFileInfo &fileInfo);

    /* copy contents from struct nc_group */
    int copyFrom(const struct nc_group &group);

    /* Convert to struct nc_group */
    int copyTo(struct nc_group &group);

    /** Add a new sub-group. */
    int def_grp(const int child_ncid, const char *name);

    /** Return the ncid of this group. */
    int inq_ncid(int *ncid);

    /** Return the ncid of a named child group. */
    int inq_ncid(const char *name, int *ncid);

    /** Return the number of groups and copy the ncids of
     *  each group in a previously allocated array of ints.
     *  If array is NULL, return the number of groups.
     */
    int inq_grps(int *numgrps, int *ncids);

    /** Copy the name of this group into a previously allocated
     *  char array.  If NULL is passed as the array, return only
     *  the length of the array.
     */
    int inq_name(int *namelen, char *name);

    /** Return the id of the groups parent.
     *
     *  @returns NC_ENOGRP if this is the root group.
     */
    int inq_grp_parent(int *parent_ncid);


    /* --------- VARIABLES --------------- */

    /** Add a new variable to this group. */
#if USE_NC_TYPE
    int def_var(const int varid, const char *name,
            const nc_type xtype, const int ndims, const int dimids[]);
#else
    int def_var(const int varid, const char *name,
            const int xtype, const int ndims, const int dimids[]);
#endif

    /** Delete a variable. */
    int del_var(const int varid);

    /** Return the number of variables and copy the varids
     *  of each variable into the previously allocated array
     *  of ints.  If the array is NULL, just return the number
     *  of variables.
     */
    int inq_vars(int *numvars, int *varids);

    /** Get info about a variable, referenced by varid. */
    int inq_var      (const int varid, char *name, nc_type *xtypep,
                         int *ndimsp, int dimids[], int *nattsp);

    /** Get the ID given a variable name. */
    int inq_varid    (const char *name, int *varidp);

    /** Get name of a variable. */
    int inq_varname  (const int varid, char *name);

    /** Get type of variable. */
    int inq_vartype  (const int varid, nc_type *xtypep);

    /** Get number of dimensions used by variable. */
    int inq_varndims (const int varid, int *ndimsp);

    /** Return the total number of values in a variable. */
    int inq_varsize(const int varid, size_t *sizep);

    /** Get dimids used by variable. */
    int inq_vardimid (const int varid, int dimids[]);

    /** Get number of attributes used by a variable. */
    int inq_varnatts (const int varid, int *nattsp);

    /* -------- DIMENSIONS -------------- */

    /** Define a new demension */
    int def_dim(const int dimid, const char *name, const size_t len);

    /** Add a new dimension. */

    /** Rename a dimension. */
    int rename_dim(const int dimid, const char *name);

    /** General inquiry */
    int inq(int *ndimsp, int *nvarsp, int *nattsp, int *unlimdimp);

    /** Find unlimited dimension IDs. */
    int inq_unlimdims(int *nunlimdimsp, int *unlimdimidsp);

    /** Return the dimension ids for all dimensions in a group, or
     *  any of its parents.  If the array is NULL, return the number
     *  of dimensions defined for all associated groups.
     */
    int inq_dims(int *numdims, int *dimids);

    /**
     * As records are written in the unlimited dimension, we have to increase
     * the length of the unlimited dimension.
     */
    int update_unlimdim_dimlen(size_t new_len);

    /** Return the name and length of a dimension referenced by dimid. */
    int inq_dim(const int dimid, char *name, size_t *lengthp);

    /** Return the name of a dimension. */
    int inq_dimname(const int dimid, char *name);

    /** Return the length of a dimension. */
    int inq_dimlen(const int dimid, size_t *lengthp);

    /** Return the dimid of the dimension named 'name' */
    int inq_dimid (const char *name, int *dimid);

    int inq_unlimdimid(int *unlimdimidp);


    /** Attributes
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
    int inq_att    (const int varid, const char *name,
            nc_type *xtypep, size_t *lenp);
    int inq_atttype(const int varid, const char *name,
            nc_type *xtypep);
#else
    int inq_att    (const int varid, const char *name,
            int *xtypep, size_t *lenp);
    int inq_atttype(const int varid, const char *name,
            int *xtypep);
#endif
    int inq_attlen  (const int varid, const char *name, size_t *lenp);
    int inq_attname(const int varid, int attnum, char *name);
    int inq_attid   (const int varid, const char *name, int *attnump);

    /** Define an attribute for this variable. */
#if USE_NC_TYPE
    int def_att(const char *name, const nc_type xtype, const size_t len);
#else
    int def_att(const char *name, const int xtype, const size_t len);
#endif

    /** Delete an attribute. */
    int del_att (const char* name);


protected:


    /** Construct a child group */
    NcGroupInfo(const int ncid, const char *name, NcGroupInfo *parent);



private:
    const int ncid;              /* id of this group */
    const NcFileInfo fileInfo;   /* points to root group (null if we are the root) */
    NcGroupInfo *parent;   /* pointer to parent (null if we are the root) */
    std::string name;            /* name of this group "/" if root */
    int unlimdimid;


public:
    std::map<int, NcVarInfo *> vars;          /* variables used by the group */
    std::map<string, NcVarInfo *> varsByName; /* variables mapped by name */
    std::map<int, NcDimInfo *> dims;          /* Dimensions used by this group */

    std::map<std::string, NcAttInfo *> atts;  /* global attributes for this group */

    std::map<int, NcGroupInfo *> children;    /* Subgroups of this group */
};


/** Global map of open groups (used by ncmpi_client and netcdf_client) */
extern map<int, NcGroupInfo *> group_map;


#endif /* NCGROUP_H_ */
