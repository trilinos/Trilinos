/*
 * ncDim.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

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
            dimid(dimid), name(name), len(len)
{
}

NcDimInfo::NcDimInfo(const nc_dim &dim) :
    dimid(dim.dimid), name(dim.name), len(dim.len)
{
}

NcDimInfo::~NcDimInfo() {

}

/**
 * Convert NcDimInfo into struct nc_dim.
 */
int NcDimInfo::copyTo(struct nc_dim &dim)
{
    dim.dimid = this->dimid;
    dim.name = strdup(this->name.c_str());
    dim.len = this->len;

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
    strcpy(name, this->name.c_str());
    return rc;
}

/** Get the length of the dimension. */
int NcDimInfo::inq_dimlen(size_t *lengthp)
{
    int rc = NC_NOERR;
    *lengthp = this->len;
    return rc;
}

/** Get the ID of the dimension. */
int NcDimInfo::inq_dimid(int *dimid)
{
    int rc = NC_NOERR;
    *dimid = this->dimid;
    return rc;
}

/** Rename a dimension. */
int NcDimInfo::rename_dim(const char *newname)
{
    int rc = NC_NOERR;
    this->name = string(newname);
    return rc;
}
