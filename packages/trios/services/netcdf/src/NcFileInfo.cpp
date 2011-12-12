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
