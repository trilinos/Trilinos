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





