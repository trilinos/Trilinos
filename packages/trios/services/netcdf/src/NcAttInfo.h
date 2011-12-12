/*
 * ncAttr.h
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#ifndef NCATTR_H_
#define NCATTR_H_

#include <string>
using namespace std;

#include "netcdf_args.h"

#undef USE_NC_TYPE
#define USE_NC_TYPE 1


/** Attributes may be associated with each netCDF variable to specify
 *  such properties as units, special values, maximum and minimum valid
 *  values, scaling factors, and offsets. Attributes for a netCDF dataset
 *  are defined when the dataset is first created, while the netCDF dataset
 *  is in define mode. Additional attributes may be added later by
 *  reentering define mode. A netCDF attribute has a netCDF variable
 *  to which it is assigned, a name, a type, a length, and a sequence
 *  of one or more values. An attribute is designated by its variable
 *  ID and name. When an attribute name is not known, it may be designated
 *  by its variable ID and number in order to determine its name, using
 *  the function nc_inq_attname.
 *
 *  The attributes associated with a variable are typically defined
 *  immediately after the variable is created, while still in define mode.
 *  The data type, length, and value of an attribute may be changed even
 *  when in data mode, as long as the changed attribute requires no more
 *  space than the attribute as originally defined.
 *
 *  It is also possible to have attributes that are not associated with any
 *  variable. These are called global attributes and are identified by using
 *  NC_GLOBAL as a variable pseudo-ID. Global attributes are usually related
 *  to the netCDF dataset as a whole and may be used for purposes such as
 *  providing a title or processing history for a netCDF dataset.
 *
 *  The NcAttribute Class does not actually contain the attribute data, it
 *  just has the attribute metadata (id, name, type, len).  Use the
 *  NcVariable::get_att methods to obtain the actual data.
 */

class NcAttInfo {

public:
#if USE_NC_TYPE
    int inq_att(char *name, nc_type *xtypep, size_t *lenp);
#else
    int inq_att(char *name, int *xtypep, size_t *lenp);
#endif
    int inq_attname(char *name);
#if USE_NC_TYPE
    int inq_atttype(nc_type *xtypep);
#else
    int inq_atttype(int *xtypep);
#endif
    int inq_attlen(size_t *attlenp);

    /**
     * Convert to struct nc_att.
     */
    int copyTo(struct nc_att &);

#if USE_NC_TYPE
    NcAttInfo(const char *name, const nc_type xtype, const size_t len);
#else
    NcAttInfo(const char *name, const int xtype, const size_t len);
#endif
    NcAttInfo(const nc_att &att);
    virtual ~NcAttInfo();

private:

    int xtype;
    string name;
    size_t len;
};

#endif /* NCATTR_H_ */
