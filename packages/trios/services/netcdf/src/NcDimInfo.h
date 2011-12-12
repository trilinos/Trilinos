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
