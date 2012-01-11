/*
 * ncDataset.h
 *
 *  Created on: Jan 22, 2009
 *      Author: raoldfi
 */

#ifndef NCFILEINFO_H_
#define NCFILEINFO_H_


#include <string>
using namespace std;

#include <Trios_nssi_client.h>


class NcFileInfo {

public:

    NcFileInfo (const char *path, const int cmode,
            const size_t initialsz, const size_t chunksize);

    NcFileInfo(const NcFileInfo &);

    virtual ~NcFileInfo();

    /** Set default creation format. */
    int set_format(const int format);

    /** Return the format of the dataset. */
    int inq_format(int *formatp);


private:

    const string path;
    const int mode;
    const size_t initialsz;
    const size_t chunksize;
    int format;

    /** Constructor for netcdf Datasets.  Only members of this class
     *  can call the constructor.
     */
    NcFileInfo();
};

#endif /* NCDATASET_H_ */
