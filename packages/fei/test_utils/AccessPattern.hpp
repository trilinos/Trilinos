#ifndef _AccessPattern_h_
#define _AccessPattern_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <cstdlib>

class AccessPattern {
 public:
  AccessPattern() : ID_(-1), numRowIDs_(0), numFieldsPerRow_(NULL),
    rowFieldIDs_(NULL), numColIDsPerRow_(0), numFieldsPerCol_(NULL),
    colFieldIDs_(NULL), interleaveStrategy_(0) {}

  ~AccessPattern()
    {
      int i;
      for(i=0; i<numRowIDs_; i++) delete [] rowFieldIDs_[i];
      for(i=0; i<numColIDsPerRow_; i++) delete [] colFieldIDs_[i];

      delete [] rowFieldIDs_;
      delete [] colFieldIDs_;
      delete [] numFieldsPerRow_;
      delete [] numFieldsPerCol_;
      numRowIDs_ = 0;
      numColIDsPerRow_ = 0;
    }

  int ID_;
  int numRowIDs_;
  int* numFieldsPerRow_;
  int** rowFieldIDs_;
  int numColIDsPerRow_;
  int* numFieldsPerCol_;
  int** colFieldIDs_;
  int interleaveStrategy_;
};

#endif // _AccessPattern_h_
