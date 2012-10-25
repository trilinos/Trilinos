#ifndef _CoefAccess_h_
#define _CoefAccess_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <cstdlib>

class CoefAccess {
 public:
  CoefAccess() : patternID_(-1), numRowIDs_(0), rowIDs_(NULL),
    numColIDsPerRow_(0), colIDs_(NULL), numRowCoefs_(0), numColCoefs_(0),
    coefs_(NULL) {}

  CoefAccess(const CoefAccess& src)
    {
      *this = src;
    }

  CoefAccess& operator=(const CoefAccess& src)
    {
      patternID_ = src.patternID_;

      numRowIDs_ = src.numRowIDs_;
      numColIDsPerRow_ = src.numColIDsPerRow_;
      numRowCoefs_ = src.numRowCoefs_;
      numColCoefs_ = src.numColCoefs_;

      if (numRowIDs_ > 0) {
	rowIDs_ = new GlobalID[numRowIDs_];
	for(int i=0; i<numRowIDs_; i++) rowIDs_[i] = src.rowIDs_[i];
      }

      if (numColIDsPerRow_ > 0 && numRowIDs_ > 0) {
	int len = numRowIDs_*numColIDsPerRow_;
	colIDs_ = new GlobalID[len];
	for(int i=0; i<len; i++) colIDs_[i] = src.colIDs_[i];
      }

      if (numRowCoefs_ > 0 && numColCoefs_ > 0) {
	int len = numRowCoefs_*numColCoefs_;
	coefs_ = new double[len];
	for(int i=0; i<len; i++) coefs_[i] = src.coefs_[i];
      }

      return(*this);
    }

  ~CoefAccess()
    {
      delete [] rowIDs_; delete [] colIDs_; delete [] coefs_;
      numRowIDs_ = 0; numColIDsPerRow_ = 0; numRowCoefs_ = 0; numColCoefs_ = 0;
    }

  int patternID_;

  int numRowIDs_;
  GlobalID* rowIDs_;

  int numColIDsPerRow_;
  GlobalID* colIDs_;

  int numRowCoefs_;
  int numColCoefs_;

  double* coefs_;
};

#endif // _CoefAccess_h_
