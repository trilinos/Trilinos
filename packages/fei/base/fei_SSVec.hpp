#ifndef _SSVec_hpp_
#define _SSVec_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_iosfwd.hpp"

#include "feiArray.hpp"

enum { SS_Constr_Default, SS_Constr_EqnBuf,
       SS_Constr_RawArrays, SS_Constr_RawArrays2, SS_Constr_RawArraysSymm };

/** SSVec stands for Super-Sparse Vector. It can hold sparse vector data
(e.g., with non-contiguous indices, etc.) and can be used in operations
with the SSMat class.

A couple of constructors are provided for wrapping an SSVec around existing
data that is either in the form of feiArrays or raw arrays. These constructors
are intended to be as light-weight as possible, so they simply keep pointers
to the data being wrapped. There is an inherent danger here, that the data may
be destroyed before the SSVec is destroyed, leaving the SSVec holding bad
pointers. USER BEWARE.
 */

class SSVec {
 public:

  /** Default constructor. */
  SSVec(int alloc_increment=32);

  /** Constructor to wrap an SSVec around existing raw array data. SSVec keeps
      pointers to this array data! In other words, don't delete the raw data
      and then try to use this SSVec.
  */
  SSVec(int numEntries, const int* indices, const double* coefs);

  /** Copy constructor. */
  SSVec(const SSVec& src);

  /** Destructor. */
  ~SSVec();

  /** assignment operator */
  SSVec& operator=(const SSVec& src);

  void logicalClear();

  int length() const;
  int size() const;

  /** Function to set the internal data. */
  int setInternalData(int numEntries, const int* indices, const double* coefs);

  int addEntry(int eqn, double coef);

  int addEntries(int numEntries, const double* coefs, const int* eqns);

  int addEntries_sortedInput(int numEntries, const double* coefs,
			     const int* eqns,
			     bool storeZeros = true);

  int putEntry(int eqn, double coef);

  int putEntries(int numEntries, const double* coefs, const int* eqns);

  feiArray<int>& indices() { return( *indices_ ); }

  feiArray<double>& coefs() { return( *coefs_ ); }

  int whichConstructor_;

  void writeToStream(FEI_OSTREAM& os);

 private:
  feiArray<int>* indices_;
  feiArray<double>* coefs_;
};

#include <fei_ostream_ops.hpp>

#endif

