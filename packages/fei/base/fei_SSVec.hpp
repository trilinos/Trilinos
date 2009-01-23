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
 */

class SSVec {
 public:

  /** Default constructor. */
  SSVec(int alloc_increment=32);

  /** Copy constructor. */
  SSVec(const SSVec& src);

  /** Destructor. */
  ~SSVec();

  /** assignment operator */
  SSVec& operator=(const SSVec& src);

  void logicalClear();

  int length() const;
  int size() const;

  int addEntry(int eqn, double coef);

  int addEntries(int numEntries, const double* coefs, const int* eqns);

  int addEntries_sortedInput(int numEntries, const double* coefs,
			     const int* eqns,
			     bool storeZeros = true);

  int putEntry(int eqn, double coef);

  int putEntries(int numEntries, const double* coefs, const int* eqns);

  feiArray<int>& indices() { return( *indices_ ); }

  const feiArray<int>& indices() const { return( *indices_ ); }

  feiArray<double>& coefs() { return( *coefs_ ); }

  const feiArray<double>& coefs() const { return( *coefs_ ); }

  int whichConstructor_;

  void writeToStream(FEI_OSTREAM& os);

 private:
  feiArray<int>* indices_;
  feiArray<double>* coefs_;
};

#include <fei_ostream_ops.hpp>

#endif

