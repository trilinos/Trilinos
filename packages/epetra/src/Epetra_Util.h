
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_UTIL_H_
#define _EPETRA_UTIL_H_
#include "Epetra_Object.h"
class Epetra_CrsMatrix;
class Epetra_MultiVector;
//! Epetra_Util:  The Epetra Util Wrapper Class.
/*! The Epetra_Util class is a collection of useful functions that cut across a broad
    set of other classes.  Specifically, sorting is provided by this class.

    Epetra_Util is a serial interface only.  This is appropriate since the standard 
    utilities are only specified for serial execution (or shared memory parallel).
*/

class Epetra_Util {
    
  public:
  //! Epetra_Util Constructor.
  /*! Builds an instance of a serial Util object.
   */
  Epetra_Util(void);


  //! Epetra_Util Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_Util instance.
  */
  Epetra_Util(const Epetra_Util& Util);

  //! Epetra_Util Destructor.
  virtual ~Epetra_Util(void);
  
  //! Epetra_Util Sort Routine (Shell sort)
  /*! 

    This function sorts a list of integer values in ascending or descending order.  Additionally it sorts any
    number of companion lists of doubles or ints.  A shell sort is used, which is fast if indices are already sorted.
    
    \param In
           SortAscending - Sort keys in ascending order if true, otherwise sort in descending order..
    \param In
           NumKeys - Number of integer values to be sorted.
    \param In/Out
           Keys - List of integers to be sorted.
    \param In
           NumDoubleCompanions - Number of lists of double precision numbers to be sorted with the key.  If set to zero,
	   DoubleCompanions is ignored and can be set to zero.
    \param In
           DoubleCompanions - DoubleCompanions[i] is a pointer to the ith list of doubles to be sorted with key.
    \param In
           NumIntCompanions - Number of lists of integers to be sorted with the key.  If set to zero, 
	   IntCompanions is ignored and can be set to zero.
    \param In
           IntCompanions - IntCompanions[i] is a pointer to the ith list of integers to be sorted with key.
	   
  */
  void Sort(bool SortAscending, int NumKeys, int * Keys, 
	    int NumDoubleCompanions,double ** DoubleCompanions, 
	    int NumIntCompanions, int ** IntCompanions) const;
  //! Epetra_Util Chop method.  Return zero if input Value is less than ChopValue
  static double Chop(const double & Value){
    if (fabs(Value) < chopVal_) return 0;
    return Value;
  };

  static const double chopVal_;
};

// Epetra_Util constructor
inline Epetra_Util::Epetra_Util(void){}
// Epetra_Util constructor
inline Epetra_Util::Epetra_Util(const Epetra_Util& Util){}
// Epetra_Util destructor
inline Epetra_Util::~Epetra_Util(){}

/** Utility function to perform a binary-search on a list of data.
    Important assumption: data is assumed to be sorted.

    @param item to be searched for
    @param list to be searched in
    @param len Length of list
    @param insertPoint Input/Output. If item is found, insertPoint is not
    referenced. If item is not found, insertPoint is set to the offset at which
    item should be inserted in list such that order (sortedness) would be
    maintained.
    @return offset Location in list at which item was found. -1 if not found.
*/
int Epetra_Util_binary_search(int item,
                              const int* list,
                              int len,
                              int& insertPoint);

/** Function to insert an item in a list, at a specified offset.
    @return error-code 0 if successful, -1 if an allocation failed or if
    input parameters seem unreasonable (offset > usedLength, offset<0, etc).

    @param item to be inserted
    @param offset location at which to insert item
    @param list array into which item is to be inserted. This array may be
           re-allocated by this function.
    @param usedLength number of items already present in list. Will be updated
          to reflect the new length.
    @param allocatedLength current allocated length of list. Will be updated
          to reflect the new allocated-length, if applicable. Re-allocation
          occurs only if usedLength==allocatedLength on entry.
    @param allocChunkSize Optional argument, defaults to 1000. Increment by
          which the array should be expanded, if re-allocation is necessary.
    @return error-code 0 if successful. -1 if allocation fails or if input
         parameters don't make sense.
 */
template<class T>
int Epetra_Util_insert(T item, int offset, T*& list,
                        int& usedLength,
                        int& allocatedLength,
                        int allocChunkSize=1000)
{
  if (offset < 0 || offset > usedLength) {
    return(-1);
  }

  if (usedLength < allocatedLength) {
    for(int i=usedLength; i>offset; --i) {
      list[i] = list[i-1];
    }
    list[offset] = item;
    ++usedLength;
    return(0);
  }

  T* newlist = new T[allocatedLength+allocChunkSize];
  if (newlist == NULL) {
    return(-1);
  }

  allocatedLength += allocChunkSize;
  int i;
  for(i=0; i<offset; ++i) {
    newlist[i] = list[i];
  }

  newlist[offset] = item;

  for(i=offset+1; i<=usedLength; ++i) {
    newlist[i] = list[i-1];
  }

  ++usedLength;
  delete [] list;
  list = newlist;
  return(0);
}

//! Harwell-Boeing data extraction routine
/*! This routine will extract data from an existing Epetra_Crs Matrix, and
    optionally from related rhs and lhs objects in a form that is compatible with
    software that requires the Harwell-Boeing data format. The matrix must be passed
    in, but the RHS and LHS arguments may be set to zero (either or both of them).
    For each of the LHS or RHS arguments, if non-trivial and contain more than one vector, the
    vectors must have strided access.  If both LHS and RHS are non-trivial, they must have the
    same number of vectors.  If the input objects are distributed, the returned matrices will 
    contain the local part of the matrix and vectors only.

    \param A (In) Epetra_CrsMatrix.
    \param LHS (In) Left hand side multivector.  Set to zero if none not available or needed.
    \param RHS (In) Right hand side multivector.  Set to zero if none not available or needed.
    \param M (Out) Local row dimension of matrix.
    \param N (Out) Local column dimension of matrix.
    \param nz (Out) Number of nonzero entries in matrix.
    \param ptr (Out) Offsets into ind and val arrays pointing to start of each row's data.
    \param ind (Out) Column indices of the matrix, in compressed form.
    \param val (Out) Matrix values, in compressed form corresponding to the ind array.
    \param Nrhs (Out) Number of right/left hand sides found (if any) in RHS and LHS.
    \param rhs (Out) Fortran-style 2D array of RHS values.
    \param ldrhs (Out) Stride between columns of rhs.
    \param lhs (Out) Fortran-style 2D array of LHS values.
    \param ldrhs (Out) Stride between columns of lhs.
*/
int Epetra_Util_ExtractHbData(Epetra_CrsMatrix * A, Epetra_MultiVector * LHS,
			      Epetra_MultiVector * RHS,
			      int & M, int & N, int & nz, int * & ptr,
			      int * & ind, double * & val, int & Nrhs,
			      double * & rhs, int & ldrhs,
			      double * & lhs, int & ldlhs);


#endif /* _EPETRA_UTIL_H_ */
