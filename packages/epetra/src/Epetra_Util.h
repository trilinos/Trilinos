/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_UTIL_H
#define EPETRA_UTIL_H

#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include <vector>
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Import;

//! Epetra_Util:  The Epetra Util Wrapper Class.
/*! The Epetra_Util class is a collection of useful functions that cut across a broad
  set of other classes.
<ul>
<li> A random number generator is provided, along with methods to set and
retrieve the random-number seed.

The random number generator is a multiplicative linear congruential generator, 
with multiplier 16807 and modulus 2^31 - 1. It is based on the algorithm described in
"Random Number Generators: Good Ones Are Hard To Find", S. K. Park and K. W. Miller, 
Communications of the ACM, vol. 31, no. 10, pp. 1192-1201.

<li> Sorting is provided by a static function on this class (i.e., it is not
necessary to construct an instance of this class to use the Sort function).

<li> A static function is provided for creating a new Epetra_Map object with
 1-to-1 ownership of entries from an existing map which may have entries that
appear on multiple processors.
</ul>

  Epetra_Util is a serial interface only.  This is appropriate since the standard 
  utilities are only specified for serial execution (or shared memory parallel).
*/
class EPETRA_LIB_DLL_EXPORT Epetra_Util {
    
  public:
  //! Epetra_Util Constructor.
  /*! Builds an instance of a serial Util object.
   */
  Epetra_Util();


  //! Epetra_Util Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_Util instance.
  */
  Epetra_Util(const Epetra_Util& Util);

  //! Epetra_Util Destructor.
  virtual ~Epetra_Util();

  //! @name Random number utilities
  //@{ 

  //! Returns a random integer on the interval (0, 2^31-1)
  unsigned int RandomInt();

  //! Returns a random double on the interval (-1.0,1.0)
  double RandomDouble();

  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  unsigned int Seed() const;

  //! Set seed for Random function.
  /*!
    \param In
    Seed - An integer on the interval [1, 2^31-2]

    \return Integer error code, set to 0 if successful.
  */
  int SetSeed(unsigned int Seed_in);

	//@}
  
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
  template<typename T>
  static void Sort(bool SortAscending, int NumKeys, T * Keys, 
		   int NumDoubleCompanions,double ** DoubleCompanions, 
		   int NumIntCompanions, int ** IntCompanions,
		   int NumLongLongCompanions, long long ** LongLongCompanions);

  static void Sort(bool SortAscending, int NumKeys, int * Keys, 
		   int NumDoubleCompanions,double ** DoubleCompanions, 
		   int NumIntCompanions, int ** IntCompanions,
		   int NumLongLongCompanions, long long ** LongLongCompanions);

  static void Sort(bool SortAscending, int NumKeys, long long * Keys, 
		   int NumDoubleCompanions,double ** DoubleCompanions, 
		   int NumIntCompanions, int ** IntCompanions,
		   int NumLongLongCompanions, long long ** LongLongCompanions);

  static void Sort(bool SortAscending, int NumKeys, int * Keys, 
		   int NumDoubleCompanions,double ** DoubleCompanions, 
		   int NumIntCompanions, int ** IntCompanions);

  //! Epetra_Util Create_Root_Map function
  /*! Function to create a new Epetra_Map object with all GIDs sent to the root processor
      which is zero by default.  All all processors will have no GIDs.  This root map can then 
      be used to create an importer or exporter that will migrate all data to the root processor.

      If root is set to -1 then the user map will be replicated completely on all processors.
  */
  static Epetra_Map Create_Root_Map(const Epetra_Map & usermap,
					int root = 0);

  //! Epetra_Util Create_OneToOne_Map function
  /*! Function to create a new Epetra_Map object with 1-to-1 ownership of
    entries from an existing map which may have entries that appear on
    multiple processors.
  */
  static Epetra_Map Create_OneToOne_Map(const Epetra_Map& usermap,
					bool high_rank_proc_owns_shared=false);

  //! Epetra_Util Create_OneToOne_Map function
  /*! Function to create a new Epetra_Map object with 1-to-1 ownership of
    entries from an existing map which may have entries that appear on
    multiple processors.
  */
  static Epetra_BlockMap Create_OneToOne_BlockMap(const Epetra_BlockMap& usermap,
						  bool high_rank_proc_owns_shared=false);


  //! Epetra_Util GetPidGidPairs function
  /*!  For each GID in the TargetMap, find who owns the GID in the SourceMap. 
    This works entirely from the Distributor and has no communication at all.  
    This routine only works if your Importer is using an Epetra_MpiDistributor under the hood.
    
    The routine returns (by reference) a std::vector of std::pair<int,int> which contains (PID,GID) pairs.
    If the use_minus_one_for_local==true, any GIDs owned by this processor get -1 instead of their PID.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static int GetPidGidPairs(const Epetra_Import & Importer,std::vector< std::pair<int,int> > & gpids, bool use_minus_one_for_local);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  static int GetPidGidPairs(const Epetra_Import & Importer,std::vector< std::pair<int,long long> > & gpids, bool use_minus_one_for_local);
#endif
  

  //! Epetra_Util GetPids function
  /*! Like GetPidGidPairs, but just gets the PIDs, ordered by the columnmap 
   */
  int GetPids(const Epetra_Import & Importer, std::vector<int> &pids, bool use_minus_one_for_local);

  //! Epetra_Util Chop method.  Return zero if input Value is less than ChopValue
  static double Chop(const double & Value);
//  {
//    if (std::abs(Value) < chopVal_) return 0;
//    return Value;
//  };

  static const double chopVal_;

 private:
	unsigned int Seed_;
};


// Epetra_Util constructor
inline Epetra_Util::Epetra_Util() : Seed_(std::rand()) {}
// Epetra_Util constructor
inline Epetra_Util::Epetra_Util(const Epetra_Util& Util) : Seed_(Util.Seed_) {}
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
template<typename T>
int Epetra_Util_binary_search(T item,
                              const T* list,
                              int len,
                              int& insertPoint);

EPETRA_LIB_DLL_EXPORT int Epetra_Util_binary_search(int item,
                              const int* list,
                              int len,
                              int& insertPoint);

EPETRA_LIB_DLL_EXPORT int Epetra_Util_binary_search(long long item,
                              const long long* list,
                              int len,
                              int& insertPoint);

template<class T>
int Epetra_Util_insert_empty_positions(T*& array, int& usedLength,
                                       int& allocatedLength,
                                       int insertOffset, int numPositions,
                                       int allocChunkSize=32)
{
  if (insertOffset < 0 || insertOffset > usedLength ||
      usedLength > allocatedLength) {
    return(-1);
  }

  if ((usedLength+numPositions) < allocatedLength) {
    for(int i=usedLength-1; i>=insertOffset; --i) {
      array[i+numPositions] = array[i];
    }
    usedLength += numPositions;
    return(0);
  }

  allocatedLength += allocChunkSize;
  //what if allocatedLength is still not large enough?
  //loop until it is large enough:
  while(allocatedLength < usedLength+numPositions) {
    allocatedLength += allocChunkSize;
  }

  T* newlist = new T[allocatedLength];

  for(int i=0; i<insertOffset; ++i) {
    newlist[i] = array[i];
  }

  for(int i=insertOffset; i<usedLength; ++i) {
    newlist[i+numPositions] = array[i];
  }

  usedLength += numPositions;
  delete [] array;
  array = newlist;
  return(0);
}

/** Function to insert an item in a list, at a specified offset.

    @param item to be inserted
    @param offset location at which to insert item
    @param list array into which item is to be inserted. This array may be
           re-allocated by this function.
    @param usedLength number of items already present in list. Will be updated
          to reflect the new length.
    @param allocatedLength current allocated length of list. Will be updated
          to reflect the new allocated-length, if applicable. Re-allocation
          occurs only if usedLength==allocatedLength on entry.
    @param allocChunkSize Optional argument, defaults to 32. Increment by
          which the array should be expanded, if re-allocation is necessary.
    @return error-code 0 if successful. -1 if input parameters don't make sense.
 */
template<class T>
int Epetra_Util_insert(T item, int offset, T*& list,
                        int& usedLength,
                        int& allocatedLength,
                        int allocChunkSize=32)
{
  int code = Epetra_Util_insert_empty_positions<T>(list, usedLength,
                                                allocatedLength, offset, 1,
                                                allocChunkSize);
  if (code != 0) {
    return(code);
  }

  list[offset] = item;

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


#endif /* EPETRA_UTIL_H */
