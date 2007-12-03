
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRA_FEVECTOR_H
#define EPETRA_FEVECTOR_H

#include <Epetra_CombineMode.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseVector;

/** Epetra Finite-Element Vector. This class inherits Epetra_MultiVector
  and thus provides all Epetra_MultiVector functionality.

  The added functionality provided by Epetra_FEVector is the ability to
  perform finite-element style vector assembly. It accepts sub-vector
  contributions, such as those that would come from element-load vectors, etc.,
  and these sub-vectors need not be owned by the local processor. In other
  words, the user can assemble overlapping data (e.g., corresponding to shared
  finite-element nodes). When the user is finished assembling their vector
  data, they then call the method Epetra_FEVector::GlobalAssemble() which
  gathers the overlapping data (all non-local data that was input on each
  processor) into the data-distribution specified by the map that the
  Epetra_FEVector is constructed with.
*/

class Epetra_FEVector : public Epetra_MultiVector {
 public:
   /** Constructor that requires a map specifying a non-overlapping
      data layout.
      @param Map Map describing a non-overlapping distribution for the
         underlying Epetra_MultiVector that this Epetra_FEVector will
         be funnelling data into.

      @param numVectors Optional argument, default value is 1. (See the
           documentation for Epetra_MultiVector for the meaning of this
          argument.

      @param ignoreNonLocalEntries Optional argument, default value is false.
             Under certain special circumstances it is desirable to have
             non-local contributions ignored rather than saving them for
             the GlobalAssemble step.
   */
   Epetra_FEVector(const Epetra_BlockMap& Map,
                   int numVectors=1,
		   bool ignoreNonLocalEntries=false);

  /** Copy constructor. */
  Epetra_FEVector(const Epetra_FEVector& source);

   /** Destructor */
   virtual ~Epetra_FEVector();

   /** Accumulate values into the vector, adding them to any values that
       already exist for the specified indices.
   */
   int SumIntoGlobalValues(int numIDs,
                           const int* GIDs, const double* values,
                           int vectorIndex=0);

   /** Accumulate values into the vector, adding them to any values that
       already exist for the specified GIDs.

       @param GIDs List of global ids. Must be the same length as the
       accompanying list of values.

       @param values List of coefficient values. Must be the same length as
       the accompanying list of GIDs.
   */
   int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
			   const Epetra_SerialDenseVector& values,
                           int vectorIndex=0);

   /** Copy values into the vector overwriting any values that already exist
        for the specified indices.
    */
   int ReplaceGlobalValues(int numIDs, const int* GIDs, const double* values,
                           int vectorIndex=0);

   /** Copy values into the vector, replacing any values that
       already exist for the specified GIDs.

       @param GIDs List of global ids. Must be the same length as the
       accompanying list of values.

       @param values List of coefficient values. Must be the same length as
       the accompanying list of GIDs.
   */
   int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
			   const Epetra_SerialDenseVector& values,
                           int vectorIndex=0);

   int SumIntoGlobalValues(int numIDs,
                           const int* GIDs,
			   const int* numValuesPerID,
			   const double* values,
                           int vectorIndex=0);

   int ReplaceGlobalValues(int numIDs, const int* GIDs,
			   const int* numValuesPerID,
			   const double* values,
                           int vectorIndex=0);

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this vector at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.
   */
   int GlobalAssemble(Epetra_CombineMode mode = Add);

   /** Set whether or not non-local data values should be ignored.
    */
   void setIgnoreNonLocalEntries(bool flag) {
     ignoreNonLocalEntries_ = flag;
   }

   Epetra_FEVector& operator=(const Epetra_FEVector& source);

 private:
  int inputValues(int numIDs,
                  const int* GIDs, const double* values,
                  bool suminto,
                  int vectorIndex);

  int inputValues(int numIDs,
                  const int* GIDs, const int* numValuesPerID,
		  const double* values,
                  bool suminto,
                  int vectorIndex);

  int inputNonlocalValue(int GID, double value, bool suminto,
                         int vectorIndex);

  int inputNonlocalValues(int GID, int numValues, const double* values,
			  bool suminto, int vectorIndex);

  void destroyNonlocalData();

  int myFirstID_;
  int myNumIDs_;

  int* nonlocalIDs_;
  int* nonlocalElementSize_;
  int numNonlocalIDs_;
  int numNonlocalIDsAlloc_;
  double** nonlocalCoefs_;
  int numNonlocalCoefs_;
  int numNonlocalCoefsAlloc_;

  bool ignoreNonLocalEntries_;
};

#endif

