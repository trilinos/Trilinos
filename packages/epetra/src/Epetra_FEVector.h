
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

#ifndef _EPETRA_FEVECTOR_H_
#define _EPETRA_FEVECTOR_H_

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>

/** Epetra Finite-Element Vector. This class inherits Epetra_MultiVector
  and thus provides all Epetra_MultiVector functionality, with one
  restriction: currently an Epetra_FEVector only has 1 internal vector.

  The added functionality provided by Epetra_FEVector is the ability to
  perform finite-element style vector assembly. It accepts sub-vector
  contributions, such as those that would come from element-load vectors, etc.,
  and these sub-vectors need not be wholly locally owned. In other words, the
  user can assemble overlapping data (e.g., corresponding to shared
  finite-element nodes). When the user is finished assembling their vector
  data, they then call the method Epetra_FEVector::GlobalAssemble() which
  gathers the overlapping data (all non-local data that was input on each
  processor) into the data-distribution specified by the map that the
  Epetra_FEVector is constructed with.

  Note: At the current time (Sept 6, 2002) the methods in this implementation
  assume that there is only 1 point associated with each map element. This 
  limitation will be removed in the near future.
*/

class Epetra_FEVector : public Epetra_MultiVector {
 public:
   /** Constructor that requires a map specifying a non-overlapping
      data layout. The methods SumIntoGlobalValues() and 
      ReplaceGlobalValues() will accept any global IDs, and GlobalAssemble()
      will move any non-local data onto the appropriate owning processors.
   */
   Epetra_FEVector(const Epetra_BlockMap& Map,
		   bool ignoreNonLocalEntries=false);

   /** Destructor */
   virtual ~Epetra_FEVector();

   /** Accumulate values into the vector, adding them to any values that
       already exist for the specified indices.
   */
   int SumIntoGlobalValues(int numIDs, const int* GIDs, const double* values);

   /** Copy values into the vector overwriting any values that already exist
        for the specified indices.
    */
   int ReplaceGlobalValues(int numIDs, const int* GIDs, const double* values);

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this vector at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.
   */
   int GlobalAssemble();

   /** Set whether or not non-local data values should be ignored.
    */
   void setIgnoreNonLocalEntries(bool flag) {
     ignoreNonLocalEntries_ = flag;
   }

   Epetra_FEVector& operator=(const Epetra_FEVector& source);

 private:
  int inputValues(int numIDs,
                  const int* GIDs, const double* values,
                  bool accumulate);

  int inputNonlocalValue(int GID, double value, bool accumulate);

  void destroyNonlocalData();

  int myFirstID_;
  int myNumIDs_;
  double* myCoefs_;

  int* nonlocalIDs_;
  int numNonlocalIDs_;
  int allocatedNonlocalLength_;
  double* nonlocalCoefs_;

  bool ignoreNonLocalEntries_;
};

#endif

