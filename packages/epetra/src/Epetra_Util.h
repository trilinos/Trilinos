
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
           NumDoubleCompanions - Number of lists of double precision numbers to be sorted with the key.
    \param In
           DoubleCompanions - DoubleCompanions[i] is a pointer to the ith list of doubles to be sorted with key.
    \param In
           NumIntCompanions - Number of lists of integers to be sorted with the key.
    \param In
           IntCompanions - IntCompanions[i] is a pointer to the ith list of integers to be sorted with key.
	   
  */
  void Sort(bool SortAscending, int NumKeys, int * Keys, int NumDoubleCompanions,double ** DoubleCompanions, 
            int NumIntCompanions, int ** IntCompanions) const;
};

// Epetra_Util constructor
inline Epetra_Util::Epetra_Util(void){}
// Epetra_Util constructor
inline Epetra_Util::Epetra_Util(const Epetra_Util& Util){}
// Epetra_Util destructor
inline Epetra_Util::~Epetra_Util(){}

#endif /* _EPETRA_UTIL_H_ */
