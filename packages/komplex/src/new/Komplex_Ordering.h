//@HEADER
/*
************************************************************************

              Komplex: Complex Linear Solver Package 
                Copyright (2002) Sandia Corporation

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

#ifndef KOMPLEX_ORDERING_H
#define KOMPLEX_ORDERING_H


class Epetra_BlockMap;
class Komplex_MultiVector;

#include "Komplex_KForms.hpp"
#include "Epetra_Vector.h"

//! Komplex_Ordering: A class for manipulating the KForm of various Komplex objects.

/*! The Komplex_Ordering class aids other Komplex classes in switching from one KForm to another with minimal amounts
  of swapping.
*/

//==========================================================================
class Komplex_Ordering {
  
public:
  
  //@{ \name Constructor/destructor.
  //! Basic Komplex_Ordering constuctor.
  /*! Creates a Komplex_Ordering object.    
    \param Map (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
  
    \warning Note that, because Epetra_LocalMap
    derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
    for all three types of Epetra map classes.
  
    \param KForm (In) The Komplex_KForms to use.
    \param IsOneObject (In) If true, this ordering is for a single object, real and imaginary values interleaved.

    \return Pointer to a Komplex_Ordering object.  
  */
  Komplex_Ordering(const Epetra_BlockMap& Map, Komplex_KForms KForm, bool IsOneObject);
  
  //! Komplex_Ordering copy constructor.
  /*! Creates a Komplex_Ordering object from a pre-existing one.    
    \param Source (In) A fully constructed Komplex_Ordering object.
  
    \return Pointer to a Komplex_Ordering object.  
  */
  Komplex_Ordering(Komplex_Ordering& Source);
  
  //! Komplex_Ordering destructor.
  virtual ~Komplex_Ordering(void);
  //@}
  
  //@{ \name Attribute access and modification functions.
  
  //! Returns the current K form.
  Komplex_KForms KForm(void);
  
  //! Switches the current K form.
  /*!
    \param NewKForm (In) The new KForms to use.

    \return Integer error code, set to 0 if successful.
  */
  int SwitchKForm(Komplex_KForms NewKForm);
  //@}
  
  //@{ \name Vector and element access functions.
  /*! P vector access function
    \param Perms (Out) Pointer to memory space that will contain the values of P.

    \return Integer error code, set to 0 if successful.
  */
  int PermutationVector(int* Perms);
  
  /*! D vector access function
    \param Scales (Out) Pointer to memory space that will contain the values of D.

    \return Integer error code, set to 0 if successful.
  */
  int ScalingVector(double* Scales); 
  
  /*! Global element in P access function
    \param GlobalRow (In) Array row to be returned.
    \param Index (Out) Integer code, 1 meaning GlobalRow is the TrueRow for a one-object object 
    or 1 meaning GlobalRow/2 in the Real object;
    -1 meaning the preceding or following row for a one-object object 
    or -1 meaning GlobalRow/2 in the Imag object.
    
    \return Integer error code
  */
  int GlobalIndex(int GlobalRow, int& Index);
  
  /*! Global element in D access function
    \param GlobalRow (In) Array row to be returned.
    \param Scalar (Out) Double address to return the value.
    
    \return Integer error code, set to 0 if successful.
  */
  int GlobalScaling(int GlobalRow, double& Scalar);
  
  /*! Local element in P access function
    \param MyRow (In) Array row to be returned.
    \param Index (Out) Integer code, 1 meaning MyRow is the TrueRow for a one-object object 
    or 1 meaning MyRow/2 in the Real object;
    -1 meaning the preceding or following row for a one-object object 
    or -1 meaning MyRow/2 in the Imag object.
    
    \return Integer error code
  */
  int MyIndex(int MyRow, int& Index);
  
  /*! Local element in D access function
    \param MyRow (In) Array row to be returned.
    \param Scalar (Out) Double address to return the value.
    
    \return Integer error code, set to 0 if successful.
  */
  int MyScaling(int MyRow, double& Scalar);
  //@}
  
  //@{ \name Expert-only unsupported modification routines. //is it really unsupported??
 
/* 
  //! Update the given Komplex_MultiVector, resetting P and D in the process. */
  /*!
    \param Source (In) Fully constructed Komplex_MultiVector that owns the \e this Komplex_Ordering.
    
    \return Integer error code, set to 0 if successful.
  */
  /*int Update(Komplex_MultiVector& Source);           
*/

  //! Reset the values of P_ and D_ to their original state and set KForm_ to NewKForm
  void Reset(Komplex_KForms NewKForm); //should it return anything???

  //@}
  
  //@{ \name Overloaded operators.
  
  //! = Operator.
  /*!
    \param Source (In) Komplex_Operator to copy.
    
    \return Komplex_Operator.
  */
  //########## Komplex_Ordering & operator = (const Komplex_Ordering& Source);
  //@}
  
protected:
  
private:
  Epetra_Vector P_;
  Epetra_Vector D_;
  Komplex_KForms KForm_;
  Komplex_KForms StoredKForm_;
  bool IsOneObject_;
  
};

#endif /* KOMPLEX_ORDERING_H */
