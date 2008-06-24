
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

// Author: Ian Karlin ikarlin@sandia.gov 05-28-2008

#ifndef EPETRA_OSKIPERMUTATION_H
#define EPETRA_OSKIPERMUTATION_H

#include "Epetra_OskiMatrix.h"
#include "Epetra_OskiMultiVector.h"
extern "C" {
#include "oski/oski.h"
}

class Epetra_OskiMultiVector;
class Epetra_OskiMatrix;

//! Epetra_OskiPermutation: A class for storing the permutation performed on a Epetra_OskiMatrix.
/*! The Epetra_OskiPermutation is one of the possible transformations that OSKI can perform
    on a matrix.  The permutation is stored with the matrix in OSKI.  Using this class a 
    Epetra_OskiPermutation can be applied to an Epetra_OskiMultiVector.
*/

//=========================================================================
class Epetra_OskiPermutation{

  public:

   //! @name Constructor/Destructor
   //@{
   
   //!Default Constructor
   Epetra_OskiPermutation();
   
   //!Copy Constructor
   Epetra_OskiPermutation (const Epetra_OskiPermutation& Source);
   //! Constructor creates an Epetra_OskiPermutation from an Epetra_OskiMatrix.
   /*! Acquires the permutation from the passed in matrix and stores it within the object.
       If RowPerm is true this is a row permutation and if RowPerm is false this is a 
       column permutation.
   */
   Epetra_OskiPermutation(bool RowPerm, const Epetra_OskiMatrix& Source);

   //! Destructor
   virtual ~Epetra_OskiPermutation();
   //@}
   
   //! @name Replace Method
   //@{
   //! Stores a permutation in the data structure.
   void ReplacePermutation(const oski_perm_t& InPerm);
   //@}


   //! @name Apply
   //@{
    //! Permutes Vector according to the Permutation.  If a transpose is desired it performs that operation.
   /*! The Vector passed into this function is a view.  It is the underlying object that is
       permuted by the function.
       \param TransA (In) If TransA = TRUE then use the transpose of the permutation and
              apply it to Vector.
       \param Vector (In/Out) The vector that is permuted by Permutation^Trans.
       \return When successful returns 0.  On error Vector is not permuted and a error
	       code is returned.
   */
   int PermuteVector(const bool TransA, Epetra_OskiMultiVector& Vector) const;
   //@}

 protected:

 private:

 const oski_perm_t* Permutation_; //this might need to be const and a pointer

};

#endif /* EPETRA_OSKIPERMUTATION_H */
