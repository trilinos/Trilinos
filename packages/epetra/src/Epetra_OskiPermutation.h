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

// Author: Ian Karlin ikarlin@sandia.gov 05-28-2008

#ifndef EPETRA_OSKIPERMUTATION_H
#define EPETRA_OSKIPERMUTATION_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



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
