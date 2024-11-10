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

#ifndef EPETRA_OSKIMULTIVECTOR_H
#define EPETRA_OSKIMULTIVECTOR_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_MultiVector.h"
extern "C" {
#include "oski/oski.h"
}

//class Epetra_MultiVector;

//! Epetra_OskiMultiVector: A class for constructing and using dense Oski multi-vectors on a single processor or a single core of a multi-processor.

/*! The Epetra_OskiMultiVector class enables the construction and use of real-valued,
  double-precision dense vectors and multi-vectors,
  in a serial environment.  The dimensions of the dense multi-vectors comes
  from the inherited Epetra_MultiVector object.  All values and data layouts are kept the
  same and the multi-vector is wrapped for use with OSKI.
*/

//==========================================================================
class Epetra_OskiMultiVector: public Epetra_MultiVector{

 public:

   //! @name Constructors/Destructor
   //@{
   //! Copy constructor
   Epetra_OskiMultiVector(const Epetra_OskiMultiVector& Source);

   //! Constructor creates and Epetra_OskiMultiVector from an Epetra_MultiVector
   /*! \param Source (In) An Epetra_MultiVector that is wrapped as an Epetra_OskiMultiVector.
       \return Pointer to an Epetra_OskiMultiVector.
       \note If the Epetra_MultiVector is not stored contigously according to the
	     BLAS standard then a deep copy is made.
   */
   Epetra_OskiMultiVector(const Epetra_MultiVector& Source);

   //! Destructor
   virtual ~Epetra_OskiMultiVector();
   //@}

   //! @name Extraction Methods
   //@{
   //! Returns true if a deep copy of the multi-vector was created by the constructor.
   bool Copy_Created() const;

   //! Returns the Oski portion of the Multi-Vector.
   oski_vecview_t Oski_View() const;

   //! Returns the Epetra portion of the Multi-Vector.
   const Epetra_MultiVector* Epetra_View() const;
   //@}

   //! @name Operators
   //@{
   //! Sets this equal to Source.
   Epetra_OskiMultiVector& operator = (const Epetra_OskiMultiVector& Source);
   //@}

 protected:

 private:
   const Epetra_MultiVector* Epetra_View_;
   oski_vecview_t Oski_View_;
   bool Copy_Created_;
};
#endif /* EPETRA_OSKIMULTIVECTOR_H */
