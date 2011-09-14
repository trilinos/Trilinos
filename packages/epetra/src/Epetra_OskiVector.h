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

#ifndef EPETRA_OSKIVECTOR_H
#define EPETRA_OSKIVECTOR_H

#include "Epetra_OskiMultiVector.h"
#include "Epetra_Vector.h"
extern "C" {
#include "oski/oski.h"
}

class Epetra_Vector;
class Epetra_OskiMultiVector;

//! Epetra_OskiVector: A class for constructing and using dense OSKI vectors on a single processor or a single core of a multi-processor.

/*! The Epetra_OskiVector class enables the construction and use of real-valued,
  double-precision dense vectors in a serial environment.  The dimensions of the dense vectors comes
  from the inherited Epetra_Vector object.  All values and data layouts are kept the
  same and the vector is wrapped for use with OSKI.
*/


//=========================================================================
class Epetra_OskiVector: public Epetra_OskiMultiVector {

  public:

   //! @name Constructors/Destructor
   //@{
   //! Copy constructor
   Epetra_OskiVector(const Epetra_OskiVector& Source);

   //! Constructor creates and Epetra_OskiVector from an Epetra_Vector
   /*! \param Source (In) An Epetra_Vector that is to be wrapped as an Epetra_OskiVector.
       \return Pointer to an Epetra_OskiVector.
   */
   Epetra_OskiVector(const Epetra_Vector& Source);

   //! Destructor
   virtual ~Epetra_OskiVector();
   //@}

   //! @name Extraction Method
   //@{
   //! Returns a view to the Epetra Object
   const Epetra_Vector* Epetra_View() const;
   //@}

   //! @name Operators
   //@{
   //! Sets this equal to Source.
   Epetra_OskiVector& operator = (const Epetra_OskiVector& Source);
   //@}

 protected:

 private:
   
    const Epetra_Vector* Epetra_View_;

};

#endif /* EPETRA_OSKIVECTOR_H */
