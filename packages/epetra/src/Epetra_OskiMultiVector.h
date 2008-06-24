
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

#ifndef EPETRA_OSKIMULTIVECTOR_H
#define EPETRA_OSKIMULTIVECTOR_H

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
