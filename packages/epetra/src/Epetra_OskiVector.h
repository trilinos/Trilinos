
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
