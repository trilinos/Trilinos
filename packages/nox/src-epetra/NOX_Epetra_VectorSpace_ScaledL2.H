// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_VECTORSPACE_SCALEDL2_H
#define NOX_EPETRA_VECTORSPACE_SCALEDL2_H

#include "NOX_Epetra_VectorSpace.H" // base class
#include "Teuchos_RefCountPtr.hpp"
#include "NOX_Epetra_Scaling.H"

// Forward declarations
class Epetra_Vector;

namespace NOX {

namespace Epetra {

  /*! \brief Concrete class for a weighted L2 vector space.
    
  This class allows users to override the inner product and norm
  used by the NOX::Epetra::Vector class.  The most frequent use of
  this class is for introducing a weighted norm throughout NOX.
  
  */
  class VectorSpaceScaledL2 : public NOX::Epetra::VectorSpace {
    
  public:
    
    //! Constructor
    VectorSpaceScaledL2(const Teuchos::RefCountPtr<NOX::Epetra::Scaling>& s,
			NOX::Epetra::Scaling::ScaleType st = 
			NOX::Epetra::Scaling::Left);
    
    //! Destructor
    virtual ~VectorSpaceScaledL2();
    
    //! Computes a scaled inner product.
    /*! Computes a scaled inner product: \f$ <Da, Db> \f$ where
      \f$D\f$ is the set of scaling vectors associated with either
      left of right scaling.
    */
    virtual double innerProduct(const Epetra_Vector& a, 
				const Epetra_Vector& b) const;
    
    //! Computes the scaled norm.
    /*! Computes the scaled norm using \f$ Da \f$ where \f$D\f$ is the
      set of scaling vectors associated with either left of right
      scaling.
    */
    virtual double norm(const Epetra_Vector& a,
			NOX::Abstract::Vector::NormType = 
			NOX::Abstract::Vector::TwoNorm) const;
    
  protected:
    
    //! Scaling vector used in the inner product
    Teuchos::RefCountPtr<NOX::Epetra::Scaling> scalingPtr;

    //! Scaling type to apply to vector space
    NOX::Epetra::Scaling::ScaleType scaleType;

    //! Temporary vector used in scaling computations
    mutable Teuchos::RefCountPtr<Epetra_Vector> tmpVectorPtr;

  };
} // namespace Epetra
} // namespace NOX

#endif
