// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_TRANSFORMER_BASE_HPP
#define THYRA_LINEAR_OP_TRANSFORMER_BASE_HPP


#include "Thyra_LinearOpBase.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Thyra {


/** \brief Base interface for transforming a LinearOpBase object. */
template<class Scalar>
class LinearOpTransformerBase
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<LinearOpTransformerBase<Scalar> >
{
public:

  /** \brief Create an uninitialized op. */
  virtual bool isCompatible(const LinearOpBase<Scalar> &op_in) const = 0;

  /** \brief Create an uninitialized op. */
  virtual RCP<LinearOpBase<Scalar> > createOutputOp() const = 0;

  /** \brief Do the transformation to a pre-created output LinearOpBase object.
   *
   * \param op_in [in] The linear operator source that will be transformed in
   * some way.  Precondition: <tt>this->isCompataible(op_in) == true</tt>.
   *
   * \param op_inout [in/out] The transformed linear operator.  This object
   * must have been created by <tt>this->createOutputOp()</tt> and may have
   * already been passed through this function before.  This allows for resuse
   * of internal structures on re-transformations.  Postcondition: On output,
   * the object <tt>op_inout</tt> will be some appropriate transformation of
   * <tt>op_in</tt>.  The exact nature of the transformation is not specified
   * in this interface.
   */
  virtual void transform(
    const LinearOpBase<Scalar> &op_in,
    const Ptr<LinearOpBase<Scalar> > &op_inout
    ) const = 0;

};


} // namespace Thyra


#endif	// THYRA_LINEAR_OP_TRANSFORMER_BASE_HPP
