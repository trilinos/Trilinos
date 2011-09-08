// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
