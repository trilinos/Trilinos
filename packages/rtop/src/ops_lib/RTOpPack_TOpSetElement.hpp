// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_SET_ELEMENT_HPP
#define RTOPPACK_TOP_SET_ELEMENT_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Teuchos_as.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation for TOpSetElement. */
template<class Scalar>
class TOpSetElementEleWiseTransformation
{
public:
  /** \brief . */
  TOpSetElementEleWiseTransformation( const Ordinal &global_i_in = -1,
    const Scalar &val_i_in = static_cast<Scalar>(0.0) )
    :global_i_(global_i_in), val_i_(val_i_in)
    {}
  /** \brief . */
  Ordinal global_i() const
    {
      return global_i_;
    }
  /** \brief . */
  void operator()( const Ordinal global_i_in, Scalar &z0 ) const
    {
      if (global_i_in == global_i_) {
        z0 = val_i_;
      }
    }
private:
  Ordinal global_i_;
  Scalar val_i_;
};


/** \brief Set the elements of a vector to: <tt>z0[i] = i+global_i+1, i=0...n-1</tt>.
 */
template<class Scalar>
class TOpSetElement
  : public TOp_0_1_CoordVariantBase<Scalar, TOpSetElementEleWiseTransformation<Scalar> >
{
public:
  /** \brief . */
  TOpSetElement(const Ordinal &global_i_in = -1,
    const Scalar &val_i_in = static_cast<Scalar>(0.0))
    {
      this->setOpNameBase("TOpSetElement");
      this->setEleWiseTransformation(
        TOpSetElementEleWiseTransformation<Scalar>(global_i_in, val_i_in));
    }
  /** \brief . */
  void initialize(const Ordinal &global_i_in, const Scalar &val_i_in)
    { 
      this->setEleWiseTransformation(
        TOpSetElementEleWiseTransformation<Scalar>(global_i_in, val_i_in));
    }
protected:
  /** \brief . */
  virtual Range1D range_impl() const
    {
      const Ordinal i = this->getEleWiseTransformation().global_i();
      return Range1D(i, i);
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_SET_ELEMENT_HPP
