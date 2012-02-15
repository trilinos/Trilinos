/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
*/

#ifndef OPTIPACK_LINE_SEARCH_POINT_EVALUATOR_BASE_HPP
#define OPTIPACK_LINE_SEARCH_POINT_EVALUATOR_BASE_HPP


#include "OptiPack_Types.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Describable.hpp"


namespace OptiPack {


/** \brief Base class interface for line search point updates.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class LineSearchPointEvaluatorBase
  : public Teuchos::Describable
{
public:

  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Compute the updated point <tt>p</tt> at <tt>alpha</tt> for a
   * linear search algorithm.
   */
  virtual void computePoint( const ScalarMag &alpha,
    const Ptr<Thyra::VectorBase<Scalar> > &p
    ) const = 0;

};


} // namespace OptiPack


#endif // OPTIPACK_LINE_SEARCH_POINT_EVALUATOR_BASE_HPP
