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
