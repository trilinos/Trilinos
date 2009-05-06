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

#ifndef THYRA_UNIVERSAL_MULTI_VECTOR_RANDOMIZER_HPP
#define THYRA_UNIVERSAL_MULTI_VECTOR_RANDOMIZER_HPP


#include "Thyra_MultiVectorRandomizerBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"


namespace Thyra {


/** \brief Univeral <tt>MultiVectorRandomizerBase</tt> subclass that is
 * compatible with all <tt>MultiVectorBase</tt> objects.
 *
 * This class simply uses <tt>randomize(-1,+1,mv)</tt> which is based on RTOp
 * and simply creates random coefficients between -1 and +1.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class UniversalMultiVectorRandomizer : public MultiVectorRandomizerBase<Scalar> {
public:

  /** \name Overridden from MultiVectorRandomizerBase */
  //@{

  /** \brief . */
  bool isCompatible( const VectorSpaceBase<Scalar> &space ) const;
  /** \brief . */
  void randomize( MultiVectorBase<Scalar> *mv );

  //@}
  
};


/** \brief Nonmember constructor.
 *
 * \relates UniversalMultiVectorRandomizer
 */
template<class Scalar>
RCP<UniversalMultiVectorRandomizer<Scalar> >
universalMultiVectorRandomizer()
{
  return Teuchos::rcp(new UniversalMultiVectorRandomizer<Scalar>());
}


// //////////////////////////////
// Definitions


template<class Scalar>
bool UniversalMultiVectorRandomizer<Scalar>::isCompatible( const VectorSpaceBase<Scalar> &space ) const
{
  return true;
}


template<class Scalar>
void UniversalMultiVectorRandomizer<Scalar>::randomize( MultiVectorBase<Scalar> *mv )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Thyra::randomize(Scalar(-ST::one()),Scalar(+ST::one()),mv);
}


} // namespace Thyra


#endif // THYRA_UNIVERSAL_MULTI_VECTOR_RANDOMIZER_HPP
