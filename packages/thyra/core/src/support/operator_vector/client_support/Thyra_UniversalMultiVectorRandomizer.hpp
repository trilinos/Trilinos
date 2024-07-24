// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  //@}

private:

  /** \name Overridded private functions */
  //@{

  /** \brief . */
  void randomizeImpl(const Ptr<MultiVectorBase<Scalar> > &mv);

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
bool UniversalMultiVectorRandomizer<Scalar>::isCompatible( const VectorSpaceBase<Scalar> &/* space */ ) const
{
  return true;
}


// Overridded private functions


template<class Scalar>
void UniversalMultiVectorRandomizer<Scalar>::randomizeImpl(
  const Ptr<MultiVectorBase<Scalar> > &mv )
{
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Thyra::randomize(as<Scalar>(-ST::one()), as<Scalar>(+ST::one()), mv);
}


} // namespace Thyra


#endif // THYRA_UNIVERSAL_MULTI_VECTOR_RANDOMIZER_HPP
