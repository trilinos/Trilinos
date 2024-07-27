// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LISTED_MULTI_VECTOR_RANDOMIZER_HPP
#define THYRA_LISTED_MULTI_VECTOR_RANDOMIZER_HPP

#include "Thyra_MultiVectorRandomizerBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"


namespace Thyra {


/** \brief <tt>MultiVectorRandomizerBase</tt> subclass that returns a revolving
 * list of preset <tt>MultiVectorBase</tt> objects.
 *
 * This class simply returns a preset list of <tt>MultiVectorBase</tt> objects
 * instead of true random multi-vectors.  This can be very useful when
 * combined with testing software.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class ListedMultiVectorRandomizer : public MultiVectorRandomizerBase<Scalar> {
public:

  /** \brief Calls <tt>this->initialize()</tt>. */
  ListedMultiVectorRandomizer(
    const Teuchos::RCP<const MultiVectorBase<Scalar> >    multiVecs[]
    ,const int                                                    numMultiVecs
    );

  /** \brief . */
  void initialize(
    const Teuchos::RCP<const MultiVectorBase<Scalar> >    multiVecs[]
    ,const int                                                    numMultiVecs
    );

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

private:

  typedef std::vector<Teuchos::RCP<const MultiVectorBase<Scalar> > > multiVecs_t;
  multiVecs_t  multiVecs_;
  int          curr_mv_i_;
  
};


// //////////////////////////////
// Definitions


template<class Scalar>
ListedMultiVectorRandomizer<Scalar>::ListedMultiVectorRandomizer(
  const Teuchos::RCP<const MultiVectorBase<Scalar> >    multiVecs[]
  ,const int                                                    numMultiVecs
  )
{
  initialize(multiVecs,numMultiVecs);
}


template<class Scalar>
void ListedMultiVectorRandomizer<Scalar>::initialize(
  const Teuchos::RCP<const MultiVectorBase<Scalar> >    multiVecs[]
  ,const int                                                    numMultiVecs
  )
{
  multiVecs_.resize(numMultiVecs);
  std::copy( multiVecs, multiVecs + numMultiVecs, multiVecs_.begin() );
  curr_mv_i_ = 0;
}


// Overridden from MultiVectorRandomizerBase


template<class Scalar>
bool ListedMultiVectorRandomizer<Scalar>::isCompatible( const VectorSpaceBase<Scalar> &space ) const
{
  return multiVecs_[curr_mv_i_]->range()->isCompatible(space);
}


// Overridded private functions


template<class Scalar>
void ListedMultiVectorRandomizer<Scalar>::randomizeImpl(
  const Ptr<MultiVectorBase<Scalar> > &mv)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( multiVecs_.size()==0 );
#endif
  const Teuchos::RCP<const MultiVectorBase<Scalar> > currMV = multiVecs_[curr_mv_i_];
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_VEC_SPACES("ListedMultiVectorRandomizer<Scalar>::randomize(mv)", *currMV->range(), *mv->range() );
  THYRA_ASSERT_VEC_SPACES("ListedMultiVectorRandomizer<Scalar>::randomize(mv)", *currMV->domain(), *mv->domain() );
#endif
  Thyra::assign( mv, *currMV );
  if( curr_mv_i_ == static_cast<int>(multiVecs_.size()) - 1 )
    curr_mv_i_ = 0;
  else
    ++curr_mv_i_;
}


} // namespace Thyra


#endif // THYRA_LISTED_MULTI_VECTOR_RANDOMIZER_HPP
