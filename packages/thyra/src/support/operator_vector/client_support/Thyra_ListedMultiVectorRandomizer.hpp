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
  /** \brief . */
  void randomize( MultiVectorBase<Scalar> *mv );

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

template<class Scalar>
void ListedMultiVectorRandomizer<Scalar>::randomize( MultiVectorBase<Scalar> *mv )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( mv==NULL );
  TEST_FOR_EXCEPT( multiVecs_.size()==0 );
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
