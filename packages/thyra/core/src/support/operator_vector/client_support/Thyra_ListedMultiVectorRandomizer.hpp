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
