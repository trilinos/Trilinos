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
