//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_EXTRACT_STATE_AND_SENS_HPP
#define RYTHMOS_EXTRACT_STATE_AND_SENS_HPP


#include "Rythmos_Types.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_ProductVectorBase.hpp"


namespace Rythmos {


/** \brief Extract the state and sensitivities from x_bar and x_bar_dot.
 *

 \verbatim

  x_bar = [ x ; s_bar ]

  x_bar_dot = [ x_dot; s_bar_dot ]

  s_bar = ..

  s_bar_dot = ...

 \endverbatim


 */
template<class Scalar>
void extractStateAndSens(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S,
  RCP<const Thyra::VectorBase<Scalar> >  *x_dot,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S_dot
  );


/** \brief Extract the state from x_bar and x_bar_dot.
 */
template<class Scalar>
void extractState(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::VectorBase<Scalar> >  *x_dot
  );


/** \brief Extract the sensitivities from x_bar and x_bar_dot.
 */
template<class Scalar>
void extractSens(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S_dot
  );


} // namespace Rythmos



//
// Implementations
//

namespace Rythmos {

namespace {


template<class Scalar>
void downcastStateAndSens(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::ProductVectorBase<Scalar> > &x_bar_pv,
  RCP<const Thyra::ProductVectorBase<Scalar> > &x_bar_dot_pv
  )
{
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_TEST_FOR_EXCEPT(is_null(x_bar));
  x_bar_pv = Thyra::productVectorBase<Scalar>(x_bar);

  TEUCHOS_TEST_FOR_EXCEPT(is_null(x_bar_dot));
  x_bar_dot_pv = Thyra::productVectorBase<Scalar>(x_bar_dot);
}


template<class Scalar>
void extractStateBlock(
  const RCP<const Thyra::ProductVectorBase<Scalar> > &x_bar_pv,
  const RCP<const Thyra::ProductVectorBase<Scalar> > &x_bar_dot_pv,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::VectorBase<Scalar> > *x_dot
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==x);
  *x = x_bar_pv->getVectorBlock(0);

  TEUCHOS_TEST_FOR_EXCEPT(0==x_dot);
  *x_dot = x_bar_dot_pv->getVectorBlock(0);
}


template<class Scalar>
void extractSensBlock(
  const RCP<const Thyra::ProductVectorBase<Scalar> > &x_bar_pv,
  const RCP<const Thyra::ProductVectorBase<Scalar> > &x_bar_dot_pv,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S_dot)
{
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_TEST_FOR_EXCEPT(0==S);
  const RCP<const Thyra::DefaultMultiVectorProductVector<Scalar> >
    s_bar = rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(
        x_bar_pv->getVectorBlock(1).assert_not_null(), true
        );
  *S = s_bar->getMultiVector();

  TEUCHOS_TEST_FOR_EXCEPT(0==S_dot);
  const RCP<const Thyra::DefaultMultiVectorProductVector<Scalar> >
    s_bar_dot = rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(
        x_bar_dot_pv->getVectorBlock(1).assert_not_null(), true
        );
  *S_dot = s_bar_dot->getMultiVector();
}


} // anonymous namespace

} // namespace Rythmos


template<class Scalar>
void Rythmos::extractStateAndSens(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S,
  RCP<const Thyra::VectorBase<Scalar> >  *x_dot,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S_dot
  )
{
  RCP<const Thyra::ProductVectorBase<Scalar> > x_bar_pv, x_bar_dot_pv;
  downcastStateAndSens(x_bar, x_bar_dot, x_bar_pv, x_bar_dot_pv);

  extractStateBlock(x_bar_pv, x_bar_dot_pv, x, x_dot);
  extractSensBlock(x_bar_pv, x_bar_dot_pv, S, S_dot);
}


template<class Scalar>
void Rythmos::extractState(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::VectorBase<Scalar> >  *x_dot
  )
{
  RCP<const Thyra::ProductVectorBase<Scalar> > x_bar_pv, x_bar_dot_pv;
  downcastStateAndSens(x_bar, x_bar_dot, x_bar_pv, x_bar_dot_pv);

  extractStateBlock(x_bar_pv, x_bar_dot_pv, x, x_dot);
}


template<class Scalar>
void Rythmos::extractSens(
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_dot,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S,
  RCP<const Thyra::MultiVectorBase<Scalar> > *S_dot
  )
{
  RCP<const Thyra::ProductVectorBase<Scalar> > x_bar_pv, x_bar_dot_pv;
  downcastStateAndSens(x_bar, x_bar_dot, x_bar_pv, x_bar_dot_pv);

  extractSensBlock(x_bar_pv, x_bar_dot_pv, S, S_dot);
}



#endif //RYTHMOS_EXTRACT_STATE_AND_SENS_HPP
