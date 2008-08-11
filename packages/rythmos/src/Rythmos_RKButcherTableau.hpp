//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_HPP


#include "Rythmos_Types.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/* \brief . */
template<class Scalar>
class RKButcherTableau {
public:
  /** \brief. */
  RKButcherTableau()
    {}
  /** \brief. */
  RKButcherTableau(
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c
    )
    : A_(A), b_(b), c_(c)
    {
      const int numStages = A.numRows();
      TEUCHOS_ASSERT_EQUALITY( A.numRows(), numStages );
      TEUCHOS_ASSERT_EQUALITY( A.numCols(), numStages );
      TEUCHOS_ASSERT_EQUALITY( b.length(), numStages );
      TEUCHOS_ASSERT_EQUALITY( c.length(), numStages );
    }
  /** \brief . */
  int numStages() const { return A_.numRows(); }
  /** \brief . */
  const Teuchos::SerialDenseMatrix<int,Scalar>& A() const { return A_; }
  /** \brief . */
  const Teuchos::SerialDenseVector<int,Scalar> b() const { return b_; }
  /** \brief . */
  const Teuchos::SerialDenseVector<int,Scalar> c() const { return c_; }
private:
  Teuchos::SerialDenseMatrix<int,Scalar> A_;
  Teuchos::SerialDenseVector<int,Scalar> b_;
  Teuchos::SerialDenseVector<int,Scalar> c_;
};


/* \brief . */
template<class Scalar>
void assembleIRKState(
  const int stageIndex,
  const Teuchos::SerialDenseMatrix<int,Scalar> &A,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::ProductVectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{

  typedef ScalarTraits<Scalar> ST;

  const int numStages = A.numRows();
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( stageIndex, 0, numStages );
  TEUCHOS_ASSERT_EQUALITY( A.numRows(), numStages );
  TEUCHOS_ASSERT_EQUALITY( A.numCols(), numStages );
  TEUCHOS_ASSERT_EQUALITY( x_stage_bar.productSpace()->numBlocks(), numStages );
  Thyra::VectorBase<Scalar>& x_out = *x_out_ptr;

  V_V( outArg(x_out), x_base );
  for ( int j = 0; j < numStages; ++j ) {
    Vp_StV( outArg(x_out), dt * A(stageIndex,j), *x_stage_bar.getVectorBlock(j) );
  }

}


/* \brief . */
template<class Scalar>
void assembleIRKSolution(
  const Teuchos::SerialDenseVector<int,Scalar> &b,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::ProductVectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{

  typedef ScalarTraits<Scalar> ST;

  const int numStages = b.length();
  TEUCHOS_ASSERT_EQUALITY( b.length(), numStages );
  TEUCHOS_ASSERT_EQUALITY( x_stage_bar.productSpace()->numBlocks(), numStages );
  Thyra::VectorBase<Scalar>& x_out = *x_out_ptr;

  V_V( outArg(x_out), x_base );
  for ( int j = 0; j < numStages; ++j ) {
    Vp_StV( outArg(x_out), dt * b(j), *x_stage_bar.getVectorBlock(j) );
  }

}



} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_HPP
