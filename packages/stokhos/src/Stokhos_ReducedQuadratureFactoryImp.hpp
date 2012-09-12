// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_TestForException.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Stokhos_SDMUtils.hpp"

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_GLPK
extern "C" {
#include "glpk.h"
}
#endif

#ifdef HAVE_STOKHOS_CLP
#include "coin/ClpSimplex.hpp"
#include "coin/CoinBuild.hpp"
#include "coin/ClpInterior.hpp"
#endif

#ifdef HAVE_STOKHOS_QPOASES
#include "qpOASES.hpp"
#endif

#ifdef HAVE_STOKHOS_DAKOTA
#include "CompressedSensing.hpp"
#endif

template <typename ordinal_type, typename value_type>
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
ReducedQuadratureFactory(
  const Teuchos::ParameterList& params_) :
  params(params_),
  reduction_method(params.get("Reduced Quadrature Method", "Q Squared")),
  solver_method(params.get("Solver Method", "TRSM")),
  eliminate_dependent_rows(params.get("Eliminate Dependent Rows", true)),
  verbose(params.get("Verbose", false)),
  reduction_tol(params.get("Reduction Tolerance", 1.0e-12)),
  objective_value(params.get("Objective Value", 0.0))
{
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::UserDefinedQuadrature<ordinal_type, value_type> >
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
createReducedQuadrature(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q2,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights) const
{
  Teuchos::RCP< Teuchos::Array<value_type> > red_weights;
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > red_points;
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > red_values;

  // Compute reduced quadrature rule
  if (reduction_method == "Q Squared") {
    if (eliminate_dependent_rows)
      reducedQuadrature_Q_Squared_CPQR(Q, F, weights, 
				       red_weights, red_points, red_values);
    else
      reducedQuadrature_Q_Squared(Q, F, weights, 
				  red_weights, red_points, red_values);
  }
  else if (reduction_method == "Q Squared2") {
    reducedQuadrature_Q_Squared_CPQR2(Q, F, weights, 
				      red_weights, red_points, red_values);
  }
  else if (reduction_method == "Q2") {
    if (eliminate_dependent_rows)
      reducedQuadrature_Q2_CPQR(Q, Q2, F, weights, 
				red_weights, red_points, red_values);
    else
      reducedQuadrature_Q2(Q, Q2, F, weights, 
			   red_weights, red_points, red_values);
  }
  else if (reduction_method == "None") {
    ordinal_type sz = Q.numCols();
    ordinal_type nqp = Q.numRows();
    ordinal_type d = F.numCols();
    red_weights =
      Teuchos::rcp(new Teuchos::Array<value_type>(weights));
    red_points =
      Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(nqp));
    red_values =
      Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(nqp));
    for (ordinal_type i=0; i<nqp; i++) {
      (*red_points)[i].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*red_points)[i][j] = F(i,j);
      (*red_values)[i].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*red_values)[i][j] = Q(i,j);
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid dimension reduction method " << reduction_method);

  // Build reduced quadrature object
  Teuchos::RCP< const Teuchos::Array<value_type> > cred_weights =
    red_weights;
  Teuchos::RCP< const Teuchos::Array< Teuchos::Array<value_type> > > cred_points = red_points;
  Teuchos::RCP< const Teuchos::Array< Teuchos::Array<value_type> > > cred_values = red_values;
  
  Teuchos::RCP<const Stokhos::UserDefinedQuadrature<ordinal_type, value_type> > 
    reduced_quad =
    Teuchos::rcp(new UserDefinedQuadrature<ordinal_type,value_type>(
		   cred_points,
		   cred_weights,
		   cred_values));

  return reduced_quad;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_Q_Squared(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q-squared matrix with all possible products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2(nqp, sz2);
  ordinal_type jdx=0;
  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type k=j; k<sz; k++) {
      for (ordinal_type i=0; i<nqp; i++)
	Q2(i,jdx) = Q(i,j)*Q(i,k);
      jdx++;
    }
  }
  TEUCHOS_ASSERT(jdx == sz2);

  // Compute e_1 = Q2^T*w which is our constraint
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz2);
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  ordinal_type ret = 
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q2, w, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Solve problem
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  underdetermined_solver(Q2, e1, u, Teuchos::TRANS, Teuchos::UNDEF_TRI);

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;

    std::cout << "reduced weights = " << u << std::endl;

    // Check residual error ||e1-Q2^T*u||
    Teuchos::SerialDenseVector<ordinal_type,value_type> err1(e1);
    ret = err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q2, u, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||e1-Q2^T*u||_infty = " << err1.normInf() << std::endl;

    // Check discrete orthogonality error ||I - Q^T*diag(u)*Q||
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> err2(sz, sz);
    err2.putScalar(0.0);
    for (ordinal_type i=0; i<sz; i++)
      err2(i,i) = 1.0;
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> WQ(nqp, sz);
    for (ordinal_type i=0; i<nqp; i++)
      for (ordinal_type j=0; j<sz; j++)
	WQ(i,j) = u[i]*Q(i,j);
    ret = err2.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, WQ, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||I-Q^T*diag(u)*Q||_infty = " << err2.normInf() << std::endl;
    print_matlab(std::cout, err2);
  }
  
  ordinal_type rank = 0;
  for (ordinal_type i=0; i<nqp; i++)
    if (std::abs(u[i]) > reduction_tol) ++rank;

  // Get reduced weights, points and values
  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (std::abs(u[i]) > reduction_tol) {
      (*reduced_weights)[idx] = u[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  
  // idx may be less than rank if we obtained zero weights in solving
  // the least-squares problem
  TEUCHOS_ASSERT(idx <= rank);
  if (idx < rank) {
    rank = idx;
    reduced_weights->resize(rank);
    reduced_points->resize(rank);
    reduced_values->resize(rank);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_Q_Squared_CPQR(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q-squared matrix with all possible products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2(nqp, sz2);
  ordinal_type jdx=0;
  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type k=j; k<sz; k++) {
      for (ordinal_type i=0; i<nqp; i++)
	Q2(i,jdx) = Q(i,j)*Q(i,k);
      jdx++;
    }
  }
  TEUCHOS_ASSERT(jdx == sz2);

  // Compute e_1 = Q2^T*w which is our constraint
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz2);
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  ordinal_type ret = 
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q2, w, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Compute QR decomposition of Q2^T:  Q2^T*P = Z*R
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2t(Q2, Teuchos::TRANS);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Z, R;
  Teuchos::Array<ordinal_type> piv;
  CPQR_Householder(Q2t, Z, R, piv);
  ordinal_type r = computeRank(R, reduction_tol);
  R.reshape(r, R.numCols());
  Z.reshape(Z.numRows(), r);

  if (verbose) {
    std::cout << "Q2 rank = " << r << std::endl;
  }

  // Compute new constraint vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1t(r);
  ret = e1t.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Z, e1, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Solve problem
  Teuchos::SerialDenseVector<ordinal_type,value_type> ut(nqp);
  underdetermined_solver(R, e1t, ut, Teuchos::NO_TRANS, Teuchos::UPPER_TRI);

  // Transform solution back to original
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[piv[i]] = ut[i];

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;

    std::cout << "reduced weights = " << u << std::endl;

    // Check residual error ||e1-Q2^T*u||
    Teuchos::SerialDenseVector<ordinal_type,value_type> err1(e1);
    ret = err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q2, u, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||e1-Q2^T*u||_infty = " << err1.normInf() << std::endl;

    // Check discrete orthogonality error ||I - Q^T*diag(u)*Q||
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> err2(sz, sz);
    err2.putScalar(0.0);
    for (ordinal_type i=0; i<sz; i++)
      err2(i,i) = 1.0;
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> WQ(nqp, sz);
    for (ordinal_type i=0; i<nqp; i++)
      for (ordinal_type j=0; j<sz; j++)
	WQ(i,j) = u[i]*Q(i,j);
    ret = err2.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, WQ, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||I-Q^T*diag(u)*Q||_infty = " << err2.normInf() << std::endl;
    print_matlab(std::cout, err2);
  }
  
  ordinal_type rank = 0;
  for (ordinal_type i=0; i<nqp; i++)
    if (std::abs(u[i]) > reduction_tol) ++rank;

  // Get reduced weights, points and values
  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (std::abs(u[i]) > reduction_tol) {
      (*reduced_weights)[idx] = u[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  
  // idx may be less than rank if we obtained zero weights in solving
  // the least-squares problem
  TEUCHOS_ASSERT(idx <= rank);
  if (idx < rank) {
    rank = idx;
    reduced_weights->resize(rank);
    reduced_points->resize(rank);
    reduced_values->resize(rank);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_Q_Squared_CPQR2(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q-squared matrix with all possible products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2(nqp, sz2);
  ordinal_type jdx=0;
  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type k=j; k<sz; k++) {
      for (ordinal_type i=0; i<nqp; i++)
	Q2(i,jdx) = Q(i,j)*Q(i,k);
      jdx++;
    }
  }
  TEUCHOS_ASSERT(jdx == sz2);

  // Compute QR decomposition of Q2:  Q2*P = Z*R
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Z, R;
  Teuchos::Array<ordinal_type> piv;
  CPQR_Householder(Q2, Z, R, piv);
  ordinal_type r = computeRank(R, reduction_tol);

  // Get the first r columns of Q2*P
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> QQ2(nqp, r);
  for (ordinal_type j=0; j<r; j++)
    for (ordinal_type i=0; i<nqp; i++)
      QQ2(i,j) = Q2(i,piv[j]);

  if (verbose) {
    std::cout << "Q2 rank = " << r << std::endl;
  }

  // Compute e_1 = QQ2^T*w which is our constraint
  Teuchos::SerialDenseVector<ordinal_type,value_type> ee1(r);
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  ordinal_type ret = 
    ee1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, QQ2, w, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Solve problem
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  underdetermined_solver(QQ2, ee1, u, Teuchos::TRANS, Teuchos::UNDEF_TRI);

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;

    std::cout << "reduced weights = " << u << std::endl;

    // Check residual error ||e1-Q2^T*u||
    Teuchos::SerialDenseVector<ordinal_type,value_type> err1(sz2);
    ret = err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q2, w, 0.0);
    TEUCHOS_ASSERT(ret == 0);
    ret = err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q2, u, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||e1-Q2^T*u||_infty = " << err1.normInf() << std::endl;

    // Check discrete orthogonality error ||I - Q^T*diag(u)*Q||
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> err2(sz, sz);
    err2.putScalar(0.0);
    for (ordinal_type i=0; i<sz; i++)
      err2(i,i) = 1.0;
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> WQ(nqp, sz);
    for (ordinal_type i=0; i<nqp; i++)
      for (ordinal_type j=0; j<sz; j++)
	WQ(i,j) = u[i]*Q(i,j);
    ret = err2.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, WQ, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||I-Q^T*diag(u)*Q||_infty = " << err2.normInf() << std::endl;
    print_matlab(std::cout, err2);
  }
  
  ordinal_type rank = 0;
  for (ordinal_type i=0; i<nqp; i++)
    if (std::abs(u[i]) > reduction_tol) ++rank;

  // Get reduced weights, points and values
  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (std::abs(u[i]) > reduction_tol) {
      (*reduced_weights)[idx] = u[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  
  // idx may be less than rank if we obtained zero weights in solving
  // the least-squares problem
  TEUCHOS_ASSERT(idx <= rank);
  if (idx < rank) {
    rank = idx;
    reduced_weights->resize(rank);
    reduced_points->resize(rank);
    reduced_values->resize(rank);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_Q2(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q2,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  
  ordinal_type sz = Q.numCols();
  ordinal_type nqp = Q.numRows();
  ordinal_type sz2 = Q2.numCols();
  ordinal_type d = F.numCols();
  
  // Compute e_1 = Q2^T*w which is our constraint
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz2);
  ordinal_type ret = 
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q2, w, 0.0);
  TEUCHOS_ASSERT(ret == 0);
  std::cout << "constraint vector = " << e1 << std::endl;
  // e1.putScalar(0.0);
  // e1[0] = 1.0;
  // ordinal_type ret;

  // Solve problem
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  underdetermined_solver(Q2, e1, u, Teuchos::TRANS, Teuchos::UNDEF_TRI);

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;
    
    std::cout << "reduced weights = " << u << std::endl;

    // Check residual error ||e1-B2^T*u||
    Teuchos::SerialDenseVector<ordinal_type,value_type> err1(e1);
    ret = err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q2, u, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||e1-Q2^T*u||_infty = " << err1.normInf() << std::endl;

    // Check discrete orthogonality error ||I - Q^T*diag(u)*Q||
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> err2(sz, sz);
    err2.putScalar(0.0);
    for (ordinal_type i=0; i<sz; i++)
      err2(i,i) = 1.0;
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> WQ(nqp, sz);
    for (ordinal_type i=0; i<nqp; i++)
      for (ordinal_type j=0; j<sz; j++)
	WQ(i,j) = u[i]*Q(i,j);
    ret = err2.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, WQ, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||I-Q^T*diag(u)*Q||_infty = " << err2.normInf() << std::endl;

    // Check discrete orthogonality error ||I - Q2^T*diag(u)*Q2||
    // Teuchos::SerialDenseMatrix<ordinal_type,value_type> err3(sz2, sz2);
    // err3.putScalar(0.0);
    // for (ordinal_type i=0; i<sz2; i++)
    //   err3(i,i) = 1.0;
    // Teuchos::SerialDenseMatrix<ordinal_type,value_type> WQ2(nqp, sz2);
    // for (ordinal_type i=0; i<nqp; i++)
    //   for (ordinal_type j=0; j<sz2; j++)
    // 	WQ2(i,j) = w[i]*Q2(i,j);
    // ret = err3.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q2, WQ2, 1.0);
    // TEUCHOS_ASSERT(ret == 0);
    // std::cout << "||I-Q2^T*diag(u)*Q2||_infty = " << err3.normInf() << std::endl;
    // print_matlab(std::cout, err3);
  }
  
  ordinal_type rank = 0;
  for (ordinal_type i=0; i<nqp; i++)
    if (std::abs(u[i]) > reduction_tol) ++rank;

  // Get reduced weights, points and values
  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (std::abs(u[i]) > reduction_tol) {
      (*reduced_weights)[idx] = u[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  
  // idx may be less than rank if we obtained zero weights in solving
  // the least-squares problem
  TEUCHOS_ASSERT(idx <= rank);
  if (idx < rank) {
    rank = idx;
    reduced_weights->resize(rank);
    reduced_points->resize(rank);
    reduced_values->resize(rank);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_Q2_CPQR(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q2,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  
  ordinal_type sz = Q.numCols();
  ordinal_type nqp = Q.numRows();
  ordinal_type sz2 = Q2.numCols();
  ordinal_type d = F.numCols();
  
  // Compute e_1 = Q2^T*w which is our constraint
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz2);
  ordinal_type ret = 
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q2, w, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Compute QR decomposition of Q2^T:  Q2^T*P = Z*R
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2t(Q2, Teuchos::TRANS);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Z, R;
  Teuchos::Array<ordinal_type> piv;
  CPQR_Householder(Q2t, Z, R, piv);
  ordinal_type r = computeRank(R, reduction_tol);
  R.reshape(r, R.numCols());
  Z.reshape(Z.numRows(), r);

  if (verbose) {
    std::cout << "Q2 rank = " << r << std::endl;
  }

  // Compute new constraint vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1t(r);
  ret = e1t.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Z, e1, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Solve problem
  Teuchos::SerialDenseVector<ordinal_type,value_type> ut(nqp);
  underdetermined_solver(R, e1t, ut, Teuchos::NO_TRANS, Teuchos::UPPER_TRI);

  // Transform solution back to original
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[piv[i]] = ut[i];

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;

    std::cout << "reduced weights = " << u << std::endl;

    // Check residual error ||e1-B2^T*u||
    Teuchos::SerialDenseVector<ordinal_type,value_type> err1(e1);
    ret = err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q2, u, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||e1-Q2^T*u||_infty = " << err1.normInf() << std::endl;

    // Check discrete orthogonality error ||I - Q^T*diag(u)*Q||
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> err2(sz, sz);
    err2.putScalar(0.0);
    for (ordinal_type i=0; i<sz; i++)
      err2(i,i) = 1.0;
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> WQ(nqp, sz);
    for (ordinal_type i=0; i<nqp; i++)
      for (ordinal_type j=0; j<sz; j++)
	WQ(i,j) = u[i]*Q(i,j);
    ret = err2.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, WQ, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||I-Q^T*diag(u)*Q||_infty = " << err2.normInf() << std::endl;
    print_matlab(std::cout, err2);
  }
  
  ordinal_type rank = 0;
  for (ordinal_type i=0; i<nqp; i++)
    if (std::abs(u[i]) > reduction_tol) ++rank;

  // Get reduced weights, points and values
  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (std::abs(u[i]) > reduction_tol) {
      (*reduced_weights)[idx] = u[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  
  // idx may be less than rank if we obtained zero weights in solving
  // the least-squares problem
  TEUCHOS_ASSERT(idx <= rank);
  if (idx < rank) {
    rank = idx;
    reduced_weights->resize(rank);
    reduced_points->resize(rank);
    reduced_values->resize(rank);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
underdetermined_solver(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
  if (solver_method == "TRSM")
    solver_TRSM(A, b, x, transa, uplo);
  else if (solver_method == "Clp")
    solver_CLP(A, b, x, transa, uplo);
  else if (solver_method == "Clp-IP")
    solver_CLP_IP(A, b, x, transa, uplo);
  else if (solver_method == "GLPK")
    solver_GLPK(A, b, x, transa, uplo);
  else if (solver_method == "qpOASES")
    solver_qpOASES(A, b, x, transa, uplo);
  else if (solver_method == "Basis Pursuit")
    solver_BasisPursuit(A, b, x, transa, uplo);
  else if (solver_method == "Orthogonal Matching Pursuit")
    solver_OrthogonalMatchingPursuit(A, b, x, transa, uplo);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Invalid solver method " << solver_method);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_TRSM(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(uplo == Teuchos::UNDEF_TRI, std::logic_error,
			     "TRSM solver requires triangular matrix!");
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_rows;
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_rows = m;
    n_cols = n;
  }
  else {
    n_rows = n;
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);
  x.putScalar(0.0);
  blas.COPY(n_rows, b.values(), 1, x.values(), 1);

  // Here we assume A is upper triangular
  blas.TRSM(Teuchos::LEFT_SIDE, uplo, transa, 
	    Teuchos::NON_UNIT_DIAG, n_rows, 1, 1.0, A.values(), A.stride(), 
	    x.values(), x.stride());
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_GLPK(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
#ifdef HAVE_STOKHOS_GLPK
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_rows;
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_rows = m;
    n_cols = n;
  }
  else {
    n_rows = n;
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);

  LPX *lp = lpx_create_prob();
  lpx_set_prob_name(lp, "Reduced Quadrature");
  lpx_set_obj_dir(lp, LPX_MIN);
  lpx_add_rows(lp, n_rows);
  lpx_add_cols(lp, n_cols);

  // Set row bounds
  for (ordinal_type i=0; i<n_rows; i++)
    lpx_set_row_bnds(lp, i+1, LPX_FX, b[i], b[i]);

  // Set columns bounds and object coefficients
  for (ordinal_type j=0; j<n_cols; j++) {
    lpx_set_col_bnds(lp, j+1, LPX_LO, 0.0, 0.0);
    lpx_set_obj_coef(lp, j+1, objective_value);
  }

  // Set constraint matrix
  // GLPK uses this silly 1-based indexing scheme which essentially
  // requires us to copy the matrix.  While copying we will also transpose
  // to make loading in GLPK easier
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA(n_rows+1,n_cols+1);
  AA.putScalar(0.0);
  Teuchos::EUplo AA_uplo = uplo;
  if (transa == Teuchos::NO_TRANS) {
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> AAA(
      Teuchos::View, AA, n_rows, n_cols, 1, 1);
    AAA.assign(A);
  }
  else {
    for (ordinal_type i=0; i<n_rows; i++)
      for (ordinal_type j=0; j<n_cols; j++)
	AA(i+1,j+1) = A(j,i);
    if (uplo == Teuchos::UPPER_TRI)
      AA_uplo = Teuchos::LOWER_TRI;
    else if (uplo == Teuchos::LOWER_TRI)
      AA_uplo = Teuchos::UPPER_TRI;
  }
  Teuchos::Array< Teuchos::Array<int> > row_indices(n_cols);
  if (AA_uplo == Teuchos::UNDEF_TRI) {
    for (ordinal_type j=0; j<n_cols; j++) {
      row_indices[j].resize(n_rows+1);
      for (ordinal_type i=0; i<n_rows; i++)
	row_indices[j][i+1] = i+1; 
      lpx_set_mat_col(lp, j+1, n_rows, &row_indices[j][0], AA[j+1]);
    }
  }
  else if (AA_uplo == Teuchos::UPPER_TRI) {
    // For AA upper-triangular, only need to include leading 
    // min(n_rows,n_cols) x n_cols block (remaining rows must be zero)
    for (ordinal_type j=0; j<n_cols; j++) {
      ordinal_type nr = n_rows;
      if (j < n_rows)
	nr = j+1;
      row_indices[j].resize(nr+1);
      for (ordinal_type i=0; i<nr; i++)
	row_indices[j][i+1] = i+1; 
      lpx_set_mat_col(lp, j+1, nr, &row_indices[j][0], AA[j+1]);
    }
  }
  else if (AA_uplo == Teuchos::LOWER_TRI) {
    // For AA lower-triangular, only need to include leading 
    // n_rows x min(n_rows,n_cols) block (remaining columns must be zero)
    for (ordinal_type j=0; j<n_cols; j++) {
      if (j < n_rows) {
	ordinal_type nr = n_rows-j;
	row_indices[j].resize(nr+1);
	for (ordinal_type i=0; i<nr; i++)
	  row_indices[j][i+1] = j+i+1; 
	lpx_set_mat_col(lp, j+1, nr, &row_indices[j][0], AA[j+1]+j);
      }
    }
  }
    
  // Write MPS-file if requested
  if (params.get("Write MPS File", false) == true)
    lpx_write_mps(lp, params.get("MPS Filename", "lp.mps").c_str());

  // Solve linear program
  lpx_simplex(lp);
  int status = lpx_get_status(lp);
  if (verbose) {
    std::cout << "glpk status = " << status << std::endl;
    double z = lpx_get_obj_val(lp);
    std::cout << "glpk objective = " << z << std::endl;
  }

  // Get solution
  for (ordinal_type i=0; i<n_cols; i++)
    x[i] = lpx_get_col_prim(lp, i+1);
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "GLPK solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_CLP(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
#ifdef HAVE_STOKHOS_CLP
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_rows;
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_rows = m;
    n_cols = n;
  }
  else {
    n_rows = n;
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);

  // Setup linear program
  ClpSimplex model;
  model.resize(n_rows, 0);

  // Set row bounds
  for (ordinal_type i=0; i<n_rows; i++) {
    model.setRowLower(i, b[i]);
    model.setRowUpper(i, b[i]);
  }

  // Set constraint matrix, columns bounds and objective coefficients
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA;
  Teuchos::EUplo AA_uplo = uplo;
  if (transa == Teuchos::NO_TRANS) {
     AA = Teuchos::SerialDenseMatrix<ordinal_type, value_type>(
       Teuchos::View, A, n_rows, n_cols);
  }
  else {
    AA.reshape(n_rows, n_cols);
    for (ordinal_type i=0; i<n_rows; i++)
      for (ordinal_type j=0; j<n_cols; j++)
	AA(i,j) = A(j,i);
    if (uplo == Teuchos::UPPER_TRI)
      AA_uplo = Teuchos::LOWER_TRI;
    else if (uplo == Teuchos::LOWER_TRI)
      AA_uplo = Teuchos::UPPER_TRI;
  }
  Teuchos::Array< Teuchos::Array<int> > row_indices(n_cols);
  CoinBuild buildObject;
  if (AA_uplo == Teuchos::UNDEF_TRI) {
    for (ordinal_type j=0; j<n_cols; j++) {
      row_indices[j].resize(n_rows);
      for (ordinal_type i=0; i<n_rows; i++)
	row_indices[j][i] = i; 
      buildObject.addColumn(
	n_rows, &row_indices[j][0], AA[j], 0.0, DBL_MAX, objective_value);
    }
  }
  else if (AA_uplo == Teuchos::UPPER_TRI) {
    // For AA upper-triangular, only need to include leading 
    // min(n_rows,n_cols) x n_cols block (remaining rows must be zero)
    for (ordinal_type j=0; j<n_cols; j++) {
      ordinal_type nr = n_rows;
      if (j < n_rows)
	nr = j+1;
      row_indices[j].resize(nr);
      for (ordinal_type i=0; i<nr; i++)
	row_indices[j][i] = i; 
      buildObject.addColumn(nr, &row_indices[j][0], AA[j], 0.0, DBL_MAX, objective_value);
    }
  }
  else if (AA_uplo == Teuchos::LOWER_TRI) {
    // For AA lower-triangular, only need to include leading 
    // n_rows x min(n_rows,n_cols) block (remaining columns must be zero)
    for (ordinal_type j=0; j<n_cols; j++) {
      if (j < n_rows) {
	ordinal_type nr = n_rows-j;
	row_indices[j].resize(nr);
	for (ordinal_type i=0; i<nr; i++)
	  row_indices[j][i] = j+i; 
	buildObject.addColumn(
	  nr, &row_indices[j][0], AA[j]+j, 0.0, DBL_MAX, objective_value);
      }
      else
	buildObject.addColumn(0, NULL, NULL, 0.0, DBL_MAX, objective_value);
    }
  }
  model.addColumns(buildObject);

  // Solve linear program
  model.primal();

  // Get solution
  double *solution = model.primalColumnSolution();
  for (ordinal_type i=0; i<n_cols; i++)
    x[i] = solution[i];
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "CLP solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_CLP_IP(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
#ifdef HAVE_STOKHOS_CLP
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_rows;
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_rows = m;
    n_cols = n;
  }
  else {
    n_rows = n;
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);

  // Setup linear program
  ClpInterior model;
  model.resize(n_rows, 0);

  // Set row bounds
  for (ordinal_type i=0; i<n_rows; i++) {
    model.setRowLower(i, b[i]);
    model.setRowUpper(i, b[i]);
  }

  // Set constraint matrix, columns bounds and objective coefficients
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA;
  Teuchos::EUplo AA_uplo = uplo;
  if (transa == Teuchos::NO_TRANS) {
     AA = Teuchos::SerialDenseMatrix<ordinal_type, value_type>(
       Teuchos::View, A, n_rows, n_cols);
  }
  else {
    AA.reshape(n_rows, n_cols);
    for (ordinal_type i=0; i<n_rows; i++)
      for (ordinal_type j=0; j<n_cols; j++)
	AA(i,j) = A(j,i);
    if (uplo == Teuchos::UPPER_TRI)
      AA_uplo = Teuchos::LOWER_TRI;
    else if (uplo == Teuchos::LOWER_TRI)
      AA_uplo = Teuchos::UPPER_TRI;
  }
  Teuchos::Array< Teuchos::Array<int> > row_indices(n_cols);
  CoinBuild buildObject;
  if (AA_uplo == Teuchos::UNDEF_TRI) {
    for (ordinal_type j=0; j<n_cols; j++) {
      row_indices[j].resize(n_rows);
      for (ordinal_type i=0; i<n_rows; i++)
	row_indices[j][i] = i; 
      buildObject.addColumn(
	n_rows, &row_indices[j][0], AA[j], 0.0, DBL_MAX, objective_value);
    }
  }
  else if (AA_uplo == Teuchos::UPPER_TRI) {
    // For AA upper-triangular, only need to include leading 
    // min(n_rows,n_cols) x n_cols block (remaining rows must be zero)
    for (ordinal_type j=0; j<n_cols; j++) {
      ordinal_type nr = n_rows;
      if (j < n_rows)
	nr = j+1;
      row_indices[j].resize(nr);
      for (ordinal_type i=0; i<nr; i++)
	row_indices[j][i] = i; 
      buildObject.addColumn(nr, &row_indices[j][0], AA[j], 0.0, DBL_MAX, objective_value);
    }
  }
  else if (AA_uplo == Teuchos::LOWER_TRI) {
    // For AA lower-triangular, only need to include leading 
    // n_rows x min(n_rows,n_cols) block (remaining columns must be zero)
    for (ordinal_type j=0; j<n_cols; j++) {
      if (j < n_rows) {
	ordinal_type nr = n_rows-j;
	row_indices[j].resize(nr);
	for (ordinal_type i=0; i<nr; i++)
	  row_indices[j][i] = j+i; 
	buildObject.addColumn(
	  nr, &row_indices[j][0], AA[j]+j, 0.0, DBL_MAX, objective_value);
      }
      else
	buildObject.addColumn(0, NULL, NULL, 0.0, DBL_MAX, objective_value);
    }
  }
  model.addColumns(buildObject);

  // Solve linear program
  model.primalDual();

  // Get solution
  double *solution = model.primalColumnSolution();
  for (ordinal_type i=0; i<n_cols; i++)
    x[i] = solution[i];
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "CLP solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_qpOASES(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
#ifdef HAVE_STOKHOS_QPOASES
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_rows;
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_rows = m;
    n_cols = n;
  }
  else {
    n_rows = n;
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);
  
  // Compute objective vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> c(n_cols);
  c.putScalar(objective_value);

  // Compute lower bounds
  Teuchos::SerialDenseVector<ordinal_type,value_type> lb(n_cols);
  lb.putScalar(0.0);

  // Setup linear program
  qpOASES::QProblem lp(n_cols, n_rows, qpOASES::HST_ZERO);

  // Solve linear program -- qpOASES expects row-wise storage
  // so transpose the matrix if necessary.  Since qpOASES uses dense
  // linear algebra, can't take advantage of upper/lower triangular info
  int nWSR = params.get("Max Working Set Recalculations", 10000);
  if (transa == Teuchos::NO_TRANS) {
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA(A, Teuchos::TRANS);
    lp.init(NULL,         // zero Hessian
	    c.values(),   // objective
	    AA.values(),  // constrait matrix
	    lb.values(),  // variable lower bounds
	    NULL,         // no upper bounds
	    b.values(),   // constraint lower bounds
	    b.values(),   // constraint upper bounds
	    nWSR,         // maximum number of working set recalculations
	    NULL          // maximum CPU time
      );
  }
  else {
    lp.init(NULL,         // zero Hessian
	    c.values(),   // objective
	    A.values(),   // constrait matrix
	    lb.values(),  // variable lower bounds
	    NULL,         // no upper bounds
	    b.values(),   // constraint lower bounds
	    b.values(),   // constraint upper bounds
	    nWSR,         // maximum number of working set recalculations
	    NULL          // maximum CPU time
      );
  }

  // Get solution
  lp.getPrimalSolution(x.values());
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "qpOASES solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_BasisPursuit(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
#ifdef HAVE_STOKHOS_DAKOTA
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_cols = n;
  }
  else {
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);

  // Setup L1 minimization problem
  Pecos::CompressedSensing CS;
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA(A, transa);
  Teuchos::SerialDenseVector<ordinal_type, value_type> bb(b);
  CS.BasisPursuit(AA, bb, x);
  
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "BasisPursuit solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
solver_OrthogonalMatchingPursuit(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  Teuchos::ETransp transa, Teuchos::EUplo uplo) const
{
#ifdef HAVE_STOKHOS_DAKOTA
  ordinal_type m = A.numRows();
  ordinal_type n = A.numCols();
  ordinal_type n_cols;
  if (transa == Teuchos::NO_TRANS) {
    n_cols = n;
  }
  else {
    n_cols = m;
  }
  if (x.length() < n_cols)
    x.shape(n_cols, 1);

  // Setup L1 minimization problem
  Pecos::CompressedSensing CS;
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA(A, transa);
  Teuchos::SerialDenseVector<ordinal_type, value_type> bb(b);
  CS.OrthogonalMatchingPursuit(AA, bb, x);
  
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "BasisPursuit solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
computeRank(
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& R,
  const value_type tol) const 
{
  // Compute rank -- since we constrain the first d+1 columns of Bp to lie
  // in the basis, the diagonal elements of Rp may not be monotonically 
  // decreasing until we get past d+1
  ordinal_type rank = 0;
  ordinal_type m = R.numRows();
  value_type r_max = std::abs(R(rank,rank));
  value_type r_min = std::abs(R(rank,rank));
  for (rank=1; rank<m; rank++) {
    if (std::abs(R(rank,rank)) > r_max)
      r_max = std::abs(R(rank,rank));
    if (std::abs(R(rank,rank)) < r_min)
      r_min = std::abs(R(rank,rank));
    if (r_min / r_max < tol)
      break;
  }  

  // Print condition number of R
  if (verbose) {
    value_type cond_r = r_max / r_min;
    std::cout << "r_max = " << r_max << std::endl;
    std::cout << "r_min = " << r_min << std::endl;
    std::cout << "Condition number of R = " << cond_r << std::endl;
  }

  return rank;
}
