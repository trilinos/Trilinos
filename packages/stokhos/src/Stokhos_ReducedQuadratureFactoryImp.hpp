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
#endif

#ifdef HAVE_STOKHOS_QPOASES
#include "qpOASES.hpp"
#endif

template <typename ordinal_type, typename value_type>
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
ReducedQuadratureFactory(
  const Teuchos::ParameterList& params_) :
  params(params_),
  reduction_method(params.get("Reduced Quadrature Method", 
			      "Column-Pivoted QR")),
  verbose(params.get("Verbose", false)),
  reduction_tol(params.get("Reduction Tolerance", 1.0e-12))
{
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::UserDefinedQuadrature<ordinal_type, value_type> >
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
createReducedQuadrature(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights) const
{
  Teuchos::RCP< Teuchos::Array<value_type> > red_weights;
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > red_points;
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > red_values;

  // Compute reduced quadrature rule
  if (reduction_method == "Column-Pivoted QR")
    reducedQuadrature_QRCP(Q, F, weights, 
			   red_weights, red_points, red_values);
  else if (reduction_method == "L1 Minimization") {
    std::string solver = params.get("LP Solver", "GLPK");
    if (solver == "GLPK")
       reducedQuadrature_GLPK(Q, F, weights, 
			      red_weights, red_points, red_values);
    else if (solver == "CLP")
       reducedQuadrature_CLP(Q, F, weights, 
			     red_weights, red_points, red_values);
    else if (solver == "qpOASES")
       reducedQuadrature_qpOASES(Q, F, weights, 
				 red_weights, red_points, red_values);
    else if (solver == "GLPK-CPQR")
       reducedQuadrature_GLPK_CPQR(Q, F, weights, 
				   red_weights, red_points, red_values);
    else if (solver == "CLP-CPQR")
       reducedQuadrature_CLP_CPQR(Q, F, weights, 
				  red_weights, red_points, red_values);
    else if (solver == "qpOASES-CPQR")
       reducedQuadrature_qpOASES_CPQR(Q, F, weights, 
				      red_weights, red_points, red_values);
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, "Invalid LP solver method " << solver);
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
reducedQuadrature_QRCP(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  //
  // Find reduced quadrature weights by applying column-pivoted QR to
  // problem Q2^T*w = e_1
  //
  ordinal_type sz = Q.numCols();
  //ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type sz2 = sz*sz;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;
  }

  // Compute Q matrix with all possible products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2t(sz2, nqp);
  ordinal_type jdx=0;
  for (ordinal_type j=0; j<sz; j++) {
    //for (ordinal_type k=j; k<sz; k++) {
    for (ordinal_type k=0; k<sz; k++) {
      for (ordinal_type i=0; i<nqp; i++)
	Q2t(jdx,i) = Q(i,j)*Q(i,k);
      jdx++;
    }
  }
  TEUCHOS_ASSERT(jdx == sz2);

  // Apply column-pivoted QR to Q2
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> QQ2, R2;
  Teuchos::Array<ordinal_type> piv2;
  Stokhos::CPQR_Householder(Q2t, QQ2, R2, piv2);
  ordinal_type rank = computeRank(R2, reduction_tol);

  if (verbose) {
    std::cout << "rank = " << rank << std::endl;

    // std::cout << "full diag(R) = [ ";
    // for (ordinal_type i=0; i<R2.numRows(); i++)
    //   std::cout << R2(i,i) << " ";
    // std::cout << "]" << std::endl;

    std::cout << "diag(R) = [ ";
    for (ordinal_type i=0; i<rank; i++)
      std::cout << R2(i,i) << " ";
    std::cout << "]" << std::endl;
  }

  // Apply c = QQ2^T*Q2^T*w
  ordinal_type n = QQ2.numCols();
  //std::cout << "n = " << n << std::endl;
  Teuchos::SerialDenseVector<ordinal_type,value_type> b(sz2), c(n);
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q2t, w, 0.0);
  c.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, QQ2, b, 0.0);

  // Compute [R11^{-1}*b1; 0] where b1 is the first r rows of b
  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
	    Teuchos::NON_UNIT_DIAG, rank, 1, 1.0, R2.values(), n, c.values(),
	    n);
  for (ordinal_type i=rank; i<n; i++)
    c[i] = 0.0;

  // Get reduced weights, points and values
  Teuchos::SerialDenseVector<ordinal_type,value_type> wt(nqp);
  wt.putScalar(0.0);
  for (ordinal_type i=0; i<n; i++)
    wt[piv2[i]] = c[i];

  if (verbose) {
    b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, Q2t, wt, 1.0);
    std::cout << "product integration error = " << b.normInf() << std::endl;

    Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz);
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, w, 0.0);
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, wt, 1.0);
    std::cout << "basis integration error = " << e1.normInf() << std::endl;
  }

  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (wt[i] != 0.0) {
      (*reduced_weights)[idx] = wt[i];
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

/*
template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_CS(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
  //
  // Find reduced quadrature weights by applying column-pivoted QR to
  // problem Q2^T*w = e_1
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  if (verbose) {
    std::cout << "sz = " << sz << std::endl;
    std::cout << "nqp = " << nqp << std::endl;
    std::cout << "sz2 = " << sz2 << std::endl;
  }

  // Compute Q matrix with all possible products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2t(sz2, nqp);
  ordinal_type jdx=0;
  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type k=j; k<sz; k++) {
      for (ordinal_type i=0; i<nqp; i++)
	Q2t(jdx,i) = Q(i,j)*Q(i,k);
      jdx++;
    }
  }
  TEUCHOS_ASSERT(jdx == sz2);

  Teuchos::SerialDenseVector<ordinal_type,value_type> b(sz2), x(nqp);
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q2t, w, 0.0);
  c.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, QQ2, b, 0.0);

  Pecos::CompressedSensing CS;
  CS.BasisPursuit(Q2t, 

  // Apply column-pivoted QR to Q2
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> QQ2, R2;
  Teuchos::Array<ordinal_type> piv2;
  Stokhos::CPQR_Householder3(Q2t, QQ2, R2, piv2);
  ordinal_type rank = computeRank(R2, reduction_tol);

  if (verbose)
    std::cout << "rank = " << rank << std::endl;

  

  // Compute [R11^{-1}*b1; 0] where b1 is the first r rows of b
  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
	    Teuchos::NON_UNIT_DIAG, rank, 1, 1.0, R2.values(), n, c.values(),
	    n);
  for (ordinal_type i=rank; i<n; i++)
    c[i] = 0.0;

  // Get reduced weights, points and values
  Teuchos::SerialDenseVector<ordinal_type,value_type> wt(nqp);
  wt.putScalar(0.0);
  for (ordinal_type i=0; i<n; i++)
    wt[piv2[i]] = c[i];

  if (verbose) {
    b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, Q2t, wt, 1.0);
    std::cout << "product integration error = " << b.normInf() << std::endl;

    Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz);
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, w, 0.0);
    e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -1.0, Q, wt, 1.0);
    std::cout << "basis integration error = " << e1.normInf() << std::endl;
  }

  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (wt[i] != 0.0) {
      (*reduced_weights)[idx] = wt[i];
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
*/

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_GLPK(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
#ifdef HAVE_STOKHOS_GLPK
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2p = Q2*R2
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q matrix with all possible products
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

  // Setup linear program
  LPX *lp = lpx_create_prob();
  lpx_set_prob_name(lp, "MonomialProj PCE Reduction");
  lpx_set_obj_dir(lp, LPX_MIN);
  lpx_add_rows(lp, sz2);
  lpx_add_cols(lp, nqp);
  for (ordinal_type i=0; i<sz2; i++)
    lpx_set_row_bnds(lp, i+1, LPX_FX, e1[i], e1[i]);

  // Set columns bounds and object coefficients
  for (ordinal_type j=0; j<nqp; j++) {
    lpx_set_col_bnds(lp, j+1, LPX_LO, 0.0, 0.0);
    lpx_set_obj_coef(lp, j+1, 1.0);
  }

  // Set constraint matrix = Q2^T
  int **cols = new int*[sz2];
  double **vals = new double*[sz2];
  for (ordinal_type i=0; i<sz2; i++) {
    cols[i] = new int[nqp+1];
    vals[i] = new double[nqp+1];
    for (ordinal_type j=0; j<nqp; j++) {
      cols[i][j+1] = j+1;
      vals[i][j+1] = Q2(j,i);
    }
    lpx_set_mat_row(lp, i+1, nqp, cols[i], vals[i]);
  }

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
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[i] = lpx_get_col_prim(lp, i+1);
  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;

  // Clean up linear program
  lpx_delete_prob(lp);
  for (ordinal_type i=0; i<sz2; i++) {
    delete [] cols[i];
    delete [] vals[i];
  }
  delete [] cols;
  delete [] vals;
  
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
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "GLPK solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_CLP(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
#ifdef HAVE_STOKHOS_CLP
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2p = Q2*R2
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q matrix with all possible products
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

  // Setup linear program
  ClpSimplex model;
  model.resize(0, nqp);

  // Set columns bounds and object coefficients
  for (ordinal_type i=0; i<nqp; i++) {
    model.setColLower(i, 0.0);
    model.setColUpper(i, DBL_MAX);
    model.setObjectiveCoefficient(i, 0.0);
  }

  // Set constraint matrix = Q2^T, set row constraints
  Teuchos::Array< Teuchos::Array<int> > cols(sz2);
  for (ordinal_type i=0; i<sz2; i++) {
    cols[i].resize(nqp);
    for (ordinal_type j=0; j<nqp; j++)
      cols[i][j] = j;
    model.addRow(nqp, &cols[i][0], Q2[i], e1[i], e1[i]);
  }
  
  // Solve linear program
  model.primal();
  double *solution = model.primalColumnSolution();
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(
    Teuchos::Copy, solution, nqp);
  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;
  
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
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "CLP solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_qpOASES(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
#ifdef HAVE_STOKHOS_QPOASES
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2p = Q2*R2
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q matrix with all possible products
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

  // Compute objective vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> c(nqp);
  c.putScalar(1.0);

  // Compute lower bounds
  Teuchos::SerialDenseVector<ordinal_type,value_type> lb(nqp);
  lb.putScalar(0.0);

  // Setup linear program
  int nWSR = 10000;
  qpOASES::QProblem lp(nqp, sz2, qpOASES::HST_ZERO);
  lp.init(NULL,         // zero Hessian
	  c.values(),   // objective
	  Q2.values(),  // constrait matrix -- qpOASES expects row-wise storage
	  lb.values(),  // variable lower bounds
	  NULL,         // no upper bounds
	  e1.values(),  // constraint lower bounds
	  e1.values(),  // constraint upper bounds
	  nWSR,         // maximum number of working set recalculations
	  NULL,         // maximum CPU time
	  w.values(),   // initial guess for primal
	  NULL,         // initial guess for dual
	  NULL,         // guessed bounds
	  NULL          // guessed constraints
    );

  // Get solution
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  lp.getPrimalSolution(u.values());
  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;
  
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
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "qpOASES solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_GLPK_CPQR(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
#ifdef HAVE_STOKHOS_GLPK
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2p = Q2*R2
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q matrix with all possible products
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

  // Compute new constraint vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1t(r);
  ret = e1t.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Z, e1, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Setup linear program
  LPX *lp = lpx_create_prob();
  lpx_set_prob_name(lp, "MonomialProj PCE Reduction");
  lpx_set_obj_dir(lp, LPX_MIN);
  lpx_add_rows(lp, r);
  lpx_add_cols(lp, nqp);
  for (ordinal_type i=0; i<r; i++)
    lpx_set_row_bnds(lp, i+1, LPX_FX, e1t[i], e1t[i]);

  // Set columns bounds and object coefficients
  for (ordinal_type j=0; j<nqp; j++) {
    lpx_set_col_bnds(lp, j+1, LPX_LO, 0.0, 0.0);
    lpx_set_obj_coef(lp, j+1, 0.0);
  }

  // Set constraint matrix = R
  int **cols = new int*[r];
  double **vals = new double*[r];
  for (ordinal_type i=0; i<r; i++) {
    cols[i] = new int[nqp+1];
    vals[i] = new double[nqp+1];
    for (ordinal_type j=0; j<nqp; j++) {
      cols[i][j+1] = j+1;
      vals[i][j+1] = R(i,j);
    }
    lpx_set_mat_row(lp, i+1, nqp, cols[i], vals[i]);
  }

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
  Teuchos::SerialDenseVector<ordinal_type,value_type> ut(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    ut[i] = lpx_get_col_prim(lp, i+1);

  // Transform solution back to original
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[piv[i]] = ut[i];

  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;

  // Clean up linear program
  lpx_delete_prob(lp);
  for (ordinal_type i=0; i<r; i++) {
    delete [] cols[i];
    delete [] vals[i];
  }
  delete [] cols;
  delete [] vals;
  
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
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "GLPK solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_CLP_CPQR(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
#ifdef HAVE_STOKHOS_CLP
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2p = Q2*R2
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q matrix with all possible products
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

  // Compute new constraint vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1t(r);
  ret = e1t.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Z, e1, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  /*
  // Setup linear program
  ClpSimplex model;
  model.resize(0, nqp);

  // Set columns bounds and object coefficients
  for (ordinal_type i=0; i<nqp; i++) {
    model.setColLower(i, 0.0);
    model.setColUpper(i, DBL_MAX);
    model.setObjectiveCoefficient(i, 0.0);
  }

  // Set constraint matrix = R.  Here we use the fact that R is upper
  // triangular
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Rt(R, Teuchos::TRANS);
  Teuchos::Array< Teuchos::Array<int> > cols(r);
  for (ordinal_type i=0; i<r; i++) {
    cols[i].resize(nqp-i);
    for (ordinal_type j=0; j<nqp-i; j++)
      cols[i][j] = i+j;
    model.addRow(nqp-i, &cols[i][0], Rt[i]+i, e1t[i], e1t[i]);
  }
  */

  // Setup linear program
  ClpSimplex model;
  model.resize(r, 0);

  // Set row bounds
  for (ordinal_type i=0; i<r; i++) {
    model.setRowLower(i, e1t[i]);
    model.setRowUpper(i, e1t[i]);
  }

  // Set constraint matrix = R.  Here we use the fact that R is upper
  // triangular
  Teuchos::Array< Teuchos::Array<int> > rows(nqp);
  CoinBuild buildObject;
  for (ordinal_type j=0; j<nqp; j++) {
    int nrows = r;
    if (j < r)
      nrows = j+1;
    rows[j].resize(nrows);
    for (ordinal_type i=0; i<nrows; i++)
      rows[j][i] = i;
    buildObject.addColumn(nrows, &rows[j][0], R[j], 0.0, DBL_MAX, 0.0);
  }
  model.addColumns(buildObject);
  
  // Solve linear program
  model.primal();
  double *solution = model.primalColumnSolution();
  Teuchos::SerialDenseVector<ordinal_type,value_type> ut(
    Teuchos::Copy, solution, nqp);

  // Transform solution back to original
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[piv[i]] = ut[i];

  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;
  
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
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "CLP solver called but not enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ReducedQuadratureFactory<ordinal_type, value_type>::
reducedQuadrature_qpOASES_CPQR(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values
  ) const
{
#ifdef HAVE_STOKHOS_QPOASES
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2p = Q2*R2
  //
  ordinal_type sz = Q.numCols();
  ordinal_type sz2 = sz*(sz+1)/2;
  ordinal_type nqp = Q.numRows();
  ordinal_type d = F.numCols();

  // Compute Q matrix with all possible products
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

  // Compute objective vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> c(nqp);
  c.putScalar(1.0);

  // Compute lower bounds
  Teuchos::SerialDenseVector<ordinal_type,value_type> lb(nqp);
  lb.putScalar(0.0);

  // Compute QR decomposition of Q2^T:  Q2^T*P = Z*R
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2t(Q2, Teuchos::TRANS);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Z, R;
  Teuchos::Array<ordinal_type> piv;
  CPQR_Householder(Q2t, Z, R, piv);
  ordinal_type r = computeRank(R, reduction_tol);
  R.reshape(r, R.numCols());
  Z.reshape(Z.numRows(), r);

  // Compute new constraint vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1t(r);
  ret = e1t.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Z, e1, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Compute new constraint matrix
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Rt(R, Teuchos::TRANS);

  // Compute new objective vector
  Teuchos::SerialDenseVector<ordinal_type,value_type> ct(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    ct[piv[i]] = c[i];

  // Compute initial guess
  Teuchos::SerialDenseVector<ordinal_type,value_type> wt(nqp);
  for (ordinal_type i=0; i<r; i++)
    wt[i] = e1t[i];
  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
	    Teuchos::NON_UNIT_DIAG, r, 1, 1.0, R.values(), r, wt.values(),
	    nqp);
  for (ordinal_type i=r; i<nqp; i++)
    wt[i] = 0.0;

  // Setup linear program
  int nWSR = 10000;
  qpOASES::QProblem lp(nqp, r, qpOASES::HST_ZERO);
  lp.init(NULL,         // zero Hessian
	  ct.values(),  // objective
	  Rt.values(),  // constrait matrix -- qpOASES expects row-wise storage
	  lb.values(),  // variable lower bounds
	  NULL,         // no upper bounds
	  e1t.values(),  // constraint lower bounds
	  e1t.values(),  // constraint upper bounds
	  nWSR,         // maximum number of working set recalculations
	  NULL,         // maximum CPU time
	  wt.values(),  // initial guess for primal
	  NULL,         // initial guess for dual
	  NULL,         // guessed bounds
	  NULL          // guessed constraints
    );

  // Get solution
  Teuchos::SerialDenseVector<ordinal_type,value_type> ut(nqp);
  lp.getPrimalSolution(ut.values());

  // Transform solution back to original
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[piv[i]] = ut[i];

  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;
  
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
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "qpOASES solver called but not enabled!");
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
