// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_KLReducedMatrixFreeOperator.hpp"
#include "Stokhos_PCEAnasaziKL.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"
#include <sstream>

Stokhos::KLReducedMatrixFreeOperator::
KLReducedMatrixFreeOperator(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
  const Teuchos::RCP<const Epetra_Map>& range_base_map_,
  const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
  const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos KL Reduced Matrix Free Operator"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_sg_map(domain_sg_map_),
  range_sg_map(range_sg_map_),
  Cijk(epetraCijk->getParallelCijk()),
  block_ops(),
  params(params_),
  useTranspose(false),
  expansion_size(sg_basis->size()),
  num_blocks(0),
  num_KL(0),
  num_KL_computed(0),
  mean(),
  block_vec_map(),
  block_vec_poly(),
  dot_products(),
  sparse_kl_coeffs(),
  kl_blocks()
{
  num_KL = params->get("Number of KL Terms", 5);
  drop_tolerance = params->get("Sparse 3 Tensor Drop Tolerance", 1e-6);
  do_error_tests = params->get("Do Error Tests", false);
}

void
Stokhos::KLReducedMatrixFreeOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops)
{
  block_ops = ops;
  num_blocks = block_ops->size();

  // Build a vector polynomial out of matrix nonzeros
  mean = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
    block_ops->getCoeffPtr(0));
  block_vec_map =
    Teuchos::rcp(new Epetra_Map(-1, mean->NumMyNonzeros(), 0,
                                domain_base_map->Comm()));
  block_vec_poly =
    Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
                   sg_basis, block_ops->map(), block_vec_map, sg_comm));

  // Setup KL blocks
  setup();
}

Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly >
Stokhos::KLReducedMatrixFreeOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly >
Stokhos::KLReducedMatrixFreeOperator::
getSGPolynomial() const
{
  return block_ops;
}

Stokhos::KLReducedMatrixFreeOperator::
~KLReducedMatrixFreeOperator()
{
}

int
Stokhos::KLReducedMatrixFreeOperator::
SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
  kl_mat_free_op->SetUseTranspose(useTranspose);
  for (int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int
Stokhos::KLReducedMatrixFreeOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return kl_mat_free_op->Apply(Input, Result);
}

int
Stokhos::KLReducedMatrixFreeOperator::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  throw "KLReducedMatrixFreeOperator::ApplyInverse not defined!";
  return -1;
}

double
Stokhos::KLReducedMatrixFreeOperator::
NormInf() const
{
  return 1.0;
}


const char*
Stokhos::KLReducedMatrixFreeOperator::
Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
Stokhos::KLReducedMatrixFreeOperator::
UseTranspose() const
{
  return useTranspose;
}

bool
Stokhos::KLReducedMatrixFreeOperator::
HasNormInf() const
{
  return false;
}

const Epetra_Comm &
Stokhos::KLReducedMatrixFreeOperator::
Comm() const
{
  return *sg_comm;
}
const Epetra_Map&
Stokhos::KLReducedMatrixFreeOperator::
OperatorDomainMap() const
{
  if (useTranspose)
    return *range_sg_map;
  return *domain_sg_map;
}

const Epetra_Map&
Stokhos::KLReducedMatrixFreeOperator::
OperatorRangeMap() const
{
  if (useTranspose)
    return *domain_sg_map;
  return *range_sg_map;
}

void
Stokhos::KLReducedMatrixFreeOperator::
setup()
{
#ifdef HAVE_STOKHOS_ANASAZI
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::KLReducedMatrixFreeOperator -- Calculation/setup of KL opeator");
#endif

  mean = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
    block_ops->getCoeffPtr(0));

  // Copy matrix coefficients into vectors
  for (int coeff=0; coeff<num_blocks; coeff++) {
    Teuchos::RCP<const Epetra_CrsMatrix> block_coeff =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>
        (block_ops->getCoeffPtr(coeff));
    int row = 0;
    for (int i=0; i<mean->NumMyRows(); i++) {
      int num_col;
      mean->NumMyRowEntries(i, num_col);
      for (int j=0; j<num_col; j++)
        (*block_vec_poly)[coeff][row++] = (*block_coeff)[i][j];
    }
  }

  int myPID = sg_comm->MyPID();

  // Compute KL expansion of solution sg_J_vec_poly
  Stokhos::PCEAnasaziKL pceKL(*block_vec_poly, num_KL);
  Teuchos::ParameterList anasazi_params = pceKL.getDefaultParams();
  bool result = pceKL.computeKL(anasazi_params);
  if (!result && myPID == 0)
    std::cout << "KL Eigensolver did not converge!" << std::endl;
  Teuchos::RCP<Epetra_MultiVector> evecs = pceKL.getEigenvectors();
  Teuchos::Array<double> evals = pceKL.getEigenvalues();
  //num_KL_computed = evecs->NumVectors();
  if (myPID == 0)
    std::cout << "num computed eigenvectors  = "
              << evecs->NumVectors() << std::endl;
  double kl_tol = params->get("KL Tolerance", 1e-6);
  num_KL_computed = 0;
  while (num_KL_computed < evals.size() &&
         std::sqrt(evals[num_KL_computed]/evals[0]) > kl_tol)
    num_KL_computed++;
  if (num_KL_computed == evals.size() && myPID == 0)
    std::cout << "Can't achieve KL tolerance " << kl_tol
              << ".  Smallest eigenvalue / largest eigenvalue = "
              << std::sqrt(evals[num_KL_computed-1]/evals[0]) << std::endl;
  if (myPID == 0)
    std::cout << "num KL eigenvectors = " << num_KL_computed << std::endl;

  // Compute dot products of Jacobian blocks and KL eigenvectors
  dot_products.resize(num_KL_computed);
  for (int rv=0; rv < num_KL_computed; rv++) {
    dot_products[rv].resize(num_blocks-1);
    for (int coeff=1; coeff < num_blocks; coeff++) {
      double dot;
      (*block_vec_poly)[coeff].Dot(*((*evecs)(rv)), &dot);
      dot_products[rv][coeff-1] = dot;
    }
  }

  // Compute KL coefficients
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  sparse_kl_coeffs =
    Teuchos::rcp(new Stokhos::Sparse3Tensor<int,double>);
  for (Cijk_type::i_iterator i_it=Cijk->i_begin();
       i_it!=Cijk->i_end(); ++i_it) {
    int i = epetraCijk->GRID(index(i_it));
    sparse_kl_coeffs->sum_term(i, i, 0, norms[i]);
  }
  Cijk_type::k_iterator l_begin = ++(Cijk->k_begin());
  Cijk_type::k_iterator l_end = Cijk->k_end();
  for (Cijk_type::k_iterator l_it=l_begin; l_it!=l_end; ++l_it) {
    int l = index(l_it);
    for (Cijk_type::kj_iterator j_it = Cijk->j_begin(l_it);
         j_it != Cijk->j_end(l_it); ++j_it) {
      int j = epetraCijk->GCID(index(j_it));
      for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
           i_it != Cijk->i_end(j_it); ++i_it) {
        int i = epetraCijk->GRID(index(i_it));
        double c = value(i_it);
        for (int k=1; k<num_KL_computed+1; k++) {
          double dp = dot_products[k-1][l-1];
          double v = dp*c;
          if (std::abs(v) > drop_tolerance)
            sparse_kl_coeffs->sum_term(i,j,k,v);
        }
      }
    }
  }
  sparse_kl_coeffs->fillComplete();

  bool save_tensor = params->get("Save Sparse 3 Tensor To File", false);
  if (save_tensor) {
    static int idx = 0;
    std::string basename = params->get("Sparse 3 Tensor Base Filename",
                                       "sparse_KL_coeffs");
    std::stringstream ss;
    ss << basename << "_" << idx++ << ".mm";
    sparse3Tensor2MatrixMarket(*sparse_kl_coeffs,
                               *(epetraCijk->getStochasticRowMap()), ss.str());
  }

  // Transform eigenvectors back to matrices
  kl_blocks.resize(num_KL_computed);
  Teuchos::RCP<Epetra_BlockMap> kl_map =
    Teuchos::rcp(new Epetra_LocalMap(num_KL_computed+1, 0,
                                     sg_comm->TimeDomainComm()));
  kl_ops =
    Teuchos::rcp(new Stokhos::EpetraOperatorOrthogPoly(
                   sg_basis, kl_map, domain_base_map, range_base_map,
                   range_sg_map, sg_comm));
  kl_ops->setCoeffPtr(0, mean);
  for (int rv=0; rv<num_KL_computed; rv++) {
    if (kl_blocks[rv] == Teuchos::null)
      kl_blocks[rv] = Teuchos::rcp(new Epetra_CrsMatrix(*mean));
    int row = 0;
    for (int i=0; i<mean->NumMyRows(); i++) {
      int num_col;
      mean->NumMyRowEntries(i, num_col);
      for (int j=0; j<num_col; j++)
        (*kl_blocks[rv])[i][j] = (*evecs)[rv][row++];
    }
    kl_ops->setCoeffPtr(rv+1, kl_blocks[rv]);
  }

  Teuchos::RCP<Stokhos::EpetraSparse3Tensor> reducedEpetraCijk =
    Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(
                   sg_basis, sparse_kl_coeffs, sg_comm,
                   epetraCijk->getStochasticRowMap(), sparse_kl_coeffs,
                   0, -1));
  reducedEpetraCijk->transformToLocal();

  // Create matrix-free op
  kl_mat_free_op = Teuchos::rcp(new Stokhos::MatrixFreeOperator(
                                  sg_comm, sg_basis, reducedEpetraCijk,
                                  domain_base_map, range_base_map,
                                  domain_sg_map, range_sg_map, params));
  kl_mat_free_op->setupOperator(kl_ops);

  // Check accuracy of KL expansion
  if (do_error_tests) {
    Teuchos::Array<double> point(sg_basis->dimension());
    for (int i=0; i<sg_basis->dimension(); i++)
      point[i] = 0.5;
    Teuchos::Array<double> basis_vals(sg_basis->size());
    sg_basis->evaluateBases(point, basis_vals);

    Epetra_Vector val(*block_vec_map);
    Epetra_Vector val_kl(*block_vec_map);
    block_vec_poly->evaluate(basis_vals, val);
    val_kl.Update(1.0, (*block_vec_poly)[0], 0.0);
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > rvs(num_KL_computed);
    Teuchos::Array<double> val_rvs(num_KL_computed);
    for (int rv=0; rv<num_KL_computed; rv++) {
      rvs[rv].reset(sg_basis);
      rvs[rv][0] = 0.0;
      for (int coeff=1; coeff<num_blocks; coeff++)
        rvs[rv][coeff] = dot_products[rv][coeff-1];
      val_rvs[rv] = rvs[rv].evaluate(point, basis_vals);
      val_kl.Update(val_rvs[rv], *((*evecs)(rv)), 1.0);
    }
    double nrm;
    val.NormInf(&nrm);
    val.Update(-1.0, val_kl, 1.0);
    double diff;
    val.NormInf(&diff);
    if (myPID == 0)
      std::cout << "Infinity norm of random field difference = " << diff/nrm
                << std::endl;

    // Check accuracy of operator
    Epetra_Vector op_input(*domain_sg_map), op_result(*range_sg_map), op_kl_result(*range_sg_map);
    op_input.PutScalar(1.0);
    Stokhos::MatrixFreeOperator op(sg_comm, sg_basis, epetraCijk,
                                   domain_base_map, range_base_map,
                                   domain_sg_map, range_sg_map, params);
    op.setupOperator(block_ops);
    op.Apply(op_input, op_result);
    this->Apply(op_input, op_kl_result);
    op_result.NormInf(&nrm);
    op_result.Update(-1.0, op_kl_result, 1.0);
    op_result.NormInf(&diff);
    if (myPID == 0)
      std::cout << "Infinity norm of operator difference = " << diff/nrm
                << std::endl;
  }

#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                     "Stokhos::KLReducedMatrixFreeOperator is available " <<
                     "only when configured with Anasazi support!")
#endif
}
