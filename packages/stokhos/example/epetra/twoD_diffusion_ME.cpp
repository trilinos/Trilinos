// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "twoD_diffusion_ME.hpp"

#include "Stokhos_KL_ExponentialRandomField.hpp"
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_Assert.hpp"
#include "Stokhos_PreconditionerFactory.hpp"

namespace {

  // Functor representing a diffusion function given by a KL expansion
  struct KL_Diffusion_Func {
    double mean;
    mutable Teuchos::Array<double> point;
    Teuchos::RCP< Stokhos::KL::ExponentialRandomField<double, Kokkos::DefaultHostExecutionSpace> > rf;

    KL_Diffusion_Func(double xyLeft, double xyRight,
                      double mean_, double std_dev,
                      double L, int num_KL) : mean(mean_), point(2)
    {
      Teuchos::ParameterList rfParams;
      rfParams.set("Number of KL Terms", num_KL);
      rfParams.set("Mean", mean);
      rfParams.set("Standard Deviation", std_dev);
      int ndim = 2;
      Teuchos::Array<double> domain_upper(ndim), domain_lower(ndim),
        correlation_length(ndim);
      for (int i=0; i<ndim; i++) {
        domain_upper[i] = xyRight;
        domain_lower[i] = xyLeft;
        correlation_length[i] = L;
      }
      rfParams.set("Domain Upper Bounds", domain_upper);
      rfParams.set("Domain Lower Bounds", domain_lower);
      rfParams.set("Correlation Lengths", correlation_length);

      rf =
        Teuchos::rcp(new Stokhos::KL::ExponentialRandomField<double, Kokkos::DefaultHostExecutionSpace>(rfParams));
    }

    ~KL_Diffusion_Func() {
    }

    double operator() (double x, double y, int k) const {
      if (k == 0)
        return mean;
      point[0] = x;
      point[1] = y;
      return rf->evaluate_eigenfunction(point, k-1);
    }
  };

  // Functor representing a diffusion function given by a KL expansion
  // where the basis is normalized
  template <typename DiffusionFunc>
  struct Normalized_KL_Diffusion_Func {
    const DiffusionFunc& func;
    int d;
    Teuchos::Array<double> psi_0, psi_1;

    Normalized_KL_Diffusion_Func(
      const DiffusionFunc& func_,
      const Stokhos::OrthogPolyBasis<int,double>& basis) :
      func(func_),
      d(basis.dimension()),
      psi_0(basis.size()),
      psi_1(basis.size())
    {
      Teuchos::Array<double> zero(d), one(d);
      for(int k=0; k<d; k++) {
        zero[k] = 0.0;
        one[k] = 1.0;
      }
      basis.evaluateBases(zero, psi_0);
      basis.evaluateBases(one, psi_1);
    }

    double operator() (double x, double y, int k) const {
      if (k == 0) {
        double val = func(x, y, 0);
        for (int i=1; i<=d; i++)
          val -= psi_0[i]/(psi_1[i]-psi_0[i])*func(x, y, i);
        val /= psi_0[0];
        return val;
      }
      else
        return 1.0/(psi_1[k]-psi_0[k])*func(x, y, k);
    }
  };

  // Functor representing a diffusion function given by a log-normal PC expansion
  template <typename DiffusionFunc>
  struct LogNormal_Diffusion_Func {
    double mean;
    const DiffusionFunc& func;
    Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis;

    LogNormal_Diffusion_Func(
      double mean_,
      const DiffusionFunc& func_,
      const Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis_)
      : mean(mean_), func(func_), prodbasis(prodbasis_) {}

    double operator() (double x, double y, int k) const {
      int d = prodbasis->dimension();
      const Teuchos::Array<double>& norms = prodbasis->norm_squared();
      const Stokhos::MultiIndex<int> multiIndex = prodbasis->term(k);
      double sum_g = 0.0, efval;
      for (int l=0; l<d; l++) {
        sum_g += std::pow(func(x,y,l+1),2);
      }
      efval = std::exp(mean + 0.5*sum_g)/norms[k];
      for (int l=0; l<d; l++) {
        efval *= std::pow(func(x,y,l+1),multiIndex[l]);
      }
      return efval;
    }
  };

}

twoD_diffusion_ME::
twoD_diffusion_ME(
  const Teuchos::RCP<const Epetra_Comm>& comm, int n, int d,
  double s, double mu,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
  bool log_normal_,
  bool eliminate_bcs_,
  const Teuchos::RCP<Teuchos::ParameterList>& precParams_) :
  mesh(n*n),
  basis(basis_),
  log_normal(log_normal_),
  eliminate_bcs(eliminate_bcs_),
  precParams(precParams_)
{
  Kokkos::initialize();

  //////////////////////////////////////////////////////////////////////////////
  // Construct the mesh.
  // The mesh is uniform and the nodes are numbered
  // LEFT to RIGHT, DOWN to UP.
  //
  // 5-6-7-8-9
  // | | | | |
  // 0-1-2-3-4
  /////////////////////////////////////////////////////////////////////////////
  double xyLeft = -.5;
  double xyRight = .5;
  h = (xyRight - xyLeft)/((double)(n-1));
  Teuchos::Array<int> global_dof_indices;
  for (int j=0; j<n; j++) {
    double y = xyLeft + j*h;
    for (int i=0; i<n; i++) {
      double x = xyLeft + i*h;
      int idx = j*n+i;
      mesh[idx].x = x;
      mesh[idx].y = y;
      if (i == 0 || i == n-1 || j == 0 || j == n-1)
        mesh[idx].boundary = true;
      if (i != 0)
        mesh[idx].left = idx-1;
      if (i != n-1)
        mesh[idx].right = idx+1;
      if (j != 0)
        mesh[idx].down = idx-n;
      if (j != n-1)
        mesh[idx].up = idx+n;
      if (!(eliminate_bcs && mesh[idx].boundary))
        global_dof_indices.push_back(idx);
    }
  }

  // Solution vector map
  int n_global_dof = global_dof_indices.size();
  int n_proc = comm->NumProc();
  int proc_id = comm->MyPID();
  int n_my_dof = n_global_dof / n_proc;
  if (proc_id == n_proc-1)
    n_my_dof += n_global_dof % n_proc;
  int *my_dof = global_dof_indices.getRawPtr() + proc_id*(n_global_dof / n_proc);
  x_map =
    Teuchos::rcp(new Epetra_Map(n_global_dof, n_my_dof, my_dof, 0, *comm));

  // Initial guess, initialized to 0.0
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->PutScalar(0.0);

  // Parameter vector map
  p_map = Teuchos::rcp(new Epetra_LocalMap(d, 0, *comm));

  // Response vector map
  g_map = Teuchos::rcp(new Epetra_LocalMap(1, 0, *comm));

  // Initial parameters
  p_init = Teuchos::rcp(new Epetra_Vector(*p_map));
  p_init->PutScalar(0.0);

  // Parameter names
  p_names = Teuchos::rcp(new Teuchos::Array<std::string>(d));
  for (int i=0;i<d;i++) {
    std::stringstream ss;
    ss << "KL Random Variable " << i+1;
    (*p_names)[i] = ss.str();
  }

  // Build Jacobian graph
  int NumMyElements = x_map->NumMyElements();
  int *MyGlobalElements = x_map->MyGlobalElements();
  graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_map, 5));
  for (int i=0; i<NumMyElements; ++i ) {

    // Center
    int global_idx = MyGlobalElements[i];
    graph->InsertGlobalIndices(global_idx, 1, &global_idx);

    if (!mesh[global_idx].boundary) {
      // Down
      if (!(eliminate_bcs && mesh[mesh[global_idx].down].boundary))
        graph->InsertGlobalIndices(global_idx, 1, &mesh[global_idx].down);

      // Left
      if (!(eliminate_bcs && mesh[mesh[global_idx].left].boundary))
        graph->InsertGlobalIndices(global_idx, 1, &mesh[global_idx].left);

      // Right
      if (!(eliminate_bcs && mesh[mesh[global_idx].right].boundary))
        graph->InsertGlobalIndices(global_idx, 1, &mesh[global_idx].right);

      // Up
      if (!(eliminate_bcs && mesh[mesh[global_idx].up].boundary))
        graph->InsertGlobalIndices(global_idx, 1, &mesh[global_idx].up);
    }
  }
  graph->FillComplete();
  graph->OptimizeStorage();

  KL_Diffusion_Func klFunc(xyLeft, xyRight, mu, s, 1.0, d);
  if (!log_normal) {
    // Fill coefficients of KL expansion of operator
    if (basis == Teuchos::null) {
      fillMatrices(klFunc, d+1);
    }
    else {
      Normalized_KL_Diffusion_Func<KL_Diffusion_Func> nklFunc(klFunc, *basis);
      fillMatrices(nklFunc, d+1);
    }
  }
  else {
    // Fill coefficients of PC expansion of operator
    int sz = basis->size();
    Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis =
      Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int, double> >(
        basis, true);
    LogNormal_Diffusion_Func<KL_Diffusion_Func> lnFunc(mu, klFunc, prodbasis);
    fillMatrices(lnFunc, sz);
  }

  // Construct deterministic operator
  A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));

  // Construct the RHS vector.
  b = Teuchos::rcp(new Epetra_Vector(*x_map));
  for( int i=0 ; i<NumMyElements; ++i ) {
    int global_idx = MyGlobalElements[i];
    if (mesh[global_idx].boundary)
      (*b)[i] = 0;
    else
      (*b)[i] = 1;
  }

  if (basis != Teuchos::null) {
    point.resize(d);
    basis_vals.resize(basis->size());
  }

  if (precParams != Teuchos::null) {
    std::string name = precParams->get("Preconditioner Type", "Ifpack");
    Teuchos::RCP<Teuchos::ParameterList> p =
      Teuchos::rcp(&(precParams->sublist("Preconditioner Parameters")), false);
    precFactory =
      Teuchos::rcp(new Stokhos::PreconditionerFactory(name, p));
  }
}

twoD_diffusion_ME::
~twoD_diffusion_ME()
{
  Kokkos::finalize();
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_x_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_f_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0,
                     std::logic_error,
                     std::endl <<
                     "Error!  twoD_diffusion_ME::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::
get_g_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0,
                     std::logic_error,
                     std::endl <<
                     "Error!  twoD_diffusion_ME::get_g_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return g_map;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
twoD_diffusion_ME::
get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0,
                     std::logic_error,
                     std::endl <<
                     "Error!  twoD_diffusion_ME::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

Teuchos::RCP<const Epetra_Vector>
twoD_diffusion_ME::
get_x_init() const
{
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
twoD_diffusion_ME::
get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0,
                     std::logic_error,
                     std::endl <<
                     "Error!  twoD_diffusion_ME::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_init;
}

Teuchos::RCP<Epetra_Operator>
twoD_diffusion_ME::
create_W() const
{
  Teuchos::RCP<Epetra_CrsMatrix> AA =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  AA->FillComplete();
  AA->OptimizeStorage();
  return AA;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
twoD_diffusion_ME::
create_WPrec() const
{
  if (precFactory != Teuchos::null) {
    Teuchos::RCP<Epetra_Operator> precOp =
      precFactory->compute(A, false);
    return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,
                                                                      true));
  }
  return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
twoD_diffusion_ME::
createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription("TwoD Diffusion Model Evaluator");

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);    // 1 parameter vector

  // Stochastic InArgs
  inArgs.setSupports(IN_ARG_x_sg,true);
  inArgs.setSupports(IN_ARG_p_sg, 0, true); // 1 SG parameter vector
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);

  // Multipoint InArgs
  inArgs.setSupports(IN_ARG_x_mp,true);
 inArgs.setSupports(IN_ARG_p_mp, 0, true); // 1 MP parameter vector

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
twoD_diffusion_ME::
createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription("TwoD Diffusion Model Evaluator");

  // Deterministic OutArgs
  outArgs.set_Np_Ng(1, 1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  if (precFactory != Teuchos::null)
    outArgs.setSupports(OUT_ARG_WPrec,true);

  // Stochastic OutArgs
  outArgs.setSupports(OUT_ARG_f_sg,true);
  outArgs.setSupports(OUT_ARG_W_sg,true);
  outArgs.setSupports(OUT_ARG_g_sg, 0, true);

  // Multipoint OutArgs
  outArgs.setSupports(OUT_ARG_f_mp,true);
  outArgs.setSupports(OUT_ARG_W_mp,true);
  outArgs.setSupports(OUT_ARG_g_mp, 0, true);

  return outArgs;
}

void
twoD_diffusion_ME::
evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{

  //
  // Determinisic calculation
  //

  // Solution vector
  Teuchos::RCP<const Epetra_Vector> det_x = inArgs.get_x();

  // Parameters
  Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(0);
  if (p == Teuchos::null)
    p = p_init;

  Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
  Teuchos::RCP<Epetra_Operator> WPrec = outArgs.get_WPrec();
  if (f != Teuchos::null || W != Teuchos::null || WPrec != Teuchos::null) {
    if (basis != Teuchos::null) {
      for (int i=0; i<point.size(); i++)
        point[i] = (*p)[i];
      basis->evaluateBases(point, basis_vals);
      A->PutScalar(0.0);
      for (int k=0;k<A_k.size();k++)
        EpetraExt::MatrixMatrix::Add((*A_k[k]), false, basis_vals[k], *A, 1.0);
    }
    else {
      *A = *(A_k[0]);
      for (int k=1;k<A_k.size();k++)
        EpetraExt::MatrixMatrix::Add((*A_k[k]), false, (*p)[k-1], *A, 1.0);
    }
    A->FillComplete();
    A->OptimizeStorage();
  }

  // Residual
  if (f != Teuchos::null) {
    Teuchos::RCP<Epetra_Vector> kx = Teuchos::rcp(new Epetra_Vector(*x_map));
    A->Apply(*det_x,*kx);
    f->Update(1.0,*kx,-1.0, *b, 0.0);
  }

  // Jacobian
  if (W != Teuchos::null) {
    Teuchos::RCP<Epetra_CrsMatrix> jac =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
    *jac = *A;
    jac->FillComplete();
    jac->OptimizeStorage();
  }

  // Preconditioner
  if (WPrec != Teuchos::null)
    precFactory->recompute(A, WPrec);

  // Responses (mean value)
  Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(0);
  if (g != Teuchos::null) {
    (det_x->MeanValue(&(*g)[0]));
    (*g)[0] *= double(det_x->GlobalLength()) / double(mesh.size());
  }

  //
  // Stochastic Galerkin calculation
  //

  // Stochastic solution vector
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();

  // Stochastic parameters
  InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(0);

  // Stochastic residual
  OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
  if (f_sg != Teuchos::null) {

    // Get stochastic expansion data
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expn =
      inArgs.get_sg_expansion();
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
    Teuchos::RCP<const Cijk_type> Cijk = expn->getTripleProduct();
    const Teuchos::Array<double>& norms = basis->norm_squared();

    if (sg_kx_vec_all.size() != basis->size()) {
      sg_kx_vec_all.resize(basis->size());
      for (int i=0;i<basis->size();i++) {
        sg_kx_vec_all[i] = Teuchos::rcp(new Epetra_Vector(*x_map));
      }
    }
    f_sg->init(0.0);

    Cijk_type::k_iterator k_begin = Cijk->k_begin();
    Cijk_type::k_iterator k_end = Cijk->k_end();
    for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int k = Stokhos::index(k_it);
      for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it);
           j_it != Cijk->j_end(k_it); ++j_it) {
        int j = Stokhos::index(j_it);
        A_k[k]->Apply((*x_sg)[j],*(sg_kx_vec_all[j]));
      }
      for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it);
           j_it != Cijk->j_end(k_it); ++j_it) {
        int j = Stokhos::index(j_it);
        for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
             i_it != Cijk->i_end(j_it); ++i_it) {
          int i = Stokhos::index(i_it);
          double c = Stokhos::value(i_it);  // C(i,j,k)
          (*f_sg)[i].Update(1.0*c/norms[i],*(sg_kx_vec_all[j]),1.0);
        }
      }
    } //End
    (*f_sg)[0].Update(-1.0,*b,1.0);
  }

  // Stochastic Jacobian
  OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
  if (W_sg != Teuchos::null) {
    for (int i=0; i<W_sg->size(); i++) {
      Teuchos::RCP<Epetra_CrsMatrix> jac =
        Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(i), true);
      *jac = *A_k[i];
      jac->FillComplete();
      jac->OptimizeStorage();
    }
  }

  // Stochastic responses
  Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > g_sg =
    outArgs.get_g_sg(0);
  if (g_sg != Teuchos::null) {
    int sz = x_sg->size();
    for (int i=0; i<sz; i++) {
      (*x_sg)[i].MeanValue(&(*g_sg)[i][0]);
      (*g_sg)[i][0] *= double((*x_sg)[i].GlobalLength()) / double(mesh.size());
    }
  }

  //
  // Multi-point calculation
  //

  // Stochastic solution vector
  mp_const_vector_t x_mp = inArgs.get_x_mp();

  // Stochastic parameters
  mp_const_vector_t p_mp = inArgs.get_p_mp(0);

  // Stochastic residual
  mp_vector_t f_mp = outArgs.get_f_mp();
  mp_operator_t W_mp = outArgs.get_W_mp();
  if (f_mp != Teuchos::null || W_mp != Teuchos::null) {
    int num_mp = x_mp->size();
    for (int i=0; i<num_mp; i++) {
      // Compute operator
      if (basis != Teuchos::null) {
        for (int k=0; k<point.size(); k++)
          point[k] = (*p_mp)[i][k];
        basis->evaluateBases(point, basis_vals);
        A->PutScalar(0.0);
        for (int k=0;k<A_k.size();k++)
          EpetraExt::MatrixMatrix::Add((*A_k[k]), false, basis_vals[k], *A,
                                       1.0);
      }
      else {
        *A = *(A_k[0]);
        for (int k=1;k<A_k.size();k++)
          EpetraExt::MatrixMatrix::Add((*A_k[k]), false, (*p_mp)[i][k-1], *A,
                                       1.0);
      }

      A->FillComplete();
      A->OptimizeStorage();

      // Compute residual
      if (f_mp != Teuchos::null) {
        A->Apply((*x_mp)[i], (*f_mp)[i]);
        (*f_mp)[i].Update(-1.0, *b, 1.0);
      }

      // Copy operator
      if (W_mp != Teuchos::null) {
        Teuchos::RCP<Epetra_CrsMatrix> jac =
          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_mp->getCoeffPtr(i),
                                                      true);
        *jac = *A;
        jac->FillComplete();
        jac->OptimizeStorage();
      }
    }
  }

  // Multipoint responses
  mp_vector_t g_mp = outArgs.get_g_mp(0);
  if (g_mp != Teuchos::null) {
    int sz = x_mp->size();
    for (int i=0; i<sz; i++) {
      (*x_mp)[i].MeanValue(&(*g_mp)[i][0]);
      (*g_mp)[i][0] *= double((*x_mp)[i].GlobalLength()) / double(mesh.size());
    }
  }

}

