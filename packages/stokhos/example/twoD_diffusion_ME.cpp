// $Id$ 
// $Source$ 
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

#include "twoD_diffusion_ME.hpp"
#include "Stokhos_Epetra.hpp"
#include "twoD_diffusion_inputs.h"
#include "twoD_diffusion_utility_functions.h"
#include <EpetraExt_MatrixMatrix.h>

#include <algorithm>
#include "Teuchos_TestForException.hpp"
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_EpetraMultiVectorOperator.hpp"

double xyLeft = -.5;
double xyRight = .5;
double meshSize;
int nel;
bool full_expn;
Teuchos::Array<double> x;
twoD_diffusion_ME::twoD_diffusion_ME(const Teuchos::RCP<Epetra_Comm>& comm, int n, int d, double s, double mu, bool full_expansion): A_k(d+1) 
{

/////////////////////////////////////////////////////////////////////////////////
//Construct the mesh.  The mesh is just the tensor of the below array with itself.
// The mesh is uniform and the nodes are numbered
// LEFT to RIGHT, DOWN to UP.
//
// 5-6-7-8-9
// | | | | |
// 0-1-2-3-4
/////////////////////////////////////////////////////////////////////////////////
//double xyLeft = -.5;
//double xyRight = .5;
double mesh_size = (xyRight - xyLeft)/((double)(n-1));
//Teuchos::Array<double> x(n);
x.resize(n);
for(int idx = 0; idx < n; idx++){
  x[idx] = xyLeft + (idx)*mesh_size;
}
int n2 = x.size()*x.size();
meshSize = x[1]-x[0];
nel = n;
full_expn = full_expansion;
////////////////////////////////////////////////////////////////////////
//Discretize the random field.
///////////////////////////////////////////////////////////////////////
sigma = s;
mean = mu;
lambda = Teuchos::Array<double>(d);
alpha = Teuchos::Array<double>(d);
omega = Teuchos::Array<double>(d);
xind = Teuchos::Array<int>(d);
yind = Teuchos::Array<int>(d);
generateExponentialRF(d, 1, lambda, alpha, omega, xind, yind);

  // Solution vector map
  x_map = Teuchos::rcp(new Epetra_Map(n*n, 0, *comm));

  // Overlapped solution vector map
  x_overlapped_map = Teuchos::rcp(new Epetra_Map(n*n, 0, *comm));

  // Importer
  importer = Teuchos::rcp(new Epetra_Import(*x_overlapped_map, *x_map));

  // Initial guess, initialized to 0.0
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->PutScalar(0.0);

  // Overlapped solution vector
  x_overlapped = Teuchos::rcp(new Epetra_Vector(*x_overlapped_map));

  // Parameter vector map
  p_map = Teuchos::rcp(new Epetra_LocalMap(d, 0, *comm));

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

  // Jacobian graph
int NumMyElements = x_map->NumMyElements();
int * MyGlobalElements = x_map->MyGlobalElements();
int * NumNz = new int[NumMyElements];

int *Indices = new int[4];
int NumEntries;
int * bcIndices = new int[NumMyElements];

  for(int i = 0; i<NumMyElements; i++){
    // MyGlobalElements[i]<x.size() ==> Boundary node on bottom edge.
    // MyGlobalElements[i]%x.size() == 0 ==> Boundary node on left edge.
    // MyGlobalElements[i]+1%x.size() == 0 ==> right edge.
    // MyGlobalElements[i] >= n - x.size() ==> top edge.

    if((MyGlobalElements[i] < static_cast<int>(x.size()) || MyGlobalElements[i]%x.size() == 0 ||
        (MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= n*n - static_cast<int>(x.size()))){
        NumNz[i] = 1;
      bcIndices[i] = 1;
    }else{
      NumNz[i] = 5;
      bcIndices[i] = 0;
    }
  }
  graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_map, NumNz));
  for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 0) {
      Indices[0] = MyGlobalElements[i]-x.size(); //Down
      Indices[1] = MyGlobalElements[i]-1;        //left
      Indices[2] = MyGlobalElements[i]+1;        //right
      Indices[3] = MyGlobalElements[i]+x.size(); //up
      NumEntries = 4;

  }
    if(bcIndices[i] == 0) graph->InsertGlobalIndices(MyGlobalElements[i], NumEntries, Indices);
    graph->InsertGlobalIndices(MyGlobalElements[i], 1, MyGlobalElements+i);
 }
  graph->FillComplete();
  graph->OptimizeStorage();
 
  std::string filename;
  std::string filename1("A");
  std::string filename2(".mm");

  double two;
  double *Values = new double[4];
  for (int k=0; k<=d; k++) {
    A_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
    for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 1 && k == 0) two = 1; //Enforce BC in mean matrix
    if (bcIndices[i] == 0) {
      Indices[0] = MyGlobalElements[i]-x.size(); //Down
      Indices[1] = MyGlobalElements[i]-1;        //left
      Indices[2] = MyGlobalElements[i]+1;        //right
      Indices[3] = MyGlobalElements[i]+x.size(); //up
      NumEntries = 4;
      Values[0] = -(1/(meshSize*meshSize))*evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k);
      Values[1] = -(1/(meshSize*meshSize))*evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()]-(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k);
      Values[2] = -(1/(meshSize*meshSize))*evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()]+(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k);
      Values[3] = -(1/(meshSize*meshSize))*evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]+(meshSize/2),k);

      two = (evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()]-(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k)
               + evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()]+(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k)
               + evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k)
               +evalEigenfunction(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]+(meshSize/2),k))
               /(meshSize*meshSize);
    }
   if(bcIndices[i] == 0) A_k[k]->ReplaceGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    if (bcIndices[i]==0 || k == 0) A_k[k]->ReplaceGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
  }
 A_k[k]->FillComplete();
 std::stringstream ss;
 ss << k ;
 std::string sstr = ss.str();
 filename = filename1 + sstr + filename2;
 //std::cout << ss.str() << std::endl;
 //EpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(), *(A_k[k]));
 }
 
  b = Teuchos::rcp(new Epetra_Vector(*x_map));
  for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 1 )
     (*b)[i] = 0;
    else 
     (*b)[i] = 1;
  }
/*  //Construct the RHS vector.
  Epetra_Vector b(*x_map,true);
  generateRHS(&RHS_function_PC, x, b,basis);
  b.Print(std::cout);
*/
  delete [] NumNz;
  delete [] Indices;
  delete [] bcIndices;
  delete [] Values;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::get_x_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::get_f_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

Teuchos::RCP<const Epetra_Map>
twoD_diffusion_ME::get_p_sg_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_sg_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
twoD_diffusion_ME::get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
twoD_diffusion_ME::get_p_sg_names(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_sg_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

Teuchos::RCP<const Epetra_Vector>
twoD_diffusion_ME::get_x_init() const
{
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
twoD_diffusion_ME::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_ME::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return p_init;
}

Teuchos::RCP<Epetra_Operator>
twoD_diffusion_ME::create_W() const
{
  Teuchos::RCP<Epetra_CrsMatrix> A = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  A->FillComplete();
  A->OptimizeStorage();
  return A;
}

EpetraExt::ModelEvaluator::InArgs
twoD_diffusion_ME::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription("TwoD Diffusion Model Evaluator");

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);    // 1 parameter vector

  // Stochastic InArgs
  inArgs.setSupports(IN_ARG_x_sg,true);
  inArgs.set_Np_sg(1); // 1 SG parameter vector
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
twoD_diffusion_ME::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription("TwoD Diffusion Model Evaluator");

  // Deterministic OutArgs
  outArgs.set_Np_Ng(1, 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  
  // Stochastic OutArgs
  outArgs.set_Np_Ng_sg(1, 0);
  outArgs.setSupports(OUT_ARG_f_sg,true);
  outArgs.setSupports(OUT_ARG_W_sg,true);

  return outArgs;
}

void 
twoD_diffusion_ME::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
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

  //
  // Stochastic Galerkin calculation
  //

  // Get stochastic expansion data
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    inArgs.get_sg_basis();
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expn = 
    inArgs.get_sg_expansion();

  Teuchos::RCP<const Stokhos::Sparse3Tensor<int, double> > Cijk;
  const Teuchos::Array<double> *norms;
  if (basis != Teuchos::null) {
    Cijk = expn->getTripleProduct(); 
    norms = &(basis->norm_squared());
  }

  // Compuet PC Coefficients of the operator B_k
  int d = basis->dimension(); 
  Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> > B_k(d+1);
  Teuchos::Array<double> zero(d), one(d);
  for(int k=0; k<d; k++) {
    zero[k] = 0.0;
    one[k] = 1.0;
  }
  for(int k=0; k<=d; k++) {
    B_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, A_k[k]->Graph()));
  }
  int sz = basis->size();
  Teuchos::Array< double > psi_0(sz), psi_1(sz);
  basis->evaluateBases(zero, psi_0);
  basis->evaluateBases(one, psi_1);
  *B_k[0] = *A_k[0];
  for(int k=1; k<=d; k++) {
    // B_0 -= a_k/(b_k-a_k)*A_k
    EpetraExt::MatrixMatrix::Add(*A_k[k], false,
   -psi_0[k]/(psi_1[k]-psi_0[k]),
                     *B_k[0], 1.0);
    // B_k = 1.0/(b_k-a_k)*A_k
    EpetraExt::MatrixMatrix::Add(*A_k[k], false, 1.0/(psi_1[k]-psi_0[k]),
                     *B_k[k], 0.0);
  }
  B_k[0]->Scale(1.0/psi_0[0]);

  // Compute eigen functions of lognormal RF 
  // from Maarten's matlab code
  Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis =
      Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int, double> >(basis, true);

  std::string filename;
  std::string filename1("C");
  std::string filename2(".mm");

  // Jacobian graph
  int NumMyElements = x_map->NumMyElements();
  int * MyGlobalElements = x_map->MyGlobalElements();
  int * NumNz = new int[NumMyElements];

  int NumEntries;
  int * bcIndices = new int[NumMyElements];
  double two;
  double *Values = new double[4];
  int *Indices = new int[4];
  int n2 = x.size()*x.size();
  meshSize = x[1]-x[0];
  for(int i = 0; i<NumMyElements; i++){
    // MyGlobalElements[i]<x.size() ==> Boundary node on bottom edge.
    // MyGlobalElements[i]%x.size() == 0 ==> Boundary node on left edge.
    // MyGlobalElements[i]+1%x.size() == 0 ==> right edge.
    // MyGlobalElements[i] >= n - x.size() ==> top edge.
    if((MyGlobalElements[i] < static_cast<int>(x.size()) || MyGlobalElements[i]%x.size() == 0 ||
        (MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= nel*nel - static_cast<int>(x.size()))){
        NumNz[i] = 1;
      bcIndices[i] = 1;
    }else{
      NumNz[i] = 5;
      bcIndices[i] = 0;
    }
  }
  Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> > C_k(sz);
  for (int k=0; k<sz; k++) {
    //C_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
    C_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, A_k[0]->Graph()));
    for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 1 && k == 0) two = 1; //Enforce BC in mean matrix
    if (bcIndices[i] == 0) {
      Indices[0] = MyGlobalElements[i]-x.size(); //Down
      Indices[1] = MyGlobalElements[i]-1;        //left
      Indices[2] = MyGlobalElements[i]+1;        //right
      Indices[3] = MyGlobalElements[i]+x.size(); //up
      NumEntries = 4;
      Values[0] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k,*prodbasis);
      Values[1] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()]-(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,*prodbasis);
      Values[2] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()]+(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,*prodbasis);
      Values[3] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]+(meshSize/2),k,*prodbasis);

      two = (evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()]-(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,*prodbasis)
               + evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()]+(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,*prodbasis)
               + evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k,*prodbasis)
               +evalPCECoefLogNormalRF(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]+(meshSize/2),k,*prodbasis))
               /(meshSize*meshSize);
/*      Values[0] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k,prodbasis);
      Values[1] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()]-(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,prodbasis);
      Values[2] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()]+(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,prodbasis);
      Values[3] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]+(meshSize/2),k,prodbasis);

      two = (evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()]-(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,prodbasis)
               + evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()]+(meshSize/2), x[(MyGlobalElements[i]%n2)/x.size()],k,prodbasis)
               + evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k,prodbasis)
               +evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]+(meshSize/2),k,prodbasis))
               /(meshSize*meshSize);*/
    }
   if(bcIndices[i] == 0) C_k[k]->ReplaceGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    if (bcIndices[i]==0 || k == 0) C_k[k]->ReplaceGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
  }
 C_k[k]->FillComplete();
 std::stringstream ss;
 ss << k ;
 std::string sstr = ss.str();
 filename = filename1 + sstr + filename2;
 //std::cout << ss.str() << std::endl;
 //EpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(), *(C_k[k]));
 }

  Teuchos::RCP<Epetra_Vector>b = Teuchos::rcp(new Epetra_Vector(*x_map));
  for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 1 )
     (*b)[i] = 0;
    else
     (*b)[i] = 1;
  }

  //Construct the RHS vector.
  //Teuchos::RCP<Epetra_Vector> b(*x_map,true);
 // generateRHS(&RHS_function_PC, x, *b,basis);
 // b->Print(std::cout);

  Teuchos::RCP<Epetra_CrsMatrix> A;
  // Residual  
  Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
  if (f != Teuchos::null || W != Teuchos::null) {
  A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  *A = *(B_k[0]);
  for (int k=1;k<B_k.size();k++) {
    EpetraExt::MatrixMatrix::Add((*B_k[k]), false, (*p)[k-1], *A, 1.0);
  }
//  EpetraExt::RowMatrixToMatrixMarketFile("A.mm", *A);
 }
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
  }

  // Stochastic solution vector
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();

  // Stochastic parameters
  InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(0);

  // Stochastic residual
  OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
  if (f_sg != Teuchos::null) {
  std::vector< Teuchos::RCP< Epetra_Vector> > sg_kx_vec_all ;
  for (int i=0;i<basis->size();i++) {
    sg_kx_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(*x_map)));
  }
  f_sg->init(0.0);

 double pc_size;
 if(!full_expn) 
  pc_size = basis->dimension()+1;
 else
  pc_size = basis->size(); 
  //for (int k=0; k<basis->dimension()+1; k++) {
  for (int k=0; k<pc_size; k++) {
    int nj = Cijk->num_j(k); 
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      //std::cout << "full_expn = " << full_expn << std::endl;
      if(!full_expn) 
        B_k[k]->Apply((*x_sg)[j],*(sg_kx_vec_all[j]));
      else
        C_k[k]->Apply((*x_sg)[j],*(sg_kx_vec_all[j]));
    }
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
      int ni = i_indices.size();
      for (int ii=0; ii<ni; ii++) {
        int i = i_indices[ii];
        double c = cijk_values[ii];  // C(i,j,k)
        (*f_sg)[i].Update(1.0*c/(*norms)[i],*(sg_kx_vec_all[j]),1.0);
      }
    }
  } //End 
  //std::cout << *f_sg <<std::endl;
  (*f_sg)[0].Update(-1.0,*b,1.0);
 }
 // Stochastic Jacobian
  OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
  if (W_sg != Teuchos::null) {
    for (int i=0; i<W_sg->size(); i++) {
      Teuchos::RCP<Epetra_CrsMatrix> jac = 
	Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(i), true);
      //*(jac) = *(B_k[i]);     
    if(!full_expn) 
       EpetraExt::MatrixMatrix::Add(*B_k[i], false, 1.0, *jac, 0.0);
     else 
       EpetraExt::MatrixMatrix::Add(*C_k[i], false, 1.0, *jac, 0.0);     
    }
    /*
    for (int i=0; i<W_sg->size(); i++) {
      Teuchos::RCP<Epetra_CrsMatrix> jac = 
	Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(i), true);
      jac->PutScalar(0.0); 
    }
    */
  }
}

// Compute PCE Coeff of LogNormal RF
double twoD_diffusion_ME::evalPCECoefLogNormalRF(double x, double y, int k,
               const Stokhos::ProductBasis<int, double>& basis) const { 
  Teuchos::Array<int> multiIndex;
  int d = basis.dimension();
  const Teuchos::Array<double> *norms;
  norms = &(basis.norm_squared());
  multiIndex = basis.getTerm(k);
  //std::cout << "multiIndex: ";
 // for (int l=0; l<d; l++)
 //   std::cout << multiIndex[l] << " ";
 // std::cout << std::endl;
  double sum_g = 0.0, efval;
  for (int l=0; l<d; l++) {
    //sum_g = sum_g + pow((sigma*evalEigenfunction(x,y,l)),2);
    sum_g = sum_g + pow((evalEigenfunction(x,y,l+1)),2);
  }
  efval = (exp(mean + 0.5*sum_g))/((*norms)[k]);
  //efval = (exp(0.5*sum_g))/((*norms)[k]);
  //std::cout << "norms[" <<k <<"] = " << (*norms)[k] << std::endl;
  for (int l=0; l<d; l++) {
    //efval *= sqrt(((*norms)[k])/(factorial(multiIndex[l])))*pow((sigma*evalEigenfunction(x,y,l)),multiIndex[l]); 
    //efval *= sqrt(((*norms)[k])/(factorial(multiIndex[l])))*pow((evalEigenfunction(x,y,l+1)),multiIndex[l]); 
    efval *= pow((evalEigenfunction(x,y,l+1)),multiIndex[l]); 
  }
  return efval;
}

// Compute PCE Coeff of LogNormal RF using qudrature.
double twoD_diffusion_ME::evalPCECoefLogNormalRFQuad(double x, double y, int k,
           const Teuchos::RCP<const Stokhos::ProductBasis<int, double> >& basis) const {
  //Stokhos::SparseGridQuadrature<int,double> quad(basis, level);
  //Stokhos::TensorProductQuadrature<int,double> quad(basis);
  //std::cout << "x y and k are " << x << " " << y << " " << k << std::endl;
  Teuchos::RCP<Stokhos::TensorProductQuadrature<int,double> >quad = 
    Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
  const Teuchos::Array<double>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<double> >& points =
        quad->getQuadPoints();
  const Teuchos::Array< Teuchos::Array<double> >& values =
        quad->getBasisAtQuadPoints();
  int nqp = weights.size();
  double efval = 0;
  for (int qp=0; qp<nqp; qp++) {
    //std::cout << "I am in qp loop" << endl;
    double klsum = 0;
    for (int l=0; l<basis->dimension(); l++) {
      klsum += evalEigenfunction(x,y,l+1)*points[qp][l];
    }
    efval += weights[qp]*(exp(mean+klsum))*values[qp][k];
    //efval += weights[qp]*(exp(klsum))*values[qp][k];
  }
  return efval;
} 
