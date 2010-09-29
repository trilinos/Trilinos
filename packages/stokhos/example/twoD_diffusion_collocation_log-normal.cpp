// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu)
// 
// ***********************************************************************
// @HEADER

#include <iostream>
#include <iomanip>
#include <fstream>

#include "Teuchos_RCP.hpp"
#include "Stokhos.hpp"
#include "AztecOO.h"
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MultiVector.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "twoD_diffusion_inputs.h"
#include "twoD_diffusion_utility_functions.h"
#include <Teuchos_Time.hpp>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_TimedEpetraOperator.hpp>
#include "EpetraExt_RowMatrixOut.h"


#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

double evalRF(double x, double y, Teuchos::Array<double> xi,
                       Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis);

// Compute PCE Coeff of LogNormal RF
double evalPCECoefLogNormalRF(double x, double y, int k,
               const Stokhos::ProductBasis<int, double>& basis) { 
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

int main(int argc, char **argv)
{

  try {
    Teuchos::Time TotalTimer("Basis Timer",false);
    TotalTimer.start();
    
    #ifdef HAVE_MPI
      MPI_Init(&argc, &argv);
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
    #else
      Epetra_SerialComm Comm;
    #endif

    //Parse the input arguments.
    //
      int n,level,d;
      if(argc < 7){
        n = 32; //Number of mesh points
        level = 5; //Polynomial degree
        d = 2;  //Terms in KL expansion
        sigma = .1;
        mean = .2;
        weightCut = 1;     // Support for distribution is +-weightCut
      }else{
        n = atoi(argv[1]);
        level = atoi(argv[2]);
        d = atoi(argv[3]);
        sigma = atof(argv[4]);
        mean = atof(argv[5]);
        weightCut = atof(argv[6]);
      }
    std::cout<< "sigma = " << sigma << " mean = " << mean << "\n";

    /////////////////////////////////////////////////////////////////////////////////
    //Construct the mesh.  The mesh is just the tensor of the below array with itself.
    // The mesh is uniform and the nodes are numbered
    // LEFT to RIGHT, DOWN to UP.
    //
    // 5-6-7-8-9
    // | | | | |
    // 0-1-2-3-4
    /////////////////////////////////////////////////////////////////////////////////
    double xyLeft = -.5;
    double xyRight = .5;
    double mesh_size = (xyRight - xyLeft)/((double)(n-1));
    Teuchos::Array<double> x(n);
    for(int idx = 0; idx < n; idx++){
      x[idx] = xyLeft + (idx)*mesh_size;
    }
    int n2 = x.size()*x.size();
    double meshSize = x[1]-x[0];

    /////////////////////////////////////////////////////////////////////////////////
    //Generate the polynomial chaos basis, the exponential random field and the 
    //sparse grid based on the zeros of the orthogonal polynomials.
    /////////////////////////////////////////////////////////////////////////////////
    const int p = level;
    const double leftEndPt = -weightCut;
    const double rightEndPt = weightCut;
    
    Teuchos::Time PolyBasisTimer("Poly Basis Timer",false);
    PolyBasisTimer.start();
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++) {
//      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p,true));
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p,true));
//      bases[i] = Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<int,double>("Beta",p,&weight,leftEndPt,rightEndPt,true));
    }
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    PolyBasisTimer.stop();
    
    //Generate the random field.
    lambda = Teuchos::Array<double>(d);
    alpha = Teuchos::Array<double>(d);
    omega = Teuchos::Array<double>(d);
    xind = Teuchos::Array<int>(d);
    yind = Teuchos::Array<int>(d);
    generateExponentialRF(d, 1,lambda, alpha, omega, xind, yind);
    
    //Construct the Sparse grid quadrature rule.
    Teuchos::Time SparseGridSetup("SP grid timer",false);
    SparseGridSetup.start();
    Stokhos::SparseGridQuadrature<int,double> my_quadrature(basis, level);
    //Stokhos::TensorProductQuadrature<int,double> my_quadrature(basis);
    Teuchos::Array< Teuchos::Array<double> > quad_points = my_quadrature.getQuadPoints();
    Teuchos::Array<double> quad_weights = my_quadrature.getQuadWeights();
    Teuchos::Array< Teuchos::Array<double> > quad_basis_values = my_quadrature.getBasisAtQuadPoints();
    SparseGridSetup.stop();


    //////////////////////////////////////////////////////////////////////
    //Generate the Deterministic stiffness matricies.
    ////////////////////////////////////////////////////////////////////// 
    Teuchos::RCP<Epetra_CrsMatrix> B_k;
    const Teuchos::RCP<const Epetra_Map> BaseMap = Teuchos::rcp(new Epetra_Map(n*n,0,Comm));
    int NumMyElements = (*BaseMap).NumMyElements();
    int * MyGlobalElements = (*BaseMap).MyGlobalElements();
    int * NumNz = new int[NumMyElements];

    double *Values = new double[4];
    int *Indices = new int[4];
    double two;
    int NumEntries;
    int * bcIndices = new int[NumMyElements];
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      if((MyGlobalElements[i] < static_cast<int>(x.size()) || MyGlobalElements[i]%x.size() == 0 || (MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= n*n - static_cast<int>(x.size()))){
      NumNz[i] = 1;
      bcIndices[i] = 1;
    }else{
      NumNz[i] = 5;
      bcIndices[i] = 0;
    }
    }
    Epetra_Vector b(*BaseMap,true);
    Epetra_Vector Eofu(*BaseMap,true);
    Epetra_Vector Eofu2(*BaseMap,true);
    Epetra_Vector varOfu(*BaseMap,true);
    Epetra_Vector EofUSquared(*BaseMap);
   
    Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> > A_k(d+1);
    Teuchos::Time AssemblyTimer("Matrix_assembly",false);
    AssemblyTimer.start();
    for(int k = 0; k<=d; k++){
      for(int i = 0; i<NumMyElements; i++){
        // MyGlobalElements[i]<x.size() ==> Boundary node on bottom edge.
        // MyGlobalElements[i]%x.size() == 0 ==> Boundary node on left edge.
        // MyGlobalElements[i]+1%x.size() == 0 ==> right edge.
        // MyGlobalElements[i] >= n - x.size() ==> top edge.

        if((MyGlobalElements[i] < static_cast<int>(x.size()) || MyGlobalElements[i]%x.size() == 0 || (MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= n*n - static_cast<int>(x.size()))){
          if(k==0){
            NumNz[i] = 1;
          } else NumNz[i] = 0;
          bcIndices[i] = 1;
        }else{
          NumNz[i] = 5;
          bcIndices[i] = 0;
        }
      }
      if(k==0) B_k = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*BaseMap,NumNz));
      A_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*BaseMap,NumNz));
  
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
        if(bcIndices[i] == 0){
          A_k[k]->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
        }
        if (bcIndices[i]==0 || k == 0) A_k[k]->InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
     }
     A_k[k]->FillComplete(); 
    }
     // Compute eigen functions of lognormal RF 
  // from Maarten's matlab code
  Teuchos::RCP<const Stokhos::ProductBasis<int, double> > prodbasis =
      Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int, double> >(basis, true);

  std::string filename;
  std::string filename1("C");
  std::string filename2(".mm");

/*  // Jacobian graph
  int NumMyElements = x_map->NumMyElements();
  int * MyGlobalElements = x_map->MyGlobalElements();
  int * NumNz = new int[NumMyElements];

  int NumEntries;
  int * bcIndices = new int[NumMyElements];
  double two;
  double *Values = new double[4];
  int *Indices = new int[4];
  int n2 = x.size()*x.size();
  meshSize = x[1]-x[0];*/
  for(int i = 0; i<NumMyElements; i++){
    // MyGlobalElements[i]<x.size() ==> Boundary node on bottom edge.
    // MyGlobalElements[i]%x.size() == 0 ==> Boundary node on left edge.
    // MyGlobalElements[i]+1%x.size() == 0 ==> right edge.
    // MyGlobalElements[i] >= n - x.size() ==> top edge.
    if((MyGlobalElements[i] < static_cast<int>(x.size()) || MyGlobalElements[i]%x.size() == 0 ||
        (MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= x.size()*x.size() - static_cast<int>(x.size()))){
        NumNz[i] = 1;
      bcIndices[i] = 1;
    }else{
      NumNz[i] = 5;
      bcIndices[i] = 0;
    }
  }
  Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> > C_k(basis->size());
  for (int k=0; k<basis->size(); k++) {
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
     /* Values[0] = -(1/(meshSize*meshSize))*evalPCECoefLogNormalRFQuad(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()]-(meshSize/2),k,prodbasis);
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
 EpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(), *(C_k[k]));
 }

/* // Teuchos::RCP<Epetra_Vector>b = Teuchos::rcp(new Epetra_Vector(*x_map));
  for( int i=0 ; i<NumMyElements; ++i ) {
    if (bcIndices[i] == 1 )
     b[i] = 0;
    else
     b[i] = 1;
  }
*/
    AssemblyTimer.stop();

    ///////////////////////////////////////////////////////////////////////
    //Construct the mean based preconditioner
    // create a parameter list for ML options
    //////////////////////////////////////////////////////////////////////
    Teuchos::Time PreconConstruct("Preconditioner_construct_time",false);
    PreconConstruct.start();
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("ML output", 10);
    MLList.set("max levels",5);
    MLList.set("increasing or decreasing","increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","ML symmetric Gauss-Seidel");
    MLList.set("smoother: sweeps",1);
    MLList.set("coarse: max size", 200);
    MLList.set("smoother: pre or post", "both");
    #ifdef HAVE_ML_AMESOS
      MLList.set("coarse: type","Amesos-KLU");
    #else
      MLList.set("coarse: type","Jacobi");
    #endif
    ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*C_k[0], MLList);
    PreconConstruct.stop();
  
    
    //////////////////////////////////////////////////////////////////////
    //Loop over sparse grid points and accumulate results for the mean
    //and variance of the solution.
    //////////////////////////////////////////////////////////////////////
    Teuchos::Array<double> xi;
    Epetra_Vector x2(*BaseMap);
    Epetra_Vector x22(*BaseMap);
    Teuchos::Time SolutionTimer("solution Timer", false);
    //Set up solver.
    const Teuchos::RCP<EpetraExt::Epetra_Timed_Operator> A = Teuchos::rcp(new EpetraExt::Epetra_Timed_Operator(B_k));
    Epetra_LinearProblem problem(&(*A),&x2,&b);
    AztecOO aztec_solver(problem);
    aztec_solver.SetPrecOperator(MLPrec);
    aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
    //aztec_solver.SetAztecOption(AZ_precond, AZ_none);
    aztec_solver.SetAztecOption(AZ_output, AZ_none);
    
    
    //LOOP OVER COLLOCATION POINTS.
    
    double coll_time = 0;
    int iters = 0;
    SolutionTimer.start();
    for( std::size_t sp_idx = 0; sp_idx < quad_points.size(); sp_idx++){
      AssemblyTimer.start();
      if(sp_idx%100 ==0) std::cout << sp_idx <<'/' << quad_points.size() << "\n";
      //xi = quad_points[sp_idx];
      B_k->PutScalar(0);
//        std::cout << "quad_points" << sp_idx << quad_points[sp_idx]<< endl;
        //std::cout << "quad_basis_values" << sp_idx << quad_basis_values[sp_idx]<< endl;
      //for(int k = 0; k<=d); k++){
      for(int k = 0; k<basis->size(); k++){
        //EpetraExt::MatrixMatrix::Add(*A_k[k],false,(k>0)?xi[k-1]:1,*B_k,1.0);
        //EpetraExt::MatrixMatrix::Add(*C_k[k],false,(k>0)?quad_basis_values[sp_idx][k-1]:1,*B_k,1.0);
        EpetraExt::MatrixMatrix::Add(*C_k[k],false,quad_basis_values[sp_idx][k],*B_k,1.0);
      }
      for(int i = 0; i<NumMyElements; i++){
         if(bcIndices[i] == 0){
          b[i] = RHS_function(x[(MyGlobalElements[i]%n2)%x.size()], x[(MyGlobalElements[i]%n2)/x.size()],quad_basis_values[sp_idx]);
         }else b[i] = 0;
      }
     /* //Teuchos::RCP<Epetra_Vector>b = Teuchos::rcp(new Epetra_Vector(*x_map));
      for( int i=0 ; i<NumMyElements; ++i ) {
        if (bcIndices[i] == 1 )
         b[i] = 0;
        else
         b[i] = 1;
      }*/
 
      B_k->FillComplete();
      AssemblyTimer.stop();
      //////////////////////////////////////////////////
      //Solve the Linear System.
      /////////////////////////////////////////////////
      {
       //TEUCHOS_FUNC_TIME_MONITOR("Total Collocation Solve Time"); 
       Teuchos::Time SolveTimer("Total Collocation Solve Time",false);
       SolveTimer.start();
       aztec_solver.Iterate(1000, 1e-12);
       SolveTimer.stop();
       //std::cout << "Solve time at qp " << sp_idx <<" is = "<< SolveTimer.totalElapsedTime(false) << std::endl;
       coll_time += SolveTimer.totalElapsedTime(false);
      }
      iters = iters + aztec_solver.NumIters();
      Eofu.Update(quad_weights[sp_idx],x2,1.0);
      x22.Multiply(1.0,x2,x2,0.0);
      Eofu2.Update(quad_weights[sp_idx], x22, 1.0);
    }
    SolutionTimer.stop();
    TotalTimer.stop();    
    std::cout << "Total Collocation Solve Time is =  " << coll_time << std::endl;
    ///////////////////////////////////////////////////////////
    //Output timings and results.
    //////////////////////////////////////////////////////////
    Teuchos::ParameterList ML_Output = MLPrec->GetOutputList();
    double precon_time = ML_Output.get("time: total apply",-2000.0);
    varOfu.Multiply(-1.0, Eofu, Eofu, 0.0);
    varOfu.Update(1.0,Eofu2,1.0);
    std::ofstream mean;
    mean.open("sp_mean.txt");
    Eofu.Print(mean);
    mean.close();

    std::ofstream var;
    var.open("sp_var.txt");
    varOfu.Print(var);
    var.close();
    
    std::ofstream time;
    time.open("sp_time.txt");
    time <<"TotalTimer= " <<TotalTimer.totalElapsedTime(false) << "\n";
    time <<"SolutionTimer= " <<SolutionTimer.totalElapsedTime(false) << "\n";
    time <<"SolveTimer= " <<coll_time << "\n";
    time <<"precon_time= " <<precon_time << "\n";
    time <<"A->ApplyTime= " <<A->ApplyTime()<< "\n";
    time <<"SparseGridSetuptime= " <<SparseGridSetup.totalElapsedTime(false) << "\n";
    time.close();

    std::ofstream dof;
    dof.open("sp_dof.txt");
    dof << quad_points.size() << "\n";
    dof.close();

    std::ofstream iters_file;
    iters_file.open("sp_iters.txt");
    iters_file << iters << "\n";
    iters_file.close();

    std::cout << "\n Polynomial Basis Cosntruction time = " << PolyBasisTimer.totalElapsedTime(false)<< "\n";
    std::cout << "Sparse Grid Setup Time = " << SparseGridSetup.totalElapsedTime(false)<< "\n";
    std::cout << " Preconsitioner Construction time = " << PreconConstruct.totalElapsedTime(false) << "\n";
    std::cout << "Stiffness Matrix Assembly time = " << AssemblyTimer.totalElapsedTime(false) << "\n";
    std::cout << "System Solution time = " << SolutionTimer.totalElapsedTime(false) << "\n";
    std::cout << "Total Running Time = " << TotalTimer.totalElapsedTime(false) << "\n";
    
  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

    
}



  

