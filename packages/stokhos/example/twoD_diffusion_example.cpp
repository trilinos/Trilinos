// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu).
// 
// ***********************************************************************
// @HEADER

// Solves a stochastic two dimensional stochastic diffusion equation in
// [-.5,.5]^2 using stochastic finite elements in the stochastic dimension
// and finite differences in the spatial dimension.  The linear solve is
// accomplished using a algebraic multigrid preconditioner based on the 
// mean stiffness matrix.  The matrix vector product is performed implicitly
// taking advantage of the tensor structure of the global stiffness matrix. 
//
// The diffusion coefficient and RHS function are defined in the file,
// twoD_diffusion_inputs.h.  Some auxillary functions are in the file,
// twoD_diffusion_utility_functions.h.
//
// At output the program writes the solution vector, the mean solution, the variance
// of the solution, and the solution realization to text files.

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

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
#include <EpetraExt_TimedEpetraOperator.hpp>
#include "Epetra_MultiVector.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include <Teuchos_Time.hpp>

#include "twoD_diffusion_inputs.h"
#include "twoD_diffusion_utility_functions.h"



#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif




int main(int argc, char **argv)
{

Teuchos::Time TotalTimer("Total Timer",false);
TotalTimer.start();

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

//Parse the input arguments.
//
int n, p, d;
//char polyType;
//double * evaluationPoint = new double(d);
if(argc < 7){
  n = 32; //Number of mesh points
  p = 5; //Polynomial degree
  d = 2;  //Terms in KL expansion
  sigma = .1;
  mean = .2;
  weightCut = 1;     // Support for distribution is +-weightCut
}else{
  n = atoi(argv[1]);
  p = atoi(argv[2]);
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
//Generate the polynomial chaos basis.
/////////////////////////////////////////////////////////////////////////////////
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
for (int i = 0; i< d; i++){
  bases[i] = Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<int,double>("beta",p,&weight,-weightCut,weightCut,true));
}

Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis =
  Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = basis->getLowOrderTripleProductTensor(d+1);



////////////////////////////////////////////////////////////////////////
//Discretize the random field.
///////////////////////////////////////////////////////////////////////
lambda = Teuchos::Array<double>(d);
alpha = Teuchos::Array<double>(d);
omega = Teuchos::Array<double>(d);
xind = Teuchos::Array<int>(d);
yind = Teuchos::Array<int>(d);
generateExponentialRF(d, 1,lambda, alpha, omega, xind, yind);


//////////////////////////////////////////////////////////////////////
//Generate the Deterministic stiffness matricies and pack them into
//the implicit MAT-VEC operator
////////////////////////////////////////////////////////////////////// 

Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> > A_k(d+1);
const Teuchos::RCP<const Epetra_Map> StochMap = Teuchos::rcp(new Epetra_Map(n*n*basis->size(),0,Comm));
const Teuchos::RCP<const Epetra_Map> BaseMap = Teuchos::rcp(new Epetra_Map(n*n,0,Comm));
int NumMyElements = (*BaseMap).NumMyElements();
int * MyGlobalElements = (*BaseMap).MyGlobalElements();
int * NumNz = new int[NumMyElements];
  
double *Values = new double[4];
int *Indices = new int[4];
double two;
int NumEntries;
int * bcIndices = new int[NumMyElements]; 

for(int k = 0; k<=d; k++){
  for(int i = 0; i<NumMyElements; i++){
    // MyGlobalElements[i]<x.size() ==> Boundary node on bottom edge.
    // MyGlobalElements[i]%x.size() == 0 ==> Boundary node on left edge.
    // MyGlobalElements[i]+1%x.size() == 0 ==> right edge.
    // MyGlobalElements[i] >= n - x.size() ==> top edge.

    if((MyGlobalElements[i] < static_cast<int>(x.size()) || MyGlobalElements[i]%x.size() == 0 ||
	(MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= n*n - static_cast<int>(x.size()))){
      if(k==0){
        NumNz[i] = 1;
      } else NumNz[i] = 0;
      bcIndices[i] = 1;
    }else{
      NumNz[i] = 5;
      bcIndices[i] = 0;
    }
  }
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
    if(bcIndices[i] == 0) A_k[k]->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    if (bcIndices[i]==0 || k == 0) A_k[k]->InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
 }
 A_k[k]->FillComplete();
}

//Construct the implicit operator.
Teuchos::RCP<Stokhos::KLMatrixFreeEpetraOp> MatFreeOp = 
   Teuchos::rcp(new Stokhos::KLMatrixFreeEpetraOp(BaseMap, StochMap, basis, Cijk, A_k)); 
EpetraExt::Epetra_Timed_Operator system(MatFreeOp);

//Construct the RHS vector.
Epetra_Vector b(*StochMap,true);
generateRHS(&RHS_function_PC, x, b,basis);

///////////////////////////////////////////////////////////////////////
//Construct the mean based preconditioner
// create a parameter list for ML options
//////////////////////////////////////////////////////////////////////
Teuchos::ParameterList MLList;
ML_Epetra::SetDefaults("SA",MLList);
MLList.set("ML output", 10);
MLList.set("max levels",5);
MLList.set("increasing or decreasing","increasing");
MLList.set("aggregation: type", "Uncoupled");
MLList.set("smoother: type","ML symmetric Gauss-Seidel");
MLList.set("smoother: sweeps",1);
MLList.set("smoother: pre or post", "both");
MLList.set("coarse: max size", 200);
#ifdef HAVE_ML_AMESOS
  MLList.set("coarse: type","Amesos-KLU");
#else
  MLList.set("coarse: type","Jacobi");
#endif
Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec = 
  Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner((*A_k[0]), MLList));
MLPrec->PrintUnused(0);
Stokhos::MeanEpetraOp MeanAMGPrecon(BaseMap, StochMap, basis->size(), MLPrec);

//////////////////////////////////////////////////
//Solve the Linear System.
/////////////////////////////////////////////////
Epetra_Vector x2(*StochMap);
Epetra_LinearProblem problem(&system, &x2, &b);
AztecOO aztec_solver(problem);
aztec_solver.SetPrecOperator(&MeanAMGPrecon);
aztec_solver.SetAztecOption(AZ_solver, AZ_cg);
//aztec_solver.SetAztecOption(AZ_precond, AZ_none);
Teuchos::Time SolutionTimer("Total Timer",false);
SolutionTimer.start();
aztec_solver.Iterate(1000, 1e-12);
SolutionTimer.stop();

//////////////////////////////////////////////////////////////////////
//Post process and output the results.
////////////////////////////////////////////////////////////////////
Epetra_Vector Eofu(*BaseMap,true);
Epetra_Vector varOfu(*BaseMap,true);
Epetra_Vector specificSol(*BaseMap,true);
computeMeanSolution(x2, Eofu, basis);
computeVarianceSolution(x2, Eofu, varOfu, basis);
TotalTimer.stop();

std::ofstream mean;
mean.open("mean.txt");
Eofu.Print(mean);
mean.close();

std::ofstream var;
var.open("var.txt");
varOfu.Print(var);
var.close();

std::ofstream time;
time.open("gal_time.txt");
Teuchos::ParameterList ML_Output = MLPrec->GetOutputList();
double precon_time = ML_Output.get("time: total apply",-2000.0);
time << TotalTimer.totalElapsedTime(false) << "\n";
time << SolutionTimer.totalElapsedTime(false) << "\n";
time << precon_time << "\n";
time << system.ApplyTime()<< "\n";
time.close();

std::ofstream dof;
dof.open("gal_dof.txt");
dof << basis->size() << "\n";
dof.close();

int iters = aztec_solver.NumIters();
std::ofstream iters_file;
iters_file.open("gal_iters.txt");
iters_file << iters << "\n";
iters_file.close();

std::cout << "\nFinished! Results written out to 'mean.txt', 'gal_time.txt', 'gal_dof.txt', and 'gal_iters.txt' '\n";



return 1;
}
 








