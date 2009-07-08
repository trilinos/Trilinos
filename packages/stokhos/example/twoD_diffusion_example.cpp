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
#include "Epetra_MultiVector.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "twoD_diffusion_inputs.h"
#include "twoD_diffusion_utility_functions.h"



#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif


int main(int argc, char **argv)
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

//Parse the input arguments.
//
int n, p, d;
char polyType;
double * evaluationPoint = new double(d);
if(argc < 8){
  std::cout<< "Usage is: Stokhos_twoD_diffusion_example.cpp <# meshPoints> <PC Degree> <PC Dimension> <sigma> <mu> <PC Type| 'r'(Rys) or 'h'(Hermite)> <Gaussian Cut (only applies to Rys)> <Evaluation Point>\n";
  std:cout << "assuming: Stokhos_twoD_diffusion_example.cpp 32 5 1 3 1 r 3 [0,0,...,0]\n";
  n = 32; //Number of mesh points
  p = 10; //Polynomial degree
  d = 1;  //Terms in KL expansion
  sigma = 3;
  mean = 1;
  polyType = 'r'; // Rys or Hermite?
  rysCut = 3;     // cutoff for truncated gaussian
}else{
  n = atoi(argv[1]);
  p = atoi(argv[2]);
  d = atoi(argv[3]);
  sigma = atof(argv[4]);
  mean = atof(argv[5]);
  polyType = (argv[6])[0];
  rysCut = atof(argv[7]);
}

//Specify where to evaluate the solution process.
for(int i = 0; i<d; i++){
  if(argc>= 9+ i){
    evaluationPoint[i] = atof(argv[8+i]);
  }else{evaluationPoint[i] = 0;}
}

double xyLeft = -.5;
double xyRight = .5;


std::cout<< "sigma = " << sigma << " mean = " << mean << "\n";

//Construct the mesh.  The mesh is just the tensor of the below array with itself.
// The mesh is uniform and the nodes are numbered
// LEFT to RIGHT, DOWN to UP.
//
// 5-6-7-8-9
// | | | | |
// 0-1-2-3-4
double mesh_size = (xyRight - xyLeft)/((double)(n-1));
std::vector<double> x(n);
for(int idx = 0; idx < n; idx++){
  x[idx] = xyLeft + (idx)*mesh_size;
}


//Generate the polynomial chaos basis.
std::cout << "Generating the polynomail chaos basis...... ";
std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
if(polyType == 'h'){
  for (int i = 0; i< d; i++){
    bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p));
  }
}

if(polyType == 'r'){
  for (int i = 0; i< d; i++){
    bases[i] = Teuchos::rcp(new Stokhos::RysBasis<int,double>(p,rysCut,true));
  }
}

Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis =
  Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = basis->getTripleProductTensor();

std::cout << "Done.\n Generating the Deterministic Stiffness Matricies...... ";

const Teuchos::RCP<const Epetra_Map> StochMap = Teuchos::rcp(new Epetra_Map(n*n*basis->size(),0,Comm));
const Teuchos::RCP<const Epetra_Map> BaseMap = Teuchos::rcp(new Epetra_Map(n*n,0,Comm));



int n2 = x.size()*x.size();
double meshSize = x[1]-x[0];
 
//Construct the deterministic stiffness matircies.  The boundary conditions are imposed by setting
//1 to the appropriate entries in the mean matrix and zero in the fluctuations.
int NumMyElements = (*BaseMap).NumMyElements();
int * MyGlobalElements = (*BaseMap).MyGlobalElements();
int * NumNz = new int[NumMyElements];


std::vector<Teuchos::RCP<Epetra_CrsMatrix> > A_k(d+1);

  
// Add  rows one-at-a-time
// Need some vectors to help
  
double *Values = new double[4];
int *Indices = new int[4];
double two;
int NumEntries;
int * bcIndices = new int[NumMyElements];  

for(int k = 0; k<=d; k++){
  for(int i = 0; i<NumMyElements; i++){
    /*
    // MyGlobalElements[i]<x.size() ==> Boundary node on bottom edge.
    // MyGlobalElements[i]%x.size() == 0 ==> Boundary node on left edge.
    // MyGlobalElements[i]+1%x.size() == 0 ==> right edge.
    // MyGlobalElements[i] >= n - x.size() ==> top edge.
    */
    if((MyGlobalElements[i] < x.size() || MyGlobalElements[i]%x.size() == 0 || (MyGlobalElements[i]+1)%x.size() == 0 || MyGlobalElements[i] >= n*n - x.size())){
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
    // Put in the diagonal entry
    if (bcIndices[i]==0 || k == 0) A_k[k]->InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
 }
 // Finish up, trasforming the matrix entries into local numbering,
 // to optimize data transfert during matrix-vector products
 A_k[k]->FillComplete();
 
}

const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> > block_ops = Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Operator>(basis));
for(int i = 0; i<=d; i++){
  block_ops->setCoeffPtr(i, A_k[i]);
}


//Construct the implicit operator.
Stokhos::KLMatrixFreeEpetraOp system(BaseMap, StochMap, basis, Cijk, A_k); 

//Generate the RHS vector
std::cout<< " Done.\n Generating the RHS Vector........ ";
double (*RHS_fn)(double,double, std::vector<double>&) = &RHS_function;
Epetra_Vector b(*StochMap,true);
generateRHS(RHS_fn, x, b,basis);

std::cout<< "Done.\n Setting up the preconditioner and solving....... \n";

//////////////////////////////////////////////////
//Solve the Linear System.
/////////////////////////////////////////////////

Epetra_Vector x2(*StochMap);
Epetra_LinearProblem problem(&system, &x2, &b);
AztecOO aztec_solver(problem);
Stokhos::StochGalerkinPrecon myPrecon(*A_k[0],basis->norm_squared(),Comm,*StochMap,*StochMap);
aztec_solver.SetPrecOperator(&myPrecon);
aztec_solver.SetAztecOption(AZ_solver, AZ_cg);
//aztec_solver.SetAztecOption(AZ_precond, AZ_none);
aztec_solver.Iterate(10000, 1e-8);


//Post process and output the results.
std::cout << "\n Done. \n Post processing results\n";
std::vector<double> xi(d);

for(int i = 0; i<d; i++){
  xi[i] = evaluationPoint[i];
}

Epetra_Vector Eofu(*BaseMap,true);
Epetra_Vector varOfu(*BaseMap,true);
Epetra_Vector specificSol(*BaseMap,true);

std::ofstream solution;
solution.open("solution.txt");
x2.Print(solution);
solution.close();

computeMeanSolution(x2, Eofu, basis);
std::ofstream mean;
mean.open("mean.txt");
Eofu.Print(mean);
mean.close();

computeVarianceSolution(x2, Eofu, varOfu, basis);
std::ofstream var;
var.open("var.txt");
varOfu.Print(var);
var.close();

generateRealization(x2, xi, specificSol, basis);
std::ofstream specific;
specific.open("specific.txt");
specificSol.Print(specific);
specific.close();

double error_norm = computeError(x2,&uExact,x, basis);
cout << "Error = " << error_norm;

std::ofstream info;
info.open("info.txt");
info << (*argv[6] == 'r')?rysCut:0;
info << "\t" << error_norm<< "\t"<<aztec_solver.NumIters();
info.close();

std::cout << "\nFinished! Results written out to 'solution.txt', 'mean.txt', 'var.txt', 'specific.txt'\n";



return 1;
}
 










