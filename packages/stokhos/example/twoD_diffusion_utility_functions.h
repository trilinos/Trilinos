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

void generateRHS(double (*rhs_function)(double, double, std::vector<double>&), std::vector<double> mesh, Epetra_Vector& RHS,Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  int N_xi = basis->size();
  int N_x = mesh.size();
  double quadOrder;
  quadOrder = 5*basis->order();
  
  Stokhos::TensorProductQuadrature<int,double> quadRule(basis,quadOrder);
  //Stokhos::TensorProductQuadrature<int,double> quadRule(basis);

  std::vector< std::vector<double> > quadPts = quadRule.getQuadPoints();
  std::vector<double> quadWeights = quadRule.getQuadWeights();
  std::vector< double > basisVals(N_xi);
  
  
  for(int quadIdx = 0; quadIdx < quadPts.size(); quadIdx++){
    basisVals = basis->evaluateBases(quadPts[quadIdx]);
    for(int stochIdx = 0; stochIdx < N_xi; stochIdx++){
      for(int ymeshIdx = 1; ymeshIdx < N_x-1; ymeshIdx++){
        for(int xmeshIdx = 1; xmeshIdx < N_x-1; xmeshIdx++){
          RHS[stochIdx*N_x*N_x + xmeshIdx + N_x*ymeshIdx] = 
           RHS[stochIdx*N_x*N_x + xmeshIdx + N_x*ymeshIdx] + basisVals[stochIdx]*rhs_function(mesh[xmeshIdx],mesh[ymeshIdx], quadPts[quadIdx])*quadWeights[quadIdx];
        }
      }
    }
  }
}

void computeMeanSolution(const Epetra_Vector& u, Epetra_Vector& Eofu, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  /*
  #ifdef HAVE_MPI
    MPI_INIT(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
  #else
    Epetra_SerialComm Comm;
  #endif
  */
  const Epetra_Comm& Comm = u.Map().Comm();

  //Roll up solution vector into a multivector containing deterministic components.
  int N_x = u.MyLength()/basis->size();
  int N_xi = basis->size();

  Eofu.PutScalar(0.0);
  
  Epetra_Map Map(N_x, 0, Comm);
  // Form x and y into block vectors.
  Epetra_MultiVector uBlock(Map,N_xi);
  
  //Get the triple product tensor.
  Teuchos::RCP< const Stokhos::Sparse3Tensor<int, double> > Cijk; 
  Cijk = basis->getTripleProductTensor();
  double val;
  int i, j;
  int n = Cijk->num_values(0);
  
  for( int l = 0; l<n; l++){
    Cijk->value(0,l,i,j,val);
    if(i==0 && j == 0) break;
  }
  for(int i = 0; i< N_x; i++){
    Eofu[i] = val*u[i];
  }

}

void computeVarianceSolution(const Epetra_Vector& u, const Epetra_Vector& Eofu, Epetra_Vector& varOfu, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  /*
  #ifdef HAVE_MPI
    MPI_INIT(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
  #else
    Epetra_SerialComm Comm;
  #endif
  */
  const Epetra_Comm& Comm = u.Map().Comm();
   
  int N_x = u.MyLength()/basis->size();
  int N_xi = basis->size();


  Epetra_Map Map(N_x, 0, Comm);
  varOfu.Multiply(-1.0, Eofu, Eofu, 0.0);
  Epetra_MultiVector uBlock(Map,N_xi);
  
  int MyLength = uBlock.MyLength();
  for( int c=0; c<N_xi ; c++){
    for( int i=0; i<MyLength; i++){
      uBlock[c][i] = (u)[c*N_x + i];
    }
  }
  
  uBlock.Multiply(1.0, uBlock, uBlock, 0.0);
  std::vector< double > norms = basis->norm_squared();
  for(int c = 0; c<N_xi; c++){
    varOfu.Update(norms[c],*uBlock(c),1.0);
  }
  
}

void generateRealization(const Epetra_Vector& u, std::vector<double>& xi, Epetra_Vector& u_xi, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  int N_xi = basis->size();
  int N_x = u.MyLength()/N_xi;
  std::vector<double> basisValues(basis->size());
  basisValues = basis->evaluateBases(xi);

  u_xi.PutScalar(0);

  //std::cout<< "N_xi = "<<N_xi<<"\n";
  for(int i = 0; i<N_xi; i++){
    for( int j = 0; j< N_x; j++){
      u_xi[j] = u_xi[j] + u[j+N_x*i]*basisValues[i];
      //std::cout << "i = "<<i<<" u_xi["<<j<<"] = " << u_xi[j] <<" " << u[j+N_x*i]*basisValues[i]<< "\n";
    }
  }
}

double computeError(const Epetra_Vector& u, double (*exact_solution)(double, double, std::vector<double>&), std::vector<double> mesh, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  /*
  #ifdef HAVE_MPI
    MPI_INIT(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
  #else
    Epetra_SerialComm Comm;
  #endif
  */
  const Epetra_Comm& Comm = u.Map().Comm();
  
  int N_xi = basis->size();
  int N_x = mesh.size();
  double quadOrder;
  quadOrder = 5*basis->order();
  Stokhos::TensorProductQuadrature<int,double> quadRule(basis, quadOrder);
  //Stokhos::TensorProductQuadrature<int,double> quadRule(basis);

  std::vector< std::vector<double> > quadPts = quadRule.getQuadPoints();
  std::vector<double> quadWeights = quadRule.getQuadWeights();
  std::vector< double > basisVals(N_xi);
  Epetra_Map EMap(N_x*N_x,0,Comm);
  Epetra_Vector localSolution(EMap,true);

  double globalError = 0;
  for(int quadIdx = 0; quadIdx < quadPts.size(); quadIdx++){
    double linfSquared = 0;
    generateRealization(u, quadPts[quadIdx], localSolution, basis);
    for(int ymeshIdx = 1; ymeshIdx < N_x-1; ymeshIdx++){
      for(int xmeshIdx = 1; xmeshIdx < N_x-1; xmeshIdx++){
        double localError = fabs(localSolution[xmeshIdx + N_x*ymeshIdx] - exact_solution(mesh[xmeshIdx],mesh[ymeshIdx],quadPts[quadIdx]));
        linfSquared = (localError>linfSquared)?localError:linfSquared;
      }
    }
    linfSquared*=linfSquared;
    globalError = globalError + quadWeights[quadIdx]*linfSquared;
  }
  return globalError;
}
