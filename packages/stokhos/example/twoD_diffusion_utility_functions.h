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
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

//If the RHS function is expanded in a KL-expansion, use this to generate the RHS vector.
void generateRHS(double (*rhs_function)(double, double, int), Teuchos::Array<double> mesh, Epetra_Vector& RHS,Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  int N_xi = basis->size();
  int N_x = mesh.size();
  int d = basis->dimension();
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = basis->getLowOrderTripleProductTensor(d);
  RHS.PutScalar(0);

  double doubleProd;
  int i,k;

  Teuchos::Array<double> one(d);
  for(int j = 0; j<d; j++)one[j] = 1;
  Teuchos::Array< double > values(N_xi);
  basis->evaluateBases(one, values);
  Teuchos::Array< double > norms = basis->norm_squared();
  for(int stochIdx = 0; stochIdx <= d; stochIdx++){
    double gamma = values[stochIdx]/(1-basis->evaluateZero(stochIdx));
    Cijk->value(0, stochIdx, i, k, doubleProd);
    if( stochIdx!=0 ){
      doubleProd /= gamma;
      if( i == k) doubleProd = doubleProd -basis->evaluateZero(stochIdx)*norms[i];
    }else{
      if( i == k ) doubleProd = norms[i];
    }
    
    for(int xmeshIdx = 1; xmeshIdx < N_x-1; xmeshIdx++){
      for(int ymeshIdx = 1; ymeshIdx < N_x-1; ymeshIdx++){
        RHS[stochIdx*N_x*N_x + xmeshIdx + N_x*ymeshIdx] = 
          RHS[stochIdx*N_x*N_x + xmeshIdx + N_x*ymeshIdx]
          +
          rhs_function(mesh[xmeshIdx],mesh[ymeshIdx],stochIdx)*doubleProd;
      }
    }
  }
}

// for an arbitary RHS function f(x,\xi) computes the RHS vector via quadrature. VERY SLOW if there are many basis functions
// and the dimension is high.
void generateRHS(double (*rhs_function)(double, double, Teuchos::Array<double>&), Teuchos::Array<double> mesh, Epetra_Vector& RHS,Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis){

  int N_xi = basis->size();
  int N_x = mesh.size();
  double quadOrder;
  quadOrder = 5*basis->order();
  std::cout << basis->order()<< "\n";
  Stokhos::SparseGridQuadrature<int,double> quadRule(basis,basis->order());
  

  Teuchos::Array< Teuchos::Array<double> > quadPts = quadRule.getQuadPoints();
  Teuchos::Array<double> quadWeights = quadRule.getQuadWeights();
  Teuchos::Array< double > basisVals(N_xi);
  
  
  for(std::size_t quadIdx = 0; quadIdx < quadPts.size(); quadIdx++){
    std::cout << quadIdx << "/" << quadPts.size()-1 << "\n";
    basis->evaluateBases(quadPts[quadIdx], basisVals);
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

//Given the computed sfem solution, computes the mean solution.
void computeMeanSolution(const Epetra_Vector& u, Epetra_Vector& Eofu, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  
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
  std::cout << "val = " << val << "\n";
  for(int i = 0; i< N_x; i++){
    Eofu[i] = val*u[i];
  }

}

//Given the computed SFEM solution and the mean solution, computes the variance via <u^2> - <u>^2
void computeVarianceSolution(const Epetra_Vector& u, const Epetra_Vector& Eofu, Epetra_Vector& varOfu, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

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
  Teuchos::Array< double > norms = basis->norm_squared();
  for(int c = 0; c<N_xi; c++){
    varOfu.Update(norms[c],*uBlock(c),1.0);
  }
  
}

//Given a value of \xi returns the SFEM approximation of the solution at \xi.
void generateRealization(const Epetra_Vector& u, Teuchos::Array<double>& xi, Epetra_Vector& u_xi, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  int N_xi = basis->size();
  int N_x = u.MyLength()/N_xi;
  Teuchos::Array<double> basisValues(basis->size());
  basis->evaluateBases(xi, basisValues);

  u_xi.PutScalar(0);

  for(int i = 0; i<N_xi; i++){
    for( int j = 0; j< N_x; j++){
      u_xi[j] = u_xi[j] + u[j+N_x*i]*basisValues[i];
    }
  }
}


//Generates the KL expansion of the field with covariance defined by C(x,y) = exp(-|x_1-x_2|*corrx - |y_1-y_2|*corrx)
void generateExponentialRF(const int orderKL, const double corrx, Teuchos::Array<double>& lambda,
                           Teuchos::Array<double>& alpha_out,Teuchos::Array<double>& omega,
                           Teuchos::Array<int>& xind, Teuchos::Array<int>& yind){
   
   
   omega.resize(orderKL);
   lambda.resize(orderKL);
   Teuchos::Array<double> alpha(orderKL);
   double lower, upper, midpoint;
   for(int i = 0; i<=ceil(((double)orderKL)/2); i++){
      
      if(i>0 && 2*i <= orderKL){
         if( (2*i-1)*PI+.000000001 > 0)lower = (2*i-1)*PI+.000000001;
         else lower = 0;
         upper = (2*i+1)*PI-.000000001;
         while(fabs(lower - upper) > 2*1e-8){
           midpoint = (lower + upper)/2;
           if((lower + corrx*tan(lower*.5)>0 && midpoint + corrx*tan(midpoint*.5)>0) || (lower + corrx*tan(lower*.5)<0 && midpoint + corrx*tan(midpoint*.5)<0)){
              lower = midpoint;
           }else{
              upper = midpoint;
           }
         }
         omega[2*i-1] = (lower + upper)/2;
         lambda[2*i-1] = 2*corrx/(omega[2*i - 1]*omega[2*i - 1] + corrx*corrx);
         alpha[2*i-1] = 1/sqrt(.5 - sin(omega[2*i - 1])/(2*omega[2*i - 1]));
      } 
      
      if(2*i + 1<= orderKL){
         if( (2*i-1)*PI+.000000001 > 0)lower = (2*i-1)*PI+.000000001;
         else lower = 0;
         upper = (2*i+1)*PI-.000000001;
         while(fabs(lower - upper) > 2*1e-8){
            midpoint = (lower + upper)/2;
            if((corrx - lower*tan(lower*.5)>0 && corrx - midpoint*tan(midpoint*.5)>0) || (corrx - lower*tan(lower*.5)<0 && corrx - midpoint*tan(midpoint*.5)<0)){
               lower = midpoint;
            }else{
               upper = midpoint;
            }
         }
         omega[2*i] = (lower + upper)/2;
         lambda[2*i] = 2*corrx/(omega[2*i]*omega[2*i] + corrx*corrx);
         alpha[2*i] = 1/sqrt(.5 + sin(omega[2*i])/(2*omega[2*i]));
      }   
   }

   Teuchos::Array<Teuchos::Array<double> > productEigs(orderKL);
   for(int i = 0; i<orderKL; i++){
     productEigs[i] = Teuchos::Array<double>(orderKL);
     for(int j = 0; j<orderKL; j++){
        productEigs[i][j] = lambda[i]*lambda[j];
     }
   }

   int maxx, maxy;
   double max;
   for(int l = 0; l<orderKL; l++){
      max = 0;
      maxx = 0;
      maxy = 0;
      for(int i = 0; i<orderKL; i++){
         
         for(int j = 0; j<orderKL; j++){
            if(productEigs[i][j] > max){
               max = productEigs[i][j];
               maxx = i;
               maxy = j;
            }
         }
      }

      alpha_out[l] = alpha[maxx]*alpha[maxy];
      lambda[l] = productEigs[maxx][maxy];
      xind[l] = maxx;
      yind[l] = maxy;
      productEigs[maxx][maxy] = 0;
   }  
   
}
