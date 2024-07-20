// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// recurrence_example
//
//  usage: 
//     recurrence_example
//
//  output:  
//     Prints the recurrence coefficients for the first 5 normalized polynomials
//     orthogonal wrt the given weight.  Follows up by printing the computed 
//     norms and outputting a 11 point gaussian quadrature rule.  
//     Demonstrate orthogonality by outputting the maximum computed 
//     |<psi_i, psi_j>| for j != i.

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"

double weightFunction(const double& x){
  if(2*fabs(x+2) < 1){
    return exp(-1/(1-pow(2*(x+2),2)));
  }else if(2*fabs(x-2)<1){
    return exp(-1/(1-pow(2*(x-2),2)));
  }else{
    return 0;
  }
}

int main(int argc, char **argv)
{
  try {
    const int p = 5;
    const double leftEndPt = -3;
    const double rightEndPt = 3;

    //Generate a basis up to order p with the support of the weight = [-5,5] using normalization.
    Teuchos::RCP<const Stokhos::DiscretizedStieltjesBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<int,double>("Beta",p,&weightFunction,leftEndPt,rightEndPt,true));
    Teuchos::Array<double> alpha;
    Teuchos::Array<double> beta;
    Teuchos::Array<double> delta;
    Teuchos::Array<double> gamma;
    Teuchos::Array<double> norm_sq;
    norm_sq = basis->norm_squared();
    basis->getRecurrenceCoefficients(alpha, beta, delta, gamma);
    
    //Fetch alpha, beta, gamma and the computed norms and print.
    for(int i = 0; i<p; i++){
      std::cout << "alpha[" << i <<"]= " << alpha[i] << "     beta[" << i <<"]= " << beta[i] << "      gamma[" << i <<"]= " << gamma[i] << "\n";
    }
   
    for(int i = 0; i<=p; i++){
      std::cout << "E(P_"<<i<<"^2) = "<< norm_sq[i] <<"\n";
    }
    
    Teuchos::Array<double> quad_points;
    Teuchos::Array<double> quad_weights;
    Teuchos::Array<Teuchos::Array< double > > quad_values;
    basis->getQuadPoints(20, quad_points, quad_weights, quad_values);
    for(int i = 0; i<quad_points.size(); i++){
      std::cout << "x_i = "<<quad_points[i]<<" w_i = "<< quad_weights[i] <<" " << i << " / " << quad_points.size()<< "\n";
    }

    double maxOffDiag = 0;
    double currentInnerProduct;
    for(int i = 0; i<=p; i++){
      for(int j = 0; j<i; j++){
        currentInnerProduct = fabs(basis->eval_inner_product(i, j));
        if(currentInnerProduct > maxOffDiag){
          maxOffDiag = currentInnerProduct;
        }
      }
    }
    std::cout<<"Maximum Off Diagonal Inner Product is: " << maxOffDiag << "\n"; 
  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

    
}
