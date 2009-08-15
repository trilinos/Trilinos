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

// -div(a(x,xi) grad u) + cu = f.
//  a is expanded in a KL expansion as a = mean + sigma\sum^M xi_i f_i(x).
//  evalEigenfunction evaluates the ith deterministic dunction in the KL
//  expansion of the diffusion coefficient.  RHS_function evaluates f.

//////////////////////////////////////////////////////
//GLOBALS
//////////////////////////////////////////////////////
double sigma, mean, weightCut;

//For the exponential random field.
std::vector<double> lambda, alpha, omega;
std::vector<int> xind, yind;


//The probability distribution of the random variables.
const double weight(const double x){
  return 1;
}


//Diffusion term is expanded in a KL expansion
//a(x,xi) = mu + sum^d xi_i f_i(x).
//This function evaluates the ith spatial component of the
//expansion.
double evalEigenfunction(double x,double y, int idx){
  
   if(idx == 0){
      return mean;
   }else{
      double omegax = omega[xind[idx-1]];
      double omegay = omega[yind[idx-1]];
      double result = 1;
      if(xind[idx-1]%2 == 1){
         result = result * sin(omegax*x);
      }else{
         result = result * cos(omegax*x);
      }

      if(yind[idx-1]%2 == 1){
         result = result * sin(omegay*y);
      }else{
         result = result * cos(omegay*y);
      }
      return sigma*result*sqrt(lambda[idx-1])*alpha[idx-1];
      
   }
 
}

//The loading term is expanded as a PC expansion,
// this function returns the value of the idx spatial
// function evaluated at x,y.
double RHS_function_PC(double x, double y, int idx){

double result;
if(idx==0)
  result = 1;
else result = 0;
return result;  

}

//Function for the righ hand side.
double RHS_function(double x, double y, std::vector<double>& xi){
   
   double result = RHS_function_PC(x,y,0);
   for(int idx = 1; idx<=xi.size(); idx++){
      result = result + RHS_function_PC(x,y,idx)*xi[idx-1];
   }
   return result;
}
