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

// -div(a(x,xi) grad u) = f.
//  a is expanded in a KL expansion as a = mean + sigma\sum^M xi_i f_i(x).
//  evalEigenfunction evaluates the ith deterministic dunction in the KL
//  expansion of the diffusion coefficient.  RHS_function evaluates f.

//////////////////////////////////////////////////////
//GLOBALS
//////////////////////////////////////////////////////
double sigma, mean, rysCut;


//Diffusion term is expanded in a KL expansion
//a(x,xi) = mu + sum^d xi_i f_i(x).
//This function evaluates the ith spatial component of the
//expansion.
double evalEigenfunction(double x,double y, int idx){

  if(idx == 0){
    return mean;
  }else{
    return sigma*(1/pow(idx*PI,2))*cos(6*PI*idx*(x*x + y*y));   
  }

  //Anisotropic diffusion
  /*
  if(idx == 0){
    return mean;
  }else{
    return sigma*(1/pow(idx*PI,2))*cos(6*PI*idx*(3*x*x + y*y));   
  }
  */
}


//Function for the righ hand side.
double RHS_function(double x, double y, std::vector<double>& xi){

int d = xi.size();
double r2 = x*x + y*y;
double z = 0;
double negExpMagxiSquared = 0;
for(int i = 0; i<d; i++) negExpMagxiSquared = negExpMagxiSquared + xi[i]*xi[i];
negExpMagxiSquared = exp(-negExpMagxiSquared);

for(int j = 0; j<d+1; j++){
  if(j == 0){
    z = z - mean*(32*(y*y - .25) + 32*(x*x - .25));
      //z = -32;
  }else{
    double pt1 = 32*(y*y - .25);
    double pt2 = 32*(x*x - .25);
    double pt3 = 2*6*32*y*y*(x*x - .25);
    double pt4 = 2*6*32*x*x*(y*y - .25);
    double pt5 = (1/(j*PI))*sin(6*PI*j*r2);
    double pt6 = (1/(j*j*PI*PI))*cos(6*PI*j*r2);
    z = z - xi[j-1]*sigma*((pt1+ pt2)*pt6 - (pt3+pt4)*pt5);
  }
}

return negExpMagxiSquared*z;

//////////////////////////////////////////////////////////
//return 1;
}


//Exact solution if known.
double uExact(double x,double y, std::vector<double>& xi){
  //double ans = 16*(x*x -.25)*(y*y -.25);
  //return ans;
  
  // u = exp(-|xi|^2) 16*(x*x -.25)*(y*y -.25)
  int d = xi.size();
  double negExpMagxiSquared = 0;
  for(int i = 0; i<d; i++) negExpMagxiSquared = negExpMagxiSquared + xi[i]*xi[i];
  negExpMagxiSquared = exp(-negExpMagxiSquared);
  double ans = 16*(x*x -.25)*(y*y -.25)*negExpMagxiSquared;
  return ans;
  
  /*
  // u = sum(xi_i^3) 16*(x*x -.25)*(y*y -.25)
  int d = xi.size();
  double MagxiSquared = 0;
  for(int i = 0; i<d; i++) MagxiSquared = MagxiSquared + xi[i]*xi[i]*xi[i];
  double ans = 16*(x*x -.25)*(y*y -.25)*MagxiSquared;
  return ans;
  */
}

