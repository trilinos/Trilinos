// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <iostream>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_Version.hpp"

int main()
{
  cout << Teuchos::Teuchos_Version() << endl << endl;

#if 0
  double* Ad = new double[16];
  double* xd = new double[4];
  double* bd = new double[4];
  float* Af = new float[16];
  float* xf = new float[4];
  float* bf = new float[4];
  int* Ai = new int[16];
  int* xi = new int[4];
  int* bi = new int[4];

  int* IPIV = new int[4];
  int info;

  int i;
  for(i = 0; i < 16; i++)
    {
      Ad[i] = 0;
      Af[i] = 0;
      Ai[i] = 0;
    }
  for(i = 0; i < 4; i++)
    {
      xd[i] = 0;
      bd[i] = 0;
      xf[i] = 0;
      bf[i] = 0;
      xi[i] = 0;
      bi[i] = 0;
    }

  Ad[0] = 1; Ad[2] = 1; Ad[5] = 1; Ad[8] = 2; Ad[9] = 1; Ad[10] = 1; Ad[14] = 2; Ad[15] = 2;
  xd[0] = -2; xd[1] = 1; xd[2] = 1; xd[3] = 1;
  bd[1] = 2; bd[2] = 1; bd[3] = 2;
  Af[0] = 1; Af[2] = 1; Af[5] = 1; Af[8] = 2; Af[9] = 1; Af[10] = 1; Af[14] = 2; Af[15] = 2;
  xf[0] = -2; xf[1] = 1; xf[2] = 1; xf[3] = 1;
  bf[1] = 2; bf[2] = 1; bf[3] = 2;
  Ai[0] = 1; Ai[2] = 1; Ai[5] = 1; Ai[8] = 2; Ai[9] = 1; Ai[10] = 1; Ai[14] = 2; Ai[15] = 2;
  xi[0] = -2; xi[1] = 1; xi[2] = 1; xi[3] = 1;
  bi[1] = 2; bi[2] = 1; bi[3] = 2;

  Teuchos::LAPACK<int,double> L;
  Teuchos::LAPACK<int,float> M;
  Teuchos::LAPACK<int,int> N;
  L.GESV(4, 1, Ad, 4, IPIV, bd, 4, &info);
  for(i = 0; i < 4; i++)
    {
      cout << bd[i] << " ";
    }
  cout << endl;
  M.GESV(4, 1, Af, 4, IPIV, bf, 4, &info);
  for(i = 0; i < 4; i++)
    {
      cout << bf[i] << " ";
    }
  cout << endl;
  N.GESV(4, 1, Ai, 4, IPIV, bi, 4, &info);
  for(i = 0; i < 4; i++)
    {
      cout << bi[i] << " ";
    }
  cout << endl;

  cout << "GEES/LAPY2" << endl;

    //   void SGEES_F77(Teuchos_fcd, Teuchos_fcd, int*, int, float*, int, int, float*, float*, float*, float*, int, float*, int, int, int*, int*);


  int* foo = new int[2];
  float* bar = new float[2];
  double* bar2 = new double[2];
  int I = 99;
  int INFO = 99;
  M.GEES('A', 'B', foo, I, bar, I, I, bar, bar, bar, bar, I, bar, I, I, I, INFO);
  cout << INFO << endl;
  L.GEES('A', 'B', foo, I, bar2, I, I, bar2, bar2, bar2, bar2, I, bar2, I, I, I, INFO);
  cout << INFO << endl;

  float xx = 9, yy = 22;
  cout << M.LAPY2(xx, yy) << endl;

#endif

  return 0;

}
