// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include <iostream>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  int numberFailedTests = 0;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  Teuchos::LAPACK<int,double> L;
  Teuchos::LAPACK<int,float> M;

  double Ad[16];
  //double xd[4];
  double bd[4];
  float Af[16];
  //float xf[4];
  float bf[4];

  int IPIV[4];
  int info;

  int i;
  for(i = 0; i < 16; i++)
    {
      Ad[i] = 0;
      Af[i] = 0;
    }
  for(i = 0; i < 4; i++)
    {
      //xd[i] = 0;
      bd[i] = 0;
      //xf[i] = 0;
      bf[i] = 0;
    }

  Ad[0] = 1; Ad[2] = 1; Ad[5] = 1; Ad[8] = 2; Ad[9] = 1; Ad[10] = 1; Ad[14] = 2; Ad[15] = 2;
  //xd[0] = -2; xd[1] = 1; xd[2] = 1; xd[3] = 1;
  bd[1] = 2; bd[2] = 1; bd[3] = 2;
  Af[0] = 1; Af[2] = 1; Af[5] = 1; Af[8] = 2; Af[9] = 1; Af[10] = 1; Af[14] = 2; Af[15] = 2;
  //xf[0] = -2; xf[1] = 1; xf[2] = 1; xf[3] = 1;
  bf[1] = 2; bf[2] = 1; bf[3] = 2;

  if (verbose) std::cout << "GESV test ... ";
  L.GESV(4, 1, Ad, 4, IPIV, bd, 4, &info);
  M.GESV(4, 1, Af, 4, IPIV, bf, 4, &info);
  for(i = 0; i < 4; i++)
    {
      if (bd[i] == bf[i]) {
        if (verbose && i==3) std::cout << "passed!" << std::endl;
      } else {
        if (verbose) std::cout << "FAILED" << std::endl;
        numberFailedTests++;	
	break;
      }
    }

  if (verbose) std::cout << "LAPY2 test ... ";
  float fx = 3, fy = 4;
  float flapy = M.LAPY2(fx, fy);
  double dx = 3, dy = 4;
  double dlapy = L.LAPY2(dx, dy);
  if ( dlapy == flapy && dlapy == 5.0 && flapy == 5.0f ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED (" << dlapy << " != " << flapy << ")" << std::endl;
    numberFailedTests++;
  }  

#if ! (defined(__INTEL_COMPILER) && defined(_WIN32) )
    
  // Check ILAENV with similarity transformation routine:  dsytrd
  // NOTE:  Do not need to put floating point specifier [s,d,c,z] before routine name, 
  //        this is handled through templating.
  if (verbose) std::cout << "ILAENV test ... ";
  int n1 = 100;
  int size = L.ILAENV(1, "sytrd", "u", n1);
  if (size > 0) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED!" << std::endl;
    numberFailedTests++;
  }

#endif

  if(numberFailedTests > 0)
    {
      if (verbose) {
        std::cout << "Number of failed tests: " << numberFailedTests << std::endl;
        std::cout << "End Result: TEST FAILED" << std::endl;
        return -1;
      }
    }
  if(numberFailedTests==0)
    std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;

}
