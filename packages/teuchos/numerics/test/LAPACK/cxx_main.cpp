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
#include <vector>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  int numberFailedTests = 0;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  using std::fabs;

  // Define some common characters
  int info=0;
  char char_N = 'N';
  char char_U = 'U';

  // Create some common typedefs
  typedef Teuchos::ScalarTraits<double> STS;
  typedef Teuchos::ScalarTraits<double>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;

  Teuchos::LAPACK<int,double> L;
  Teuchos::LAPACK<int,float> M;

  const int n_gesv = 4;
  std::vector<double> Ad(n_gesv*n_gesv,0.0);
  std::vector<double> bd(n_gesv,0.0);
  std::vector<float> Af(n_gesv*n_gesv,0.0);
  std::vector<float> bf(n_gesv,0.0);
  int IPIV[n_gesv];

  Ad[0] = 1; Ad[2] = 1; Ad[5] = 1; Ad[8] = 2; Ad[9] = 1; Ad[10] = 1; Ad[14] = 2; Ad[15] = 2;
  bd[1] = 2; bd[2] = 1; bd[3] = 2;
  Af[0] = 1; Af[2] = 1; Af[5] = 1; Af[8] = 2; Af[9] = 1; Af[10] = 1; Af[14] = 2; Af[15] = 2;
  bf[1] = 2; bf[2] = 1; bf[3] = 2;

  if (verbose) std::cout << "GESV test ... ";
  L.GESV(n_gesv, 1, &Ad[0], n_gesv, IPIV, &bd[0], n_gesv, &info);
  M.GESV(n_gesv, 1, &Af[0], n_gesv, IPIV, &bf[0], n_gesv, &info);
  for(int i = 0; i < 4; i++)
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

  if (verbose) std::cout << "LAMCH test ... ";

  char char_E = 'E';
  double d_eps = L.LAMCH( char_E );
  float f_eps = M.LAMCH( char_E );
  if (verbose)
    std::cout << "[ Double-precision eps = " << d_eps << ", single-precision eps = " << f_eps << " ] passed!" << std::endl;

  if (verbose) std::cout << "POTRF test ... ";

  int n_potrf = 5;
  std::vector<double> diag_a(n_potrf*n_potrf, 0.0);
  for (int i=0; i<n_potrf; i++)
    diag_a[i*n_potrf + i] = (i+1)*(i+1);
  L.POTRF(char_U, n_potrf, &diag_a[0], n_potrf, &info);

  if (info != 0)
  {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }
  else
  { 
    for (int i=0; i<n_potrf; i++)
    {
      if ( diag_a[i*n_potrf + i] == (i+1) ) 
      {
        if (verbose && i==(n_potrf-1)) std::cout << "passed!" << std::endl;
      } 
      else
      {
        if (verbose) std::cout << "FAILED" << std::endl;
        numberFailedTests++;
       break;
      }
    }
  }

  if (verbose) std::cout << "POCON test ... ";
  
  double anorm = (n_potrf*n_potrf), rcond;
  std::vector<double> work(3*n_potrf);
  std::vector<int> iwork(n_potrf);
  
  L.POCON(char_U, n_potrf, &diag_a[0], n_potrf, anorm, &rcond, &work[0], &iwork[0], &info);
  if (info != 0 || (rcond != 1.0/anorm))
  {
    numberFailedTests++;
    if (verbose) std::cout << "FAILED" << std::endl;
  }
  else
  { 
    if (verbose) std::cout << "passed!" << std::endl;
  } 

  if (verbose) std::cout << "POTRI test ... ";
  std::vector<double> diag_a_trtri(diag_a); // Save a copy for TRTRI test 
   
  L.POTRI(char_U, n_potrf, &diag_a[0], n_potrf, &info);

  if (info != 0 || (diag_a[n_potrf+1] != 1.0/4.0))
  {
    numberFailedTests++;
    if (verbose) std::cout << "FAILED" << std::endl;
  }
  else
  { 
    if (verbose) std::cout << "passed!" << std::endl;
  } 
 
  if (verbose) std::cout << "TRTRI test ... ";
  
  int n_trtri = n_potrf;
  L.TRTRI( char_U, char_N, n_trtri, &diag_a_trtri[0], n_trtri, &info );
  for (int i=0; i<n_trtri; i++)
  {
    if ( diag_a_trtri[i*n_trtri + i] == 1.0/(i+1) ) 
    {
      if (verbose && i==(n_trtri-1)) std::cout << "passed!" << std::endl;
    } 
    else
    {
      if (verbose) std::cout << "FAILED" << std::endl;
      numberFailedTests++;
      break;
    }
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
 
  if (verbose) std::cout << "STEQR test ... ";

#ifndef TEUCHOSNUMERICS_DISABLE_STEQR_TEST

  const int n_steqr = 10;
  std::vector<MagnitudeType> diagonal(n_steqr);
  std::vector<MagnitudeType> subdiagonal(n_steqr-1);

  for (int i=0; i < n_steqr; ++i) {
    diagonal[i] = n_steqr - i;
    if (i < n_steqr-1)
      subdiagonal[i] = STM::eps() * (i+1);
  }

  std::vector<double> scalar_dummy(1,0.0);
  std::vector<MagnitudeType> mag_dummy(4*n_steqr,0.0);

  L.STEQR (char_N, n_steqr, &diagonal[0], &subdiagonal[0],
           &scalar_dummy[0], n_steqr, &mag_dummy[0], &info);

  if (info != 0)
  {
    if (verbose)  std::cout << "STEQR: compute symmetric tridiagonal eigenvalues: "
		  << "LAPACK's _STEQR failed with info = "
		  << info;

      numberFailedTests++;
  }

  MagnitudeType lambda_min = diagonal[0];
  MagnitudeType lambda_max = diagonal[n_steqr-1];
  MagnitudeType exp_lambda_min = STM::one();
  MagnitudeType exp_lambda_max = STM::one()*n_steqr;

  if ((fabs(lambda_min-exp_lambda_min)<1e-12) && (fabs(lambda_max-exp_lambda_max)<1e-12))
  {
    if (verbose) std::cout << "passed!" << std::endl;
  } 
  else 
  {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }

#else // TEUCHOSNUMERICS_DISABLE_STEQR_TEST

  if (verbose) std::cout << "SKIPPED!\n";

#endif // TEUCHOSNUMERICS_DISABLE_STEQR_TEST

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
