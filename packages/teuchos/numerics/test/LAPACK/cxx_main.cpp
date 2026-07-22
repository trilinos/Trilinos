// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <vector>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"


template <typename T>
int lapackTest( bool verbose );

template <typename T>
struct specializedLAPACK 
{
  // Specialized test for real-valued vs. complex valued types
  static int test(bool verbose);

  // Specializations for mixed complex-real operations
  template<class R>
  void add( const R& a, const T& b, T& result ){ result = a + b; }

  template<class R>
  void multiply( const R& a, const T& b, T& result ){ result = a * b; }

  template<class R>
  void divide( const T& a, const R& b, T& result ){ result = a / b; }
};

#ifdef HAVE_TEUCHOS_COMPLEX

// Partial specialization for std::complex numbers templated on real type T
template <typename T>
struct specializedLAPACK< std::complex<T> >
{
  // Specialized test for real-valued vs. complex valued types
  static int test(bool verbose);

  // Specializations for mixed complex-real operations
  template<class R>
  void add( const R& a, const std::complex<T>& b, std::complex<T>& result )
  { std::complex<T> tmp( a, 0 ); result = b + tmp; } 

  template<class R>
  void multiply( const R& a, const std::complex<T>& b, std::complex<T>& result )
  { std::complex<T> tmp( a, 0 ); result = tmp * b; }

  template<class R>
  void divide( const std::complex<T>& a, const R& b, std::complex<T>& result )
  { std::complex<T> tmp( b, 0 ); result = a / tmp; }
};

#endif

// Main test
int main(int argc, char* argv[])
{
  int numberFailedTests = 0;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  using std::fabs;

#ifdef HAVE_TEUCHOS_INST_FLOAT
  if (verbose)
    std::cout << std::endl << "LAPACK test for float" << std::endl; 
  numberFailedTests += lapackTest<float>(verbose);
#endif
  if (verbose)
    std::cout << std::endl << "LAPACK test for double" << std::endl; 
  numberFailedTests += lapackTest<double>(verbose);

#ifdef HAVE_TEUCHOS_COMPLEX
#ifdef HAVE_TEUCHOS_INST_COMPLEX_FLOAT
  if (verbose)
    std::cout << std::endl << "LAPACK test for std::complex<float>" << std::endl; 
  numberFailedTests += lapackTest<std::complex<float> >(verbose);
#endif
#ifdef HAVE_TEUCHOS_INST_COMPLEX_DOUBLE
  if (verbose)
    std::cout << std::endl << "LAPACK test for std::complex<double>" << std::endl; 
  numberFailedTests += lapackTest<std::complex<double> >(verbose);
#endif
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

// Common test for all four types: float, double, std::complex<float>, std::complex<double>
// Calls the specialized test for types whose interfaces are different or undefined
template <typename T>
int lapackTest( bool verbose )
{ 
  int numberFailedTests = 0;

  // Define some common characters
  int info=0;
  char char_G = 'G';
  char char_N = 'N';
  char char_U = 'U';

  // Create some common typedefs
  typedef Teuchos::ScalarTraits<T> STS;
  typedef typename STS::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;

  T one = STS::one();
  MagnitudeType m_one = STM::one();
  T zero = STS::zero();

  Teuchos::LAPACK<int,T> L;
  specializedLAPACK<T> sL; 

  const int n_gesv = 4;
  std::vector<T> Ad(n_gesv*n_gesv,zero);
  std::vector<T> bd(n_gesv,zero);
  int IPIV[n_gesv];

  Ad[0] = 1; Ad[2] = 1; Ad[5] = 1; Ad[8] = 2; Ad[9] = 1; Ad[10] = 1; Ad[14] = 2; Ad[15] = 2;
  bd[1] = 2; bd[2] = 1; bd[3] = 2;

  if (verbose) std::cout << "LASCL test ... ";
  L.LASCL(char_G, 1, 1, m_one, m_one, n_gesv, n_gesv, &Ad[0], n_gesv, &info);
  if ( !info ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }

  if (verbose) std::cout << "GESV test ... ";
  L.GESV(n_gesv, 1, &Ad[0], n_gesv, IPIV, &bd[0], n_gesv, &info);
  if ( !info ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED" << std::endl;
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

  // Create a simple diagonal linear system
  std::vector<T> Ad2_sub(n_gesv-1, zero), b2(n_gesv, one);
  std::vector<MagnitudeType> Ad2(n_gesv, m_one);

  if (verbose) std::cout << "PTTRF test ... ";
  L.PTTRF(n_gesv, &Ad2[0], &Ad2_sub[0], &info);
  if ( !info ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }

  int n_potrf = 5;
  std::vector<T> diag_a(n_potrf*n_potrf, zero);
  for (int i=0; i<n_potrf; i++)
  {
    T tmp = zero;
    sL.add( i, one, tmp );
    diag_a[i*n_potrf + i] = tmp*tmp;
  }

  if (verbose) std::cout << "POTRF test ... ";
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
      T tmp = zero;
      sL.add( i, one, tmp );
      if ( diag_a[i*n_potrf + i] == tmp )
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

  if (verbose) std::cout << "POTRI test ... ";
  std::vector<T> diag_a_trtri(diag_a); // Save a copy for TRTRI test

  L.POTRI(char_U, n_potrf, &diag_a[0], n_potrf, &info);

  T tmp = zero;
  sL.multiply( 1.0/4.0, one, tmp );
  if ( info != 0 || (diag_a[n_potrf+1] != tmp) )
  {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }
  else
    if (verbose) std::cout << "passed!" << std::endl;

  if (verbose) std::cout << "TRTRI test ... ";

  int n_trtri = n_potrf;
  L.TRTRI( char_U, char_N, n_trtri, &diag_a_trtri[0], n_trtri, &info );
  for (int i=0; i<n_trtri; i++)
  {
    tmp = zero;
    sL.divide( one, i+1.0, tmp );
    if ( info != 0 )
    {
      numberFailedTests++;
      break;
    }
    else if ( diag_a_trtri[i*n_trtri + i] == tmp )
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

#ifndef TEUCHOSNUMERICS_DISABLE_STEQR_TEST

  if (verbose) std::cout << "STEQR test ... ";

  const int n_steqr = 10;
  std::vector<MagnitudeType> diagonal(n_steqr);
  std::vector<MagnitudeType> subdiagonal(n_steqr-1);

  for (int i=0; i < n_steqr; ++i) {
    diagonal[i] = n_steqr - i;
    if (i < n_steqr-1)
      subdiagonal[i] = STM::eps() * (i+1);
  }

  std::vector<T> scalar_dummy(1,0.0);
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

#endif // TEUCHOSNUMERICS_DISABLE_STEQR_TEST

  numberFailedTests += specializedLAPACK<T>::test( verbose ); 

  return numberFailedTests; 
}

template<class T>
int specializedLAPACK<T>::test(bool verbose)
{
  // Create some common typedefs
  typedef Teuchos::ScalarTraits<T> STS;
  typedef typename STS::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;

  T one = STS::one();
  MagnitudeType m_one = STM::one();
  T zero = STS::zero();

  char char_E = 'E';
  char char_U = 'U';

  int info=0;
  int numberFailedTests = 0;

  Teuchos::LAPACK<int,T> L;

  if (verbose) std::cout << "LAPY2 test ... ";
  T x = 3*one, y = 4*one;
  T lapy = L.LAPY2(x, y);
  if ( lapy == 5*one ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED ( " << lapy << " != 5 )" << std::endl;
    numberFailedTests++;
  }

  if (verbose) std::cout << "LAMCH test ... ";

  T st_eps = L.LAMCH( char_E );
  if (verbose)
    std::cout << "[ eps = " << st_eps << " ] passed!" << std::endl;

  // Create a simple diagonal linear system
  const int n = 4;
  std::vector<T> Ad2_sub(n-1, zero), b2(n, one);
  std::vector<MagnitudeType> Ad2(n, m_one);

  if (verbose) std::cout << "PTTRS test ... ";
  L.PTTRS(n, 1, &Ad2[0], &Ad2_sub[0], &b2[0], n, &info);
  if ( !info ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }

  if (verbose) std::cout << "POCON test ... ";

  std::vector<T> diag_a(n*n);
  for (int i=0; i<n; i++)
  {
    diag_a[i*n + i] = one;
  }
  MagnitudeType rcond, anorm = m_one;
  std::vector<T> work(3*n);
  std::vector<int> iwork(n);

  L.POCON(char_U, n, &diag_a[0], n, anorm, &rcond, &work[0], &iwork[0], &info);
  if (info != 0 || (rcond != m_one))
  {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }
  else
    if (verbose) std::cout << "passed!" << std::endl;


  return numberFailedTests;
}

#ifdef HAVE_TEUCHOS_COMPLEX

template<class T>
int specializedLAPACK<std::complex<T> >::test( bool verbose )
{
  // Create some common typedefs
  typedef Teuchos::ScalarTraits<std::complex<T> > STS;
  typedef typename STS::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;

  std::complex<T> one = STS::one();
  MagnitudeType m_one = STM::one();
  std::complex<T> zero = STS::zero();

  char char_L = 'L';
  char char_U = 'U';

  int info=0;
  int numberFailedTests = 0;

  Teuchos::LAPACK<int,std::complex<T> > L;

  // Create a simple diagonal linear system
  const int n = 4;
  std::vector<std::complex<T> > Ad2_sub(n-1, zero), b2(n, one);
  std::vector<MagnitudeType> Ad2(n, m_one);

  if (verbose) std::cout << "PTTRS test ... ";
  L.PTTRS(char_L, n, 1, &Ad2[0], &Ad2_sub[0], &b2[0], n, &info);
  if ( !info ) {
    if (verbose) std::cout << "passed!" << std::endl;
  } else {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }

  if (verbose) std::cout << "POCON test ... ";

  std::vector<std::complex<T> > diag_a(n*n);
  for (int i=0; i<n; i++)
  {
    diag_a[i*n + i] = one;
  }
  MagnitudeType rcond, anorm = m_one;
  std::vector<std::complex<T> > work(2*n);
  std::vector<MagnitudeType> rwork(n);
  std::vector<int> iwork(n);

  L.POCON(char_U, n, &diag_a[0], n, anorm, &rcond, &work[0], &rwork[0], &info);
  if (info != 0 || (rcond != m_one))
  {
    if (verbose) std::cout << "FAILED" << std::endl;
    numberFailedTests++;
  }
  else
    if (verbose) std::cout << "passed!" << std::endl;
  
return numberFailedTests;
}

#endif
