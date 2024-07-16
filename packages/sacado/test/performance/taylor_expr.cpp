// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_Random.hpp"
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// ADOL-C includes
#ifdef HAVE_ADOLC
#ifdef PACKAGE
#undef PACKAGE
#endif
#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif
#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif
#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif
#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif
#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif
#ifdef VERSION
#undef VERSION
#endif
#include "adolc/adouble.h"
#include "adolc/interfaces.h"
#include "adolc/taping.h"
#endif

using std::cout;
using std::endl;

// A simple performance test that computes a Taylor expansion of a simple
// expression using many variants of Taylor polynomial classes.

template <typename T>
inline void
func(const T& x1, const T& x2, T& y) {
  y = x1*x2 + sin(x1)/x2;
}

template <typename TaylorType>
double
do_time(int degree, int nloop)
{
  TaylorType x1, x2, y;
  Sacado::Random<double> urand(0.0, 1.0);

  x1 = TaylorType(degree, urand.number());
  x2 = TaylorType(degree, urand.number());
  y = 0.0;
  for (int j=0; j<=degree; j++) {
    x1.fastAccessCoeff(j) = urand.number();
    x2.fastAccessCoeff(j) = urand.number();
  }
  
  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (int j=0; j<nloop; j++) {
    func(x1, x2, y);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

#ifdef HAVE_ADOLC
double
do_time_adolc(int degree, int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  double **X, **Y;

  X = new double*[2];
  X[0] = new double[degree+1];
  X[1] = new double[degree+1];
  Y = new double*[1];
  Y[0] = new double[degree+1]; 
  for (int j=0; j<=degree; j++) {
    X[0][j] = urand.number();
    X[1][j] = urand.number();
    Y[0][j] = 0.0;
  }

  trace_on(0);
  adouble x1, x2, y;
  x1 <<= X[0][0];
  x2 <<= X[1][0];
  func(x1, x2, y);
  y >>= Y[0][0];
  trace_off();
  
  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (int j=0; j<nloop; j++) {
    forward(0,1,2,degree,0,X,Y);
  }
  timer.stop();

  delete [] X[1];
  delete [] X[0];
  delete [] X;

  delete [] Y[0];
  delete [] Y;

  return timer.totalElapsedTime() / nloop;
}

double
do_time_adolc_retape(int degree, int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  double **X, **Y;

  X = new double*[2];
  X[0] = new double[degree+1];
  X[1] = new double[degree+1];
  Y = new double*[1];
  Y[0] = new double[degree+1]; 
  for (int j=0; j<=degree; j++) {
    X[0][j] = urand.number();
    X[1][j] = urand.number();
    Y[0][j] = 0.0;
  }
  adouble x1, x2, y;
  
  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (int j=0; j<nloop; j++) {
    trace_on(0);
    x1 <<= X[0][0];
    x2 <<= X[1][0];
    func(x1, x2, y);
    y >>= Y[0][0];
    trace_off();
    forward(0,1,2,degree,0,X,Y);
  }
  timer.stop();

  delete [] X[1];
  delete [] X[0];
  delete [] X;

  delete [] Y[0];
  delete [] Y;

  return timer.totalElapsedTime() / nloop;
}
#endif

int main(int argc, char* argv[]) {
  int ierr = 0;

  try {
    double t;
    int p = 2;
    int w = p+7;

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a single multiplication operation");
    int degree = 10;
    clp.setOption("degree", &degree, "Polynomial degree");
    int nloop = 1000000;
    clp.setOption("nloop", &nloop, "Number of loops");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "Times (sec) for degree = " << degree
	      << " nloop =  " << nloop << ":  " << std::endl;

    t = do_time< Sacado::Tay::Taylor<double> >(degree, nloop);
    std::cout << "Taylor:          " << std::setw(w) << t << std::endl;

    t = do_time< Sacado::Tay::CacheTaylor<double> >(degree, nloop);
    std::cout << "CacheTaylor:     " << std::setw(w) << t << std::endl;

#ifdef HAVE_ADOLC
    t = do_time_adolc(degree, nloop);
    std::cout << "ADOL-C:          " << std::setw(w) << t << std::endl;
    t = do_time_adolc_retape(degree, nloop);
    std::cout << "ADOL-C (retape): " << std::setw(w) << t << std::endl;
#endif
    
  }
  catch (std::exception& e) {
    cout << e.what() << endl;
    ierr = 1;
  }
  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

  return ierr;
}
