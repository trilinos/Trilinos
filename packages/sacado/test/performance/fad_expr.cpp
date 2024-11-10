// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"
#include "Sacado_Fad_SimpleFad.hpp"

#include "Fad/fad.h"
#include "TinyFadET/tfad.h"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// A simple performance test that computes the derivative of a simple
// expression using many variants of Fad.

void FAD::error(const char *msg) {
  std::cout << msg << std::endl;
}

template <typename T>
inline void
func1(const T& x1, const T& x2, T& y) {
  y = (x1*x2 + sin(x1)/x2);
}

inline void
func1_and_deriv(int n, double x1, double x2, double* x1dot, double* x2dot,
		double& y, double* ydot) {
  double s = sin(x1);
  double c = cos(x1);
  double t = s/x2;
  double t1 = x2 + c/x2;
  double t2 = x1 - t/x2;
  y = x1*x2 + t;
  for (int i=0; i<10; i++)
    ydot[i] = t1*x1dot[i] + t2*x2dot[i];
}

template <typename FadType>
double
do_time(int nderiv, int nloop)
{
  FadType x1, x2, y;
  Sacado::Random<double> urand(0.0, 1.0);

  x1 = FadType(nderiv,  urand.number());
  x2 = FadType(nderiv,  urand.number());
  y = 0.0;
  for (int j=0; j<nderiv; j++) {
    x1.fastAccessDx(j) = urand.number();
    x2.fastAccessDx(j) = urand.number();
  }

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (int j=0; j<nloop; j++) {
    func1(x1, x2, y);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

double
do_time_analytic(int nderiv, int nloop)
{
  double x1, x2, y;
  double *x1dot, *x2dot, *ydot;
  Sacado::Random<double> urand(0.0, 1.0);

  x1 = urand.number();
  x2 = urand.number();
  y = 0.0;
  x1dot = new double[nderiv];
  x2dot = new double[nderiv];
  ydot = new double[nderiv];
  for (int j=0; j<nderiv; j++) {
    x1dot[j] = urand.number();
    x2dot[j] = urand.number();
  }

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (int j=0; j<nloop; j++) {
    func1_and_deriv(nderiv, x1, x2, x1dot, x2dot, y, ydot);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

int main(int argc, char* argv[]) {
  int ierr = 0;

  try {
    double t, ta;
    int p = 2;
    int w = p+7;

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a single multiplication operation");
    int nderiv = 10;
    clp.setOption("nderiv", &nderiv, "Number of derivative components");
    int nloop = 1000000;
    clp.setOption("nloop", &nloop, "Number of loops");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "Times (sec) for nderiv = " << nderiv
	      << " nloop =  " << nloop << ":  " << std::endl;

    ta = do_time_analytic(nderiv, nloop);
    std::cout << "Analytic:  " << std::setw(w) << ta << std::endl;

    t = do_time< Sacado::Fad::SimpleFad<double> >(nderiv, nloop);
    std::cout << "SimpleFad: " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< FAD::TFad<10,double> >(nderiv, nloop);
    std::cout << "TFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< FAD::Fad<double> >(nderiv, nloop);
    std::cout << "Fad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::SFad<double,10> >(nderiv, nloop);
    std::cout << "SFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::SLFad<double,10> >(nderiv, nloop);
    std::cout << "SLFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::DFad<double> >(nderiv, nloop);
    std::cout << "DFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::ELRFad::SFad<double,10> >(nderiv, nloop);
    std::cout << "ELRSFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::ELRFad::SLFad<double,10> >(nderiv, nloop);
    std::cout << "ELRSLFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::ELRFad::DFad<double> >(nderiv, nloop);
    std::cout << "ELRDFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::CacheFad::DFad<double> >(nderiv, nloop);
    std::cout << "CacheFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::DVFad<double> >(nderiv, nloop);
    std::cout << "DVFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  return ierr;
}
