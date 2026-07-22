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
#include "Teuchos_BLAS.hpp"
#include "Sacado_Fad_BLAS.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// A performance test that compares the cost of differentiating BLAS routines
// with Fad

double
do_time_teuchos_double_gemm(unsigned int m, unsigned int n, unsigned int k,
                            unsigned int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  Teuchos::BLAS<int,double> blas;

  std::vector<double> A(m*k), B(k*n), C(m*n);
  for (unsigned int j=0; j<k; j++)
    for (unsigned int i=0; i<m; i++)
      A[i+j*m] = urand.number();
  for (unsigned int j=0; j<n; j++)
    for (unsigned int i=0; i<k; i++)
    B[i+j*k] = urand.number();
  for (unsigned int j=0; j<n; j++)
    for (unsigned int i=0; i<m; i++)
      C[i+j*m] = urand.number();
  double alpha = urand.number();
  double beta = urand.number();

  Teuchos::Time timer("Teuchos Double GEMM", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, k, alpha, &A[0], m,
              &B[0], k, beta, &C[0], m);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

double
do_time_teuchos_double_gemv(unsigned int m, unsigned int n, unsigned int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  Teuchos::BLAS<int,double> blas;

  std::vector<double> A(m*n), B(n), C(m);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++)
      A[i+j*m] = urand.number();
    B[j] = urand.number();
  }
  for (unsigned int i=0; i<m; i++)
    C[i] = urand.number();
  double alpha = urand.number();
  double beta = urand.number();

  Teuchos::Time timer("Teuchos Double GEMV", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, beta, &C[0], 1);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

double
do_time_teuchos_double_dot(unsigned int m, unsigned int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  Teuchos::BLAS<int,double> blas;

  std::vector<double> X(m), Y(m);
  for (unsigned int i=0; i<m; i++) {
    X[i] = urand.number();
    Y[i] = urand.number();
  }

  Teuchos::Time timer("Teuchos Double DOT", false);
  timer.start(true);
  double z = 0.0;
  for (unsigned int j=0; j<nloop; j++) {
    z += blas.DOT(m, &X[0], 1, &Y[0], 1);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

template <typename FadType>
double
do_time_teuchos_fad_gemm(unsigned int m, unsigned int n, unsigned int k,
                         unsigned int ndot, unsigned int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  Teuchos::BLAS<int,FadType> blas;

  std::vector<FadType> A(m*k), B(k*n), C(m*n);
  for (unsigned int j=0; j<k; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int l=0; l<ndot; l++)
        A[i+j*m].fastAccessDx(l) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<k; i++) {
      B[i+j*k] = FadType(ndot, urand.number());
      for (unsigned int l=0; l<ndot; l++)
        B[i+j*k].fastAccessDx(l) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      C[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int l=0; l<ndot; l++)
        C[i+j*m].fastAccessDx(l) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int l=0; l<ndot; l++) {
    alpha.fastAccessDx(l) = urand.number();
    beta.fastAccessDx(l) = urand.number();
  }

  Teuchos::Time timer("Teuchos Fad GEMM", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, k, alpha, &A[0], m,
              &B[0], k, beta, &C[0], m);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

template <typename FadType>
double
do_time_teuchos_fad_gemv(unsigned int m, unsigned int n, unsigned int ndot,
                         unsigned int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  Teuchos::BLAS<int,FadType> blas;

  std::vector<FadType> A(m*n), B(n), C(m);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      //A[i+j*m] = urand.number();
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = urand.number();
    }
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<m; i++) {
    C[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      C[i].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  Teuchos::Time timer("Teuchos Fad GEMV", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, beta, &C[0], 1);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

template <typename FadType>
double
do_time_teuchos_fad_dot(unsigned int m, unsigned int ndot, unsigned int nloop)
{
  Sacado::Random<double> urand(0.0, 1.0);
  Teuchos::BLAS<int,FadType> blas;

  std::vector<FadType> X(m), Y(m);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, urand.number());
    Y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
      Y[i].fastAccessDx(k) = urand.number();
    }
  }

  Teuchos::Time timer("Teuchos Fad DOT", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    FadType z = blas.DOT(m, &X[0], 1, &Y[0], 1);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

template <typename FadType>
double
do_time_sacado_fad_gemm(unsigned int m, unsigned int n, unsigned int k,
                        unsigned int ndot, unsigned int nloop, bool use_dynamic)
{
  Sacado::Random<double> urand(0.0, 1.0);
  unsigned int sz = (m*k+k*n+m*n)*(1+ndot);
  Teuchos::BLAS<int,FadType> blas(false,use_dynamic,sz);

  Sacado::Fad::Vector<unsigned int, FadType> A(m*k,ndot), B(k*n,ndot), C
    (m*n,ndot);
  for (unsigned int j=0; j<k; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int l=0; l<ndot; l++)
        A[i+j*m].fastAccessDx(l) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<k; i++) {
      B[i+j*k] = FadType(ndot, urand.number());
      for (unsigned int l=0; l<ndot; l++)
        B[i+j*k].fastAccessDx(l) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      C[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int l=0; l<ndot; l++)
        C[i+j*m].fastAccessDx(l) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int l=0; l<ndot; l++) {
    alpha.fastAccessDx(l) = urand.number();
    beta.fastAccessDx(l) = urand.number();
  }

  Teuchos::Time timer("Teuchos Fad GEMM", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, k, alpha, &A[0], m,
              &B[0], k, beta, &C[0], m);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

template <typename FadType>
double
do_time_sacado_fad_gemv(unsigned int m, unsigned int n, unsigned int ndot,
                        unsigned int nloop, bool use_dynamic)
{
  Sacado::Random<double> urand(0.0, 1.0);
  unsigned int sz = m*n*(1+ndot) + 2*n*(1+ndot);
  Teuchos::BLAS<int,FadType> blas(false,use_dynamic,sz);

  Sacado::Fad::Vector<unsigned int, FadType> A(m*n,ndot), B(n,ndot), C(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      //A[i+j*m] = urand.number();
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = urand.number();
    }
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<m; i++) {
    C[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      C[i].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  Teuchos::Time timer("Teuchos Fad GEMV", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, beta, &C[0], 1);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

template <typename FadType>
double
do_time_sacado_fad_dot(unsigned int m, unsigned int ndot,
                       unsigned int nloop, bool use_dynamic)
{
  Sacado::Random<double> urand(0.0, 1.0);
  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> blas(false,use_dynamic,sz);

  Sacado::Fad::Vector<unsigned int, FadType> X(m,ndot), Y(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, urand.number());
    Y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
      Y[i].fastAccessDx(k) = urand.number();
    }
  }

  Teuchos::Time timer("Teuchos Fad DOT", false);
  timer.start(true);
  for (unsigned int j=0; j<nloop; j++) {
    FadType z = blas.DOT(m, &X[0], 1, &Y[0], 1);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

int main(int argc, char* argv[]) {
  int ierr = 0;

  try {
    double t, tb;
    int p = 2;
    int w = p+7;

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of differentiating BLAS routines using Fad");
    int m = 10;
    clp.setOption("m", &m, "Number of rows");
    int n = 10;
    clp.setOption("n", &n, "Number of columns");
    int k = 10;
    clp.setOption("k", &k, "Number of columns for GEMM");
    int ndot = 10;
    clp.setOption("ndot", &ndot, "Number of derivative components");
    int nloop = 100000;
    clp.setOption("nloop", &nloop, "Number of loops");
    int dynamic = 1;
    clp.setOption("dynamic", &dynamic, "Use dynamic allocation");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;
    bool use_dynamic = (dynamic != 0);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "Times (sec) for m = " << m << ", n = " << n
              << ", ndot = " << ndot << ", nloop =  " << nloop
              << ", dynamic = " << use_dynamic << ":  "
              << std::endl;

    tb = do_time_teuchos_double_gemm(m,n,k,nloop);
    std::cout << "GEMM:                 " << std::setw(w) << tb << std::endl;

    t = do_time_sacado_fad_gemm< Sacado::Fad::DVFad<double> >(m,n,k,ndot,nloop,use_dynamic);
    std::cout << "Sacado DVFad GEMM:    " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    t = do_time_sacado_fad_gemm< Sacado::Fad::DFad<double> >(m,n,k,ndot,nloop,use_dynamic);
    std::cout << "Sacado DFad GEMM:     " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    t = do_time_teuchos_fad_gemm< Sacado::Fad::DFad<double> >(m,n,k,ndot,nloop);
    std::cout << "Teuchos DFad GEMM:    " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    // t = do_time_teuchos_fad_gemm< Sacado::ELRFad::DFad<double> >(m,n,k,ndot,nloop);
    // std::cout << "Teuchos ELRDFad GEMM:  " << std::setw(w) << t << "\t"
    //        << std::setw(w) << t/tb << std::endl;

    t = do_time_teuchos_fad_gemm< Sacado::Fad::DVFad<double> >(m,n,k,ndot,nloop);
    std::cout << "Teuchos DVFad GEMM:   " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    std::cout << std::endl;

    tb = do_time_teuchos_double_gemv(m,n,nloop);
    std::cout << "GEMV:                 " << std::setw(w) << tb << std::endl;

    t = do_time_sacado_fad_gemv< Sacado::Fad::DVFad<double> >(m,n,ndot,nloop*10,use_dynamic);
    std::cout << "Sacado DVFad GEMV:    " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    t = do_time_sacado_fad_gemv< Sacado::Fad::DFad<double> >(m,n,ndot,nloop*10,use_dynamic);
    std::cout << "Sacado DFad GEMV:     " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    t = do_time_teuchos_fad_gemv< Sacado::Fad::DFad<double> >(m,n,ndot,nloop*10);
    std::cout << "Teuchos DFad GEMV:    " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    // t = do_time_teuchos_fad_gemv< Sacado::ELRFad::DFad<double> >(m,n,ndot,nloop*10);
    // std::cout << "Teuchos ELRDFad GEMV:  " << std::setw(w) << t << "\t"
    //        << std::setw(w) << t/tb << std::endl;

    t = do_time_teuchos_fad_gemv< Sacado::Fad::DVFad<double> >(m,n,ndot,nloop*10);
    std::cout << "Teuchos DVFad GEMV:   " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    std::cout << std::endl;

    tb = do_time_teuchos_double_dot(m,nloop*100);
    std::cout << "DOT:                  " << std::setw(w) << tb << std::endl;

    t = do_time_sacado_fad_dot< Sacado::Fad::DVFad<double> >(m,ndot,nloop*100,use_dynamic);
    std::cout << "Sacado DVFad DOT:     " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    t = do_time_sacado_fad_dot< Sacado::Fad::DFad<double> >(m,ndot,nloop*100,use_dynamic);
    std::cout << "Sacado DFad DOT:      " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    t = do_time_teuchos_fad_dot< Sacado::Fad::DFad<double> >(m,ndot,nloop*100);
    std::cout << "Teuchos DFad DOT:     " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

    // t = do_time_teuchos_fad_dot< Sacado::ELRFad::DFad<double> >(m,ndot,nloop*100);
    // std::cout << "Teuchos ELRDFad DOT:  " << std::setw(w) << t << "\t"
    //        << std::setw(w) << t/tb << std::endl;

    t = do_time_teuchos_fad_dot< Sacado::Fad::DVFad<double> >(m,ndot,nloop*100);
    std::cout << "Teuchos DVFad DOT:    " << std::setw(w) << t << "\t"
              << std::setw(w) << t/tb << std::endl;

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
