// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Sacado_Random.hpp"
#include "Sacado.hpp"

#include "Fad/fad.h"
#include "TinyFadET/tfad.h"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// A simple performance test that computes the derivative of a simple
// expression using many variants of Fad.

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;

void FAD::error(const char *msg) {
  std::cout << msg << std::endl;
}

namespace {
  double xi[3], xj[3], pa[4], f[3], delr[3];
}

template <typename T>
inline T
vec3_distsq(const T xi[], const double xj[]) {
  T delr0 = xi[0]-xj[0];
  T delr1 = xi[1]-xj[1];
  T delr2 = xi[2]-xj[2];
  return delr0*delr0 + delr1*delr1 + delr2*delr2;
}

template <typename T>
inline T
vec3_distsq(const T xi[], const double xj[], T delr[]) {
  delr[0] = xi[0]-xj[0];
  delr[1] = xi[1]-xj[1];
  delr[2] = xi[2]-xj[2];
  return delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];
}

template <typename T>
inline void
lj(const T xi[], const double xj[], T& energy) {
  T delr2 = vec3_distsq(xi,xj);
  T delr_2 = 1.0/delr2;
  T delr_6 = delr_2*delr_2*delr_2;
  energy = (pa[1]*delr_6 - pa[2])*delr_6 - pa[3];
}

inline void
lj_and_grad(const double xi[], const double xj[], double& energy,
	    double f[]) {
  double delr2 = vec3_distsq(xi,xj,delr);
  double delr_2 = 1.0/delr2;
  double delr_6 = delr_2*delr_2*delr_2;
  energy = (pa[1]*delr_6 - pa[2])*delr_6 - pa[3];
  double tmp = (-12.0*pa[1]*delr_6 - 6.0*pa[2])*delr_6*delr_2;
  f[0] = delr[0]*tmp;
  f[1] = delr[1]*tmp;
  f[2] = delr[2]*tmp;
}

template <typename FadType>
double
do_time(int nloop)
{
  Teuchos::Time timer("lj", false);
  FadType xi_fad[3], energy;

  for (int i=0; i<3; i++) {
    xi_fad[i] = FadType(3, i, xi[i]);
  }
  
  timer.start(true);
  for (int j=0; j<nloop; j++) {

    lj(xi_fad, xj, energy);

    for (int i=0; i<3; i++)
      f[i] += -energy.fastAccessDx(i);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

double
do_time_analytic(int nloop)
{
  Teuchos::Time timer("lj", false);
  double energy, ff[3];

  timer.start(true);
  for (int j=0; j<nloop; j++) {

    lj_and_grad(xi, xj, energy, ff);

    for (int i=0; i<3; i++)
      f[i] += -ff[i];

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
    int nloop = 1000000;
    clp.setOption("nloop", &nloop, "Number of loops");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    // Memory pool & manager
    Sacado::Fad::MemPoolManager<double> poolManager(3);
    Sacado::Fad::MemPool* pool = poolManager.getMemoryPool(3);
    Sacado::Fad::DMFad<double>::setDefaultPool(pool);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "Times (sec) nloop =  " << nloop << ":  " << std::endl;

    Sacado::Random urand(0.0, 1.0);
    for (int i=0; i<3; i++) {
      xi[i] = urand.number();
      xj[i] = urand.number();
      pa[i] = urand.number();
    }
    pa[3] = urand.number();

    ta = do_time_analytic(nloop);
    std::cout << "Analytic:  " << std::setw(w) << ta << std::endl;

    t = do_time< FAD::TFad<3,double> >(nloop);
    std::cout << "TFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< FAD::Fad<double> >(nloop);
    std::cout << "Fad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::SFad<double,3> >(nloop);
    std::cout << "SFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::SLFad<double,3> >(nloop);
    std::cout << "SLFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;
    
    t = do_time< Sacado::Fad::DFad<double> >(nloop);
    std::cout << "DFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::DMFad<double> >(nloop);
    std::cout << "DMFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl; 

    t = do_time< Sacado::ELRFad::SFad<double,3> >(nloop);
    std::cout << "ELRSFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::ELRFad::SLFad<double,3> >(nloop);
    std::cout << "ELRSLFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::ELRFad::DFad<double> >(nloop);
    std::cout << "ELRDFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;
    
    t = do_time< Sacado::CacheFad::DFad<double> >(nloop);
    std::cout << "CacheFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;

    t = do_time< Sacado::Fad::DVFad<double> >(nloop);
    std::cout << "DVFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << std::endl;
    
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
