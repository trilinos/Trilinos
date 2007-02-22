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

//#include "Sacado_Random.hpp"
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_Fad_DMFad.hpp"
#include "Sacado_Fad_MemPoolManager.hpp"
#include "Sacado_Fad_SFad.hpp"
#include "Sacado_Fad_SLFad.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Fad/fad.h"
#include "TinyFadET/tfad.h"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;

void FAD::error(char *msg) {
  std::cout << msg << std::endl;
}

template <typename T>
inline void
mult1(const T& x1, const T& x2, T& y) {
  y = x1*x2;
}

template <typename FadType>
double
do_time(int nderiv, int nloop)
{
  FadType x1, x2, y;
//   Sacado::Random urand(0.0, 1.0);

//   x1 = FadType(nderiv,  urand.number());
//   x2 = FadType(nderiv,  urand.number());
  x1 = FadType(nderiv,  2.0);
  x2 = FadType(nderiv,  3.0);
  y = 0.0;
  for (int j=0; j<nderiv; j++) {
//     x1.fastAccessDx(j) = urand.number();
//     x2.fastAccessDx(j) = urand.number();
    x1.fastAccessDx(j) = 2.0 + static_cast<double>(j)/10.0;
    x2.fastAccessDx(j) = 3.0 + static_cast<double>(j)/10.0;;
  }
  
  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (int j=0; j<nloop; j++) {
    mult1(x1, x2, y);
  }
  timer.stop();

  return timer.totalElapsedTime() / nloop;
}

int main(int argc, char* argv[]) {
  int ierr = 0;

  try {
    double t;
    int p = 2;
    int w = p+7;

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a single multiplication operation");
    int nderiv = 10;
    clp.setOption("nderiv", &nderiv, "Number of derivative components");
    int nloop = 10000000;
    clp.setOption("nloop", &nloop, "Number of loops");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    // Memory pool & manager
    Sacado::Fad::MemPoolManager<double> poolManager(10);
    Sacado::Fad::MemPool* pool = poolManager.getMemoryPool(nderiv);
    Sacado::Fad::DMFad<double>::setDefaultPool(pool);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "Times (sec) for nderiv = " << nderiv 
	      << " nloop =  " << nloop << ":  " << std::endl;

    t = do_time< FAD::Fad<double> >(nderiv, nloop);
    std::cout << "Fad:       " << std::setw(w) << t << std::endl;

    t = do_time< FAD::TFad<10,double> >(nderiv, nloop);
    std::cout << "TFad:      " << std::setw(w) << t << std::endl;
    
    t = do_time< Sacado::Fad::DFad<double> >(nderiv, nloop);
    std::cout << "DFad:      " << std::setw(w) << t << std::endl;

   t = do_time< Sacado::Fad::DMFad<double> >(nderiv, nloop);
    std::cout << "DMFad:     " << std::setw(w) << t << std::endl; 

    t = do_time< Sacado::Fad::SFad<double,10> >(nderiv, nloop);
    std::cout << "SFad:      " << std::setw(w) << t << std::endl;

    t = do_time< Sacado::Fad::SLFad<double,10> >(nderiv, nloop);
    std::cout << "SLFad:     " << std::setw(w) << t << std::endl;
    
    t = do_time< Sacado::CacheFad::DFad<double> >(nderiv, nloop);
    std::cout << "CacheFad:  " << std::setw(w) << t << std::endl;
    
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
