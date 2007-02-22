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
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Fad/fad.h"
#include "Teuchos_Time.hpp"

void FAD::error(char *msg) {
  std::cout << msg << std::endl;
}

template <typename T>
inline void
mult1(const T& x1, const T& x2, T& y) {
  y = x1*x2;
}

template <typename T>
inline void
mult2(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1;
}

template <typename T>
inline void
mult3(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1*x2;
}

template <typename T>
inline void
mult4(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1*x2*x1;
}

template <typename T>
inline void
mult5(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1*x2*x1*x2;
}

template <typename T>
inline void
mult10(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1;
}

template <typename T>
inline void
mult15(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2;
}

template <typename T>
inline void
mult20(const T& x1, const T& x2, T& y) {
  y = x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1*x2*x1;
}

template <typename FadType>
void
do_times(const std::string& name)
{
  const int nfunc = 8;
  const int nderiv = 10;
  int deriv_dim[nderiv];
  int nloop[nderiv];
  double times[nfunc][nderiv];
  int p = 1;
  int w = p+7;

  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << name << " Times (sec): " << std::endl;
  std::cout << std::setw(5) << "deriv" << " "
	    << std::setw(w) << "mult1" << " "
	    << std::setw(w) << "mult2" << " "
	    << std::setw(w) << "mult3" << " "
	    << std::setw(w) << "mult4" << " "
	    << std::setw(w) << "mult5" << " "
	    << std::setw(w) << "mult10" << " "
	    << std::setw(w) << "mult15" << " "
	    << std::setw(w) << "mult20" << std::endl;
  std::cout << "===== ";
  for (int i=0; i<nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;
  
  for (int i=0; i<5; i++)
    deriv_dim[i] = i;
  for (int i=5; i<nderiv; i++)
    deriv_dim[i] = 5*(i-4);
  for (int i=0; i<nderiv; i++)
    nloop[i] = static_cast<int>(1000000.0/(deriv_dim[i]+1));

  FadType x1, x2, y;
  Sacado::Random urand(0.0, 1.0);
  for (int i=0; i<nderiv; i++) {
    std::cout << std::setw(5) << deriv_dim[i] << " ";

    x1 = FadType(deriv_dim[i],  urand.number());
    x2 = FadType(deriv_dim[i],  urand.number());
    y = 0.0;
    for (int j=0; j<deriv_dim[i]; j++) {
      x1.fastAccessDx(j) = urand.number();
      x2.fastAccessDx(j) = urand.number();
    }

    Teuchos::Time timer("mult", false);

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult1(x1, x2, y);
    timer.stop();
    times[0][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[0][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult2(x1, x2, y);
    timer.stop();
    times[1][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[1][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult3(x1, x2, y);
    timer.stop();
    times[2][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[2][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult4(x1, x2, y);
    timer.stop();
    times[3][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[3][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult5(x1, x2, y);
    timer.stop();
    times[4][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[4][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult10(x1, x2, y);
    timer.stop();
    times[5][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[5][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult15(x1, x2, y);
    timer.stop();
    times[6][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[6][i] << " ";

    timer.start(true);
    for (int j=0; j<nloop[i]; j++)
      mult20(x1, x2, y);
    timer.stop();
    times[7][i] = timer.totalElapsedTime() / nloop[i];
    y = 0.0;
    std::cout << std::setw(w) << times[7][i] << " ";

    std::cout << endl;
  }
}

int main() {
  do_times< FAD::Fad<double> >("FAD::Fad");
  do_times< Sacado::Fad::DFad<double> >("Sacado::Fad::DFad");
  do_times< Sacado::CacheFad::DFad<double> >("Sacado::CacheFad::DFad");

  return 0;
}
