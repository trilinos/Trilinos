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

#include "Teuchos_Time.hpp"
#include "Teuchos_Array.hpp"
#include <fstream>

#include "fad_expr_funcs.hpp"

// A simple performance test that computes the derivative of expressions of
// various depths.

template <typename T>
void do_times_mult(const T x[], int nloop, Teuchos::Array<double>& times)
{
  times.resize(ExprFuncs<T>::nfunc);
    
  T y = 0.0;
  Teuchos::Time timer("mult", false);
  int i = 0;

  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult1(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;

  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult2(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult3(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult4(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult5(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult10(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult15(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::mult20(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
}

template <typename T>
void do_times_add(const T x[], int nloop, Teuchos::Array<double>& times)
{
  times.resize(ExprFuncs<T>::nfunc);
    
  T y = 0.0;
  Teuchos::Time timer("add", false);
  int i = 0;

  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add1(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;

  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add2(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add3(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add4(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add5(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add10(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add15(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::add20(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
}

template <typename T>
void do_times_nest(const T x[], int nloop, Teuchos::Array<double>& times)
{
  times.resize(ExprFuncs<T>::nfunc);
    
  T y = 0.0;
  Teuchos::Time timer("nest", false);
  int i = 0;

  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest1(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;

  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest2(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest3(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest4(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest5(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest10(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest15(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
  
  timer.start(true);
  for (int j=0; j<nloop; j++)
    ExprFuncs<T>::nest20(x, y);
  timer.stop();
  times[i++] = timer.totalElapsedTime() / nloop;
  y = 0.0;
}

template <typename FadType>
void print_times(const std::string& screen_name, const std::string& file_name)
{
  const int nderiv = 10;
  int deriv_dim[nderiv] = { 0, 1, 3, 5, 10, 15, 20, 30, 40, 50 };
  Sacado::Random<double> urand(0.0, 1.0);
  int p = 1;
  int w = p+7;
  std::ofstream file(file_name.c_str(), std::ios::out);

  {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << screen_name << " Relative times (time/(func_time*nderiv)): " 
	    << std::endl;
  std::cout << std::setw(5) << "deriv" << " ";
  for (int i=0; i<ExprFuncs<FadType>::nfunc; i++)
    std::cout << std::setw(w) << ExprFuncs<FadType>::mult_names[i] << " ";
  std::cout << std::endl;
  std::cout << "===== ";
  for (int i=0; i<ExprFuncs<FadType>::nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;

  // Get function evaluation times
  double x[ExprFuncs<double>::nx_max];
  for (int i=0; i<ExprFuncs<double>::nx_max; i++)
    x[i] = urand.number();
  int nloop_func = 10000000;
  Teuchos::Array<double> times_func;
  do_times_mult(x, nloop_func, times_func);
  
  // Get times for each derivative dimension
  for (int i=0; i<nderiv; i++) {
    FadType fx[ExprFuncs<FadType>::nx_max];
    for (int k=0; k<ExprFuncs<FadType>::nx_max; k++) {
      fx[k] = FadType(deriv_dim[i],  urand.number());
      for (int j=0; j<deriv_dim[i]; j++) {
	fx[k].fastAccessDx(j) = urand.number();
      }
    }
    
    int nloop = static_cast<int>(1000000.0/(deriv_dim[i]+1));
    Teuchos::Array<double> times;
    do_times_mult(fx, nloop, times);
     
    // Print times
    int d = deriv_dim[i];
    if (d == 0)
      d = 1;
    std::cout << std::setw(5) << deriv_dim[i] << " ";
    file << deriv_dim[i] << " ";
    for (int j=0; j<times.size(); j++) {
      double rel_time = times[j]/(times_func[j]*d);
      std::cout << std::setw(w) << rel_time << " ";
      file << rel_time << " ";
    }
    std::cout << std::endl;
    file << std::endl;
  }
  }

  {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << screen_name << " Relative times (time/(func_time*nderiv)): " 
	    << std::endl;
  std::cout << std::setw(5) << "deriv" << " ";
  for (int i=0; i<ExprFuncs<FadType>::nfunc; i++)
    std::cout << std::setw(w) << ExprFuncs<FadType>::add_names[i] << " ";
  std::cout << std::endl;
  std::cout << "===== ";
  for (int i=0; i<ExprFuncs<FadType>::nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;

  // Get function evaluation times
  double x[ExprFuncs<double>::nx_max];
  for (int i=0; i<ExprFuncs<double>::nx_max; i++)
    x[i] = urand.number();
  int nloop_func = 10000000;
  Teuchos::Array<double> times_func;
  do_times_add(x, nloop_func, times_func);
  
  // Get times for each derivative dimension
  for (int i=0; i<nderiv; i++) {
    FadType fx[ExprFuncs<FadType>::nx_max];
    for (int k=0; k<ExprFuncs<FadType>::nx_max; k++) {
      fx[k] = FadType(deriv_dim[i],  urand.number());
      for (int j=0; j<deriv_dim[i]; j++) {
	fx[k].fastAccessDx(j) = urand.number();
      }
    }
    
    int nloop = static_cast<int>(1000000.0/(deriv_dim[i]+1));
    Teuchos::Array<double> times;
    do_times_add(fx, nloop, times);
     
    // Print times
    int d = deriv_dim[i];
    if (d == 0)
      d = 1;
    std::cout << std::setw(5) << deriv_dim[i] << " ";
    file << deriv_dim[i] << " ";
    for (int j=0; j<times.size(); j++) {
      double rel_time = times[j]/(times_func[j]*d);
      std::cout << std::setw(w) << rel_time << " ";
      file << rel_time << " ";
    }
    std::cout << std::endl;
    file << std::endl;
  }
  }

  {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << screen_name << " Relative times (time/(func_time*nderiv)): " 
	    << std::endl;
  std::cout << std::setw(5) << "deriv" << " ";
  for (int i=0; i<ExprFuncs<FadType>::nfunc; i++)
    std::cout << std::setw(w) << ExprFuncs<FadType>::nest_names[i] << " ";
  std::cout << std::endl;
  std::cout << "===== ";
  for (int i=0; i<ExprFuncs<FadType>::nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;

  // Get function evaluation times
  double x[ExprFuncs<double>::nx_max];
  for (int i=0; i<ExprFuncs<double>::nx_max; i++)
    x[i] = urand.number();
  int nloop_func = 10000000;
  Teuchos::Array<double> times_func;
  do_times_nest(x, nloop_func, times_func);
  
  // Get times for each derivative dimension
  for (int i=0; i<nderiv; i++) {
    FadType fx[ExprFuncs<FadType>::nx_max];
    for (int k=0; k<ExprFuncs<FadType>::nx_max; k++) {
      fx[k] = FadType(deriv_dim[i],  urand.number());
      for (int j=0; j<deriv_dim[i]; j++) {
	fx[k].fastAccessDx(j) = urand.number();
      }
    }
    
    int nloop = static_cast<int>(1000000.0/(deriv_dim[i]+1));
    Teuchos::Array<double> times;
    do_times_nest(fx, nloop, times);
     
    // Print times
    int d = deriv_dim[i];
    if (d == 0)
      d = 1;
    std::cout << std::setw(5) << deriv_dim[i] << " ";
    file << deriv_dim[i] << " ";
    for (int j=0; j<times.size(); j++) {
      double rel_time = times[j]/(times_func[j]*d);
      std::cout << std::setw(w) << rel_time << " ";
      file << rel_time << " ";
    }
    std::cout << std::endl;
    file << std::endl;
  }
  }
}

int main() {
  print_times< Sacado::Fad::DFad<double> >(
    "Sacado::Fad::DFad", "fad_expr_depth_dfad.txt");
  print_times< Sacado::ELRFad::DFad<double> >(
    "Sacado::ELRFad::DFad", "fad_expr_depth_elr_dfad.txt");
  print_times< Sacado::CacheFad::DFad<double> >(
    "Sacado::CacheFad::DFad", "fad_expr_depth_cache_dfad.txt");
  print_times< Sacado::ELRCacheFad::DFad<double> >(
    "Sacado::ELRCacheFad::DFad", "fad_expr_depth_elr_cache_dfad.txt");

  return 0;
}
