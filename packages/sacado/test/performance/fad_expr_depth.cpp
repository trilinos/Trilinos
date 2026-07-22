// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_Random.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_Array.hpp"
#include <fstream>

#include "fad_expr_funcs.hpp"

// A simple performance test that computes the derivative of expressions of
// various depths.

template <typename T, typename F>
double do_time(const T x[], int nloop, const F& f)
{
  T y = 0.0;
  Teuchos::Time timer("F", false);

  timer.start(true);
  for (int j=0; j<nloop; j++)
    f(x, y);
  timer.stop();
  return timer.totalElapsedTime() / nloop;
}

template <typename T, template <typename, int> class F>
void do_times(const T x[], int nloop, Teuchos::Array<double>& times)
{
  int i = 0;
  times[i++] = do_time(x, nloop, F<T,1>());
  times[i++] = do_time(x, nloop, F<T,2>());
  times[i++] = do_time(x, nloop, F<T,3>());
  times[i++] = do_time(x, nloop, F<T,4>());
  times[i++] = do_time(x, nloop, F<T,5>());
  times[i++] = do_time(x, nloop, F<T,10>());
  times[i++] = do_time(x, nloop, F<T,15>());
  times[i++] = do_time(x, nloop, F<T,20>());
}

#ifdef HAVE_ADOLC
template <typename F>
double do_time_adolc(double *x, double **seed, int d, int nloop,
                     int tag, const F& f)
{
  Teuchos::Time timer("F", false);
  int n = F::n;
  trace_on(tag);
  adouble *x_ad = new adouble[n];
  for (int i=0; i<n; i++)
    x_ad[i] <<= x[i];
  adouble y_ad = 0.0;
  f(x_ad, y_ad);
  double y;
  y_ad >>= y;
  delete [] x_ad;
  trace_off();

  double **jac = new double*[1];
  jac[0] = new double[d];
  timer.start(true);
  for (int j=0; j<nloop; j++)
    forward(tag, 1, n, d, x, seed, &y, jac);
  timer.stop();

  delete [] jac[0];
  delete [] jac;

  return timer.totalElapsedTime() / nloop;
}

template <template <typename,int> class F>
void do_times_adolc(double *x, double **seed, int d, int nloop,
                    int& tag, Teuchos::Array<double>& times)
{
  int i = 0;
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,1>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,2>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,3>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,4>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,5>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,10>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,15>());
  times[i++] = do_time_adolc(x, seed, d, nloop, tag++, F<adouble,20>());
}
#endif

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
  for (int i=0; i<ExprFuncs::nfunc; i++)
    std::cout << std::setw(w) << ExprFuncs::mult_names[i] << " ";
  std::cout << std::endl;
  std::cout << "===== ";
  for (int i=0; i<ExprFuncs::nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;

  // Get function evaluation times
  double x[ExprFuncs::nx_max];
  for (int i=0; i<ExprFuncs::nx_max; i++)
    x[i] = urand.number();
  int nloop_func = 10000000;
  Teuchos::Array<double> times_func(ExprFuncs::nfunc);
  do_times<double,ExprFuncs::mult>(x, nloop_func, times_func);

  // Get times for each derivative dimension
  for (int i=0; i<nderiv; i++) {
    FadType fx[ExprFuncs::nx_max];
    for (int k=0; k<ExprFuncs::nx_max; k++) {
      fx[k] = FadType(deriv_dim[i],  urand.number());
      for (int j=0; j<deriv_dim[i]; j++) {
        fx[k].fastAccessDx(j) = urand.number();
      }
    }

    //int nloop = static_cast<int>(1000000.0/(deriv_dim[i]+1));
    int nloop = static_cast<int>(100000.0);
    Teuchos::Array<double> times(ExprFuncs::nfunc);
    do_times<FadType,ExprFuncs::mult>(fx, nloop, times);

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
  for (int i=0; i<ExprFuncs::nfunc; i++)
    std::cout << std::setw(w) << ExprFuncs::add_names[i] << " ";
  std::cout << std::endl;
  std::cout << "===== ";
  for (int i=0; i<ExprFuncs::nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;

  // Get function evaluation times
  double x[ExprFuncs::nx_max];
  for (int i=0; i<ExprFuncs::nx_max; i++)
    x[i] = urand.number();
  int nloop_func = 10000000;
  Teuchos::Array<double> times_func(ExprFuncs::nfunc);
  do_times<double,ExprFuncs::add>(x, nloop_func, times_func);

  // Get times for each derivative dimension
  for (int i=0; i<nderiv; i++) {
    FadType fx[ExprFuncs::nx_max];
    for (int k=0; k<ExprFuncs::nx_max; k++) {
      fx[k] = FadType(deriv_dim[i],  urand.number());
      for (int j=0; j<deriv_dim[i]; j++) {
        fx[k].fastAccessDx(j) = urand.number();
      }
    }

    //int nloop = static_cast<int>(1000000.0/(deriv_dim[i]+1));
    int nloop = static_cast<int>(100000.0);
    Teuchos::Array<double> times(ExprFuncs::nfunc);
    do_times<FadType,ExprFuncs::add>(fx, nloop, times);

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
  for (int i=0; i<ExprFuncs::nfunc; i++)
    std::cout << std::setw(w) << ExprFuncs::nest_names[i] << " ";
  std::cout << std::endl;
  std::cout << "===== ";
  for (int i=0; i<ExprFuncs::nfunc; i++) {
    for (int j=0; j<w; j++)
      std::cout << '=';
    std::cout << " ";
  }
  std::cout << std::endl;

  // Get function evaluation times
  double x[ExprFuncs::nx_max];
  for (int i=0; i<ExprFuncs::nx_max; i++)
    x[i] = urand.number();
  int nloop_func = 10000000;
  Teuchos::Array<double> times_func(ExprFuncs::nfunc);
  do_times<double,ExprFuncs::nest>(x, nloop_func, times_func);

  // Get times for each derivative dimension
  for (int i=0; i<nderiv; i++) {
    FadType fx[ExprFuncs::nx_max];
    for (int k=0; k<ExprFuncs::nx_max; k++) {
      fx[k] = FadType(deriv_dim[i],  urand.number());
      for (int j=0; j<deriv_dim[i]; j++) {
        fx[k].fastAccessDx(j) = urand.number();
      }
    }

    //int nloop = static_cast<int>(1000000.0/(deriv_dim[i]+1));
    int nloop = static_cast<int>(100000.0);
    Teuchos::Array<double> times(ExprFuncs::nfunc);
    do_times<FadType,ExprFuncs::nest>(fx, nloop, times);

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

#ifdef HAVE_ADOLC
void print_times_adolc(const std::string& screen_name,
                       const std::string& file_name)
{
  const int nderiv = 10;
  int deriv_dim[nderiv] = { 0, 1, 3, 5, 10, 15, 20, 30, 40, 50 };
  const int deriv_max = 50;
  Sacado::Random<double> urand(0.0, 1.0);
  int p = 1;
  int w = p+7;
  std::ofstream file(file_name.c_str(), std::ios::out);
  int tag = 0;

  double **seed = new double*[ExprFuncs::nx_max];
  for (int i=0; i<ExprFuncs::nx_max; i++) {
    seed[i] = new double[deriv_max];
    for (int j=0; j<deriv_max; j++)
      seed[i][j] = urand.number();
  }

  {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << screen_name << " Relative times (time/(func_time*nderiv)): "
              << std::endl;
    std::cout << std::setw(5) << "deriv" << " ";
    for (int i=0; i<ExprFuncs::nfunc; i++)
      std::cout << std::setw(w) << ExprFuncs::mult_names[i] << " ";
    std::cout << std::endl;
    std::cout << "===== ";
    for (int i=0; i<ExprFuncs::nfunc; i++) {
      for (int j=0; j<w; j++)
        std::cout << '=';
      std::cout << " ";
    }
    std::cout << std::endl;

    // Get function evaluation times
    double x[ExprFuncs::nx_max];
    for (int i=0; i<ExprFuncs::nx_max; i++)
      x[i] = urand.number();
    int nloop_func = 10000000;
    Teuchos::Array<double> times_func(ExprFuncs::nfunc);
    do_times<double,ExprFuncs::mult>(x, nloop_func, times_func);

    // Get times for each derivative dimension
    for (int i=0; i<nderiv; i++) {
      //int nloop = static_cast<int>(100000.0/(deriv_dim[i]+1));
      int nloop = static_cast<int>(10000.0);
      Teuchos::Array<double> times(ExprFuncs::nfunc);
      do_times_adolc<ExprFuncs::mult>(x, seed, deriv_dim[i], nloop, tag, times);

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
    for (int i=0; i<ExprFuncs::nfunc; i++)
      std::cout << std::setw(w) << ExprFuncs::add_names[i] << " ";
    std::cout << std::endl;
    std::cout << "===== ";
    for (int i=0; i<ExprFuncs::nfunc; i++) {
      for (int j=0; j<w; j++)
        std::cout << '=';
      std::cout << " ";
    }
    std::cout << std::endl;

    // Get function evaluation times
    double x[ExprFuncs::nx_max];
    for (int i=0; i<ExprFuncs::nx_max; i++)
      x[i] = urand.number();
    int nloop_func = 10000000;
    Teuchos::Array<double> times_func(ExprFuncs::nfunc);
    do_times<double,ExprFuncs::add>(x, nloop_func, times_func);

    // Get times for each derivative dimension
     for (int i=0; i<nderiv; i++) {
       //int nloop = static_cast<int>(100000.0/(deriv_dim[i]+1));
       int nloop = static_cast<int>(10000.0);
       Teuchos::Array<double> times(ExprFuncs::nfunc);
       do_times_adolc<ExprFuncs::add>(x, seed, deriv_dim[i], nloop, tag, times);

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
    for (int i=0; i<ExprFuncs::nfunc; i++)
      std::cout << std::setw(w) << ExprFuncs::nest_names[i] << " ";
    std::cout << std::endl;
    std::cout << "===== ";
    for (int i=0; i<ExprFuncs::nfunc; i++) {
      for (int j=0; j<w; j++)
        std::cout << '=';
      std::cout << " ";
    }
    std::cout << std::endl;

    // Get function evaluation times
    double x[ExprFuncs::nx_max];
    for (int i=0; i<ExprFuncs::nx_max; i++)
      x[i] = urand.number();
    int nloop_func = 10000000;
    Teuchos::Array<double> times_func(ExprFuncs::nfunc);
    do_times<double,ExprFuncs::nest>(x, nloop_func, times_func);

    // Get times for each derivative dimension
    for (int i=0; i<nderiv; i++) {
      //int nloop = static_cast<int>(10000.0/(deriv_dim[i]+1));
      int nloop = static_cast<int>(1000.0);
      Teuchos::Array<double> times(ExprFuncs::nfunc);
      do_times_adolc<ExprFuncs::nest>(x, seed, deriv_dim[i], nloop, tag, times);

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

  delete [] seed;
}
#endif

int main() {
  print_times< Sacado::Fad::DFad<double> >(
    "Sacado::Fad::DFad", "fad_expr_depth_dfad.txt");
  print_times< Sacado::ELRFad::DFad<double> >(
    "Sacado::ELRFad::DFad", "fad_expr_depth_elr_dfad.txt");
  print_times< Sacado::CacheFad::DFad<double> >(
    "Sacado::CacheFad::DFad", "fad_expr_depth_cache_dfad.txt");
  print_times< Sacado::ELRCacheFad::DFad<double> >(
    "Sacado::ELRCacheFad::DFad", "fad_expr_depth_elr_cache_dfad.txt");
  print_times< Sacado::Fad::SimpleFad<double> >(
    "Sacado::Fad::SimpleFad", "fad_expr_depth_simple_fad.txt");
#ifdef HAVE_ADOLC
  print_times_adolc("ADOL-C", "fad_expr_depth_adolc.txt");
#endif

  return 0;
}
