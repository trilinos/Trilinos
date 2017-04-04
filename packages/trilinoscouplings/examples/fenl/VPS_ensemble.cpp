#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include "VPS_ensemble.hpp"

EnsembleVPS::EnsembleVPS(const size_t dim,
                         const size_t num_samples,
                         const size_t ensemble_size,
                         const int rank)
{
  // **********************
  // PARAMETERS
  // **********************

  _num_dim = dim; // Dimensionality
  _num_fun = 1; // DON'T CHANGE
  _num_itr = 1; // DON'T CHANGE

  _num_samples = num_samples;      // Sample budget
  _ensemble_size = ensemble_size;     // Ensemble Size
  _sampling = 0; // MC
  _surrogate_order = 1;

  _MCratio = 5;

  _maxMisses = 100;
  _shrinkage_rho = 0.5;

  _proc_rank = rank;

  _eval_test_function = SmoothHerbie;

  initiate_random_number_generator(1234567890);

  // ----------------------------
  // SET DOMAIN BOUNDS and VOLUME
  // ----------------------------
  _xmin = new double[_num_dim];
  _xmax = new double[_num_dim];
  _Dvol = 1.0;
  _diag = 0.0;

  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    _xmin[idim] = -1.0; _xmax[idim] = 1.0;
    _Dvol = _Dvol * (_xmax[idim] - _xmin[idim]);
    double dx = _xmax[idim] - _xmin[idim];
    _diag += dx * dx;
  }
  _diag = sqrt(_diag);

  // ----------------------------
  // INITIATE SAMPLES & FUNCTIONS
  // ----------------------------
  _x = new double*[_num_samples];
  _f = new double*[_num_samples];
  _g = new double*[_num_samples];

  for (size_t iSample = 0; iSample < _num_samples; iSample++)
  {
    _x[iSample] = new double[_num_dim];
    for (size_t idim = 0; idim < _num_dim; idim++) _x[iSample][idim] = 0.0;

    _f[iSample] = new double[_num_fun];
    for (size_t ifun = 0; ifun < _num_fun; ifun++) _f[iSample][ifun] = 0.0;

    _g[iSample] = new double[_num_itr];
    for (size_t iitr = 0; iitr < _num_itr; iitr++) _g[iSample][iitr] = 0.0;
  }

  // -----------------
  // INITIATE Ensemble
  // -----------------
  _xens = new double*[_ensemble_size];
  _fsur = new double*[_ensemble_size];
  _fens = new double[_ensemble_size];
  _gens = new size_t[_ensemble_size];

  for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
  {
    _xens[iSample] = new double[_num_dim];
    for (size_t idim = 0; idim < _num_dim; idim++) _xens[iSample][idim] = 0.0;

    _fsur[iSample] = new double[_num_fun];
    for (size_t ifun = 0; ifun < _num_fun; ifun++) _fsur[iSample][ifun] = 0.0;

    _fens[iSample] = 0.0;
    _gens[iSample] = 0;
  }
}

EnsembleVPS::~EnsembleVPS()
{
  delete[] _xmin;
  delete[] _xmax;

  for (size_t iSample = 0; iSample < _num_samples; iSample++) delete[] _x[iSample];
  delete[] _x;

  for (size_t iSample = 0; iSample < _num_samples; iSample++) delete[] _f[iSample];
  delete[] _f;

  for (size_t iSample = 0; iSample < _num_samples; iSample++) delete[] _g[iSample];
  delete[] _g;

  for (size_t iSample = 0; iSample < _ensemble_size; iSample++) delete[] _xens[iSample];
  delete[] _xens;

  for (size_t iSample = 0; iSample < _ensemble_size; iSample++) delete[] _fsur[iSample];
  delete[] _fsur;

  delete [] _fens;
  delete [] _gens;
}

///////////////////////////////////////////////////////////////////////

int EnsembleVPS::generate_MC_points(size_t num_dim, size_t Npoints, double** &x) {
  double* dart = new double[num_dim];
  size_t num_points = 0;
  while (num_points < Npoints)
  {
    for (size_t idim = 0; idim < num_dim; idim++)
    {
      double u = generate_a_random_number();
      dart[idim] = _xmin[idim] + u * (_xmax[idim] - _xmin[idim]);
    }
    x[num_points] = dart;
    dart = new double[num_dim];
    num_points++;
  }
  delete[] dart;
  return 0;
}

///////////////////////////////////////////////////////////////////////

int EnsembleVPS::generate_WS_points(size_t num_dim, size_t Npoints, double** &x, size_t maxMisses, double rho) {
  double r = _diag;
  size_t num_points = 0;

  double* dart = new double[num_dim];
  for (size_t idim = 0; idim < num_dim; idim++)
    dart[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);

  x[0] = dart;
  num_points++;

  // ------------------------------
  // Loop until you collect Npoints
  while (num_points < Npoints)
  {
    bool point_placed(false);
    size_t numMisses = 0;

    // ------------------------------
    // Keep trying to find one as long as you don't exceed maxMisses
    while (numMisses < maxMisses)
    {
      for (size_t idim = 0; idim < num_dim; idim++)
        dart[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);

      // ------------------------------
      // check against all existing spheres
      bool dart_in_sphere(false);

      for (size_t ip = 0; ip < num_points; ip++)
      {
        double dstsq = 0;
        for (size_t idim = 0; idim < num_dim; idim++)
          dstsq += (dart[idim] - x[ip][idim])*(dart[idim] - x[ip][idim]);

        if ((dstsq - (r*r)) < 1E-10)
        {
          dart_in_sphere = true;
          break;
        }
      }
      // ------------------------------
      if (dart_in_sphere == false)
      {
        point_placed = true;
        x[num_points] = dart;
        dart = new double[num_dim];
        num_points++;
        break;
      }
      else
      {
        dart = new double[num_dim];
        numMisses++;
      }
      // ------------------------------
    }

    if (point_placed == false)
    {
      r *= rho;
    }
    // ------------------------------
  }
  delete[] dart;
  return 0;
}

///////////////////////////////////////////////////////////////////////

int EnsembleVPS::pick_MD_points(size_t numIn, double** x, double* g, size_t numOut, double** &y) {
  std::vector<double> vg(g, g + numIn);
  std::vector<int> index(vg.size(), 0);
  for (size_t i = 0; i < index.size(); i++) index[i] = i;
  //std::cout << "Unsorted: " << std::endl;
  //for (size_t i = 0; i < vg.size(); ++i) std::cout << index[i] << ": " << vg[i] << std::endl;
  //std::cout << std::endl;
  sort(index.begin(), index.end(), [&](const int& a, const int& b) {return (vg[a] < vg[b]); });
  //std::cout << "Sorted: " << std::endl;
  //for (size_t i = 0; i < vg.size(); ++i) std::cout << index[i] << ": " << vg[index[i]] << std::endl;

  double* gwindow = new double[numOut];
  for (size_t ig = 0; ig < numOut; ig++)
    gwindow[ig] = g[index[ig]];

  size_t i_min = 0;
  double v_min = g_dispersion(gwindow, numOut);

  //std::cout << "Window " << i_min << ", variance = " << v_min << std::endl;

  // Loop to find window with min variance
  double v_window;
  for (size_t iw = 1; iw < (numIn - numOut + 1); iw++)
  {
    gwindow = new double[numOut];
    for (size_t ig = 0; ig < numOut; ig++)
      gwindow[ig] = g[index[iw + ig]];

    v_window = g_dispersion(gwindow, numOut);

    //std::cout << "Window " << iw << ", variance = " << v_window << std::endl;

    if (v_window < v_min)
    {
      //std::cout << "Update .... " << std::endl;
      i_min = iw;
      v_min = v_window;
    }
  }

  //std::cout << "Min Variance = " << v_min << std::endl;
  // clear memory
  delete[] gwindow;

  //std::cout << "window_minV: " << i_min << std::endl;
  for (size_t iy = 0; iy < numOut; iy++)
  {
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      y[iy][idim] = x[index[i_min + iy]][idim];
    }
    //std::cout << index[i_min + iy] << std::endl;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////

double EnsembleVPS::g_dispersion(double *x, size_t num_x) {
  double xmean(0.0), xstdv(0.0);
  for (size_t i = 0; i < num_x; i++) xmean += (x[i] / num_x);
  for (size_t i = 0; i < num_x; i++) xstdv += ((x[i] - xmean)*(x[i] - xmean) / (num_x - 1));
  return sqrt(xstdv) / xmean;
}

///////////////////////////////////////////////////////////////////////

void EnsembleVPS::print_data() {
  // uses data stored in _x[][] and _f[][]
  for (size_t iSample = 0; iSample < _num_samples; iSample++)
  {
    std::cout << "<P " << iSample << ">: ";
    for (size_t idim = 0; idim < _num_dim; idim++)
      std::cout << _x[iSample][idim] << " ";
    std::cout << ", f = " << _f[iSample][0];
    std::cout << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////

void EnsembleVPS::save_samples_and_evals() {
  // uses data stored in _x[][], _f[][], and _g[]
  std::ofstream txt_samples;
  txt_samples.open("x_f_g.txt");

  for (size_t iSample = 0; iSample < _num_samples; iSample++)
  {
    for (size_t idim = 0; idim < _num_dim; idim++)
      txt_samples << _x[iSample][idim] << " ";
    txt_samples << _f[iSample][0] << " ";
    txt_samples << _g[iSample][0];
    txt_samples << std::endl;
  }
  txt_samples.close();
}

///////////////////////////////////////////////////////////////////////

void EnsembleVPS::save_ei_sur() {
  // uses data stored in _ei_sur[]
  std::ofstream txt_error;
  txt_error.open("ei_sur.txt");

  for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++)
    txt_error << _ensemble_size*(iEnsemble + 1) << " " << _ei_sur[iEnsemble] << std::endl;

  txt_error.close();
}

///////////////////////////////////////////////////////////////////////

void EnsembleVPS::save_di() {
  // uses data stored in _di[]
  std::ofstream txt_error;
  txt_error.open("di.txt");

  for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++)
    txt_error << _ensemble_size*(iEnsemble + 1) << " " << _di[iEnsemble] << std::endl;

  txt_error.close();
}

///////////////////////////////////////////////////////////////////////

void EnsembleVPS::save_Ri() {
  // uses data stored in _Ri[]
  std::ofstream txt_error;
  txt_error.open("Ri.txt");

  for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++)
    txt_error << _ensemble_size*(iEnsemble + 1) << " " << _Ri[iEnsemble] << std::endl;

  txt_error.close();
}

///////////////////////////////////////////////////////////////////////

void EnsembleVPS::initiate_random_number_generator(unsigned long x) {
  //assert(sizeof (double) >= 54) ;

  cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
  size_t i;
  size_t qlen = indx = sizeof Q / sizeof Q[0];
  for (i = 0; i < qlen; i++) Q[i] = 0;

  //double c = 0.0; zc = 0.0;     /* current CSWB and SWB `borrow` */
  c = 0.0; // ETP
  zc = 0.0; // ETP
  zx = 5212886298506819.0 / 9007199254740992.0;   /* SWB seed1 */
  zy = 2020898595989513.0 / 9007199254740992.0;   /* SWB seed2 */

  size_t j;
  double s, t;     /* Choose 32 bits for x, 32 for y */
  if (x == 0) x = 123456789; /* default seeds */
  unsigned long y = 362436069; /* default seeds */

  /* Next, seed each Q[i], one bit at a time, */
  for (i = 0; i < qlen; i++)
  { /* using 9th bit from Cong+Xorshift */
    s = 0.0;
    t = 1.0;
    for (j = 0; j < 52; j++)
    {
      t = 0.5 * t; /* make t=.5/2^j */
      x = 69069 * x + 123;
      y ^= (y << 13);
      y ^= (y >> 17);
      y ^= (y << 5);
      if (((x + y) >> 23) & 1) s = s + t; /* change bit of s, maybe */
    }        /* end j loop */
    Q[i] = s;
  } /* end i seed loop, Now generate 10^9 dUNI's: */
}

///////////////////////////////////////////////////////////////////////

double EnsembleVPS::generate_a_random_number() {
  /* Takes 14 nanosecs, Intel Q6600,2.40GHz */
  int i, j;
  double t; /* t: first temp, then next CSWB value */
  /* First get zy as next lag-2 SWB */
  t = zx - zy - zc;
  zx = zy;
  if (t < 0)
  {
    zy = t + 1.0;
    zc = cc;
  }
  else
  {
    zy = t;
    zc = 0.0;
  }
  /* Then get t as the next lag-1220 CSWB value */
  if (indx < 1220)
    t = Q[indx++];
  else
  { /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
    for (i = 0; i < 1220; i++)
    {
      j = (i < 30) ? i + 1190 : i - 30;
      t = Q[j] - Q[i] + c; /* Get next CSWB element */
      if (t > 0)
      {
        t = t - cc;
        c = cc;
      }
      else
      {
        t = t - cc + 1.0;
        c = 0.0;
      }
      Q[i] = t;
    }        /* end i loop */
    indx = 1;
    t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
  } /* end else segment; return t-zy mod 1 */

  return ((t < zy) ? 1.0 + (t - zy) : t - zy);
}

///////////////////////////////////////////////////////////////////////

double EnsembleVPS::f_test(double* x)
{
  if (_eval_test_function == SmoothHerbie)
  {
    double fval = 1.0;
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double xm = x[idim] - 1.0;
      double xp = x[idim] + 1.0;
      double wherb = exp(-xm * xm) + exp(-0.8 * xp * xp);
      fval *= wherb;
    }
    fval = -fval;
    return fval;
  }
  else if (_eval_test_function == Herbie)
  {
    double fval = 1.0;
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double xm = x[idim] - 1.0;
      double xp = x[idim] + 1.0;
      double wherb = exp(-xm * xm) + exp(-0.8 * xp * xp) - 0.05 * sin(8 * (x[idim] + 0.1));
      fval *= wherb;
    }
    fval = -fval;
    return fval;
  }
  else if (_eval_test_function == Cone)
  {
    double fval = 0.0;
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double xm = x[idim];
      fval += xm * xm;
    }
    fval = sqrt(fval);
    return fval;
  }
  else if (_eval_test_function == Cross)
  {
    double fval = 1.0;
    const double pi = 3.14159265358979324;
    double dpow = 1.0 / _num_dim;
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      fval *= 0.5 * (1.0 + std::cos(2.0 * pi * x[idim]));
    }
    fval = std::pow(fval, dpow);
    return fval;
  }
  else if (_eval_test_function == Trench)
  {
    // step function in sphere
    double h = 0.0;
    double fval = 0.0;
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = x[idim];
      h += dx * dx;
    }
    h = std::sqrt(h);
    if ((h > 1.0) || (h < 0.5))
    {
      fval = 1.0;
    }
    return fval;
  }
  return 0.0;
}

///////////////////////////////////////////////////////////////////////

double EnsembleVPS::g_test(double* x)
{
  double* a = new double[_num_dim];
  double* u = new double[_num_dim];
  double xterm;
  double fval = 1.0;
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    a[idim] = 5.0 / (_xmax[idim] - _xmin[idim]);
    u[idim] = 0.5 * (_xmax[idim] - _xmin[idim]) + _xmin[idim];
    xterm = (a[idim] * a[idim]) * ((x[idim] - u[idim])*(x[idim] - u[idim]));
    fval += xterm;
  }
  //fval = exp((-1) * fval);
  fval *= 100.0;
  return fval;
}
