#include <iostream>
#include "VPS.hpp"

class EnsembleVPS {
private:
  const double PI = 3.141592653589793238462643383279502;

  size_t _num_dim; // problem dimensions
  double** _x;     // data points
  double** _xens;  // One ensemble
  double* _fens;
  size_t* _gens;
  double* _xmin;   //
  double* _xmax;   //
  double  _Dvol;   // domain volume
  double _diag;    // domain diagonal

  size_t _num_fun; // number of potential EVAL functions
  double** _f;     // evaluations
  double** _fsur;  // One ensemble
  double* _fmin;   //
  double* _fmax;   //

  size_t _num_itr; // number of potential ITER functions
  double** _g;     // iterations
  double* _gmin;   //
  double* _gmax;   //

  size_t _num_samples;     // Budget
  size_t _ensemble_size;   // Ensemble size
  size_t _sampling;
  size_t _surrogate_order;

  size_t _MCratio;        // MC candidates ratio

  size_t _maxMisses;       // Well-spacedness parameter
  double _shrinkage_rho;   // Well-spacedness parameter

  double* _ei_sur;          // Interpolation error
  double* _di;             // iterations dispersion
  double* _Ri;             // Grouping efficiency
  double _Rt;              // Total grouping efficiency

  double Q[1220];          // variables for Random number generator
  int indx;
  double cc, c, zc, zx, zy;
  size_t qlen;

  int _proc_rank;

  // Test functions
  enum testfunction { SmoothHerbie, Herbie, Cone, Cross, Trench };
  testfunction _eval_test_function;

public:

  EnsembleVPS(const size_t dim,
              const size_t num_samples,
              const size_t ensemble_size,
              const int rank);

  ~EnsembleVPS();

  template< class Model >
  void run(const Model& model)
  {
    // **************
    // INITIATE ERROR
    // **************
    _ei_sur = new double[_num_samples / _ensemble_size];
    for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++) _ei_sur[iEnsemble] = 0.0;

    _di = new double[_num_samples / _ensemble_size];
    for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++) _di[iEnsemble] = 0.0;

    _Ri = new double[_num_samples / _ensemble_size];
    for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++) _Ri[iEnsemble] = 0.0;

    // ********************************* //
    _Rt = 0.0;
    // ********************************* //

    double Rt_num = 0.0;
    double Rt_den = 0.0;

    // **********************
    // ENSEMBLE SAMPLING
    // **********************
    size_t iEnsemble = 0;
    while (iEnsemble < (_num_samples / _ensemble_size))
    {
      // --------------------
      // FIRST ENSEMBLE
      // --------------------
      if (iEnsemble == 0)
      {
        // populates _xens
        generate_WS_points(_num_dim, _ensemble_size, _xens, _maxMisses, _shrinkage_rho);
      }

      // Evaluate response and iterations at each sample point
      model(_ensemble_size, _num_dim, _xens, _fens, _gens);
      // for (size_t iSample = 0; iSample < _ensemble_size; iSample++) {
      //   _fens[iSample] = f_test(_xens[iSample]);
      //   _gens[iSample] = g_test(_xens[iSample]);
      // }

      // -----------------------------------------------
      // Eval previous surrogate before updating points
      // -----------------------------------------------
      double esurr(0.0);
      for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
        esurr += (_fens[iSample] - _fsur[iSample][0]) * (_fens[iSample] - _fsur[iSample][0]);

      _ei_sur[iEnsemble] += (esurr / _ensemble_size);

      // ----------------------
      // Update _x, _f, and _g
      // ----------------------
      for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
      {
        for (size_t idim = 0; idim < _num_dim; idim++)
          _x[iSample + (iEnsemble*_ensemble_size)][idim] = _xens[iSample][idim];

        //
        _f[iSample + (iEnsemble*_ensemble_size)][0] = _fens[iSample];
        _g[iSample + (iEnsemble*_ensemble_size)][0] = _gens[iSample];
      }

      // ----------------------------------------
      // VPS surrogate(s) construction
      // ----------------------------------------
      VPS f_surrogate;
      VPS g_surrogate;

      f_surrogate.build_surrogate(_num_dim, _xmin, _xmax, 1, VPS::Regression, VPS::monomials, _surrogate_order, (iEnsemble + 1)*_ensemble_size, (iEnsemble + 1)*_ensemble_size, _x, _f, 0, 0);
      g_surrogate.build_surrogate(_num_dim, _xmin, _xmax, 1, VPS::Regression, VPS::monomials, 0, (iEnsemble + 1)*_ensemble_size, (iEnsemble + 1)*_ensemble_size, _x, _g, 0, 0);

      // -------------------------
      // ADAPTATION: Start collecting Candidates
      // -------------------------

      size_t _numMC = _MCratio * _ensemble_size;
      double** yMC = new double*[_numMC];
      double*  yMC_ErrorEst = new double[_numMC];

      double r = _diag;
      double* dart = new double[_num_dim];
      double* p_facet = new double[_num_dim];
      double p_errest;

      // --------------------------
      // Find maximum potential error
      // --------------------------
      double facetMaxErrEst = 0.0;
      size_t iMCp(0);
      while (iMCp < 10000)
      {
        for (size_t idim = 0; idim < _num_dim; idim++)
          dart[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);
        f_surrogate.suggest_new_sample(dart, p_facet, p_errest);
        if (p_errest > facetMaxErrEst) facetMaxErrEst = p_errest;
        iMCp++;
      }
      // --------------------------
      // Now loop to collect points
      // --------------------------
      size_t num_points = 0;
      while (num_points < _numMC)
      {
        size_t numMisses = 0;

        while (numMisses < _maxMisses)
        {
          for (size_t idim = 0; idim < _num_dim; idim++)
            dart[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);

          f_surrogate.suggest_new_sample(dart, p_facet, p_errest);
          // ------------------
          bool miss(false);
          bool conflict(false);
          // ------------------
          for (size_t ip = 0; ip < num_points; ip++)
          {
            double dstsq = 0;
            for (size_t idim = 0; idim < _num_dim; idim++)
              dstsq += (p_facet[idim] - yMC[ip][idim])*(p_facet[idim] - yMC[ip][idim]);

            if ((dstsq - (r * r)) < 1E-10)
            {
              if (p_errest > yMC_ErrorEst[ip])
                conflict = true; // Conflict = point in a sphere and has a higher error
              else
              {
                miss = true;    // Miss = point in a sphere, but its error is lower
                break;
              }
            }
          }
          // ------------------
          if ((miss) || (p_errest < 0.1 * facetMaxErrEst))
          {
            numMisses++;
            continue;    // don't add points, just skip and try a new dart
          }
          // --------------------------
          // Add point
          yMC[num_points] = p_facet;
          yMC_ErrorEst[num_points] = p_errest;
          num_points++;
          numMisses = 0;
          // --------------------------
          if (conflict)
          {
            // point has higher error
            for (size_t ip = 0; ip < num_points; ip++)
            {
              double dstsq = 0.0;
              for (size_t idim = 0; idim < _num_dim; idim++)
                dstsq += (p_facet[idim] - yMC[ip][idim])*(p_facet[idim] - yMC[ip][idim]);

              if ((dstsq - (r * r)) > 1E-10)
                continue;

              yMC[ip] = yMC[num_points - 1];
              yMC_ErrorEst[ip] = yMC_ErrorEst[num_points - 1];
              num_points--;
            }
          }
          dart = new double[_num_dim];
          p_facet = new double[_num_dim];
          break;
        }

        // -----------------------------------------
        // if you don't have enough points yet, then you need to open up space

        if ((numMisses == _maxMisses) && (num_points < _numMC))
        {
          r *= _shrinkage_rho;
        }
      }

      delete[] dart;
      delete[] p_facet;
      // ----------------------

      // -----------------------------
      //pick minimum dispersion points
      // -----------------------------
      double* gsMC = new double[1];
      double* g_estimMC = new double[_numMC];

      size_t ic = 0;
      while (ic < _numMC)
      {
        // Evaluate surrogate of iterations at MC points
        g_surrogate.evaluate_surrogate(yMC[ic], gsMC);
        g_estimMC[ic] = gsMC[0];
        //g_estimMC[ic] = g_test(yMC[ic]);
        ic++;
      }

      double** yMD = new double*[_ensemble_size];    // MD points

      for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
      {
        yMD[iSample] = new double[_num_dim];
        for (size_t idim = 0; idim < _num_dim; idim++)
          yMD[iSample][idim] = 0.0;
      }
      //populates yMD
      pick_MD_points(_numMC, yMC, g_estimMC, _ensemble_size, yMD);


      // ----------------------
      // ----------------------
      // ----------------------
      // PASS ENSEMBLE BACK
      // ----------------------
      // ----------------------
      // ----------------------
      double* fs = new double[1];
      for (size_t ix = 0; ix < _ensemble_size; ix++)
      {
        for (size_t idim = 0; idim < _num_dim; idim++)
          _xens[ix][idim] = yMD[ix][idim];

        f_surrogate.evaluate_surrogate(_xens[ix], fs);
        _fsur[ix][0] = fs[0];
      }
      delete[] fs;

      // ----------------

      // clear MC/MD memory
      for (size_t iSample = 0; iSample < _numMC; iSample++)
        delete[] yMC[iSample];
      delete[] yMC;
      delete[] yMC_ErrorEst;
      for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
        delete[] yMD[iSample];
      delete[] yMD;

      // ----------------------------------------
      // Quantify points' iterations dispersion
      // ----------------------------------------
      double* g_estim = new double[_ensemble_size];
      double* gs = new double[1];

      double max_g_estim = _gens[0];
      double mean_g_estim = 0.0;

      for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
      {
        g_estim[iSample] = _gens[iSample];

        // find max and mean
        if (g_estim[iSample] > max_g_estim)
          max_g_estim = g_estim[iSample];
        mean_g_estim += (g_estim[iSample] / _ensemble_size);
      }

      _di[iEnsemble] += g_dispersion(g_estim, _ensemble_size);
      _Ri[iEnsemble] += max_g_estim / mean_g_estim;

      Rt_num += max_g_estim;
      Rt_den += mean_g_estim;

      delete[] g_estim;
      delete[] gs;

      // ----------------------------------------
      iEnsemble++;
    }

    _Rt = Rt_num / Rt_den;

    // -----------------------------------------------------
    // Average interpolation error over experiments and save
    // -----------------------------------------------------
    if (_proc_rank == 0) {
      std::cout << "ei_surrogate" << std::endl;
      std::cout << "=============" << std::endl;
      for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++)
      {
        std::cout << _ei_sur[iEnsemble] << std::endl;
      }
      std::cout << std::endl;
      save_ei_sur();
    }

    // -----------------------------------------------------
    // Average d over experiments and save
    // -----------------------------------------------------
    if (_proc_rank == 0) {
      std::cout << "d_ensemble" << std::endl;
      std::cout << "==========" << std::endl;
      for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++)
      {
        std::cout << _di[iEnsemble] << std::endl;
      }
      std::cout << std::endl;
      save_di();
    }

    // -----------------------------------------------------------------------------------------
    // Compute R_{ensemble} = avg_{experiments} [[ (max)_{ensemble} / (avg)_{ensemble} ]]
    // -----------------------------------------------------------------------------------------
    if (_proc_rank == 0) {
      std::cout << "R_ensemble" << std::endl;
      std::cout << "==========" << std::endl;
      for (size_t iEnsemble = 0; iEnsemble < (_num_samples / _ensemble_size); iEnsemble++)
      {
        std::cout << _Ri[iEnsemble] << std::endl;
      }
      std::cout << std::endl;
      save_Ri();
    }

    // -----------------------------------------------------------------------------------------
    // Compute R_total = avg_{experiments} [[ sum_{ensembles} (max) / sum_{ensembles} (avg) ]]
    // -----------------------------------------------------------------------------------------
    double R_total = _Rt;
    if (_proc_rank == 0) {
      std::cout << "R_total" << std::endl;
      std::cout << "=======" << std::endl;
      std::cout << R_total << std::endl;
    }

    // ---------------------
    // Clear metrics memory
    // ---------------------
    delete _ei_sur;
    delete _di;
    delete _Ri;
  }

  int generate_MC_points(size_t num_dim, size_t Npoints, double** &x);

  int generate_WS_points(size_t num_dim, size_t Npoints, double** &x, size_t maxMisses, double rho);

  int pick_MD_points(size_t numIn, double** x, double* g, size_t numOut, double** &y);

  double g_dispersion(double *x, size_t num_x);

  void print_data();

  void save_samples_and_evals();

  void save_ei_sur();

  void save_di();

  void save_Ri();

  void initiate_random_number_generator(unsigned long x);

  double generate_a_random_number();

  double f_test(double* x);

  double g_test(double* x);

};
