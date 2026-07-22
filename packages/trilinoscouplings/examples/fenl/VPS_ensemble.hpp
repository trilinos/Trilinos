// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include "VPS.hpp"

class EnsembleVPS {
private:
  size_t _num_dim; // problem dimensions
  size_t _num_fun; // number of potential EVAL functions
  size_t _num_samples;     // Budget
  size_t _ensemble_size;   // Ensemble size
  size_t _num_ensembles;
  size_t _surrogate_order;
  int _proc_rank;

  double* _xmin;
  double* _xmax;

  double** _x;
  double** _f;

  double** _xens;
  double** _fens;

public:

  EnsembleVPS(const size_t dim,
              const size_t num_samples,
              const size_t ensemble_size,
              const int rank);

  ~EnsembleVPS();

  template< class Model >
  void run(const Model& model, double& mean, double& stddev)
  {
    double R_num = 0.0, R_den = 0.0;

    for (size_t iensemble = 0; iensemble < _num_ensembles; iensemble++)
    {
      size_t num_prior_samples = iensemble * _ensemble_size;

      if (iensemble == 0) {
        VPS vps;
        vps.get_initial_well_spaced_samples(_num_dim, _xmin, _xmax,
                                            _ensemble_size, _xens);
      }
      else {
        VPS vps;
        vps.build_surrogate(_num_dim, _xmin, _xmax, _num_fun,
                            VPS::Regression, VPS::monomials, _surrogate_order,
                            num_prior_samples, num_prior_samples, _x, _f, 0, 0);
        vps.get_ensemble(_ensemble_size, _xens, _proc_rank);
      }

      model(_ensemble_size, _num_dim, _xens, _fens);

      double sum_iter(0.0), max_iter(0.0);
      for (size_t i = 0; i < _ensemble_size; i++) {
        sum_iter += _fens[i][1];
        if (_fens[i][1] > max_iter)
          max_iter = _fens[i][1];

        for (size_t idim = 0; idim < _num_dim; idim++)
          _x[num_prior_samples + i][idim] = _xens[i][idim];
        for (size_t ifunc = 0; ifunc < _num_fun; ifunc++)
          _f[num_prior_samples + i][ifunc] = _fens[i][ifunc];
      }

      R_num += max_iter * _ensemble_size;
      R_den += sum_iter;
      double ensemble_R = max_iter * _ensemble_size / sum_iter;
      double R = R_num / R_den;
      if (_proc_rank == 0)
        std::cout << "Ensemble Actual R = " << ensemble_R
                  << ", Ensembles R = " << R
                  << std::endl << std::endl;
    }

    // Compute surrogate mean, standard deviation
    VPS vps;
    vps.build_surrogate(_num_dim, _xmin, _xmax, _num_fun,
                        VPS::Regression, VPS::monomials, _surrogate_order,
                        _num_samples, _num_samples, _x, _f, 0, 0);
    vps.get_stats(1000000, 0, mean, stddev);
    stddev = std::sqrt(stddev);

    double R = R_num / R_den;
    if (_proc_rank == 0)
      std::cout << std::endl << "R_total = " << R << std::endl;
  }

};
