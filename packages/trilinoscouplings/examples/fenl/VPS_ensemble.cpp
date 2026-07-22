// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  _num_fun = 2; // DON'T CHANGE
  _num_samples = num_samples;      // Sample budget
  _ensemble_size = ensemble_size;     // Ensemble Size

  // The chosen ensemble size may not evenly divide the number of samples
  _num_ensembles = (_num_samples+_ensemble_size-1) / _ensemble_size;
  _num_samples = _ensemble_size * _num_ensembles;

  _surrogate_order = 2;
  _proc_rank = rank;

  // ----------------------------
  // SET DOMAIN BOUNDS and VOLUME
  // ----------------------------
  _xmin = new double[_num_dim];
  _xmax = new double[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    _xmin[idim] = -1.0;
    _xmax[idim] = 1.0;
  }

  // ----------------------------
  // INITIATE SAMPLES & FUNCTIONS
  // ----------------------------
  _x = new double*[_num_samples];
  _f = new double*[_num_samples];
  for (size_t iSample = 0; iSample < _num_samples; iSample++)
  {
    _x[iSample] = new double[_num_dim];
    _f[iSample] = new double[_num_fun];
  }

  // -----------------
  // INITIATE Ensemble
  // -----------------
  _xens = new double*[_ensemble_size];
  _fens = new double*[_ensemble_size];
  for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
  {
    _xens[iSample] = new double[_num_dim];
    _fens[iSample] = new double[_num_fun];
  }
}

EnsembleVPS::~EnsembleVPS()
{
  delete[] _xmin;
  delete[] _xmax;

  for (size_t iSample = 0; iSample < _num_samples; iSample++)
    delete[] _x[iSample];
  delete[] _x;

  for (size_t iSample = 0; iSample < _num_samples; iSample++)
    delete[] _f[iSample];
  delete[] _f;

  for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
    delete[] _xens[iSample];
  delete[] _xens;

  for (size_t iSample = 0; iSample < _ensemble_size; iSample++)
    delete[] _fens[iSample];
  delete[] _fens;
}
