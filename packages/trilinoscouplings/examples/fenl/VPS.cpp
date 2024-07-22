// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 2016/09/23
// Automatic Discontinuity Detection using Spline Extrapolation

#include "VPS.hpp"

VPS::VPS()
{
  /////////// kd-tree  Variables /////////////////////////////////////////////////
  _tree_left = 0; _tree_right = 0;

  /////////// Delaunay Graph  Variables //////////////////////////////////////////
  _seed_box = 0; _seed_rf = 0; _seed_rc = 0; _seed_neighbors = 0;

  _seed_disc_neighbors = 0; // ETP

  /////////// Surrogate  Variables //////////////////////////////////////////
  _p = 0; _basis_coef = 0; _basis_index = 0;

  initiate_random_number_generator(1234567890);

  size_t seed(size_t(time(0)));
  seed = 1234567890;
  initiate_random_number_generator(seed);
  //std::cout << "RNG seed = " << seed << std::endl;
  _num_vs = 0;
}

VPS::~VPS()
{
  clear_memory();
}

int VPS::clear_memory()
{
  /////////// kd-tree  Variables /////////////////////////////////////////////////
  _tree_max_height = 0;
  if (_tree_left != 0){ delete[] _tree_left; _tree_left = 0; }
  if (_tree_right != 0){ delete[] _tree_right; _tree_right = 0; }

  /////////// Delaunay Graph  Variables //////////////////////////////////////////
  if (_seed_box != 0)
  {
    for (size_t i = 0; i < _num_samples; i++) delete[] _seed_box[i];
    delete[] _seed_box;
    _seed_box = 0;
  }
  if (_seed_rf != 0){ delete[] _seed_rf; _seed_rf = 0; }
  if (_seed_rc != 0){ delete[] _seed_rc; _seed_rc = 0; }

  if (_seed_neighbors != 0)
  {
    for (size_t i = 0; i < _num_samples; i++) delete[] _seed_neighbors[i];
    delete[] _seed_neighbors;
    _seed_neighbors = 0;
  }
  if (_seed_disc_neighbors != 0)
  {
    for (size_t iseed = 0; iseed < _num_samples; iseed++)
    {
      if (_seed_disc_neighbors[iseed] != 0)
      {
        for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
        {
          if (_seed_disc_neighbors[iseed][ifunc] != 0) delete[] _seed_disc_neighbors[iseed][ifunc];
        }
        delete[] _seed_disc_neighbors[iseed];
      }

    }
    delete[] _seed_disc_neighbors;
    _seed_disc_neighbors = 0;
  }

  /////////// Surrogate  Variables //////////////////////////////////////////
  if (_p != 0)
  {
    for (size_t ibasis = 0; ibasis < _num_basis_pool; ibasis++) delete[] _p[ibasis];
    delete[] _p;
    _p = 0;
  }

  if (_basis_coef != 0)
  {
    for (size_t i = 0; i < _num_samples; i++)
    {
      for (size_t j = 0; j < _num_functions; j++) delete[] _basis_coef[i][j];
      delete[] _basis_coef[i];
    }
    delete[] _basis_coef;
    _basis_coef = 0;
  }

  if (_basis_index != 0)
  {
    for (size_t i = 0; i < _num_samples; i++)
    {
      for (size_t j = 0; j < _num_functions; j++) delete[] _basis_index[i][j];
      delete[] _basis_index[i];
    }
    delete[] _basis_index;
    _basis_index = 0;
  }

  return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// VPS Public Methods               //////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int VPS::get_initial_well_spaced_samples(size_t num_dim, double* xmin, double* xmax, size_t num_samples, double** x)
{
  double r_sq(0.0);
  for (size_t idim = 0; idim < num_dim; idim++)
  {
    double dx = xmax[idim] - xmin[idim];
    r_sq += dx * dx;
  }

  size_t isample(0);
  size_t num_successive_misses(0), max_num_successive_misses(100);
  while (isample < num_samples)
  {
    for (size_t idim = 0; idim < num_dim; idim++) x[isample][idim] = xmin[idim] + generate_a_random_number() * (xmax[idim] - xmin[idim]);
    bool miss(false);
    for (size_t i = 0; i < isample; i++)
    {
      double dst_sq(0.0);
      for (size_t idim = 0; idim < num_dim; idim++)
      {
        double dx = x[i][idim] - x[isample][idim];
        dst_sq += dx * dx;
      }
      if (dst_sq < r_sq)
      {
        miss = true; break;
      }
    }
    if (miss) num_successive_misses++;
    else
    {
      num_successive_misses = 0;
      isample++;
    }

    if (num_successive_misses == max_num_successive_misses)
    {
      num_successive_misses = 0;
      r_sq *= 0.81;
    }
  }
  return 0;
}

int VPS::build_surrogate(size_t num_dim,
                         double* xmin, double* xmax,
                         size_t num_functions,
                         surrogate_method method,
                         surrogate_basis basis, size_t desired_order,
                         size_t budget, size_t num_samples,
                         double** x,
                         double** f,
                         double*** g,
                         double**** h
  )
{
  clear_memory();

  _num_dim = num_dim; _xmin = xmin; _xmax = xmax;
  _num_functions = num_functions; _method = method; _basis = basis;
  _desired_order = desired_order; _budget = budget; _num_samples = num_samples;
  _x = x; _f = f; _g = g; _h = h;

  find_domain_diagonal();

  // 1. Build kd-tree
  kd_tree_build_balanced();

  // 2. Construct Delaunay Graph
  if (desired_order > 0) construct_delaunay_graph();

  //plot_delaunay_graph("graph_nodisc.ps");

  // 3. Detect Discontinuites
  //detect_discontinuities();

  //plot_delaunay_graph("graph_disc.ps");

  // 3. Construct local surrogates
  construct_local_surrogates();


  std::string name_1;
  std::string name_2;

  if (true)
  {
    std::stringstream sstm;
    sstm << "vps_qoi_" << _num_samples << ".ps";
    name_1 = sstm.str();
  }
  if (true)
  {
    std::stringstream sstm;
    sstm << "vps_nit_" << _num_samples << ".ps";
    name_2 = sstm.str();
  }

  if (_num_samples > 1190 && _num_samples < 2100)
  {
    //plot_vps_surrogate("vps_2d_qoi.ps", 0, 21, true);
    //plot_vps_surrogate("vps_2d_nit.ps", 1, 21, true);
    // plot_vps_frames(name_1, 0, 11, 11, 21, false);
    // plot_vps_frames(name_2, 1, 11, 11, 21, false);
  }

  //plot_vps_frames("vps_surrogate_r.ps", 0, 11, 11, 21, false);
  //plot_vps_frames("vps_surrogate_g.ps", 1, 11, 11, 21, false);
  //plot_vps_frames("vps_surrogate_b.ps", 2, 5, 5, 21, false);
  return 0;
}

int VPS::evaluate_surrogate(double* x, double* fs)
{
  double closest_distance = _diag; size_t closest_seed;
  get_closest_seed_tree(x, closest_seed, closest_distance);
  evaluate_surrogate(closest_seed, x, fs);
  return 0;
}

int VPS::evaluate_surrogate(size_t cell_index, double* x, double* fs)
{
  for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
  {
    fs[ifunc] = 0.0;

    if (_desired_order == 0)
    {
      fs[ifunc] = _f[cell_index][ifunc];
    }
    else
    {
      for (size_t ibasis = 0; ibasis < _num_basis; ibasis++)
      {
        size_t basis_index = _basis_index[cell_index][ifunc][ibasis];
        fs[ifunc] += _basis_coef[cell_index][ifunc][ibasis] * evaluate_basis_function(x, cell_index, basis_index);
      }
    }
  }
  return 0;
}

int VPS::suggest_new_sample(double* x, double &r, double &err_est)
{
  double* xmc = new double[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++) xmc[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);

  size_t closest_seed; double dst(DBL_MAX);
  get_closest_seed_tree(xmc, closest_seed, dst);

  if (_num_samples < 5 * _num_dim)
  {
    // find distance to closest boundary
    double h_dim; double h_shortest(DBL_MAX);

    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      h_dim = fabs(xmc[idim] - _xmin[idim]);
      if (h_dim < h_shortest) h_shortest = h_dim;

      h_dim = fabs(xmc[idim] - _xmax[idim]);
      if (h_dim < h_shortest) h_shortest = h_dim;
    }

    // error estimate = min(distance to closest seed, distance to closest boundary)
    r = dst;
    if (h_shortest < dst) dst = h_shortest;

    for (size_t idim = 0; idim < _num_dim; idim++) x[idim] = xmc[idim];
    err_est = dst;
  }
  else
  {
    size_t closest_seed;
    double dst(DBL_MAX);
    get_closest_seed_tree(xmc, closest_seed, dst);
    sample_voronoi_vertex(closest_seed, _xmin, _xmax, _diag, x);

    r = 0.0;
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = _x[closest_seed][idim] - x[idim];
      r += dx * dx;
    }
    r = sqrt(r);

    size_t capacity = 10;
    size_t* neighbor_seeds = new size_t[capacity];
    size_t num_neighbor_seeds = 0;
    get_seeds_in_sphere_tree(x, r + 1E-10, num_neighbor_seeds, neighbor_seeds, capacity);

    if (num_neighbor_seeds == 1)
    {
      // switch to neighbors of closest seed
      num_neighbor_seeds = _seed_neighbors[closest_seed][1] + 1;
      neighbor_seeds[0] = closest_seed;
      for (size_t j = 1; j < num_neighbor_seeds; j++)
      {
        neighbor_seeds[j] = _seed_neighbors[closest_seed][2 + j];
      }
    }

    double* fneighbor = new double[num_neighbor_seeds];
    double* fs = new double[_num_functions];
    for (size_t i = 0; i < num_neighbor_seeds; i++)
    {
      size_t seed = neighbor_seeds[i];
      evaluate_surrogate(seed, x, fs);
      fneighbor[i] = fs[0];
    }

    err_est = 0.0;
    for (size_t i = 0; i < num_neighbor_seeds; i++)
    {
      for (size_t j = i + 1; j < num_neighbor_seeds; j++)
      {
        double df = fabs(fneighbor[i] - fneighbor[j]);
        if (df > err_est) err_est = df;
      }
    }
    err_est *= pow(r, _num_dim);

    delete[] fneighbor;
    delete[] fs;
    delete[] neighbor_seeds;
  }

  delete[]xmc;
  return 0;
}

int VPS::get_ensemble(size_t num_ensemble_points, double** ensemble_points,
                      int proc_rank)
{
  double* ensemble_err = new double[num_ensemble_points];
  double* ensemble_nit = new double[num_ensemble_points];
  double* ensemble_r = new double[num_ensemble_points];

  // collect points based on WS, Error Estimate, and min variance in nit:
  double* x = new double[_num_dim]; double r, err_est;

  double* fs = new double[_num_functions];

  size_t num_points(0);
  size_t num_misses(0), max_num_misses(100);

  double nit_av(DBL_MAX), nit_max(0.0), err_av(0.0), err_max(0.0);

  bool improve_error_only(true);
  while (num_misses < max_num_misses)
  {
    suggest_new_sample(x, r, err_est);

    // make sure this point is outside of ensemble spheres
    bool conflict(false);
    for (size_t ipoint = 0; ipoint < num_points; ipoint++)
    {
      double dstsq(0.0);
      for (size_t idim = 0; idim < _num_dim; idim++)
      {
        double dx = x[idim] - ensemble_points[ipoint][idim];
        dstsq += dx * dx;
      }
      if (dstsq < r * r || dstsq < ensemble_r[ipoint] * ensemble_r[ipoint])
      {
        conflict = true; break;
      }
    }

    if (conflict) { num_misses++; continue; } // Violation of Wellspasedness condition

    evaluate_surrogate(x, fs);

    size_t iloc(num_ensemble_points);
    if (num_points < num_ensemble_points) iloc = num_points;

    if (iloc == num_ensemble_points)
    {
      // try replacing all prior points
      for (size_t jloc = 0; jloc < num_points; jloc++)
      {
        // attempt same/better efficiency and better accuracy
        double new_nit_av(0.0), new_nit_max(0.0), new_err_av(0.0), new_err_max(0.0);
        for (size_t ipoint = 0; ipoint < num_points; ipoint++)
        {
          double nit = ensemble_nit[ipoint]; double err = ensemble_err[ipoint];
          if (ipoint == jloc)
          {
            nit = fs[1]; err = err_est;
          }
          new_nit_av += nit;
          new_err_av += err;

          if (nit > new_nit_max) new_nit_max = nit;
          if (err > new_err_max) new_err_max = err;
        }
        new_nit_av /= num_points;
        new_err_av /= num_points;

        if (improve_error_only)
        {
          if (new_err_av > err_av)
          {
            iloc = jloc;
            break;
          }
        }
        else
        {
          double oldR = nit_max / nit_av;
          double newR = new_nit_max / new_nit_av;
          if (newR < oldR && new_err_av > 0.8 * err_av) // better effeciency at least same accuracy
          {
            iloc = jloc; break;
          }
        }
      }
    }

    if (iloc < num_ensemble_points)
    {
      // Add point to ensemble
      num_misses = 0;
      for (size_t idim = 0; idim < _num_dim; idim++) ensemble_points[iloc][idim] = x[idim];
      ensemble_nit[iloc] = fs[1];
      ensemble_err[iloc] = err_est;
      ensemble_r[iloc] = r;
      if (iloc == num_points) num_points++;

      // update error and nit bounds
      nit_av = 0.0; nit_max = 0.0;
      if (improve_error_only)
      {
        err_av = 0.0; err_max = 0.0; // update error metrics
      }

      for (size_t ipoint = 0; ipoint < num_points; ipoint++)
      {
        if (ensemble_nit[ipoint] > nit_max) nit_max = ensemble_nit[ipoint];
        nit_av += ensemble_nit[ipoint];

        if (improve_error_only)
        {
          if (ensemble_err[ipoint] > err_max) err_max = ensemble_err[ipoint];
          err_av += ensemble_err[ipoint];
        }
      }
      nit_av /= num_points;
      if (improve_error_only) err_av /= num_points;
    }
    else num_misses++;

    if (num_misses == max_num_misses && improve_error_only)
    {
      num_misses = 0; improve_error_only = false;
      max_num_misses = 1000;
    }
  }

  delete[] x; delete[] fs;

  delete[] ensemble_err; delete[] ensemble_nit; delete[] ensemble_r;

  if (proc_rank == 0)
    std::cout << "Number of samples = " << _num_samples << " , Ensemble expected R = " << nit_max / nit_av << std::endl;
  return 0;
}

int VPS::get_stats(size_t num_mc_points, size_t function_index, double& mean, double& var)
{
  double* xs = new double[_num_dim];
  double* fs = new double[_num_functions];

  mean = 0.0;
  for (size_t imc = 0; imc < num_mc_points; imc++)
  {
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      xs[idim] = _xmin[idim] + generate_a_random_number()*(_xmax[idim] - _xmin[idim]);
    }
    evaluate_surrogate(xs, fs);
    mean += fs[function_index];
  }
  mean /= num_mc_points;

  var = 0.0;
  for (size_t imc = 0; imc < num_mc_points; imc++)
  {
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      xs[idim] = _xmin[idim] + generate_a_random_number()*(_xmax[idim] - _xmin[idim]);
    }
    evaluate_surrogate(xs, fs);
    var += (fs[function_index] - mean) * (fs[function_index] - mean);
  }
  var /= num_mc_points;

  delete[] xs; delete[] fs;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RNG               /////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VPS::initiate_random_number_generator(unsigned long x)
{
  //assert(sizeof (double) >= 54) ;

  cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
  size_t i;
  size_t qlen = indx = sizeof Q / sizeof Q[0];
  for (i = 0; i < qlen; i++) Q[i] = 0;

  c = 0.0;     /* current CSWB */
  zc = 0.0;       /* current SWB `borrow` */
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

double VPS::generate_a_random_number()
{
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

int VPS::sample_uniformly_from_unit_sphere(double* dart, size_t num_dim)
{
  size_t idim = 0;
  while (true)
  {
    double u1 = generate_a_random_number();
    double u2 = generate_a_random_number();
    double r = sqrt(-2 * log(u1));
    double theta = 2 * PI * u2;
    double n1 = r * cos(theta);
    double n2 = r * sin(theta);
    dart[idim] = n1; idim++; if (idim == num_dim) break;
    dart[idim] = n2; idim++; if (idim == num_dim) break;
  }
  double norm(0.0);
  for (idim = 0; idim < num_dim; idim++) norm += dart[idim] * dart[idim];

  norm = 1.0 / sqrt(norm);
  for (idim = 0; idim < num_dim; idim++) dart[idim] *= norm;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// kd-tree  Methods  /////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int VPS::kd_tree_build_balanced()
{
//#pragma region Build Balanced kd-tree:
  //std::cout << "\n    * Build Balanced kd-tree .......... ";
  //clock_t start_time, end_time; double cpu_time;
  //start_time = clock();

  if (_tree_left == 0) _tree_left = new size_t[_budget];
  if (_tree_right == 0) _tree_right = new size_t[_budget];

  size_t* tree_nodes_sorted = new size_t[_num_samples];
  for (size_t i = 0; i < _num_samples; i++) tree_nodes_sorted[i] = i;

  for (size_t iseed = 0; iseed < _budget; iseed++)
  {
    _tree_left[iseed] = iseed; _tree_right[iseed] = iseed;
  }
  _tree_origin = _budget;

  size_t target_pos = _num_samples / 2; _tree_max_height = 0;
  kd_tree_balance_quicksort(target_pos, 0, _num_samples - 1, 0, tree_nodes_sorted);

  delete[] tree_nodes_sorted;

  //end_time = clock();
  //cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
  //std::cout << cpu_time << " seconds." << std::endl;
  return 0;
//#pragma endregion
}

int VPS::kd_tree_balance_quicksort(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
//#pragma region kd tree balance:
  kd_tree_quicksort_adjust_target_position(target_pos, left, right, active_dim, tree_nodes_sorted);

  // target position is correct .. add to tree
  if (_tree_origin == _budget)      _tree_origin = tree_nodes_sorted[target_pos];
  else                              kd_tree_add_point(tree_nodes_sorted[target_pos]);

  /* recursion */
  active_dim++;
  if (active_dim == _num_dim) active_dim = 0;

  if (target_pos > left + 1)  kd_tree_balance_quicksort((left + target_pos - 1) / 2, left, target_pos - 1, active_dim, tree_nodes_sorted);
  else if (left < target_pos) kd_tree_add_point(tree_nodes_sorted[left]);

  if (target_pos + 1 < right)  kd_tree_balance_quicksort((target_pos + 1 + right) / 2, target_pos + 1, right, active_dim, tree_nodes_sorted);
  else if (right > target_pos) kd_tree_add_point(tree_nodes_sorted[right]);

  return 0;
//#pragma endregion
}

int VPS::kd_tree_quicksort_adjust_target_position(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
//#pragma region kd tree Quick sort pivot:
  size_t i = left, j = right;

  size_t pivot_seed = tree_nodes_sorted[(left + right) / 2];
  double pivot = _x[pivot_seed][active_dim];

  /* partition */
  while (i <= j)
  {
    while (_x[tree_nodes_sorted[i]][active_dim] < pivot)
      i++;
    while (_x[tree_nodes_sorted[j]][active_dim] > pivot)
      j--;
    if (i <= j)
    {
      size_t tmp_index = tree_nodes_sorted[i];
      tree_nodes_sorted[i] = tree_nodes_sorted[j];
      tree_nodes_sorted[j] = tmp_index;

      i++;
      if (j > 0) j--;
    }
  };

  /* recursion */

  if (j > 0 && left < j && left <= target_pos && j >= target_pos)
    kd_tree_quicksort_adjust_target_position(target_pos, left, j, active_dim, tree_nodes_sorted);
  if (i < right && i <= target_pos && right >= target_pos)
    kd_tree_quicksort_adjust_target_position(target_pos, i, right, active_dim, tree_nodes_sorted);

  return 0;
//#pragma endregion
}

int VPS::kd_tree_add_point(size_t seed_index)
{
//#pragma region kd tree add point:
  // insert sphere into tree
  size_t parent_index(_tree_origin); size_t d_index(0);
  size_t branch_height(1);
  while (true)
  {
    if (_x[seed_index][d_index] >  _x[parent_index][d_index])
    {
      if (_tree_right[parent_index] == parent_index)
      {
        _tree_right[parent_index] = seed_index;
        branch_height++;
        break;
      }
      else
      {
        parent_index = _tree_right[parent_index];
        branch_height++;
      }
    }
    else
    {
      if (_tree_left[parent_index] == parent_index)
      {
        _tree_left[parent_index] = seed_index;
        branch_height++;
        break;
      }
      else
      {
        parent_index = _tree_left[parent_index];
        branch_height++;
      }
    }
    d_index++;
    if (d_index == _num_dim) d_index = 0;
  }
  if (branch_height > _tree_max_height) _tree_max_height = branch_height;
  return 0;
//#pragma endregion
}

int VPS::kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
                                     double r,                                                     // neighborhood radius
                                     size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
                                     size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
                                     size_t &capacity                                              // Size of points in sphere array
  )
{
//#pragma region kd tree neighbor search:
  if (d_index == _num_dim) d_index = 0;

  double dst_sq(0.0);
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    double dx = _x[node_index][idim] - x[idim];
    dst_sq += dx * dx;
  }
  if (dst_sq < r * r + DST_TOL) add_entry(node_index, num_points_in_sphere, points_in_sphere, capacity);

  bool check_right(false), check_left(false);
  double neighbor_min(x[d_index] - r), neighbor_max(x[d_index] + r);

  if (_tree_right[node_index] != node_index && neighbor_max >  _x[node_index][d_index]) check_right = true;
  if (_tree_left[node_index] != node_index && neighbor_min < _x[node_index][d_index]) check_left = true;

  if (check_right) kd_tree_get_seeds_in_sphere(x, r, d_index + 1, _tree_right[node_index], num_points_in_sphere, points_in_sphere, capacity);
  if (check_left)  kd_tree_get_seeds_in_sphere(x, r, d_index + 1, _tree_left[node_index], num_points_in_sphere, points_in_sphere, capacity);

  return 0;
//#pragma endregion
}

int VPS::kd_tree_get_closest_seed(size_t seed_index,                                  // seed index
                                  size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
                                  size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
  )
{
//#pragma region kd tree closest neighbor search:
  if (d_index == _num_dim) d_index = 0;

  if (seed_index != node_index)
  {
    double dst_sq(0.0);
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = _x[seed_index][idim] - _x[node_index][idim];
      dst_sq += dx * dx;
    }

    if (dst_sq < closest_distance * closest_distance)
    {
      // add to neighbors
      closest_seed = node_index;
      closest_distance = sqrt(dst_sq);
    }
  }

  bool check_right(false), check_left(false);

  double neighbor_min, neighbor_max;
  if (closest_distance == DBL_MAX){ check_right = true; check_left = true; }
  else
  {
    neighbor_min = _x[seed_index][d_index] - closest_distance; neighbor_max = _x[seed_index][d_index] + closest_distance;
    if (_tree_right[node_index] != node_index && neighbor_max >  _x[node_index][d_index]) check_right = true;
    if (_tree_left[node_index] != node_index && neighbor_min < _x[node_index][d_index]) check_left = true;
  }


  if (check_right && check_left)
  {
    // check the half that x blongs to first
    if (_x[seed_index][d_index] > _x[node_index][d_index])
    {
      kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
      neighbor_min = _x[seed_index][d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
      if (neighbor_min < _x[node_index][d_index])
      {
        kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
      }
    }
    else
    {
      kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
      neighbor_max = _x[seed_index][d_index] + closest_distance;
      if (neighbor_max > _x[node_index][d_index])
      {
        kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
      }
    }
  }
  else if (check_right) kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
  else if (check_left)  kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

  return 0;
//#pragma endregion
}

int VPS::kd_tree_get_closest_seed(double* x,                                          // point location
                                  size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
                                  size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
  )
{
//#pragma region kd tree closest neighbor search:
  if (d_index == _num_dim) d_index = 0;


  double dst_sq(0.0);
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    double dx = x[idim] - _x[node_index][idim];
    dst_sq += dx * dx;
  }

  if (dst_sq < closest_distance * closest_distance)
  {
    // add to neighbors
    closest_seed = node_index;
    closest_distance = sqrt(dst_sq);
  }

  bool check_right(false), check_left(false);

  double neighbor_min, neighbor_max;
  if (closest_distance == DBL_MAX){ check_right = true; check_left = true; }
  else
  {
    neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
    if (_tree_right[node_index] != node_index && neighbor_max > _x[node_index][d_index]) check_right = true;
    if (_tree_left[node_index] != node_index  && neighbor_min < _x[node_index][d_index]) check_left = true;
  }

  if (check_right && check_left)
  {
    // check the half that x blongs to first
    if (x[d_index] > _x[node_index][d_index])
    {
      kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
      neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
      if (neighbor_min < _x[node_index][d_index])
      {
        kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
      }

    }
    else
    {
      kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
      neighbor_max = x[d_index] + closest_distance;
      if (neighbor_max > _x[node_index][d_index])
      {
        kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
      }
    }
  }
  else if (check_right) kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
  else if (check_left)  kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

  return 0;
//#pragma endregion
}

int VPS::kd_tree_get_closest_seed(double* x,                                          // a point in space
                                  size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
                                  size_t num_exculded_seeds, size_t* exculded_seeds,  // Number of excluded seeds and their indices
                                  size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
  )
{
//#pragma region kd tree closest neighbor search:
  if (d_index == _num_dim) d_index = 0;


  if (!find_brute(node_index, exculded_seeds, num_exculded_seeds))
  {
    double dst_sq(0.0);
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = x[idim] - _x[node_index][idim];
      dst_sq += dx * dx;
    }

    if (dst_sq < closest_distance * closest_distance)
    {
      // add to neighbors
      closest_seed = node_index;
      closest_distance = sqrt(dst_sq);
    }
  }

  bool check_right(false), check_left(false);

  double neighbor_min, neighbor_max;
  if (closest_distance == DBL_MAX){ check_right = true; check_left = true; }
  else
  {
    neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
    if (_tree_right[node_index] != node_index && neighbor_max > _x[node_index][d_index]) check_right = true;
    if (_tree_left[node_index] != node_index && neighbor_min < _x[node_index][d_index]) check_left = true;
  }

  if (check_right && check_left)
  {
    // check the half that x blongs to first
    if (x[d_index] > _x[node_index][d_index])
    {
      kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
      neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
      if (neighbor_min < _x[node_index][d_index])
      {
        kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
      }

    }
    else
    {
      kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
      neighbor_max = x[d_index] + closest_distance;
      if (neighbor_max > _x[node_index][d_index])
      {
        kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
      }
    }
  }
  else if (check_right) kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
  else if (check_left)  kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);

  return 0;
//#pragma endregion
}

int VPS::get_closest_seed_tree(size_t seed_index, size_t &closest_seed, double &closest_distance)
{
//#pragma region Closest Neighbor Search using tree:

  closest_seed = _budget;
  kd_tree_get_closest_seed(seed_index, 0, _tree_origin, closest_seed, closest_distance);

  return 0;
//#pragma endregion
}

int VPS::get_closest_seed_tree(double* x, size_t &closest_seed, double &closest_distance)
{
//#pragma region Closest Neighbor Search using tree:
  closest_seed = _budget;
  kd_tree_get_closest_seed(x, 0, _tree_origin, closest_seed, closest_distance);
  return 0;
//#pragma endregion
}

int VPS::get_closest_seed_tree(double* x, size_t num_exculded_seeds, size_t* exculded_seeds, size_t &closest_seed, double& closest_distance)
{
//#pragma region Closest Neighbor Search using a tree:

  closest_seed = _budget;
  kd_tree_get_closest_seed(x, 0, _tree_origin, num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
  return 0;

//#pragma endregion
}

int VPS::get_seeds_in_sphere_tree(double* x, double r, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity)
{
//#pragma region tree neighbor search:
  num_points_in_sphere = 0;
  kd_tree_get_seeds_in_sphere(x, r, 0, _tree_origin, num_points_in_sphere, points_in_sphere, capacity);

  return 0;
//#pragma endregion
}

int VPS::get_closest_seed_brute(double* x, size_t num_significant_neighbors, size_t* significant_neighbors, size_t num_exculded_seeds, size_t* exculded_seeds, size_t &closest_seed)
{
//#pragma region Brute Force closest Neighbor search over significant non-excluded seeds:
  // Brute Force to retrieve Significant closest Seed

  double closest_distance_sq = DBL_MAX; closest_seed = _budget;
  for (size_t isig_seed = 0; isig_seed < num_significant_neighbors; isig_seed++)
  {
    size_t iseed = significant_neighbors[isig_seed];

    bool excluded(false);
    for (size_t iex = 0; iex < num_exculded_seeds; iex++)
    {
      if (exculded_seeds[iex] == iseed){ excluded = true; break; }
    }
    if (excluded) continue;

    double dst_sq(0.0);
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = _x[iseed][idim] - x[idim];
      dst_sq += dx * dx;
    }
    if (dst_sq < closest_distance_sq)
    {
      closest_distance_sq = dst_sq;
      closest_seed = iseed;
    }
  }
  return 0;
//#pragma endregion
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// General Methods      //////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int VPS::quicksort(size_t* I, size_t left, size_t right)
{
//#pragma region Quick Sort:
  size_t i = left, j = right;
  size_t pivot = I[(left + right) / 2];

  /* partition */
  while (i <= j)
  {
    while (I[i] < pivot)
      i++;
    while (I[j] > pivot)
      j--;
    if (i <= j)
    {
      size_t tmp = I[i];
      I[i] = I[j];
      I[j] = tmp;

      i++;
      if (j > 0) j--;
    }
  };

  /* recursion */

  if (j > 0 && left < j)
    quicksort(I, left, j);
  if (i < right)
    quicksort(I, i, right);

  return 0;
//#pragma endregion
}

int VPS::quicksort(double* V, size_t* I, size_t left, size_t right)
{
//#pragma region Quick Sort:
  size_t i = left, j = right;
  double pivot = V[(left + right) / 2];

  /* partition */
  while (i <= j)
  {
    while (V[i] < pivot)
      i++;
    while (V[j] > pivot)
      j--;
    if (i <= j)
    {
      double tmpd = V[i]; V[i] = V[j]; V[j] = tmpd;
      size_t tmpi = I[i]; I[i] = I[j]; I[j] = tmpi;
      i++;
      if (j > 0) j--;
    }
  };

  /* recursion */
  if (j > 0 && left < j)
    quicksort(V, I, left, j);

  if (i < right)
    quicksort(V, I, i, right);

  return 0;
//#pragma endregion
}

bool VPS::find_binary(size_t entry, size_t* I, size_t left, size_t right)
{
//#pragma region find using Binary search:
  size_t ipivot = (left + right) / 2;
  size_t pivot = I[ipivot];
  if (pivot == entry) return true;
  if (left == right)  return false;

  if (pivot < entry && ipivot + 1 <= right) return find_binary(entry, I, ipivot + 1, right);
  if (pivot > entry && left + 1 <= ipivot) return find_binary(entry, I, left, ipivot - 1);

  return false;
//#pragma endregion
}

bool VPS::find_brute(size_t entry, size_t* I, size_t num_entries)
{
//#pragma region find using Brutal search:
  for (size_t i = 0; i < num_entries; i++) if (I[i] == entry) return true;

  return false;
//#pragma endregion
}


int VPS::add_entry(size_t entry, size_t &num_entries, size_t* &I, size_t &capacity)
{
//#pragma region add an entry to a list:
  // add to neighbors
  I[num_entries] = entry;
  num_entries++;
  if (num_entries == capacity)
  {
    capacity *= 2;
    size_t* tmp = new size_t[capacity];
    for (size_t i = 0; i < num_entries; i++) tmp[i] = I[i];
    delete[] I;
    I = tmp;
  }
  return 0;
//#pragma endregion
}



size_t VPS::retrieve_num_permutations(size_t num_dim, size_t upper_bound, bool force_sum_constraint, size_t sum_constraint)
{
//#pragma region Retrieve Number of Permutation:
  size_t* t = new size_t[num_dim];
  for (size_t idim = 0; idim < num_dim; idim++) t[idim] = 0;

  // count output
  size_t num_perm = 0;
  size_t k_dim(num_dim - 1);
  while (true)
  {
    while (t[k_dim] <= upper_bound)
    {
      bool valid(true);

      if (force_sum_constraint)
      {
        size_t s_const(0);
        for (size_t idim = 0; idim < num_dim; idim++) s_const += t[idim];

        if (s_const > upper_bound)
        {
          valid = false;
        }
      }
      if (valid) num_perm++;

      t[k_dim]++; // move to the next enumeration
    }


    size_t kk_dim(k_dim - 1);

    bool done(false);
    while (true)
    {
      t[kk_dim]++;
      if (t[kk_dim] > upper_bound)
      {
        t[kk_dim] = 0;
        if (kk_dim == 0)
        {
          done = true;
          break;
        }
        kk_dim--;
      }
      else break;
    }
    if (done) break;
    t[k_dim] = 0;
  }
  delete[] t;
  return num_perm;
//#pragma endregion
}


int VPS::retrieve_permutations(size_t &num_perm, size_t** &perm, size_t num_dim, size_t upper_bound, bool force_sum_constraint, size_t sum_constraint)
{
//#pragma region Retrieve Permutation:
  size_t* t = new size_t[num_dim];
  num_perm = retrieve_num_permutations(num_dim, upper_bound, force_sum_constraint, sum_constraint);

  //std::cout<< "Number of permutations = " << m << std::endl;

  perm = new size_t*[num_perm];
  for (size_t i = 0; i < num_perm; i++)
  {
    perm[i] = new size_t[num_dim];
    for (size_t idim = 0; idim < num_dim; idim++) perm[i][idim] = 0;
  }

  // Taxi Counter to evaluate the surrogate
  for (size_t idim = 0; idim < num_dim; idim++) t[idim] = 0;
  size_t k_dim = num_dim - 1;
  num_perm = 0; // index of alpha

  while (true)
  {
    while (t[k_dim] <= upper_bound)
    {
      bool valid(true);

      if (force_sum_constraint)
      {
        size_t s_const(0);
        for (size_t idim = 0; idim < num_dim; idim++) s_const += t[idim];

        if (s_const > upper_bound)
        {
          valid = false;
        }
      }
      if (valid)
      {
        // store t in t_alpha of counter
        for (size_t idim = 0; idim < num_dim; idim++) perm[num_perm][idim] = t[idim];

        //for (size_t idim = 0; idim < num_dim; idim++) std::cout << perm[m][idim] << " ";
        //std::cout << std::endl;

        num_perm++;
      }

      t[k_dim]++; // move to the next enumeration
    }


    size_t kk_dim(k_dim - 1);

    bool done(false);
    while (true)
    {
      t[kk_dim]++;
      if (t[kk_dim] > upper_bound)
      {
        t[kk_dim] = 0;
        if (kk_dim == 0)
        {
          done = true;
          break;
        }
        kk_dim--;
      }
      else break;
    }
    if (done) break;
    t[k_dim] = 0;
  }

  // re-order perm:
  size_t num_basis = num_perm;
  //std::cout << "*** VPS:: num_basis = " << num_basis << std::endl;
  if (num_basis >= num_dim + 1)
  {
    // reorder first order terms
    for (size_t idim = 0; idim < num_dim; idim++)
    {
      // search fot the corresponding basis and bring it to its proper location
      for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
      {
        size_t sum(0), sig_dim(0);
        for (size_t jdim = 0; jdim < num_dim; jdim++)
        {
          sum += perm[ibasis][jdim];
          if (perm[ibasis][jdim] == 1) sig_dim = jdim;
        }
        if (sum == 1 && sig_dim == idim)
        {
          // move this basis to location idim + 1
          size_t* tmp = perm[ibasis];
          perm[ibasis] = perm[idim + 1];
          perm[idim + 1] = tmp;
          break;
        }
      }
    }
  }

  if (num_basis >= num_dim * (num_dim + 1) / 2 + num_dim + 1)
  {
    // reorder second order terms
    size_t iloc(num_dim + 1);
    for (size_t idim = 0; idim < num_dim; idim++)
    {
      for (size_t jdim = idim; jdim < num_dim; jdim++)
      {
        // search fot the corresponding basis and bring it to its proper location
        for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
        {
          size_t sum(0), sig_dim_i(num_dim), sig_dim_j(num_dim);
          for (size_t kdim = 0; kdim < num_dim; kdim++)
          {
            sum += perm[ibasis][kdim];
            if (perm[ibasis][kdim] > 0 && sig_dim_i == num_dim)
            {
              sig_dim_i = kdim;
              sig_dim_j = sig_dim_i;
            }
            else if (perm[ibasis][kdim] > 0 && sig_dim_i < num_dim) sig_dim_j = kdim;
          }
          if (sum == 2 && sig_dim_i == idim && sig_dim_j == jdim)
          {
            // move this basis to location iloc
            size_t* tmp = perm[ibasis];
            perm[ibasis] = perm[iloc];
            perm[iloc] = tmp; iloc++;
            break;
          }
        }
      }
    }
  }

  // make sure that the permutation order increases
  for (size_t i = 0; i < num_perm; i++)
  {
    size_t sum_i = 0;
    for (size_t idim = 0; idim < num_dim; idim++) sum_i += perm[i][idim];
    for (size_t j = i + 1; j < num_perm; j++)
    {
      size_t sum_j = 0;
      for (size_t idim = 0; idim < num_dim; idim++) sum_j += perm[j][idim];

      if (sum_j < sum_i)
      {
        size_t* tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
        sum_i = sum_j;
      }
    }
  }

  /*
    for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
    {
    for (size_t idim = 0; idim < num_dim; idim++) std::cout << perm[ibasis][idim] << " ";
    std::cout << std::endl;
    }
  */

  delete[] t;

  return 0;
//#pragma endregion
}


void VPS::plot_graph(std::string file_name, size_t num_points, double* x, size_t num_functions, double** f)
{
//#pragma region Plot Graph:
  std::fstream file(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  double xmin(x[0]), ymin(f[0][0]), xmax(xmin), ymax(ymin);
  for (size_t ipoint = 0; ipoint < num_points; ipoint++)
  {
    if (x[ipoint] < xmin) xmin = x[ipoint];
    if (x[ipoint] > xmax) xmax = x[ipoint];

    for (size_t ifunc = 0; ifunc < num_functions; ifunc++)
    {
      if (f[ifunc][ipoint] < ymin)
        ymin = f[ifunc][ipoint];
      if (f[ifunc][ipoint] > ymax)
        ymax = f[ifunc][ipoint];
    }
  }

  if (xmax < xmin + 0.01 * (xmax - xmin)) xmax = xmin + 0.01 * (xmax - xmin);
  if (ymax < ymin + 0.01 * (ymax - ymin)) ymax = ymin + 0.01 * (ymax - ymin);

  double Lx(xmax - xmin);
  double Ly(ymax - ymin);

  double scale_x, scale_y;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  shift_x = 1.0 - xmin * scale_x;
  shift_y = 1.0 - ymin * scale_y;
  file << shift_x << " " << shift_y << " translate" << std::endl;

  file << "/Courier findfont" << std::endl;
  file << "0.12 scalefont" << std::endl;
  file << "setfont" << std::endl;

  for (size_t ifunc = 0; ifunc < num_functions; ifunc++)
  {
    double gs = ifunc * 1.0 / num_functions;

    double r, g, b;

    if (gs < 0.25)     r = 1.0;
    //else if (gs < 0.5) r = 2.0 - 4.0 * gs;
    else if (gs < 0.5) r = 1.0 - 16.0 * (gs - 0.25) * (gs - 0.25);
    else               r = 0.0;

    double go(0.25), gn(1.0 - go);
    if (gs < go)      g = gs / go;
    else if (gs < gn) g = 1.0;
    else              g = 1.0 / (1.0 - gn) - gs / (1.0 - gn);


    if (gs < 0.5)       b = 0.0;
    else if (gs < 0.75) b = 1.0 - 16.0 * (gs - 0.75) * (gs - 0.75);
    else                b = 1.0;

    for (size_t ipoint = 0; ipoint < num_points - 1; ipoint++)
    {
      size_t jpoint = ipoint + 1;

      file << "newpath" << std::endl;
      file << x[ipoint] * scale_x << " " << f[ifunc][ipoint] * scale_y << " moveto" << std::endl;
      file << x[jpoint] * scale_x << " " << f[ifunc][jpoint] * scale_y << " lineto" << std::endl;
      file << " closepath" << std::endl;
      file << " gsave" << std::endl;
      file << " grestore" << std::endl;
      file << " " << r << " " << g << " " << b << " setrgbcolor" << std::endl;
      file << " 0.01 setlinewidth" << std::endl;
      file << " stroke" << std::endl;
    }
  }

  size_t num_steps(10);
  double dx((xmax - xmin) / num_steps);
  double dy((ymax - ymin) / num_steps);

  // Vertical Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin;
    double yy = ymin + istep * dy;

    std::stringstream ss;
    ss << trunc((yy + 0.0001) * 1000) / 1000;
    std::string str = ss.str();

    file << "newpath " << (xx - dx) * scale_x << " " << yy * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xmin * scale_x << " " << yy * scale_y << " moveto" << std::endl;
    file << xmax * scale_x << " " << yy * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  // Horizontal Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin + istep * dx;
    double yy = ymin;

    std::stringstream ss;
    ss << trunc((xx + 0.001) * 100) / 100;
    std::string str = ss.str();

    file << "newpath " << xx * scale_x << " " << (yy - 0.2 * dy) * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xx * scale_x << " " << ymin * scale_y << " moveto" << std::endl;
    file << xx * scale_x << " " << ymax * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

//#pragma endregion
}

void VPS::plot_polynomial(std::string file_name, size_t num_basis, double* c, double xmin, double xmax, size_t num_points, double* px, double* py)
{
//#pragma region Plot Polynomial:
  std::fstream file(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  size_t num_p_points = 101;
  double dx = (xmax - xmin) / (num_p_points - 1);

  double ymin(DBL_MAX), ymax(-DBL_MAX);
  for (size_t ipoint = 0; ipoint < num_p_points - 1; ipoint++)
  {
    //size_t jpoint = ipoint + 1; // ETP

    double x = xmin + ipoint * dx;
    double f = 0.0;
    for (size_t ibasis = 0; ibasis < num_basis; ibasis++) f += c[ibasis] * pow(x, ibasis);
    if (f < ymin) ymin = f;
    if (f > ymax) ymax = f;
  }

  double Lx(xmax - xmin);
  double Ly(ymax - ymin);

  if (Ly < 1E-10){ymin -= 0.5; ymax += 0.5; Ly = 1.0;}

  double scale_x, scale_y;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  shift_x = 1.0 - xmin * scale_x;
  shift_y = 1.0 - ymin * scale_y;
  file << shift_x << " " << shift_y << " translate" << std::endl;

  file << "/Courier findfont" << std::endl;
  file << "0.12 scalefont" << std::endl;
  file << "setfont" << std::endl;


  size_t num_steps(10);
  double dxg((xmax - xmin) / num_steps);
  double dyg((ymax - ymin) / num_steps);

  // Vertical Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin;
    double yy = ymin + istep * dyg;

    std::stringstream ss;
    ss << trunc((yy + 0.0001) * 1000) / 1000;
    std::string str = ss.str();

    file << "newpath " << (xx - dxg) * scale_x << " " << yy * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xmin * scale_x << " " << yy * scale_y << " moveto" << std::endl;
    file << xmax * scale_x << " " << yy * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  // Horizontal Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin + istep * dxg;
    double yy = ymin;

    std::stringstream ss;
    ss << trunc((xx + 0.001) * 100) / 100;
    std::string str = ss.str();

    file << "newpath " << xx * scale_x << " " << (yy - 0.2 * dyg) * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xx * scale_x << " " << ymin * scale_y << " moveto" << std::endl;
    file << xx * scale_x << " " << ymax * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }


  for (size_t ipoint = 0; ipoint < num_p_points - 1; ipoint++)
  {
    //size_t jpoint = ipoint + 1; // ETP

    double x = xmin + ipoint * dx;
    double f = 0.0;
    for (size_t ibasis = 0; ibasis < num_basis; ibasis++) f += c[ibasis] * pow(x, ibasis);

    double xp = x + dx;
    double fp = 0.0;
    for (size_t ibasis = 0; ibasis < num_basis; ibasis++) fp += c[ibasis] * pow(xp, ibasis);

    file << "newpath" << std::endl;
    file << x * scale_x << " " << f * scale_y << " moveto" << std::endl;
    file << xp * scale_x << " " << fp * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 1.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  for (size_t ipoint = 0; ipoint < num_points; ipoint++)
  {
    file << "newpath" << std::endl;
    file << px[ipoint] * scale_x << " " << py[ipoint] * scale_y << " " << 0.005 * scale_x << " 0 360 arc" << std::endl;
    file << "closepath" << std::endl;
    file << "gsave" << std::endl;
    file << "0 0 1 setrgbcolor" << std::endl; // Inactive seeds
    file << "fill" << std::endl;
    file << "grestore" << std::endl;
    file << "0.0 setlinewidth" << std::endl;
    file << "0 0 0 setrgbcolor" << std::endl; // discs borders
    file << "stroke " << std::endl;
  }

//#pragma endregion
}

void VPS::plot_piecewise_polynomial(std::string file_name, size_t num_pieces, double* pmin, double* pmax, size_t* num_basis, double** c, size_t num_points, double* px, double* py)
{
//#pragma region Plot Polynomial:
  std::fstream file(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  double xmin = pmin[0]; double xmax = pmax[num_pieces - 1];

  double ymin(DBL_MAX), ymax(-DBL_MAX);
  for (size_t ipiece = 0; ipiece < num_pieces; ipiece++)
  {
    size_t num_steps(50);
    double ds = (pmax[ipiece] - pmin[ipiece]) / num_steps;
    for (size_t istep = 0; istep <= num_steps; istep++)
    {
      double x = pmin[ipiece] + istep * ds;
      double f = 0.0;
      for (size_t ibasis = 0; ibasis < num_basis[ipiece]; ibasis++) f += c[ipiece][ibasis] * pow(x, ibasis);
      if (f < ymin) ymin = f;
      if (f > ymax) ymax = f;
    }
  }

  double Lx(xmax - xmin);
  double Ly(ymax - ymin);

  if (Ly < 1E-10){ ymin -= 0.5; ymax += 0.5; Ly = 1.0; }

  double scale_x, scale_y;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  shift_x = 1.0 - xmin * scale_x;
  shift_y = 1.0 - ymin * scale_y;
  file << shift_x << " " << shift_y << " translate" << std::endl;

  file << "/Courier findfont" << std::endl;
  file << "0.12 scalefont" << std::endl;
  file << "setfont" << std::endl;


  size_t num_steps(10);
  double dxg((xmax - xmin) / num_steps);
  double dyg((ymax - ymin) / num_steps);

  // Vertical Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin;
    double yy = ymin + istep * dyg;

    std::stringstream ss;
    ss << trunc((yy + 0.0001) * 1000) / 1000;
    std::string str = ss.str();

    file << "newpath " << (xx - dxg) * scale_x << " " << yy * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xmin * scale_x << " " << yy * scale_y << " moveto" << std::endl;
    file << xmax * scale_x << " " << yy * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  // Horizontal Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin + istep * dxg;
    double yy = ymin;

    std::stringstream ss;
    ss << trunc((xx + 0.001) * 100) / 100;
    std::string str = ss.str();

    file << "newpath " << xx * scale_x << " " << (yy - 0.2 * dyg) * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xx * scale_x << " " << ymin * scale_y << " moveto" << std::endl;
    file << xx * scale_x << " " << ymax * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  for (size_t ipiece = 0; ipiece < num_pieces; ipiece++)
  {
    size_t num_steps(50);
    double ds = (pmax[ipiece] - pmin[ipiece]) / num_steps;
    for (size_t istep = 0; istep < num_steps; istep++)
    {
      double x = pmin[ipiece] + istep * ds;
      double f = 0.0;
      for (size_t ibasis = 0; ibasis < num_basis[ipiece]; ibasis++) f += c[ipiece][ibasis] * pow(x, ibasis);

      double xp = x + ds;
      double fp = 0.0;
      for (size_t ibasis = 0; ibasis < num_basis[ipiece]; ibasis++) fp += c[ipiece][ibasis] * pow(xp, ibasis);

      file << "newpath" << std::endl;
      file << x * scale_x << " " << f * scale_y << " moveto" << std::endl;
      file << xp * scale_x << " " << fp * scale_y << " lineto" << std::endl;
      file << " closepath" << std::endl;
      file << " gsave" << std::endl;
      file << " grestore" << std::endl;
      file << " " << 1.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
      file << " 0.01 setlinewidth" << std::endl;
      file << " stroke" << std::endl;
    }
  }

  for (size_t ipoint = 0; ipoint < num_points; ipoint++)
  {
    file << "newpath" << std::endl;
    file << px[ipoint] * scale_x << " " << py[ipoint] * scale_y << " " << 0.005 * scale_x << " 0 360 arc" << std::endl;
    file << "closepath" << std::endl;
    file << "gsave" << std::endl;
    file << "0 0 1 setrgbcolor" << std::endl; // Inactive seeds
    file << "fill" << std::endl;
    file << "grestore" << std::endl;
    file << "0.0 setlinewidth" << std::endl;
    file << "0 0 0 setrgbcolor" << std::endl; // discs borders
    file << "stroke " << std::endl;
  }

//#pragma endregion
}

void VPS::plot_FourierExpansion(std::string file_name, size_t num_basis, double xmin, double xmax, double* a, double* b, size_t num_points, double* px, double* py, size_t* num_cbasis, double** c)
{
//#pragma region Plot Fourier Expansion:
  std::fstream file(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  size_t num_p_points = 501;
  double dx = (xmax - xmin) / (num_p_points - 1);


  double DX = xmax - xmin;

  double ymin(DBL_MAX), ymax(-DBL_MAX);
  for (size_t ipoint = 0; ipoint < num_p_points; ipoint++)
  {
    // size_t jpoint = ipoint + 1; // ETP

    double x = xmin + ipoint * dx;
    double f = 0.5 * a[0];
    for (size_t ibasis = 1; ibasis < num_basis; ibasis++) f += a[ibasis] * cos(2.0 * PI * ibasis * x / DX) + b[ibasis] * sin(2.0 * PI * ibasis * x / DX);
    if (f < ymin) ymin = f;
    if (f > ymax) ymax = f;
  }

  double Lx(xmax - xmin);
  double Ly(ymax - ymin);

  if (Ly < 1E-10){ ymin -= 0.5; ymax += 0.5; Ly = 1.0; }

  double scale_x, scale_y;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  shift_x = 1.0 - xmin * scale_x;
  shift_y = 1.0 - ymin * scale_y;
  file << shift_x << " " << shift_y << " translate" << std::endl;

  file << "/Courier findfont" << std::endl;
  file << "0.12 scalefont" << std::endl;
  file << "setfont" << std::endl;


  size_t num_steps(10);
  double dxg((xmax - xmin) / num_steps);
  double dyg((ymax - ymin) / num_steps);

  // Vertical Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin;
    double yy = ymin + istep * dyg;

    std::stringstream ss;
    ss << trunc((yy + 0.0001) * 1000) / 1000;
    std::string str = ss.str();

    file << "newpath " << (xx - dxg) * scale_x << " " << yy * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xmin * scale_x << " " << yy * scale_y << " moveto" << std::endl;
    file << xmax * scale_x << " " << yy * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  // Horizontal Grid
  for (size_t istep = 0; istep <= num_steps; istep++)
  {
    double xx = xmin + istep * dxg;
    double yy = ymin;

    std::stringstream ss;
    ss << trunc((xx + 0.001) * 100) / 100;
    std::string str = ss.str();

    file << "newpath " << xx * scale_x << " " << (yy - 0.2 * dyg) * scale_y << " moveto (" << str << ") show" << std::endl;

    file << "newpath" << std::endl;
    file << xx * scale_x << " " << ymin * scale_y << " moveto" << std::endl;
    file << xx * scale_x << " " << ymax * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }


  for (size_t ipoint = 1; ipoint < num_p_points; ipoint++)
  {
    double xm = xmin + (ipoint - 1) * dx;
    double fm = 0.5 * a[0];
    for (size_t ibasis = 1; ibasis < num_basis; ibasis++) fm += a[ibasis] * cos(2.0 * PI * ibasis * xm / DX) + b[ibasis] * sin(2.0 * PI * ibasis * xm / DX);

    double xx = xm + dx;
    double ff = 0.5 * a[0];
    for (size_t ibasis = 1; ibasis < num_basis; ibasis++) ff += a[ibasis] * cos(2.0 * PI * ibasis * xx / DX) + b[ibasis] * sin(2.0 * PI * ibasis * xx / DX);

    file << "newpath" << std::endl;
    file << xm * scale_x << " " << fm * scale_y << " moveto" << std::endl;
    file << xx * scale_x << " " << ff * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 1.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;

    xm = xmin + (ipoint - 1) * dx;
    fm = 0.5 * a[0];
    for (size_t ibasis = 1; ibasis < num_basis/2; ibasis++) fm += a[ibasis] * cos(2.0 * PI * ibasis * xm / DX) + b[ibasis] * sin(2.0 * PI * ibasis * xm / DX);

    xx = xm + dx;
    ff = 0.5 * a[0];
    for (size_t ibasis = 1; ibasis < num_basis/2; ibasis++) ff += a[ibasis] * cos(2.0 * PI * ibasis * xx / DX) + b[ibasis] * sin(2.0 * PI * ibasis * xx / DX);

    file << "newpath" << std::endl;
    file << xm * scale_x << " " << fm * scale_y << " moveto" << std::endl;
    file << xx * scale_x << " " << ff * scale_y << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " gsave" << std::endl;
    file << " grestore" << std::endl;
    file << " " << 0.0 << " " << 0.0 << " " << 1.0 << " setrgbcolor" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
  }

  size_t num_pieces = num_points - 1;
  for (size_t ipiece = 0; ipiece < num_pieces; ipiece++)
  {
    size_t num_steps(50);
    double pmin = px[ipiece]; double pmax = px[ipiece + 1];
    double ds = (pmax - pmin) / num_steps;
    for (size_t istep = 0; istep < num_steps; istep++)
    {
      double x = pmin + istep * ds;
      double f = 0.0;
      for (size_t ibasis = 0; ibasis < num_cbasis[ipiece]; ibasis++) f += c[ipiece][ibasis] * pow(x, ibasis);

      double xp = x + ds;
      double fp = 0.0;
      for (size_t ibasis = 0; ibasis < num_cbasis[ipiece]; ibasis++) fp += c[ipiece][ibasis] * pow(xp, ibasis);

      file << "newpath" << std::endl;
      file << x * scale_x << " " << f * scale_y << " moveto" << std::endl;
      file << xp * scale_x << " " << fp * scale_y << " lineto" << std::endl;
      file << " closepath" << std::endl;
      file << " gsave" << std::endl;
      file << " grestore" << std::endl;
      file << " " << 0.0 << " " << 1.0 << " " << 0.0 << " setrgbcolor" << std::endl;
      file << " 0.01 setlinewidth" << std::endl;
      file << " stroke" << std::endl;
    }
  }

  for (size_t ipoint = 0; ipoint < num_points; ipoint++)
  {
    file << "newpath" << std::endl;
    file << px[ipoint] * scale_x << " " << py[ipoint] * scale_y << " " << 0.005 * scale_x << " 0 360 arc" << std::endl;
    file << "closepath" << std::endl;
    file << "gsave" << std::endl;
    file << "0 0 1 setrgbcolor" << std::endl; // Inactive seeds
    file << "fill" << std::endl;
    file << "grestore" << std::endl;
    file << "0.0 setlinewidth" << std::endl;
    file << "0 0 0 setrgbcolor" << std::endl; // discs borders
    file << "stroke " << std::endl;
  }

//#pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// VPS General Methods      //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int VPS::find_domain_diagonal()
{
//#pragma region Calculate Domain Diagonal
  _diag = 0.0;
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    double dx = _xmax[idim] - _xmin[idim];
    _diag += dx * dx;
  }
  _diag = sqrt(_diag);
  return 0;
//#pragma endregion
}

void VPS::plot_delaunay_graph(const std::string outFile)
{
//#pragma region Plot a Delaunay Mesh:
  std::fstream file(outFile.c_str(), std::ios::out);

  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  double Lx(_xmax[0] - _xmin[0]), Ly(_xmax[1] - _xmin[1]);

  double scale_x, scale_y, scale;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  if (scale_x < scale_y)
  {
    scale = scale_x;
    shift_x = 1.0 - _xmin[0] * scale;
    shift_y = 0.5 * (11.0 - Ly * scale) - _xmin[1] * scale;
  }
  else
  {
    scale = scale_y;
    shift_x = 0.5 * (8.5 - Lx * scale) - _xmin[0] * scale;
    shift_y = 1.0 - _xmin[1] * scale;
  }
  file << shift_x << " " << shift_y << " translate" << std::endl;


  file << "/Courier findfont" << std::endl;
  file << "0.12 scalefont" << std::endl;
  file << "setfont" << std::endl;

  for (size_t i = 0; i < _num_samples; i++)
  {
    for (size_t j = 0; j < _seed_neighbors[i][1]; j++)
    {
      size_t neighbor = _seed_neighbors[i][2 + j];

      file << "newpath" << std::endl;
      // draw a line between isample and neighbor
      file << _x[i][0] * scale << "  " << _x[i][1] * scale << "  ";
      file << "moveto" << std::endl;
      file << _x[neighbor][0] * scale << "  " << _x[neighbor][1] * scale << "  ";
      file << "lineto" << std::endl;
      file << "0.0 setlinewidth" << std::endl;
      file << "0 0 0 setrgbcolor" << std::endl; // Element borders
      file << "stroke " << std::endl;
    }
  }

  for (size_t i = 0; i < _num_samples; i++)
  {
    file << "newpath" << std::endl;
    file << _x[i][0] * scale << " " << _x[i][1] * scale << " " << 0.005 * scale << " 0 360 arc" << std::endl;
    file << "closepath" << std::endl;
    file << "gsave" << std::endl;
    file << "0 0 0 setrgbcolor" << std::endl; // Inactive seeds
    file << "fill" << std::endl;
    file << "grestore" << std::endl;
    file << "0.0 setlinewidth" << std::endl;
    file << "0 0 0 setrgbcolor" << std::endl; // discs borders
    file << "stroke " << std::endl;

    std::stringstream ss;
    ss << i;
    std::string str = ss.str();

    file << "newpath " << _x[i][0] * scale << " " << _x[i][1] * scale << " moveto (" << str << ") show" << std::endl;
  }

  for (size_t i = 0; i < _num_vs; i++)
  {
    if (_vs[i] == 0) continue;

    file << "newpath" << std::endl;
    file << _vs[i][0] * scale << " " << _vs[i][1] * scale << " " << 0.005 * scale << " 0 360 arc" << std::endl;
    file << "closepath" << std::endl;
    file << "gsave" << std::endl;
    file << "1 0 0 setrgbcolor" << std::endl; // Inactive seeds
    file << "fill" << std::endl;
    file << "grestore" << std::endl;
    file << "0.0 setlinewidth" << std::endl;
    file << "1 0 0 setrgbcolor" << std::endl; // discs borders
    file << "stroke " << std::endl;

    std::stringstream ss;
    ss << i;
    std::string str = ss.str();

    file << "newpath " << _vs[i][0] * scale << " " << _vs[i][1] * scale << " moveto (" << str << ") show" << std::endl;
  }

  // draw bounding box
  file << "newpath" << std::endl;
  file << _xmin[0] * scale << "  " << _xmin[1] * scale << "  ";
  file << "moveto" << std::endl;
  file << _xmax[0] * scale << "  " << _xmin[1] * scale << "  ";
  file << "lineto" << std::endl;
  file << _xmax[0] * scale << "  " << _xmax[1] * scale << "  ";
  file << "lineto" << std::endl;
  file << _xmin[0] * scale << "  " << _xmax[1] * scale << "  ";
  file << "lineto" << std::endl;
  file << _xmin[0] * scale << "  " << _xmin[1] * scale << "  ";
  file << "lineto" << std::endl;
  file << "0.02 setlinewidth" << std::endl;
  file << "0 0 1 setrgbcolor" << std::endl; // Element borders
  file << "stroke " << std::endl;

  file << "showpage" << std::endl;
//#pragma endregion
}

void VPS::plot_vps_surrogate(std::string file_name, size_t function_index, size_t num_contours, bool plot_graph)
{
//#pragma region Plot Solid Isocontours:
  //std::cout << ".: VPS Debug Mode :. Plotting ps files .... " << std::endl;

  std::fstream file(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  double xmin(_xmin[0]);
  double ymin(_xmin[1]);
  double Lx(_xmax[0] - _xmin[0]);
  double Ly(_xmax[1] - _xmin[0]);

  double scale_x, scale_y, scale;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  if (scale_x < scale_y)
  {
    scale = scale_x;
    shift_x = 1.0 - xmin * scale;
    shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
  }
  else
  {
    scale = scale_y;
    shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
    shift_y = 1.0 - ymin * scale;
  }
  file << shift_x << " " << shift_y << " translate" << std::endl;

  file << "/redseg      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " 1 0 0 setrgbcolor" << std::endl;
  file << " 0.01 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/greenseg      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 1 0 setrgbcolor" << std::endl;
  file << " 0.01 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;


  file << "/blueseg      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 1 setrgbcolor" << std::endl;
  file << " 0.02 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/blackquad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.02 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/circ    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/blackfcirc    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.0 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/redfcirc    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 1 0 0 setrgbcolor" << std::endl;
  file << " fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.0 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/bluefcirc    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0 0 1 setrgbcolor" << std::endl;
  file << " fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.0 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/greenfcirc    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0 1 0 setrgbcolor" << std::endl;
  file << " fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.0 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 1.0 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << "} def" << std::endl;

  file << "/quad_bold      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.01 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/Courier findfont" << std::endl;
  file << "0.12 scalefont" << std::endl;
  file << "setfont" << std::endl;

  std::vector<double> poly_x;
  std::vector<double> poly_y;

  size_t num_cells(100);
  double* xx = new double[_num_dim];
  double* f = new double[_num_functions];
  double sx = (_xmax[0] - _xmin[0]) / num_cells;
  double sy = (_xmax[1] - _xmin[1]) / num_cells;


  double fmin, fmax;

  for (size_t i = 0; i < num_cells; i++)
  {
    xx[0] = _xmin[0] + i * sx;
    for (size_t j = 0; j < num_cells; j++)
    {
      xx[1] = _xmin[1] + j * sy;
      evaluate_surrogate(xx, f);
      if (i == 0)
      {
        fmin = f[function_index];
        fmax = f[function_index];
        continue;
      }
      if (f[function_index] < fmin) fmin = f[function_index];
      if (f[function_index] > fmax) fmax = f[function_index];
    }
  }

  std::vector<double> contours;
  contours.push_back(fmin - 2 * (fmax - fmin));
  for (size_t i = 0; i < num_contours; i++) contours.push_back(fmin + (1.0 / num_contours) * i * (fmax - fmin));
  contours.push_back(fmax + 2 * (fmax - fmin));

  for (size_t i = 0; i < num_cells; i++)
  {
    double xo = _xmin[0] + i * sx;
    for (size_t j = 0; j < num_cells; j++)
    {
      double fo(0.0), f1(0.0), f2(0.0), f3(0.0);

      double yo = _xmin[1] + j * sy;
      xx[0] = xo; xx[1] = yo;
      evaluate_surrogate(xx, f);
      fo = f[function_index];

      xx[0] = xo + sx; xx[1] = yo;
      evaluate_surrogate(xx, f);
      f1 = f[function_index];

      xx[0] = xo + sx; xx[1] = yo + sy;
      evaluate_surrogate(xx, f);
      f2 = f[function_index];

      xx[0] = xo; xx[1] = yo + sy;
      evaluate_surrogate(xx, f);
      f3 = f[function_index];


      size_t num_isocontours = contours.size();
      for (size_t icont = 0; icont < num_isocontours; icont++)
      {
        double contour = contours[icont];
        double contour_m = -1000.00;
        if (icont > 0) contour_m = contours[icont - 1];

        //std::cout<< "contour_m = " << contour_m << " , contour = " << contour << std::endl;

        poly_x.clear(); poly_y.clear();

        // moving right
        if (fo >= contour_m - 1E-10 && fo < contour + 1E-10)
        {
          poly_x.push_back(xo);
          poly_y.push_back(yo);
          if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
          {
            double h = sx * (contour - fo) / (f1 - fo);
            poly_x.push_back(xo + h);
            poly_y.push_back(yo);
          }
          else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
          {
            double h = sx * (contour_m - fo) / (f1 - fo);
            poly_x.push_back(xo + h);
            poly_y.push_back(yo);
          }
        }
        else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
        {
          double hm = sx * (contour_m - fo) / (f1 - fo);
          double h = hm;
          if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
          {
            h = sx * (contour - fo) / (f1 - fo);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo + hm);
          poly_y.push_back(yo);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo + h);
            poly_y.push_back(yo);
          }
        }
        else if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
        {
          double h = sx * (contour - fo) / (f1 - fo);
          poly_x.push_back(xo + h);
          poly_y.push_back(yo);
        }

        // moving up
        if (f1 >= contour_m - 1E-10 && f1 < contour + 1E-10)
        {
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo);
          if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
          {
            double h = sy * (contour - f1) / (f2 - f1);
            poly_x.push_back(xo + sx);
            poly_y.push_back(yo + h);
          }
          else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
          {
            double h = sy * (contour_m - f1) / (f2 - f1);
            poly_x.push_back(xo + sx);
            poly_y.push_back(yo + h);
          }

        }
        else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
        {
          double hm = sy * (contour_m - f1) / (f2 - f1);
          double h = hm;
          if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
          {
            h = sy * (contour - f1) / (f2 - f1);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo + hm);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo + sx);
            poly_y.push_back(yo + h);
          }
        }
        else if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
        {
          double h = sy * (contour - f1) / (f2 - f1);
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo + h);
        }

        // moving left
        if (f2 >= contour_m - 1E-10 && f2 < contour + 1E-10)
        {
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo + sy);
          if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
          {
            double h = sx * (contour - f2) / (f3 - f2);
            poly_x.push_back(xo + sx - h);
            poly_y.push_back(yo + sy);
          }
          else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
          {
            double h = sx * (contour_m - f2) / (f3 - f2);
            poly_x.push_back(xo + sx - h);
            poly_y.push_back(yo + sy);
          }
        }
        else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
        {
          double hm = sx * (contour_m - f2) / (f3 - f2);
          double h = hm;
          if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
          {
            h = sx * (contour - f2) / (f3 - f2);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo + sx - hm);
          poly_y.push_back(yo + sy);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo + sx - h);
            poly_y.push_back(yo + sy);
          }
        }
        else if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
        {
          double h = sx * (contour - f2) / (f3 - f2);
          poly_x.push_back(xo + sx - h);
          poly_y.push_back(yo + sy);
        }

        // moving down
        if (f3 >= contour_m - 1E-10 && f3 < contour + 1E-10)
        {
          poly_x.push_back(xo);
          poly_y.push_back(yo + sy);
          if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
          {
            double h = sy * (contour - f3) / (fo - f3);
            poly_x.push_back(xo);
            poly_y.push_back(yo + sy - h);
          }
          else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
          {
            double h = sy * (contour_m - f3) / (fo - f3);
            poly_x.push_back(xo);
            poly_y.push_back(yo + sy - h);
          }
        }
        else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
        {
          double hm = sy * (contour_m - f3) / (fo - f3);
          double h = hm;
          if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
          {
            h = sy * (contour - f3) / (fo - f3);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo);
          poly_y.push_back(yo + sy - hm);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo);
            poly_y.push_back(yo + sy - h);
          }
        }
        else if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
        {
          double h = sy * (contour - f3) / (fo - f3);
          poly_x.push_back(xo);
          poly_y.push_back(yo + sy - h);
        }


        size_t num_corners(poly_x.size());
        if (num_corners > 1)
        {
          double gs = 1.0 - icont * 1.0 / num_isocontours;
          file << "newpath" << std::endl;
          file << poly_x[0] * scale << " " << poly_y[0] * scale << " moveto" << std::endl;
          //std::cout<< "*** x = " <<  poly_x[0] << ", y = " << poly_y[0] << std::endl;
          for (size_t icorner = 1; icorner < num_corners; icorner++)
          {
            file << poly_x[icorner] * scale << " " << poly_y[icorner] * scale << " lineto" << std::endl;
            //std::cout << "*** x = " <<  poly_x[icorner] << ", y = " << poly_y[icorner] << std::endl;
          }
          //std::cout << std::endl;

          file << "closepath" << std::endl;
          file << "gsave" << std::endl;
          file << "grestore" << std::endl;

          double r, g, b;

          if (gs < 0.25)     r = 1.0;
          //else if (gs < 0.5) r = 2.0 - 4.0 * gs;
          else if (gs < 0.5) r = 1.0 - 16.0 * (gs - 0.25) * (gs - 0.25);
          else               r = 0.0;

          double go(0.25), gn(1.0 - go);
          if (gs < go)      g = gs / go;
          else if (gs < gn) g = 1.0;
          else              g = 1.0 / (1.0 - gn) - gs / (1.0 - gn);


          if (gs < 0.5)       b = 0.0;
          else if (gs < 0.75) b = 1.0 - 16.0 * (gs - 0.75) * (gs - 0.75);
          else                b = 1.0;

          file << r << " " << g << " " << b << " setrgbcolor" << std::endl;

          file << " fill" << std::endl;
        }
      }
    }
  }
  delete[] xx; delete[] f;

  if (plot_graph)
  {
    for (size_t i = 0; i < _num_samples; i++)
    {
      for (size_t j = 0; j < _seed_neighbors[i][1]; j++)
      {

        if (_seed_disc_neighbors != 0 && _seed_disc_neighbors[i] != 0)
        {
          if (_seed_disc_neighbors[i][function_index][j]) continue;
        }

        size_t neighbor = _seed_neighbors[i][2 + j];
        file << "newpath" << std::endl;
        // draw a line between isample and neighbor
        file << _x[i][0] * scale << "  " << _x[i][1] * scale << "  ";
        file << "moveto" << std::endl;
        file << _x[neighbor][0] * scale << "  " << _x[neighbor][1] * scale << "  ";
        file << "lineto" << std::endl;
        file << "0.0 setlinewidth" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl; // Element borders
        file << "stroke " << std::endl;
      }
    }

    for (size_t i = 0; i < _num_samples; i++)
    {
      file << "newpath" << std::endl;
      file << _x[i][0] * scale << " " << _x[i][1] * scale << " " << 0.005 * scale << " 0 360 arc" << std::endl;
      file << "closepath" << std::endl;
      file << "gsave" << std::endl;
      file << "0 0 0 setrgbcolor" << std::endl; // Inactive seeds
      file << "fill" << std::endl;
      file << "grestore" << std::endl;
      file << "0.0 setlinewidth" << std::endl;
      file << "0 0 0 setrgbcolor" << std::endl; // discs borders
      file << "stroke " << std::endl;

      std::stringstream ss;
      ss << i;
      std::string str = ss.str();

      file << "newpath " << _x[i][0] * scale << " " << _x[i][1] * scale << " moveto (" << str << ") show" << std::endl;
    }

    for (size_t i = 0; i < _num_vs; i++)
    {
      if (_vs[i] == 0) continue;

      file << "newpath" << std::endl;
      file << _vs[i][0] * scale << " " << _vs[i][1] * scale << " " << 0.005 * scale << " 0 360 arc" << std::endl;
      file << "closepath" << std::endl;
      file << "gsave" << std::endl;
      file << "1 0 0 setrgbcolor" << std::endl; // Inactive seeds
      file << "fill" << std::endl;
      file << "grestore" << std::endl;
      file << "0.0 setlinewidth" << std::endl;
      file << "1 0 0 setrgbcolor" << std::endl; // discs borders
      file << "stroke " << std::endl;

      std::stringstream ss;
      ss << i;
      std::string str = ss.str();

      file << "newpath " << _vs[i][0] * scale << " " << _vs[i][1] * scale << " moveto (" << str << ") show" << std::endl;
    }
  }

  double DX = _xmax[0] - _xmin[0];
  double DY = _xmax[1] - _xmin[1];

  // plot domain boundaries
  file << (_xmin[0] - DX) * scale << "  " << _xmin[1] * scale << "  ";
  file << (_xmax[0] + DX) * scale << "  " << _xmin[1] * scale << "  ";
  file << (_xmax[0] + DX) * scale << "  " << (_xmin[1] - DY) * scale << "  ";
  file << (_xmin[0] - DX) * scale << "  " << (_xmin[1] - DY) * scale << "  ";
  file << "quad_white" << std::endl;

  file << (_xmin[0] - DX) * scale << "  " << _xmax[1] * scale << "  ";
  file << (_xmax[0] + DX) * scale << "  " << _xmax[1] * scale << "  ";
  file << (_xmax[0] + DX) * scale << "  " << (_xmax[1] + DY) * scale << "  ";
  file << (_xmin[0] - DX) * scale << "  " << (_xmax[1] + DY) * scale << "  ";
  file << "quad_white" << std::endl;

  file << _xmax[0] * scale << "  " << (_xmin[1] - DY) * scale << "  ";
  file << (_xmax[0] + DX) * scale << "  " << (_xmin[1] - DY) * scale << "  ";
  file << (_xmax[0] + DX) * scale << "  " << (_xmax[1] + DY) * scale << "  ";
  file << _xmax[0] * scale << "  " << (_xmax[1] + DY) * scale << "  ";
  file << "quad_white" << std::endl;

  file << (_xmin[0] - DX) * scale << "  " << (_xmin[1] - DY) * scale << "  ";
  file << _xmin[0] * scale << "  " << (_xmin[1] - DY) * scale << "  ";
  file << _xmin[0] * scale << "  " << (_xmax[1] + DY) * scale << "  ";
  file << (_xmin[0] - DX) * scale << "  " << (_xmax[1] + DY) * scale << "  ";
  file << "quad_white" << std::endl;


  // plot domain boundaries
  file << _xmin[0] * scale << "  " << _xmin[1] * scale << "  ";
  file << _xmax[0] * scale << "  " << _xmin[1] * scale << "  ";
  file << _xmax[0] * scale << "  " << _xmax[1] * scale << "  ";
  file << _xmin[0] * scale << "  " << _xmax[1] * scale << "  ";
  file << "quad_bold" << std::endl;

  file << "showpage" << std::endl;

//#pragma endregion
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Delaunay Graph Methods      ///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int VPS::construct_delaunay_graph()
{
//#pragma region Construct Delaunay Graph:

  if (_seed_neighbors == 0)
  {
    _seed_neighbors = new size_t*[_budget];
    _seed_box = new double*[_budget];
    _seed_rf = new double[_budget];
    _seed_rc = new double[_budget];
    for (size_t i = 0; i < _budget; i++)
    {
      _seed_neighbors[i] = 0;
      _seed_box[i] = 0;
      _seed_rf[i] = 0; // ETP
      _seed_rc[i] = 0; // ETP
    }
  }

  for (size_t i = 0; i < _num_samples; i++)
  {
    update_delaunay_graph(i);
  }

  if (_seed_disc_neighbors == 0)
  {
    _seed_disc_neighbors = new bool**[_budget];

    for (size_t i = 0; i < _budget; i++) _seed_disc_neighbors[i] = 0;
  }

  for (size_t iseed = 0; iseed < _num_samples; iseed++)
  {
    size_t num_neighbors; get_num_seed_neighbors(iseed, num_neighbors);

    if (num_neighbors == 0) continue;

    if (_seed_disc_neighbors[iseed] != 0)
    {
      for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) delete[] _seed_disc_neighbors[iseed][ifunc];
      delete[] _seed_disc_neighbors[iseed];
    }
    _seed_disc_neighbors[iseed] = new bool*[_num_functions];
    for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
    {
      _seed_disc_neighbors[iseed][ifunc] = new bool[num_neighbors];
      for (size_t i = 0; i < num_neighbors; i++) _seed_disc_neighbors[iseed][ifunc][i] = false;
    }
  }

  return 0;

//#pragma endregion
}

int VPS::update_delaunay_graph(size_t seed_index)
{
//#pragma region Update Delaunay Graph of a given Seed:
  size_t max_num_misses(5);
  if (_num_samples > 19000) std::cout << "Del graph:" << seed_index << std::endl;
  if (_seed_neighbors[seed_index] != 0)
  {
    while (_seed_neighbors[seed_index][1] != 0) disconnect_seeds(seed_index, _seed_neighbors[seed_index][2]);
    delete[] _seed_neighbors[seed_index];
    _seed_neighbors[seed_index] = 0;
  }

  if (_seed_box[seed_index] == 0) _seed_box[seed_index] = new double[2 * _num_dim];

  _seed_rc[seed_index] = 0;
  _seed_rf[seed_index] = _diag;
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    _seed_box[seed_index][idim] = _x[seed_index][idim];
    _seed_box[seed_index][_num_dim + idim] = _x[seed_index][idim];
  }

  double* v = new double[_num_dim];

  size_t num_misses(0), old_neighbors(0);
  get_num_seed_neighbors(seed_index, old_neighbors);
  while (num_misses < max_num_misses)
  {
    sample_voronoi_vertex(seed_index, _xmin, _xmax, _diag, v);

    // update Coverage radius
    double rc(0.0);
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = v[idim] - _x[seed_index][idim];
      rc += dx * dx;
    }
    rc = sqrt(rc);
    if (rc > _seed_rc[seed_index]) _seed_rc[seed_index] = rc;

    // update bounding box for Voronoi Cell
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      if (v[idim] < _seed_box[seed_index][idim]) _seed_box[seed_index][idim] = v[idim];
      if (v[idim] > _seed_box[seed_index][_num_dim + idim]) _seed_box[seed_index][_num_dim + idim] = v[idim];
    }

    // Collect all seeds that are equidistant from xend
    size_t num_neighbors(0);
    size_t neighbors_capacity(100);
    size_t* neighbors = new size_t[neighbors_capacity];

    get_seeds_in_sphere_tree(v, rc + 1E-10, num_neighbors, neighbors, neighbors_capacity);

    // Sorting ensures
    quicksort(neighbors, 0, num_neighbors - 1);
    for (size_t i = 0; i < num_neighbors; i++)
    {
      for (size_t j = i + 1; j < num_neighbors; j++)
      {
        connect_seeds(neighbors[i], neighbors[j]);
      }
    }

    delete[] neighbors;

    size_t new_neighbors;
    get_num_seed_neighbors(seed_index, new_neighbors);

    if (new_neighbors == old_neighbors) num_misses++;
    else
    {
      old_neighbors = new_neighbors;
      num_misses = 0;
    }
  }
  delete[] v;
  return 0;
//#pragma endregion
}

int VPS::sample_voronoi_vertex(double* v)
{
//#pragma region Sample Voronoi Vertex:
  double* dart = new double[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++) dart[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);

  double closest_distance = _diag; size_t closest_seed;
  get_closest_seed_tree(dart, closest_seed, closest_distance);

  sample_voronoi_vertex(closest_seed, _xmin, _xmax, _diag, v);

  delete[] dart;
  return 0;
//#pragma endregion
}

int VPS::connect_seeds(size_t seed_i, size_t seed_j)
{
//#pragma region Connect Seeds:

  if (seed_i == seed_j) return 1;

  double dst(0.0);
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    double dx = _x[seed_i][idim] - _x[seed_j][idim];
    dst += dx * dx;
  }
  dst = sqrt(dst);
  if (dst < _seed_rf[seed_i])  _seed_rf[seed_i] = dst;
  if (dst < _seed_rf[seed_j])  _seed_rf[seed_j] = dst;

  if (_seed_neighbors[seed_i] == 0)
  {
    _seed_rf[seed_i] = dst;
    _seed_neighbors[seed_i] = new size_t[12];
    _seed_neighbors[seed_i][0] = 12;
    _seed_neighbors[seed_i][1] = 1;
    _seed_neighbors[seed_i][2] = seed_j;
    if (seed_i < seed_j) connect_seeds(seed_j, seed_i);
    return 0;
  }
  for (size_t i = 0; i < _seed_neighbors[seed_i][1]; i++)
  {
    if (_seed_neighbors[seed_i][2 + i] == seed_j) return 0; // two seeds are already connected
  }

  size_t num_neighbors = _seed_neighbors[seed_i][1];
  _seed_neighbors[seed_i][num_neighbors + 2] = seed_j;
  _seed_neighbors[seed_i][1]++;

  if (_seed_neighbors[seed_i][1] + 2 == _seed_neighbors[seed_i][0])
  {
    size_t new_capacity = _seed_neighbors[seed_i][0] * 2;
    size_t* tmp = new size_t[new_capacity];
    tmp[0] = new_capacity; tmp[1] = _seed_neighbors[seed_i][1];
    for (size_t i = 0; i < _seed_neighbors[seed_i][1]; i++) tmp[2 + i] = _seed_neighbors[seed_i][2 + i];
    delete[] _seed_neighbors[seed_i];
    _seed_neighbors[seed_i] = tmp;
  }
  if (seed_i < seed_j) connect_seeds(seed_j, seed_i);
  return 0;
//#pragma endregion
}

int VPS::disconnect_seeds(size_t seed_i, size_t seed_j)
{
//#pragma region Disconnect Seeds:

  if (seed_i == seed_j) return 1;

  if (_seed_neighbors[seed_i] == 0) return 1;

  for (size_t i = 0; i < _seed_neighbors[seed_i][1]; i++)
  {
    if (_seed_neighbors[seed_i][2 + i] == seed_j)
    {
      _seed_neighbors[seed_i][2 + i] = _seed_neighbors[seed_i][1 + _seed_neighbors[seed_i][1]];
      _seed_neighbors[seed_i][1]--;
      break;
    }
  }

  for (size_t i = 0; i < _seed_neighbors[seed_j][1]; i++)
  {
    if (_seed_neighbors[seed_j][2 + i] == seed_i)
    {
      _seed_neighbors[seed_j][2 + i] = _seed_neighbors[seed_j][1 + _seed_neighbors[seed_j][1]];
      _seed_neighbors[seed_j][1]--;
      break;
    }
  }
  return 0;
//#pragma endregion
}

int VPS::get_num_seed_neighbors(size_t seed_index, size_t &num_seed_neighbors)
{
//#pragma region Get Number od Delaunay Edges:
  if (_seed_neighbors[seed_index] == 0) num_seed_neighbors = 0;
  else                                  num_seed_neighbors = _seed_neighbors[seed_index][1];
  return 0;
//#pragma endregion
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Least Square QR Solver Methods           //////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int VPS::LeastSquare_QR(size_t num_data_points, double xmin, double xmax, double* x, double* f, size_t num_basis, double* c)
{
//#pragma region Least Square QR graph fit:
  double** A = new double*[num_data_points];

  for (size_t idatapoint = 0; idatapoint < num_data_points; idatapoint++)
  {
    A[idatapoint] = new double[num_basis];
    for (size_t ibasis = 0; ibasis < num_basis; ibasis++) A[idatapoint][ibasis] = pow(x[idatapoint], ibasis);
  }

  LS_QR_Solver(num_data_points, num_basis, A, f, c);

  for (size_t idatapoint = 0; idatapoint < num_data_points; idatapoint++) delete[] A[idatapoint];
  delete[] A;

  return 0;
//#pragma endregion
}

int VPS::NaturalCubicSplineInterpolation(size_t num_data_points, double* x, double* y, double* co, double* c1, double* c2, double* c3)
{
//#pragma region Natural Cubic Spline Interpolation:
  // Second Derivative at start points of each quadratic Curve
  size_t n = num_data_points - 1;

  if (n == 1)
  {
    // linear interpolation
    co[0] = y[0] - x[0] * (y[1] - y[0]) / (x[1] - x[0]);
    c1[0] = (y[1] - y[0]) / (x[1] - x[0]);
    c2[0] = 0.0;
    c3[0] = 0.0;
    return 0;
  }

  double* G = new double[n + 1];
  G[0] = 0.0; G[n] = 0.0;

  double* LD = new double[n - 1];
  double* DD = new double[n - 1];
  double* UD = new double[n - 1];
  double* b = new double[n - 1];
  double dm = x[1] - x[0]; double dp = x[2] - x[1];
  double Dm = y[1] - y[0]; double Dp = y[2] - y[1];
  DD[0] = 2 * (dm + dp); UD[0] = dp; b[0] = 6 * (Dp / dp - Dm / dm);

  dp = x[n] - x[n - 1]; dm = x[n - 1] - x[n - 2];
  Dp = y[n] - y[n - 1]; Dm = y[n - 1] - y[n - 2];
  LD[n - 2] = dm; DD[n - 2] = 2 * (dm + dp); b[n - 2] = 6 * (Dp / dp - Dm / dm);
  if (n > 2)
  {
    for (size_t i = 1; i < n - 2; i++)
    {
      dp = x[i + 2] - x[i + 1]; dm = x[i + 1] - x[i];
      Dp = y[i + 2] - y[i + 1]; Dm = y[i + 1] - y[i];
      LD[i] = dm;
      DD[i] = 2 * (dm + dp);
      UD[i] = dp;
      b[i] = 6 * (Dp / dp - Dm / dm);
    }
  }

  double* kk = new double[n - 1];

  solve_LU_Doolittle(n - 1, LD, DD, UD, b, kk);
  for (size_t i = 0; i < n - 1; i++) G[i + 1] = kk[i];

  for (size_t i = 0; i < num_data_points - 1; i++)
  {
    dp = x[i + 1] - x[i]; Dp = y[i + 1] - y[i];
    double ki = (Dp / dp) - ((2.0 / 3.0) * dp * 0.5 * G[i]) - (1.0 / 3.0) * dp * 0.5 * G[i + 1];
    double di = 0.5 * (G[i + 1] - G[i]) / (3.0 * dp);

    co[i] = y[i] - ki * x[i] + 0.5 * G[i] * pow(x[i], 2) - di * pow(x[i], 3);

    c1[i] = ki - G[i] * x[i] + 3 * di * pow(x[i], 2);;

    c2[i] = 0.5 * G[i] - 3 * di * x[i];

    c3[i] = di;

    //yp[j] = y[i] + ki * (xp[j] - x[i]) + 0.5 * G[i] * pow((xp[j] - x[i]), 2) + di * pow((xp[j] - x[i]), 3);
  }
  delete[] G; delete[] LD; delete[] DD; delete[] UD; delete[] b; delete[] kk;
  return 0;
//#pragma endregion
}

int VPS::FourierExpansion(size_t num_data_points, double xmin, double xmax, double* x, double* f, double** c, size_t num_basis, double* a, double* b)
{
//#pragma region Fourier Expansion of a function:

  //size_t num_pieces(num_data_points - 1); // ETP

  double P = xmax - xmin;

  for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
  {
    a[ibasis] = 0.0; b[ibasis] = 0.0;

    size_t num_pieces(num_data_points - 1);
    for (size_t ipiece = 0; ipiece < num_pieces; ipiece++)
    {
      double xm = x[ipiece]; double xp = x[ipiece + 1];

      if (ibasis == 0)
      {
        a[ibasis] += c[ipiece][0] * (xp - xm);
        a[ibasis] += c[ipiece][1] * (pow(xp, 2) - pow(xm, 2)) / 2;
        a[ibasis] += c[ipiece][2] * (pow(xp, 3) - pow(xm, 3)) / 3;
        a[ibasis] += c[ipiece][3] * (pow(xp, 4) - pow(xm, 4)) / 4;
        continue;
      }

      double tm = 2.0 * PI * ibasis * xm / P; double tp = 2.0 * PI * ibasis * xp / P;

      a[ibasis] += c[ipiece][0] * P * (sin(tp) - sin(tm)) / (2.0 * PI * ibasis);

      a[ibasis] += c[ipiece][1] * pow(P, 2) * (cos(tp) - cos(tm)) / (4 * pow(PI, 2) * pow(ibasis, 2));
      a[ibasis] += c[ipiece][1] * P * (xp * sin(tp) - xm * sin(tm)) / (2.0 * PI * ibasis);

      a[ibasis] += c[ipiece][2] * pow(P, 2) * (xp * cos(tp) - xm * cos(tm)) / (2.0 * pow(PI, 2) * pow(ibasis, 2));
      a[ibasis] -= c[ipiece][2] * pow(P, 3) * (sin(tp) - sin(tm)) / (4.0 * pow(PI, 3) * pow(ibasis, 3));
      a[ibasis] += c[ipiece][2] * P * (pow(xp, 2) * sin(tp) - pow(xm, 2) * sin(tm)) / (2.0 * PI * ibasis);

      a[ibasis] += c[ipiece][3] * P * (pow(xp, 3) * sin(tp) - pow(xm, 3) * sin(tm)) / (2.0 * PI * ibasis);
      a[ibasis] -= c[ipiece][3] * 3 * pow(P, 3) * (xp * sin(tp) - xm * sin(tm)) / (4.0 * pow(PI, 3) * pow(ibasis, 3));
      a[ibasis] -= c[ipiece][3] * 3 * pow(P, 4) * (cos(tp) - cos(tm)) / (8.0 * pow(PI, 4) * pow(ibasis, 4));
      a[ibasis] += c[ipiece][3] * 3 * pow(P, 2) * (pow(xp, 2) * cos(tp) - pow(xm, 2) * cos(tm)) / (4.0 * pow(PI, 2) * pow(ibasis, 2));



      b[ibasis] -= c[ipiece][0] * P * (cos(tp) - cos(tm)) / (2.0 * PI * ibasis);

      b[ibasis] += c[ipiece][1] * pow(P, 2) * (sin(tp) - sin(tm)) / (4 * pow(PI, 2) * pow(ibasis, 2));
      b[ibasis] -= c[ipiece][1] * P * (xp * cos(tp) - xm * cos(tm)) / (2.0 * PI * ibasis);


      b[ibasis] += c[ipiece][2] * pow(P, 3) * (cos(tp) - cos(tm)) / (4.0 * pow(PI, 3) * pow(ibasis, 3));
      b[ibasis] -= c[ipiece][2] * P * (pow(xp, 2) * cos(tp) - pow(xm, 2) * cos(tm)) / (2.0 * PI * ibasis);
      b[ibasis] += c[ipiece][2] * pow(P, 2) * (xp * sin(tp) - xm * sin(tm)) / (2.0 * pow(PI, 2) * pow(ibasis, 2));

      b[ibasis] += c[ipiece][3] * 3 * pow(P, 3) * (xp * cos(tp) - xm * cos(tm)) / (4.0 * pow(PI, 3) * pow(ibasis, 3));
      b[ibasis] -= c[ipiece][3] * P * (pow(xp, 3) * cos(tp) - pow(xm, 3) * cos(tm)) / (2.0 * PI * ibasis);
      b[ibasis] -= c[ipiece][3] * 3 * pow(P, 4) * (sin(tp) - sin(tm)) / (8.0 * pow(PI, 4) * pow(ibasis, 4));
      b[ibasis] += c[ipiece][3] * 3 * pow(P, 2) * (pow(xp, 2) * sin(tp) - pow(xm, 2) * sin(tm)) / (4.0 * pow(PI, 2) * pow(ibasis, 2));

    }
    a[ibasis] *= (2.0 / P); b[ibasis] *= (2.0 / P);
  }
  return 0;
//#pragma endregion
}


int VPS::solve_LU_Doolittle(size_t num_eq, double* LD, double* DD, double* UD, double* b, double* x)
{
//#pragma region Solve LU Doolittle:
  if (num_eq == 0) return 1;

  // L and U Matrices
  double* L = new double[num_eq];
  double** U = new double*[num_eq];
  for (size_t i = 0; i < num_eq; i++) U[i] = new double[2];

  U[0][0] = DD[0];
  for (size_t i = 0; i < num_eq - 1; i++) U[i][1] = UD[i];

  for (size_t i = 1; i < num_eq; i++)
  {
    L[i] = LD[i] / U[i - 1][0];
    U[i][0] = DD[i] - L[i] * U[i - 1][1];
  }

  // Solution Vectors:
  double* y = new double[num_eq];
  y[0] = b[0];
  for (size_t i = 1; i < num_eq; i++)
  {
    y[i] = b[i] - L[i] * y[i - 1];
  }
  // Back step substitution:
  x[num_eq - 1] = y[num_eq - 1] / U[num_eq - 1][0];
  if (num_eq > 1)
  {
    size_t i(num_eq - 2);
    while (true)
    {
      x[i] = (y[i] - U[i][1] * x[i + 1]) / U[i][0];
      if (i == 0) break;
      i--;
    }
  }
  for (size_t i = 0; i < num_eq; i++) delete[] U[i];
  delete[] L; delete[] U; delete[] y;
  return 0;
//#pragma endregion
}

void VPS::LS_QR_Solver(size_t nrow, size_t ncol, double** A, double* b, double* x)
{
//#pragma region Least Square QR Solver:
  double** Q = new double*[nrow];
  for (size_t irow = 0; irow < nrow; irow++) Q[irow] = new double[ncol];

  double** R = new double*[ncol]; // A square matrix
  for (size_t irow = 0; irow < ncol; irow++) R[irow] = new double[ncol];

  QR_Factorization(nrow, ncol, A, Q, R);

  double* bmod = new double[ncol];
  for (size_t irow = 0; irow < ncol; irow++)
  {
    bmod[irow] = 0.0;
    for (size_t icol = 0; icol < nrow; icol++) bmod[irow] += Q[icol][irow] * b[icol];
  }

  // backward solve
  size_t jrow = ncol - 1;
  while (true)
  {
    if (fabs(R[jrow][jrow]) > 1E-8)
    {
      double sum(0.0);
      for (size_t jcol = ncol - 1; jcol > jrow; jcol--) sum += R[jrow][jcol] * x[jcol];
      x[jrow] = (bmod[jrow] - sum) / R[jrow][jrow];
    }
    else
    {
      // This column is in the Null Space of A
      x[jrow] = 0.0;
    }

    if (jrow == 0) break;
    jrow--;
  }

  for (size_t irow = 0; irow < nrow; irow++) delete[] Q[irow];
  delete[] Q;

  for (size_t irow = 0; irow < ncol; irow++) delete[] R[irow];
  delete[] R;

  delete[] bmod;
//#pragma endregion
}

void VPS::GramSchmidt(size_t nrow, size_t ncol, double** A, double** Q)
{
//#pragma region GramSchmidt Process:

  // norw >= ncol (for LS problems)
  bool* null_basis = new bool[ncol];
  for (size_t ibasis = 0; ibasis < ncol; ibasis++) null_basis[ibasis] = false;

  double** ebasis = new double*[ncol];
  for (size_t ibasis = 0; ibasis < ncol; ibasis++) ebasis[ibasis] = new double[nrow];

  double* v = new double[nrow];

  size_t num_basis(0);

  for (size_t icol = 0; icol < ncol; icol++)
  {
    for (size_t irow = 0; irow < nrow; irow++) ebasis[num_basis][irow] = A[irow][icol];

    // normalize
    double norm(0.0);
    for (size_t irow = 0; irow < nrow; irow++) norm += ebasis[num_basis][irow] * ebasis[num_basis][irow];
    norm = sqrt(norm);

    if (norm > 1E-8)
    {
      for (size_t irow = 0; irow < nrow; irow++)  ebasis[num_basis][irow] /= norm;
    }
    else
    {
      null_basis[num_basis] = true;
      num_basis++;
      continue;
    }

    orthonormalize_vector(num_basis, nrow, ebasis, ebasis[num_basis], norm);

    if (norm <= 1E-8)
    {
      null_basis[num_basis] = true;
    }
    num_basis++;
  }

  // Replacing Null basis with basis from Null Space of A
  for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
  {
    if (!null_basis[ibasis]) continue;
    for (size_t idim = 0; idim < nrow; idim++)
    {
      for (size_t jdim = 0; jdim < nrow; jdim++) v[jdim] = 0.0;
      v[idim] = 1.0;
      double norm;
      orthonormalize_vector(num_basis, nrow, ebasis, v, norm);

      if (norm <= 1E-8) continue;

      for (size_t jdim = 0; jdim < nrow; jdim++) ebasis[ibasis][jdim] = v[jdim];
      break;
    }
  }

  MAT_Transpose(num_basis, nrow, ebasis, Q);

  for (size_t ibasis = 0; ibasis < ncol; ibasis++) delete[] ebasis[ibasis];
  delete[] ebasis; delete[] null_basis; delete[] v;

//#pragma endregion
}


void VPS::orthonormalize_vector(size_t nbasis, size_t ndim, double** e_basis, double* v, double &norm)
{
//#pragma region orthonomalize Vector:
  for (size_t ibasis = 0; ibasis < nbasis; ibasis++)
  {
    // project
    double dot(0.0);
    for (size_t idim = 0; idim < ndim; idim++) dot += v[idim] * e_basis[ibasis][idim];
    for (size_t idim = 0; idim < ndim; idim++) v[idim] -= dot * e_basis[ibasis][idim];

    // normalize
    norm = 0.0;
    for (size_t idim = 0; idim < ndim; idim++) norm += v[idim] * v[idim];
    norm = sqrt(norm);

    if (norm > 1E-8)
    {
      for (size_t idim = 0; idim < ndim; idim++) v[idim] /= norm;
    }
    else
    {
      norm = 0.0;
      return;
    }
  }
//#pragma endregion
}


void VPS::QR_Factorization(size_t nrow, size_t ncol, double** A, double** Q, double** R)
{
//#pragma region QR Factorization:
  // norw >= ncol (for LS problems)
  GramSchmidt(nrow, ncol, A, Q);
  double** QT = new double*[ncol];
  for (size_t irow = 0; irow < ncol; irow++) QT[irow] = new double[nrow];

  MAT_Transpose(nrow, ncol, Q, QT);

  MAT_MAT(ncol, nrow, QT, ncol, A, R);

  for (size_t irow = 0; irow < ncol; irow++) delete[] QT[irow];
  delete[] QT;
//#pragma endregion
}

void VPS::MAT_MAT(size_t nrow, size_t ncol, double** A, size_t mcol, double** B, double** C)
{
//#pragma region Matrix Matrix Multiplication:
  for (size_t irow = 0; irow < nrow; irow++)
  {
    for (size_t jcol = 0; jcol < mcol; jcol++)
    {
      C[irow][jcol] = 0.0;
      for (size_t icol = 0; icol < ncol; icol++)
      {
        C[irow][jcol] += A[irow][icol] * B[icol][jcol];
      }
    }
  }
//#pragma endregion
}

void VPS::MAT_Transpose(size_t nrow, size_t ncol, double** A, double** AT)
{
//#pragma region Matrix Transpose:
  for (size_t irow = 0; irow < nrow; irow++)
  {
    for (size_t icol = 0; icol < ncol; icol++)
    {
      AT[icol][irow] = A[irow][icol];
    }
  }
//#pragma endregion
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Surrogate Methods           ///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int VPS::construct_local_surrogates()
{
  init_vps_containers();
  for (size_t isample = 0; isample < _num_samples; isample++)
  {
    if (_num_samples > 19000) std::cout << "Building Surrogate:" << isample << std::endl;
    if (_method == Regression) retrieve_weights_regression(isample);
  }
  return 0;
}


int VPS::init_vps_containers()
{
//#pragma region init_vps_containers:
  if (_p == 0)
  {
    _num_basis = retrieve_num_permutations(_num_dim, _desired_order, true, _desired_order);
    _num_basis_pool = _num_basis; size_t order = _desired_order;
    while (_num_basis_pool < 10 * _budget)
    {
      order++;
      _num_basis_pool = retrieve_num_permutations(_num_dim, order, true, order);
    }
    retrieve_permutations(_num_basis_pool, _p, _num_dim, order, true, order);
  }

  if (_basis_coef == 0)
  {
    // init Surrogate container
    _basis_coef = new double**[_budget];
    _basis_index = new size_t**[_budget];
    for (size_t i = 0; i < _budget; i++)
    {
      _basis_coef[i] = 0; _basis_index[i] = 0;
    }
  }
  return 0;
//#pragma endregion
}

int VPS::add_neighbor_layer(size_t cell_index, size_t function_index, size_t &num_neighbors, size_t* &neighbors, size_t &neighbors_capacity)
{
//#pragma region Add neighbor Layer:

  if (_seed_neighbors == 0 || _seed_neighbors[cell_index] == 0) return 0;

  size_t num_new_neighbors(0), new_neighbors_capacity(100);
  size_t* new_neighbors = new size_t[new_neighbors_capacity];

  if (num_neighbors == 1)
  {
    // add delaunay neighbors
    for (size_t j = 0; j < _seed_neighbors[cell_index][1]; j++)
    {
      size_t potential_neighbor = _seed_neighbors[cell_index][2 + j];

      bool disc(false);
      if (_seed_disc_neighbors != 0 && _seed_disc_neighbors[cell_index] != 0)
      {
        disc = _seed_disc_neighbors[cell_index][function_index][j];
      }
      if (disc) continue;
      add_entry(potential_neighbor, num_new_neighbors, new_neighbors, new_neighbors_capacity);
    }
  }
  else
  {
    // add neighbors of existing neighbors
    for (size_t ineighbor = 0; ineighbor < num_neighbors; ineighbor++)
    {
      size_t neighbor = neighbors[ineighbor];
      for (size_t j = 0; j < _seed_neighbors[neighbor][1]; j++)
      {
        size_t potential_neighbor = _seed_neighbors[neighbor][2 + j];

        bool disc(false);
        if (_seed_disc_neighbors != 0 && _seed_disc_neighbors[neighbor] != 0)
        {
          disc = _seed_disc_neighbors[neighbor][function_index][j];
        }
        if (disc) continue;

        if (potential_neighbor == cell_index) continue;
        if (find_brute(potential_neighbor, neighbors, num_neighbors)) continue;
        if (find_brute(potential_neighbor, new_neighbors, num_new_neighbors)) continue;

        add_entry(potential_neighbor, num_new_neighbors, new_neighbors, new_neighbors_capacity);
      }
    }
  }

  for (size_t inew = 0; inew < num_new_neighbors; inew++)
  {
    add_entry(new_neighbors[inew], num_neighbors, neighbors, neighbors_capacity);
  }
  delete[] new_neighbors;
  return 0;
//#pragma endregion
}

int VPS::retrieve_weights_regression(size_t cell_index)
{
//#pragma region Retrieve Weight Using Regression:

  if (_desired_order == 0)
  {
    _basis_coef[cell_index] = new double*[_num_functions];
    _basis_index[cell_index] = new size_t*[_num_functions];
    for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
    {
      _basis_coef[cell_index][ifunc] = new double[_num_basis];
      _basis_index[cell_index][ifunc] = new size_t[_num_basis];
      _basis_coef[cell_index][ifunc][0] = _f[cell_index][ifunc];
      _basis_index[cell_index][ifunc][0] = 0;
    }
    return 0;
  }

  if (_basis_coef[cell_index] == 0)
  {
    _basis_coef[cell_index] = new double*[_num_functions];
    _basis_index[cell_index] = new size_t*[_num_functions];
    for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
    {
      _basis_coef[cell_index][ifunc] = new double[_num_basis];
      _basis_index[cell_index][ifunc] = new size_t[_num_basis];
    }
  }

  for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
  {
    for (size_t ibasis = 0; ibasis < _num_basis; ibasis++)
    {
      _basis_coef[cell_index][ifunc][ibasis] = 0.0;
      _basis_index[cell_index][ifunc][ibasis] = ibasis;
    }
  }



  for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
  {
    // collect neighbors
    size_t num_data_points(1), data_points_capacity(10);
    size_t* data_points = new size_t[data_points_capacity];
    data_points[0] = cell_index;
    while (_desired_order > 0 && num_data_points < 2 * _num_basis)
    {
      size_t old_num = num_data_points;
      add_neighbor_layer(cell_index, ifunc, num_data_points, data_points, data_points_capacity);
      if (old_num == num_data_points) break; // no more neighbors to add
    }

    size_t num_basis(_num_basis), order(_desired_order);
    while (num_data_points < num_basis)
    {
      order--;
      num_basis = retrieve_num_permutations(_num_dim, order, true, order);
    }

    // Select basis that is less than or equal to order
    size_t nb(0);
    for (size_t ibasis = 0; ibasis < _num_basis; ibasis++)
    {
      size_t sum(0);
      for (size_t idim = 0; idim < _num_dim; idim++) sum += _p[ibasis][idim];
      if (sum > order) continue;

      for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) _basis_index[cell_index][ifunc][nb] = ibasis;
      nb++;
    }




    double** A = new double*[num_data_points];
    double* b = new double[num_data_points];
    for (size_t ipoint = 0; ipoint < num_data_points; ipoint++)
    {
      size_t data_point_index = data_points[ipoint];

      A[ipoint] = new double[num_basis];

      for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
      {
        size_t basis_index = _basis_index[cell_index][ifunc][ibasis];
        A[ipoint][ibasis] = evaluate_basis_function(_x[data_point_index], cell_index, basis_index);
      }
      b[ipoint] = _f[data_point_index][ifunc];
    }

    bool constrained_regression(true);
    if (constrained_regression)
    {
      for (size_t ipoint = 0; ipoint < num_data_points; ipoint++)
      {
        size_t data_point_index = data_points[ipoint];
        if (data_point_index == cell_index) continue;
        A[ipoint][0] = 0.0;
        b[ipoint] -= _f[cell_index][ifunc];
      }
    }

    LS_QR_Solver(num_data_points, /*_*/num_basis, A, b, _basis_coef[cell_index][ifunc]); // ETP


    delete[] data_points;

    for (size_t ipoint = 0; ipoint < num_data_points; ipoint++) delete[] A[ipoint];
    delete[] A; delete[] b;

  }
  return 0;
//#pragma endregion
}

double VPS::evaluate_basis_function(double* x, size_t cell_index, size_t ibasis)
{
//#pragma region Evaluate Basis Function:
  if (_basis == monomials || _basis == legendre || _basis == chebyshev)
  {
    double f_basis(1.0);
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      f_basis *= evaluate_one_dimensional_basis_function(x[idim], _x[cell_index][idim], _seed_box[cell_index][idim], _seed_box[cell_index][_num_dim + idim], _p[ibasis][idim]);
    }
    return f_basis;
  }
  else if (_basis == radial)
  {
    /*
      size_t basis_index = icell;

      double dst_sq = 0.0;
      for (size_t idim = 0; idim < _num_dim; idim++)
      {
      double dx = x[idim] - _sample_basis[icell][ibasis][idim];
      dst_sq += dx * dx;
      }
      h = _seed_rc[cell_index];
      double r_basis = 4.0 * h;
      double r_sq = r_basis * r_basis;

      return exp(-dst_sq / r_sq);
    */
  }
  return 0.0;
//#pragma endregion
}

double VPS::evaluate_one_dimensional_basis_function(double x, double xo, double xmin, double xmax, size_t ibasis)
{
//#pragma region Evaluate One Dimensional Basis Function:
  if (_basis == monomials)
  {
    double fbasis = 0.0;
    x -= xo; // shift origin
    double dxmin = xo - xmin;
    double dxmax = xmax - xo;
    if (dxmin > dxmax) x /= dxmin;
    else               x /= dxmax;

    fbasis = pow(x, ibasis);
    return fbasis;
  }
  else if (_basis == legendre)
  {
    double fbasis = 0.0;

    x -= xo; // shift origin
    double dxmin = xo - xmin;
    double dxmax = xmax - xo;
    if (dxmin > dxmax) x /= dxmin;
    else               x /= dxmax;

    if (ibasis == 0) fbasis = 1.0;
    else if (ibasis == 1) fbasis = x;
    else
    {
      double fm = evaluate_one_dimensional_basis_function(x, 0.0, -1.0, 1.0, ibasis - 1);
      double fmm = evaluate_one_dimensional_basis_function(x, 0.0, -1.0, 1.0, ibasis - 2);
      fbasis = ((2 * ibasis - 1) * x * fm - (ibasis - 1) * fmm) / ibasis;
    }
    return fbasis;
  }
  else if (_basis == chebyshev)
  {
    double fbasis = 0.0;
    x -= xo; // shift origin
    double dxmin = xo - xmin;
    double dxmax = xmax - xo;
    if (dxmin > dxmax) x /= dxmin;
    else               x /= dxmax;

    if (ibasis == 0) fbasis = 1.0;
    else if (ibasis == 1) fbasis = x;
    else
    {
      double fm = evaluate_one_dimensional_basis_function(x, 0.0, -1.0, 1.0, ibasis - 1);
      double fmm = evaluate_one_dimensional_basis_function(x, 0.0, -1.0, 1.0, ibasis - 2);
      fbasis = 2 * x * fm - fmm;
    }
    return fbasis;
  }
  return 0.0;
//#pragma endregion
}


int VPS::detect_discontinuities()
{
//#pragma region Detect Discontinuities:
  for (size_t iseed = 0; iseed < _num_samples; iseed++)
  {
    size_t num_neighbors; get_num_seed_neighbors(iseed, num_neighbors);

    if (num_neighbors == 0) return 0;

    for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
    {
      for (size_t i = 0; i < num_neighbors; i++)
      {
        if (_seed_disc_neighbors[iseed][ifunc][i]) continue;
        size_t neighbor_i = _seed_neighbors[iseed][2 + i];

        for (size_t j = i + 1; j < num_neighbors; j++)
        {
          if (_seed_disc_neighbors[iseed][ifunc][j]) continue;
          size_t neighbor_j = _seed_neighbors[iseed][2 + j];

          // Asser that neighbor_i and neighbor_j lie on opposite direction with regard to iseed

          double hi(0.0), hj(0.0), dot(0.0);
          for (size_t idim = 0; idim < _num_dim; idim++)
          {
            double dxi = _x[neighbor_i][idim] - _x[iseed][idim];
            double dxj = _x[neighbor_j][idim] - _x[iseed][idim];
            hi += dxi * dxi;
            hj += dxj * dxj;
            dot += dxi * dxj;
          }
          hi = sqrt(hi); hj = sqrt(hj);
          dot /= hi; dot /= hj;
          if (dot > -0.5) continue;

          double dfi = _f[neighbor_i][ifunc] - _f[iseed][ifunc];
          double dfj = _f[neighbor_j][ifunc] - _f[iseed][ifunc];

          double norm_hfi = sqrt(hi * hi + dfi * dfi);
          double norm_hfj = sqrt(hj * hj + dfj * dfj);
          double dot_hf = -hi * hj + dfi * dfj;

          dot_hf /= norm_hfi; dot_hf /= norm_hfj;

          if (dot_hf < -0.5) continue;

          double si = fabs(dfi) / hi;
          double sj = fabs(dfj) / hj;

          if (si > sj)
          {
            _seed_disc_neighbors[iseed][ifunc][i] = true;
            size_t num_ii_neighbors; get_num_seed_neighbors(neighbor_i, num_ii_neighbors);
            for (size_t ii = 0; ii < num_ii_neighbors; ii++)
            {
              if (_seed_neighbors[neighbor_i][2 + ii] == iseed) _seed_disc_neighbors[neighbor_i][ifunc][ii] = true;
            }
          }
          else
          {
            _seed_disc_neighbors[iseed][ifunc][j] = true;
            size_t num_jj_neighbors; get_num_seed_neighbors(neighbor_j, num_jj_neighbors);
            for (size_t jj = 0; jj < num_jj_neighbors; jj++)
            {
              if (_seed_neighbors[neighbor_j][2 + jj] == iseed) _seed_disc_neighbors[neighbor_j][ifunc][jj] = true;
            }
          }
        }
      }
    }
  }
  return 0;
//#pragma endregion
}

int VPS::sample_voronoi_facet(size_t seed_index, double* xmin, double* xmax, double diag, double* v)
{

//#pragma region Sample A Voronoi facet bounding Cell with seed_index:

  size_t num_basis(0);
  double** basis = new double*[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++)
  {

    basis[idim] = new double[_num_dim];

    for (size_t jdim = 0; jdim < _num_dim; jdim++) basis[idim][jdim] = 0.0;

    basis[idim][idim] = 1.0;

  }

  double* xst = new double[_num_dim];
  double* xend = new double[_num_dim];
  double* dart = new double[_num_dim];
  size_t* simplex = new size_t[_num_dim + 1];
  double* vect = new double[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++) xst[idim] = _x[seed_index][idim];

  double closest_distance = 0.0;
  // add first point
  size_t num_simplex_seeds(0);
  simplex[num_simplex_seeds] = seed_index; num_simplex_seeds++;

  size_t num_misses(0); size_t vertex_dim(0);

  while (num_simplex_seeds <= _num_dim - vertex_dim)
  {
    // Throw a random spoke
    sample_uniformly_from_unit_sphere(dart, _num_dim - num_basis);

    for (size_t idim = 0; idim < _num_dim; idim++) xend[idim] = xst[idim];

    for (size_t ibasis = num_basis; ibasis < _num_dim; ibasis++)
    {
      for (size_t idim = 0; idim < _num_dim; idim++) xend[idim] += diag * dart[ibasis - num_basis] * basis[ibasis][idim];
    }

    closest_distance = 0.0;

    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = xend[idim] - _x[seed_index][idim];
      closest_distance += dx * dx;
    }
    closest_distance = sqrt(closest_distance);

    bool done(true);
    size_t neighbor(seed_index);
    while (true)
    {
      // get closest point to xend
      size_t closest_seed(seed_index);

      double excluded_closest_distance = closest_distance + 1E-10;

      get_closest_seed_tree(xend, num_simplex_seeds, simplex, closest_seed, excluded_closest_distance);

      bool vertex_out_of_bounding_box(false);

      if (neighbor == closest_seed)
      {
        for (size_t idim = 0; idim < _num_dim; idim++)
        {
          if (xend[idim] > xmin[idim] - 1E-10 && xend[idim] < xmax[idim] + 1E-10) continue;
          vertex_out_of_bounding_box = true;
          break;
        }
      }

      if (closest_seed == _num_samples || (neighbor == closest_seed && vertex_out_of_bounding_box))
      {
        // reset sampling
        num_simplex_seeds = 1;
        for (size_t idim = 0; idim < _num_dim; idim++) xst[idim] = _x[seed_index][idim];

        num_basis = 0;
        for (size_t idim = 0; idim < _num_dim; idim++)
        {
          for (size_t jdim = 0; jdim < _num_dim; jdim++) basis[idim][jdim] = 0.0;
          basis[idim][idim] = 1.0;
        }
        num_misses++;

        if (num_misses == 100)
        {
          vertex_dim++; num_misses = 0;
        }
        done = false;
        break;
      }

      if (neighbor == closest_seed || fabs(excluded_closest_distance - closest_distance) < 1E-10)
      {
//#pragma region A voronoi neighbor apply recursion:
        simplex[num_simplex_seeds] = neighbor; num_simplex_seeds++;

        // update spoke start
        for (size_t idim = 0; idim < _num_dim; idim++) xst[idim] = xend[idim];

        // update basis
        for (size_t idim = 0; idim < _num_dim; idim++) vect[idim] = _x[closest_seed][idim] - _x[seed_index][idim];

        double norm;
        get_normal_component(_num_dim, num_basis, basis, vect, norm);
        for (size_t idim = 0; idim < _num_dim; idim++) basis[num_basis][idim] = vect[idim];
        num_basis++;

        // update remaining basis
        for (size_t ibasis = num_basis; ibasis < _num_dim; ibasis++)
        {
          for (size_t idim = 0; idim < _num_dim; idim++)
          {
            for (size_t jdim = 0; jdim < _num_dim; jdim++) vect[jdim] = 0.0;
            vect[idim] = 1.0;
            get_normal_component(_num_dim, num_basis, basis, vect, norm);
            if (norm > 0.1) break;
          }

          for (size_t idim = 0; idim < _num_dim; idim++) basis[ibasis][idim] = vect[idim];
        }
        break;
//#pragma endregion
      }
      else
      {
        neighbor = closest_seed;
        trim_spoke(_num_dim, xst, xend, _x[seed_index], _x[neighbor]);
        closest_distance = 0.0;
        for (size_t idim = 0; idim < _num_dim; idim++)
        {
          double dx = xend[idim] - _x[seed_index][idim];
          closest_distance += dx * dx;
        }
        closest_distance = sqrt(closest_distance);
      }
    }
    if (done) break;
  }

  for (size_t idim = 0; idim < _num_dim; idim++) v[idim] = xend[idim];
  for (size_t idim = 0; idim < _num_dim; idim++) delete[] basis[idim];
  delete[] basis;

  delete[] xst; delete[] xend; delete[] dart; delete[] simplex; delete[] vect;
  return 0;
//#pragma endregion
}

int VPS::sample_voronoi_vertex(size_t seed_index, double* xmin, double* xmax, double diag, double* v)
{

//#pragma region Sample A Voronoi Vertex bounding Cell with seed_index:

  size_t num_basis(0);
  double** basis = new double*[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++)
  {

    basis[idim] = new double[_num_dim];

    for (size_t jdim = 0; jdim < _num_dim; jdim++) basis[idim][jdim] = 0.0;

    basis[idim][idim] = 1.0;

  }

  double* xst = new double[_num_dim];
  double* xend = new double[_num_dim];
  double* dart = new double[_num_dim];
  size_t* simplex = new size_t[_num_dim + 1];
  double* vect = new double[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++) xst[idim] = _x[seed_index][idim];

  double closest_distance = 0.0;
  // add first point
  size_t num_simplex_seeds(0);
  simplex[num_simplex_seeds] = seed_index; num_simplex_seeds++;

  size_t num_misses(0); size_t vertex_dim(0);

  while (num_simplex_seeds <= _num_dim - vertex_dim)
  {
    // Throw a random spoke
    sample_uniformly_from_unit_sphere(dart, _num_dim - num_basis);

    for (size_t idim = 0; idim < _num_dim; idim++) xend[idim] = xst[idim];

    for (size_t ibasis = num_basis; ibasis < _num_dim; ibasis++)
    {
      for (size_t idim = 0; idim < _num_dim; idim++) xend[idim] += diag * dart[ibasis - num_basis] * basis[ibasis][idim];
    }

    closest_distance = 0.0;

    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      double dx = xend[idim] - _x[seed_index][idim];
      closest_distance += dx * dx;
    }
    closest_distance = sqrt(closest_distance);

    size_t neighbor(seed_index);
    while (true)
    {
      // get closest point to xend
      size_t closest_seed(seed_index);

      double excluded_closest_distance = closest_distance + 1E-10;

      get_closest_seed_tree(xend, num_simplex_seeds, simplex, closest_seed, excluded_closest_distance);

      bool vertex_out_of_bounding_box(false);

      if (neighbor == closest_seed)
      {
        for (size_t idim = 0; idim < _num_dim; idim++)
        {
          if (xend[idim] > xmin[idim] - 1E-10 && xend[idim] < xmax[idim] + 1E-10) continue;
          vertex_out_of_bounding_box = true;
          break;
        }
      }

      if (closest_seed == _num_samples || (neighbor == closest_seed && vertex_out_of_bounding_box))
      {
        // reset sampling
        num_simplex_seeds = 1;
        for (size_t idim = 0; idim < _num_dim; idim++) xst[idim] = _x[seed_index][idim];

        num_basis = 0;
        for (size_t idim = 0; idim < _num_dim; idim++)
        {
          for (size_t jdim = 0; jdim < _num_dim; jdim++) basis[idim][jdim] = 0.0;
          basis[idim][idim] = 1.0;
        }
        num_misses++;

        if (num_misses == 100)
        {
          vertex_dim++; num_misses = 0;
        }
        break;
      }

      if (neighbor == closest_seed || fabs(excluded_closest_distance - closest_distance) < 1E-10)
      {
//#pragma region A voronoi neighbor apply recursion:
        simplex[num_simplex_seeds] = neighbor; num_simplex_seeds++;

        // update spoke start
        for (size_t idim = 0; idim < _num_dim; idim++) xst[idim] = xend[idim];

        // update basis
        for (size_t idim = 0; idim < _num_dim; idim++) vect[idim] = _x[closest_seed][idim] - _x[seed_index][idim];

        double norm;
        get_normal_component(_num_dim, num_basis, basis, vect, norm);
        for (size_t idim = 0; idim < _num_dim; idim++) basis[num_basis][idim] = vect[idim];
        num_basis++;

        // update remaining basis
        for (size_t ibasis = num_basis; ibasis < _num_dim; ibasis++)
        {
          for (size_t idim = 0; idim < _num_dim; idim++)
          {
            for (size_t jdim = 0; jdim < _num_dim; jdim++) vect[jdim] = 0.0;
            vect[idim] = 1.0;
            get_normal_component(_num_dim, num_basis, basis, vect, norm);
            if (norm > 0.1) break;
          }

          for (size_t idim = 0; idim < _num_dim; idim++) basis[ibasis][idim] = vect[idim];
        }
        break;
//#pragma endregion
      }
      else
      {
        neighbor = closest_seed;
        trim_spoke(_num_dim, xst, xend, _x[seed_index], _x[neighbor]);
        closest_distance = 0.0;
        for (size_t idim = 0; idim < _num_dim; idim++)
        {
          double dx = xend[idim] - _x[seed_index][idim];
          closest_distance += dx * dx;
        }
        closest_distance = sqrt(closest_distance);
      }
    }
  }

  for (size_t idim = 0; idim < _num_dim; idim++) v[idim] = xend[idim];
  for (size_t idim = 0; idim < _num_dim; idim++) delete[] basis[idim];
  delete[] basis;

  delete[] xst; delete[] xend; delete[] dart; delete[] simplex; delete[] vect;
  return 0;
//#pragma endregion
}

bool VPS::trim_spoke(size_t num_dim, double* xst, double* xend, double* p, double* q)
{

//#pragma region Trim a Spoke:
  double* nH = new double[num_dim];
  double* xH = new double[num_dim];

  // trim spoke using Voronoi hyperplane between ipoint and iclosest
  //double norm(0.0);
  for (size_t idim = 0; idim < num_dim; idim++)
  {
    nH[idim] = q[idim] - p[idim];
    xH[idim] = 0.5 * (q[idim] + p[idim]);
  }

  double dotv(0.0), dote(0.0);
  for (size_t idim = 0; idim < num_dim; idim++)
  {
    double dxv = xH[idim] - xst[idim];
    double dxe = xend[idim] - xst[idim];
    dotv += dxv * nH[idim];
    dote += dxe * nH[idim];
  }
  delete[] nH; delete[] xH;
  bool trimmed = false;

  if (fabs(dote) > 1E-10)
  {
    double u = dotv / dote;
    if (u < 0.0) u = 0.0;
    if (u > -1.0E-10 && u < 1.0 - 1.0E-10)
    {
      for (size_t idim = 0; idim < num_dim; idim++) xend[idim] = xst[idim] + u * (xend[idim] - xst[idim]);
      trimmed = true;
    }
  }
  return trimmed;
//#pragma endregion
}

int VPS::get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm)
{
//#pragma region Get Normal component to some basis:
  double* comp = new double[num_basis];

  // project point to current basis
  for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
  {
    comp[ibasis] = 0.0;
    for (size_t idim = 0; idim < num_dim; idim++) comp[ibasis] += vect[idim] * basis[ibasis][idim];
  }
  // get vector component orthogonal to current basis
  for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
  {
    for (size_t idim = 0; idim < num_dim; idim++) vect[idim] -= comp[ibasis] * basis[ibasis][idim];
  }
  delete[] comp;

  norm = 0.0;
  for (size_t idim = 0; idim < num_dim; idim++) norm += vect[idim] * vect[idim];

  if (fabs(norm) < 1E-10) return 1;

  norm = 1.0 / sqrt(norm);

  for (size_t idim = 0; idim < num_dim; idim++) vect[idim] *= norm;

  return 0;

//#pragma endregion
}


void VPS::plot_vps_frames(std::string file_name, size_t function_index, size_t nx, size_t ny, size_t num_contours, bool plot_graph)
{
//    #pragma region Plot Frames:
  double* xsec = new double[_num_dim];
  for (size_t idim = 0; idim < _num_dim; idim++) xsec[idim] = 0.5 * (_xmin[idim] + _xmax[idim]);
  size_t dim_i(0), dim_j(1), dim_k(2), dim_l(3);

  std::fstream ps_file;
  double scale;
  create_ps_file(file_name, nx, ny, ps_file, scale);

  double fmin(DBL_MAX), fmax(-DBL_MAX);

  size_t num_misses(0), max_misses(100);
  double* xx = new double[_num_dim];
  double* ff = new double[_num_functions];
  while (num_misses < max_misses)
  {
    for (size_t idim = 0; idim < _num_dim; idim++)
    {
      xx[idim] = _xmin[idim] + generate_a_random_number() * (_xmax[idim] - _xmin[idim]);
    }
    evaluate_surrogate(xx, ff);
    if (ff[function_index] < fmin)
    {
      fmin = ff[function_index]; num_misses = 0;
    }
    else num_misses++;
    if (ff[function_index] > fmax)
    {
      fmax = ff[function_index];
      num_misses = 0;
    }
    else num_misses++;
  }
  delete[] xx;
  delete[] ff;

  std::vector<double> contours;
  contours.push_back(fmin - 2 * (fmax - fmin));
  for (size_t i = 0; i < num_contours; i++) contours.push_back(fmin + (1.0 / num_contours) * i * (fmax - fmin));
  contours.push_back(fmax + 2 * (fmax - fmin));

  for (size_t i = 0; i < nx; i++)
  {
    xsec[dim_k] = _xmin[dim_k] + i * 1.0 / (nx - 1) * (_xmax[dim_k] - _xmin[dim_k]);
    for (size_t j = 0; j < ny; j++)
    {
      xsec[dim_l] = _xmin[dim_l] + j * 1.0 / (ny - 1) * (_xmax[dim_l] - _xmin[dim_l]);
      plot_vps_surrogate_frame(ps_file, scale, function_index, xsec, dim_i, dim_j, i, j, nx, ny, false, contours);
    }
  }
  ps_file << "showpage" << std::endl;

  delete[] xsec;
//#pragma endregion
}

void VPS::create_ps_file(std::string file_name, size_t nx, size_t ny, std::fstream &file, double &scale)
{
//#pragma region Create PS file:
  //std::cout << ".: VPS Debug Mode :. Plotting ps files .... " << std::endl;

  file.open(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

  double xmin(_xmin[0]);
  double ymin(_xmin[1]);
  double Lx(_xmax[0] - _xmin[0]);
  double Ly(_xmax[1] - _xmin[0]);

  Lx *= (nx + 1); Ly *= (ny + 1);

  double scale_x, scale_y;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  if (scale_x < scale_y)
  {
    scale = scale_x;
    shift_x = 1.0 - xmin * scale;
    shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
  }
  else
  {
    scale = scale_y;
    shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
    shift_y = 1.0 - ymin * scale;
  }
  file << shift_x << " " << shift_y << " translate" << std::endl;

  file << "/Courier findfont" << std::endl;
  file << "0.05 scalefont" << std::endl;
  file << "setfont" << std::endl;

//#pragma endregion
}

void VPS::plot_vps_surrogate_frame(std::fstream &file, double scale, size_t function_index, double* xsec, size_t dim_i, size_t dim_j,
                                   size_t frame_i, size_t frame_j, size_t frame_ni, size_t frame_nj, bool plot_graph, std::vector<double> &contours)
{
//#pragma region Plot Solid Isocontours Of a Single Frame:

  double dx = (_xmax[dim_i] - _xmin[dim_i]);
  double DX = frame_i * dx *(1.0 + 1.0 / (frame_ni - 1));

  double dy = (_xmax[dim_j] - _xmin[dim_j]);
  double DY = frame_j * dy *(1.0 + 1.0 / (frame_nj - 1));

  std::vector<double> poly_x;
  std::vector<double> poly_y;

  size_t num_cells(100);
  double* xx = new double[_num_dim];
  double* f = new double[_num_functions];
  double sx = (_xmax[dim_i] - _xmin[dim_i]) / num_cells;
  double sy = (_xmax[dim_j] - _xmin[dim_j]) / num_cells;

  for (size_t idim = 0; idim < _num_dim; idim++) xx[idim] = xsec[idim];

  for (size_t i = 0; i < num_cells; i++)
  {
    double xo = _xmin[dim_i] + i * sx;
    for (size_t j = 0; j < num_cells; j++)
    {
      double fo(0.0), f1(0.0), f2(0.0), f3(0.0);

      double yo = _xmin[dim_j] + j * sy;
      xx[dim_i] = xo; xx[dim_j] = yo;
      evaluate_surrogate(xx, f);
      fo = f[function_index];

      xx[dim_i] = xo + sx; xx[dim_j] = yo;
      evaluate_surrogate(xx, f);
      f1 = f[function_index];

      xx[dim_i] = xo + sx; xx[dim_j] = yo + sy;
      evaluate_surrogate(xx, f);
      f2 = f[function_index];

      xx[dim_i] = xo; xx[dim_j] = yo + sy;
      evaluate_surrogate(xx, f);
      f3 = f[function_index];

      size_t num_isocontours = contours.size();
      for (size_t icont = 0; icont < num_isocontours; icont++)
      {
        double contour = contours[icont];
        double contour_m = -1000.00;
        if (icont > 0) contour_m = contours[icont - 1];

        //std::cout<< "contour_m = " << contour_m << " , contour = " << contour << std::endl;

        poly_x.clear(); poly_y.clear();

        // moving right
        if (fo >= contour_m - 1E-10 && fo < contour + 1E-10)
        {
          poly_x.push_back(xo);
          poly_y.push_back(yo);
          if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
          {
            double h = sx * (contour - fo) / (f1 - fo);
            poly_x.push_back(xo + h);
            poly_y.push_back(yo);
          }
          else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
          {
            double h = sx * (contour_m - fo) / (f1 - fo);
            poly_x.push_back(xo + h);
            poly_y.push_back(yo);
          }
        }
        else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
        {
          double hm = sx * (contour_m - fo) / (f1 - fo);
          double h = hm;
          if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
          {
            h = sx * (contour - fo) / (f1 - fo);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo + hm);
          poly_y.push_back(yo);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo + h);
            poly_y.push_back(yo);
          }
        }
        else if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
        {
          double h = sx * (contour - fo) / (f1 - fo);
          poly_x.push_back(xo + h);
          poly_y.push_back(yo);
        }

        // moving up
        if (f1 >= contour_m - 1E-10 && f1 < contour + 1E-10)
        {
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo);
          if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
          {
            double h = sy * (contour - f1) / (f2 - f1);
            poly_x.push_back(xo + sx);
            poly_y.push_back(yo + h);
          }
          else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
          {
            double h = sy * (contour_m - f1) / (f2 - f1);
            poly_x.push_back(xo + sx);
            poly_y.push_back(yo + h);
          }

        }
        else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
        {
          double hm = sy * (contour_m - f1) / (f2 - f1);
          double h = hm;
          if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
          {
            h = sy * (contour - f1) / (f2 - f1);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo + hm);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo + sx);
            poly_y.push_back(yo + h);
          }
        }
        else if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
        {
          double h = sy * (contour - f1) / (f2 - f1);
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo + h);
        }

        // moving left
        if (f2 >= contour_m - 1E-10 && f2 < contour + 1E-10)
        {
          poly_x.push_back(xo + sx);
          poly_y.push_back(yo + sy);
          if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
          {
            double h = sx * (contour - f2) / (f3 - f2);
            poly_x.push_back(xo + sx - h);
            poly_y.push_back(yo + sy);
          }
          else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
          {
            double h = sx * (contour_m - f2) / (f3 - f2);
            poly_x.push_back(xo + sx - h);
            poly_y.push_back(yo + sy);
          }
        }
        else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
        {
          double hm = sx * (contour_m - f2) / (f3 - f2);
          double h = hm;
          if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
          {
            h = sx * (contour - f2) / (f3 - f2);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo + sx - hm);
          poly_y.push_back(yo + sy);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo + sx - h);
            poly_y.push_back(yo + sy);
          }
        }
        else if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
        {
          double h = sx * (contour - f2) / (f3 - f2);
          poly_x.push_back(xo + sx - h);
          poly_y.push_back(yo + sy);
        }

        // moving down
        if (f3 >= contour_m - 1E-10 && f3 < contour + 1E-10)
        {
          poly_x.push_back(xo);
          poly_y.push_back(yo + sy);
          if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
          {
            double h = sy * (contour - f3) / (fo - f3);
            poly_x.push_back(xo);
            poly_y.push_back(yo + sy - h);
          }
          else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
          {
            double h = sy * (contour_m - f3) / (fo - f3);
            poly_x.push_back(xo);
            poly_y.push_back(yo + sy - h);
          }
        }
        else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
        {
          double hm = sy * (contour_m - f3) / (fo - f3);
          double h = hm;
          if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
          {
            h = sy * (contour - f3) / (fo - f3);
          }
          if (h < hm)
          {
            double tmp = h; h = hm; hm = tmp;
          }
          poly_x.push_back(xo);
          poly_y.push_back(yo + sy - hm);

          if (h - hm > 1E-10)
          {
            poly_x.push_back(xo);
            poly_y.push_back(yo + sy - h);
          }
        }
        else if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
        {
          double h = sy * (contour - f3) / (fo - f3);
          poly_x.push_back(xo);
          poly_y.push_back(yo + sy - h);
        }


        size_t num_corners(poly_x.size());

        // shift ploy
        for (size_t icorner = 0; icorner < num_corners; icorner++)
        {
          poly_x[icorner] += DX; poly_y[icorner] += DY;
        }

        if (num_corners > 1)
        {
          double gs = 1.0 - icont * 1.0 / num_isocontours;
          file << "newpath" << std::endl;
          file << poly_x[0] * scale << " " << poly_y[0] * scale << " moveto" << std::endl;
          //std::cout<< "*** x = " <<  poly_x[0] << ", y = " << poly_y[0] << std::endl;
          for (size_t icorner = 1; icorner < num_corners; icorner++)
          {
            file << poly_x[icorner] * scale << " " << poly_y[icorner] * scale << " lineto" << std::endl;
            //std::cout << "*** x = " <<  poly_x[icorner] << ", y = " << poly_y[icorner] << std::endl;
          }
          //std::cout << std::endl;

          file << "closepath" << std::endl;
          file << "gsave" << std::endl;
          file << "grestore" << std::endl;

          double r, g, b;

          if (gs < 0.25)     r = 1.0;
          //else if (gs < 0.5) r = 2.0 - 4.0 * gs;
          else if (gs < 0.5) r = 1.0 - 16.0 * (gs - 0.25) * (gs - 0.25);
          else               r = 0.0;

          double go(0.25), gn(1.0 - go);
          if (gs < go)      g = gs / go;
          else if (gs < gn) g = 1.0;
          else              g = 1.0 / (1.0 - gn) - gs / (1.0 - gn);


          if (gs < 0.5)       b = 0.0;
          else if (gs < 0.75) b = 1.0 - 16.0 * (gs - 0.75) * (gs - 0.75);
          else                b = 1.0;

          file << r << " " << g << " " << b << " setrgbcolor" << std::endl;

          file << " fill" << std::endl;
        }
      }
    }
  }
  delete[] xx; delete[] f;
}
