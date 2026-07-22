// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 2016/09/23

#ifndef _TRILINOS_COUPLINGS_FENL_VPS_H_
#define _TRILINOS_COUPLINGS_FENL_VPS_H_

#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include <time.h>

#include <cfloat> // ETP

class VPS
{
public:

  const double DST_TOL = 1E-10;
  const double PI = 3.141592653589793;

  VPS();

  ~VPS();

  int clear_memory();


  enum surrogate_method{ Interpolation, Regression, Integration };

  enum surrogate_basis{ monomials, legendre, chebyshev, radial };

  int get_initial_well_spaced_samples(size_t num_dim, double* xmin, double* xmax, size_t num_samples, double** x);

  int build_surrogate(size_t num_dim,
                      double* xmin, double* xmax,
                      size_t num_functions,
                      surrogate_method method,
                      surrogate_basis basis, size_t desired_order,
                      size_t budget, size_t num_samples,
                      double** x,
                      double** f,
                      double*** g,
                      double**** h);

  int evaluate_surrogate(double* x, double* fs);

  int evaluate_surrogate(size_t cell_index, double* x, double* fs);

  int suggest_new_sample(double* x, double &r, double &err_est);
 
  int get_ensemble(size_t num_ensemble_points, double** ensemble_points,
                   int proc_rank = 0);

  int get_stats(size_t num_mc_points, size_t function_index, double& mean, double& var);

private:

  int sample_voronoi_vertex(size_t seed_index, double* xmin, double* xmax, double diag, double* v);

  int sample_voronoi_facet(size_t seed_index, double* xmin, double* xmax, double diag, double* v);

  int get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm);

  bool trim_spoke(size_t num_dim, double* xst, double* xend, double* p, double* q);



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// RNG Methods       /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void initiate_random_number_generator(unsigned long x);
  double generate_a_random_number();

  int sample_uniformly_from_unit_sphere(double* dart, size_t num_dim);

  ////////////////////////////////////////////////////////////////////////////////
  /////////// RNG  Variables /////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  double Q[1220];
  int indx;
  double cc;
  double c; /* current CSWB */
  double zc;      /* current SWB `borrow` */
  double zx;      /* SWB seed1 */
  double zy;      /* SWB seed2 */
  size_t qlen;/* length of Q array */

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// kd-tree Methods       /////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int kd_tree_build_balanced();

  int kd_tree_balance_quicksort(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

  int kd_tree_quicksort_adjust_target_position(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

  int kd_tree_add_point(size_t point_index);

  int kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
                                  double r,                                                     // Sphere radius
                                  size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
                                  size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
                                  size_t &capacity                                              // Size of points in sphere array
    );

  int kd_tree_get_closest_seed(size_t seed_index,                                  // seed index
                               size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
                               size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
    );

  int kd_tree_get_closest_seed(double* x,                                          // a point in space
                               size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
                               size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
    );

  int kd_tree_get_closest_seed(double* x,                                          // a point in space
                               size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
                               size_t num_exculded_seeds, size_t* exculded_seeds,  // Number of excluded seeds and their indices
                               size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
    );

  int get_closest_seed_tree(size_t seed_index, size_t &closest_seed, double &closest_distance);

  int get_closest_seed_tree(double* x, size_t &closest_seed, double &closest_distance);

  int get_closest_seed_tree(double* x, size_t num_exculded_seeds, size_t* exculded_seeds, size_t &closest_seed, double &closest_distance);

  int get_seeds_in_sphere_tree(double* x, double r, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);

  int get_closest_seed_brute(double* x, size_t num_significant_neighbors, size_t* significant_neighbors, size_t num_exculded_seeds, size_t* exculded_seeds, size_t &closest_seed);


  ////////////////////////////////////////////////////////////////////////////////
  /////////// kd-tree  Variables /////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  size_t            _tree_origin;                // index of tree root
  size_t            _tree_max_height;            // Max height of a tree branch
  size_t*           _tree_left;                  // left pointer of a tree node
  size_t*           _tree_right;                 // right pointer of a tree node

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// General Methods      //////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int quicksort(size_t* I, size_t left, size_t right);

  int quicksort(double* V, size_t* I, size_t left, size_t right);

  bool find_binary(size_t entry, size_t* I, size_t left, size_t right);

  bool find_brute(size_t entry, size_t* I, size_t num_entries);

  int add_entry(size_t entry, size_t &num_entries, size_t* &I, size_t &capacity);

  size_t retrieve_num_permutations(size_t num_dim, size_t upper_bound, bool force_sum_constraint, size_t sum_constraint);

  int retrieve_permutations(size_t &num_perm, size_t** &perm, size_t num_dim, size_t upper_bound, bool force_sum_constraint, size_t sum_constraint);

  void plot_graph(std::string file_name, size_t num_points, double* x, size_t num_functions, double** f);

  void plot_polynomial(std::string file_name, size_t num_basis, double* c, double xmin, double xmax, size_t num_points, double* px, double* py);

  void plot_piecewise_polynomial(std::string file_name, size_t num_pieces, double* pmin, double* pmax, size_t* num_basis, double** c, size_t num_points, double* px, double* py);

  void plot_FourierExpansion(std::string file_name, size_t num_basis, double xmin, double xmax, double* a, double* b, size_t num_points, double* px, double* py, size_t* num_cbasis, double** c);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// VPS General Methods      //////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int find_domain_diagonal();

  // Plot Delaunay Mesh in a postscript
  void plot_delaunay_graph(const std::string outFile);

  // Plot 2d Function
  void plot_vps_surrogate(std::string file_name, size_t function_index, size_t num_contours, bool plot_graph);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// Delaunay Graph Methods      ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int construct_delaunay_graph();

  int update_delaunay_graph(size_t seed_index);

  int sample_voronoi_vertex(double* v);

  int connect_seeds(size_t seed_i, size_t seed_j);

  int disconnect_seeds(size_t seed_i, size_t seed_j);

  int get_num_seed_neighbors(size_t seed_index, size_t &num_seed_neighbors);

  ////////////////////////////////////////////////////////////////////////////////
  /////////// Delaunay Graph  Variables //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  size_t**           _seed_neighbors;                             // indices of delaunay neighbors around a given seed
  bool***            _seed_disc_neighbors;                        // boolean to indicate discontinuous neighbors
  double**           _seed_box;                                   // A box around every seed
  double*            _seed_rf;                                    // Seed disk-free Radius
  double*            _seed_rc;                                    // Seed coverage Radius

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// Fitting Methods                          //////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Data Fit
  int LeastSquare_QR(size_t num_data_points, double xmin, double xmax, double* x, double* f, size_t num_basis, double* c);

  int NaturalCubicSplineInterpolation(size_t num_data_points, double* x, double* y, double* co, double* c1, double* c2, double* c3);

  int FourierExpansion(size_t num_data_points, double xmin, double xmax, double* x, double* f, double** c, size_t num_basis, double* a, double* b);

  // Linear Solvers

  int solve_LU_Doolittle(size_t num_eq, double* LD, double* DD, double* UD, double* b, double* x);

  void LS_QR_Solver(size_t nrow, size_t ncol, double** A, double* b, double* x);

  void GramSchmidt(size_t nrow, size_t ncol, double** A, double** Q);

  void orthonormalize_vector(size_t nbasis, size_t ndim, double** e_basis, double* v, double &norm);

  void QR_Factorization(size_t nrow, size_t ncol, double** A, double** Q, double** R);

  void MAT_MAT(size_t nrow, size_t ncol, double** A, size_t mcol, double** B, double** C);

  void MAT_Transpose(size_t nrow, size_t ncol, double** A, double** AT);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// Surrogate Methods           ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int construct_local_surrogates();

  int init_vps_containers();

  int add_neighbor_layer(size_t cell_index, size_t function_index, size_t &num_neighbors, size_t* &neighbors, size_t &neighbors_capacity);

  int retrieve_weights_regression(size_t cell_index);

  double evaluate_basis_function(double* x, size_t cell_index, size_t ibasis);

  double evaluate_one_dimensional_basis_function(double x, double xo, double xmin, double xmax, size_t ibasis);

  int detect_discontinuities();

  void plot_vps_frames(std::string file_name, size_t function_index, size_t nx, size_t ny, size_t num_contours, bool plot_graph);

  void create_ps_file(std::string file_name, size_t nx, size_t ny, std::fstream &file, double &scale);

  void plot_vps_surrogate_frame(std::fstream &file, double scale, size_t function_index, double* xsec, size_t dim_i, size_t dim_j,
	                            size_t frame_i, size_t frame_j, size_t frame_ni, size_t frame_nj, bool plot_graph, std::vector<double> &contours);


  ////////////////////////////////////////////////////////////////////////////////
  /////////// Surrogate Variables       //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  surrogate_method _method;
  surrogate_basis  _basis;

  size_t             _pool_order;                                 // order of basis pool
  size_t             _num_basis_pool;                             // number of global basis
  size_t**           _p;                                          // global basis indicator

  size_t             _num_basis;                                  // number of global basis
  double***          _basis_coef;                                 // basis coeficients
  size_t***          _basis_index;                                // basis index



  size_t _num_dim;
  size_t _num_functions;
  size_t _desired_order;
  size_t _budget;
  size_t _num_samples;

  double             _diag;                                       // domain diagonal

  // pointers
  double*            _xmin;                                       // domain lower bound
  double*            _xmax;                                       // domain upper bound
  double**           _x;                                          // sample locations
  double**           _f;                                          // function values
  double***          _g;                                          // function gradients
  double****         _h;                                          // functions hessians


  size_t            _num_vs;
  double**          _vs;


};

#endif
