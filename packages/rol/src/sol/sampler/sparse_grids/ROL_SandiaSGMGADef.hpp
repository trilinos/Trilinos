// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


namespace ROL {

//****************************************************************************80

double *SandiaSGMGA::sgmga_aniso_balance 
(
  double alpha_max,
  int dim_num, 
  double level_weight[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_ANISO_BALANCE "balances" an anisotropic weight vector.
//
//  Discussion:
//
//    The entries in LEVEL_WEIGHT are essentially arbitrary nonnegative numbers.
//
//    The ratio between two entries indicates their relative importance.
//    For instance,
//
//      LEVEL_WEIGHT(1) / LEVEL_WEIGHT(2) = 10
//
//    means that variable 2 is 10 times more important than variable 1.
//    Here, being 10 times more important means that we will generate 10 levels
//    of sparse grid in direction 2 as we generate 1 level in direction 1.
//
//    Under this interpretation, a ratio of 10 already indicates an extreme 
//    imbalanace in the variables, since 10 sparse grid levels in the second
//    variable corresponds roughly to approximating x^1 only, and 
//    all of y^1 through y^10.  A ratio higher than this seems unreasonable.
//
//    Therefore, this function tries to take a somewhat arbitrary level weight
//    vector, and produce a "balanced" level weight vector with the properties
//    that the mininum entry is 1 (representing the item of most importance)
//    and the maximum entry is ALPHA_MAX.  A reasonable value of ALPHA_MAX
//    might be 10 or even 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, double ALPHA_MAX, the maximum legal value of 
//    LEVEL_WEIGHT, after all entries have been divided by the minimum 
//    nonzero entry.  1 <= ALPHA_MAX.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.  
//    The values must be positive.  
//
//    Output, double SGMGA_ANISO_BALANCE[DIM_NUM], the balanced 
//    anisotropic weights.  The smallest nonzero entry is 1.0 and 
//    no entry is greater than ALPHA_MAX.
//
{
  int dim;
  double *level_weight2;
  double level_weight_min;
  int nonzero_num;

  if ( alpha_max < 1.0 )
  {
    std::cerr << "\n";
    std::cerr << "SGMGA_ANISO_BALANCE - Fatal error!\n";
    std::cerr << "  ALPHA_MAX < 1.0\n";
    std::exit ( 1 );
  }
//
//  Find the smallest nonzero entry.
//
  level_weight_min = webbur->r8_huge ( );
  nonzero_num = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      if ( level_weight[dim] < level_weight_min )
      {
        level_weight_min = level_weight[dim];
        nonzero_num = nonzero_num + 1;
      }
    }
  }

  if ( nonzero_num == 0 )
  {
    std::cerr << "\n";
    std::cerr << "SGMGA_ANISO_BALANCE - Fatal error!\n";
    std::cerr << "  Could not find a positive entry in LEVEL_WEIGHT.\n";
    std::exit ( 1 );
  }
//
//  Rescale so the smallest nonzero entry is 1.
//
  level_weight2 = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    level_weight2[dim] = level_weight[dim] / level_weight_min;
  }
//
//  Set the maximum entry to no more than ALPHA_MAX.
//
  for ( dim = 0; dim < dim_num; dim++ )
  {
    level_weight2[dim] = webbur->r8_min ( alpha_max, level_weight2[dim] );
  }
  return level_weight2;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_aniso_normalize 
( 
  int option,
  int dim_num, 
  double level_weight[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_ANISO_NORMALIZE normalizes the anisotropic weight vector.
//
//  Discussion:
//
//    It is convenient for the user to initialize the anisotropic weight
//    vector with any set of positive values.  These values are to be used
//    as coefficients of the 1D levels, to evaluate an expression which 
//    determines which 1D levels will be included in a given rule.
//
//    This means that a relatively LARGE coefficient forces the corresponding 
//    level to be relatively SMALL.  This is perhaps the opposite of what
//    a user might expect.  If a user wishes to use an importance vector,
//    so that a relatively large importance should correspond to more levels,
//    and hence more points, in that dimension, then the function
//    SGMGA_IMPORTANCE_TO_ANISO should be called first!
//
//    Since the weights only represent the relative importance of the
//    components, they may be multiplied by any (positive) scale factor.
//    Nonetheless, it may be convenient to choose a particular normalization
//    for the weights.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int OPTION, the normalization option.
//    0, no scaling is applied.
//    1, the weights are scaled so that the minimum nonzero entry is 1.
//    2, the weights are scaled so that they sum to DIM_NUM.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input/output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
//    weights.  The input values must be strictly positive.  
//    On output, these have been normalized.
//
{
  int dim;
  int found;
  double level_weight_min;
  double level_weight_sum;
//
//  Option 0, no normalization.
//
  if ( option == 0 )
  {
  }
//
//  Option 1, the minimum nonzero entry is 1.
//
  else if ( option == 1 )
  {
    level_weight_min = webbur->r8_huge ( );
    found = 0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( 0.0 < level_weight[dim] )
      {
        if ( level_weight[dim] < level_weight_min )
        {
          level_weight_min = level_weight[dim];
          found = found + 1;
        }
      }
    }

    if ( found == 0 )
    {
      std::cerr << "\n";
      std::cerr << "SGMGA_ANISO_NORMALIZE - Fatal error!\n";
      std::cerr << "  Could not find a positive entry in LEVEL_WEIGHT.\n";
      std::exit ( 1 );
    }

    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_weight[dim] = level_weight[dim] / level_weight_min;
    }
  }
//
//  Option 2, rescale so sum of weights is DIM_NUM.
//
  else if ( option == 2 )
  {
    level_weight_sum = webbur->r8vec_sum ( dim_num, level_weight );

    if ( level_weight_sum <= 0.0 )
    {
      std::cerr << "\n";
      std::cerr << "SGMGA_ANISO_NORMALIZE - Fatal error!\n";
      std::cerr << "  Sum of level weights is not positive.\n";
      std::exit ( 1 );
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_weight[dim] = ( ( double ) ( dim_num ) * level_weight[dim] )
        / level_weight_sum;
    }
  }

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_importance_to_aniso 
( 
  int dim_num, 
  double importance[], 
  double level_weight[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_IMPORTANCE_TO_ANISO: importance to anisotropic weight.
//
//  Discussion:
//
//    To specify the anisotropy of a multidimensional problem, the user is
//    allowed to specify an "importance vector".  This vector can contain
//    any set of positive values.  These values represent the relative
//    importance of each dimension.  These values, with a suitable normalization,
//    will be used to evaluate a constraint of the following form:
//
//      QMIN < Level(1) / Importance(1) + Level(2) / Importance(2) + ...
//             Level(N) / Importance(N) <= QMAX
//
//    and a set of levels that satisfies this constraint will then be included
//    in a given anistotropic sparse grid rule.  Thus, increasing the
//    importance value of a particular dimension allows larger level values
//    in that dimension to satisfy the constraint.
//
//    The program actually works with coefficients LEVEL_WEIGHT that are
//    the inverse of the importance vector entries, with a suitable
//    normalization.  This function is supplied to convert between the
//    more natural "importance vector" and the internally useful 
//    "level_weight" vector.
//
//    This function converts the importance vector to an unnormalized 
//    anisotropy weight vector.
//
//    Note that some (but not all) of the IMPORTANCE vector entries may be zero.
//    This indicates that the corresponding dimension is of "zero" or
//    rather "minimal" importance.  In such a case, only a one-point quadrature
//    rule will be applied for that dimension, no matter what sparse grid
//    level is requested for the overall problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double IMPORTANCE[DIM_NUM], the importance vector.
//    All entries must be nonnegative, and at least one must be positive.
//
//    Output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
//    weights.
//
{
  int dim;
  int found;
  double level_weight_norm;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( importance[dim] < 0.0 )
    {
      std::cerr << "\n";
      std::cerr << "SGMGA_IMPORTANCE_TO_ANISO - Fatal error!\n";
      std::cerr << "  Some IMPORTANCE entries are not positive.\n";
      std::exit ( 1 );
    }
  }

  found = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < importance[dim] )
    {
      level_weight[dim] = 1.0 / importance[dim];
      found = found + 1;
    }
    else
    {
      level_weight[dim] = 0.0;
    }
  }

  if ( found == 0 )
  {
    std::cerr << "\n";
    std::cerr << "SGMGA_IMPORTANCE_TO_ANISO - Fatal error!\n";
    std::cerr << "  No importance entry is positive.\n";
    std::exit ( 1 );
  }

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_index 
(
  int dim_num,
  double level_weight[],
  int level_max, 
  int point_num,
  int point_total_num,
  int sparse_unique_index[],
  int growth,
  int ( SandiaRules::*gw_compute_order[] ) ( int level, int growth ),
  int sparse_order[],
  int sparse_index[]
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_INDEX indexes an SGMGA grid.
//
//  Discussion:
//
//    For each "unique" point in the sparse grid, we return its INDEX and ORDER.
//
//    That is, for the I-th unique point P, we determine the product grid which
//    first generated this point, and  and we return in SPARSE_ORDER the orders 
//    of the 1D rules in that grid, and  and in SPARSE_INDEX the component 
//    indexes in those rules that generated this specific point.
//
//    For instance, say P was first generated by a rule which was a 3D product
//    of a 9th order CC rule and  and a 15th order GL rule, and  and that to 
//    generate P, we used the 7-th point of the CC rule and  and the 3rh point 
//    of the GL rule.  Then the SPARSE_ORDER information would be (9,15) and
//    the SPARSE_INDEX information would be (7,3).  This, combined with the 
//    information in RULE, is enough to regenerate the value of P.
//
//    The user must preallocate space for the output arrays SPARSE_ORDER and
//    SPARSE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, int POINT_TOTAL_NUM, the total number of points in the grid.
//
//    Input, int SPARSE_UNIQUE_INDEX[POINT_TOTAL_NUM], associates each
//    point in the grid with its unique representative.
//
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
//    Output, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, 
//    for each point, the order of the 1D rules used in the grid that 
//    generated it.
//
//    Output, int SPARSE_INDEX[DIM_NUM*POINT_NUM)] lists, for 
//    each point, its index in each of the 1D rules in the grid that generated 
//    it.  The indices are 1-based.
//
{
  double coef;
  int dim;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  bool more_grids;
  bool more_points;
  int *order_1d;
  int point;
  int point_count;
  int *point_index;
  int point_unique;
  double q_max;
  double q_min;

//
//  Special cases.
//
  if ( level_max < 0 )
  {
    return;
  }

  if ( level_max == 0 )
  {
    point = 0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_order[dim+point*dim_num] = 1;
      sparse_index[dim+point*dim_num] = 1;
    }
    return;
  }
//
//  Initialize the INDEX and ORDER arrays to -1 to help catch errors.
//
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_order[dim+point*dim_num] = -1;
      sparse_index[dim+point*dim_num] = -1;
    }
  }

  point_count = 0;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  order_1d = new int[dim_num];
  point_index = new int[dim_num];
//
//  Initialization for SGMGA_VCN_ORDERED.
//
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - webbur->r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int )( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = false;
//
//  Seek all vectors LEVEL_1D which satisfy the constraint:
//
//    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
//      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
//      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
//
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, 
      level_1d, q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
//
//  Compute the combinatorial coefficient.
//
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d, 
      q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
//
//  Transform each 1D level to a corresponding 1D order.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      order_1d[dim] = ((*webbur).*gw_compute_order[dim])(level_1d[dim],growth);
    }
//
//  The inner loop generates a POINT of the GRID of the LEVEL.
//
    more_points = false;

    for ( ; ; )
    {
      webbur->vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

      if ( !more_points )
      {
        break;
      }
      point_unique = sparse_unique_index[point_count];
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_order[dim+point_unique*dim_num] = order_1d[dim];
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_index[dim+point_unique*dim_num] = point_index[dim];
      }
      point_count = point_count + 1;
    }
  }

  delete [] level_1d;
  delete [] level_1d_max;
  delete [] order_1d;
  delete [] point_index;

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_point 
( 
  int dim_num,
  double level_weight[],
  int level_max, 
  void ( SandiaRules2::*gw_compute_points[] ) ( int order, int dim, double x[] ),
  int point_num,
  int sparse_order[],
  int sparse_index[],
  int growth,
  int ( SandiaRules::*gw_compute_order[] ) ( int level, int growth ),
  double sparse_point[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_POINT computes the points of an SGMGA rule.
//
//  Discussion:
//
//    The sparse grid is the logical sum of low degree product rules.
//
//    Each product rule is the product of 1D factor rules.
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//    * the quadrature rules.
//    * the number of points.
//
//    The user must preallocate space for the output array SPARSE_POINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int LEVEL_MAX, controls the size of the final sparse grid.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int dim, double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, int POINT_NUM, the number of points in the grid,
//    as determined by SGMGA_SIZE.
//
//    Input, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, for each point,
//    the order of the 1D rules used in the grid that generated it.
//
//    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM], lists, for each point,
//    its index in each of the 1D rules in the grid that generated it.
//    The indices are 1-based.
//
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
//    Output, double SPARSE_POINT[DIM_NUM*POINT_NUM], the points.
//
{
  int dim;
  int level;
  int *level_1d_max;
  double level_weight_min_pos;
  int order;
  int point;
  double *points;
  double q_max;

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_point[dim+point*dim_num] = - webbur->r8_huge ( );
    }
  }
//
//  Compute the point coordinates.
//
  level_1d_max = new int[dim_num];
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }

    for ( level = 0; level <= level_1d_max[dim]; level++ )
    {
      order = ((*webbur).*gw_compute_order[dim]) ( level, growth );

      points = new double[order];

      ((*webbur2).*gw_compute_points[dim]) ( order, dim, points );

      for ( point = 0; point < point_num; point++ )
      {
        if ( sparse_order[dim+point*dim_num] == order )
        {
          sparse_point[dim+point*dim_num] = 
            points[sparse_index[dim+point*dim_num]-1];
        }
      }
      delete [] points;
    }
  }
//
//  Check to see if we missed any points.
//
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( sparse_point[dim+point*dim_num] == - webbur-> r8_huge ( ) )
      {
        std::cerr << "\n";
        std::cerr << "SGMGA_POINT - Fatal error!\n";
        std::cerr << "  At least one point component was not assigned.\n";
        std::cerr << "  POINT = " << point << "\n";
        std::cerr << "  DIM = " << dim << "\n";
        std::cerr << "  SPARSE_ORDER(DIM,POINT) = " 
          << sparse_order[dim+point*dim_num] << "\n";
        std::cerr << "  LEVEL_WEIGHT(DIM) = " << level_weight[dim] << "\n";
        std::exit ( 1 );
      }
    }
  }

  delete [] level_1d_max;

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_product_weight 
(
  int dim_num,
  int order_1d[],
  int order_nd, 
  void ( SandiaRules2::*gw_compute_weights[] ) ( int order, int dim, double w[] ),
  double weight_nd[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_PRODUCT_WEIGHT computes the weights of a mixed product rule.
//
//  Discussion:
//
//    This routine computes the weights for a quadrature rule which is
//    a product of 1D rules of varying order and kind.
//
//    The user must preallocate space for the output array WEIGHT_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
//
//    Input, int ORDER_ND, the order of the product rule.
//
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int dim, double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Output, double WEIGHT_ND[ORDER_ND], the product rule weights.
//
{
  int dim;
  int i;
  double *weight_1d;
  
  for ( i = 0; i < order_nd; i++ )
  {
    weight_nd[i] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    weight_1d = new double[order_1d[dim]];

    ((*webbur2).*gw_compute_weights[dim]) ( order_1d[dim], dim, weight_1d );

    webbur->r8vec_direct_product2 ( dim, order_1d[dim], weight_1d, dim_num, 
      order_nd, weight_nd );

    delete [] weight_1d;
  }
  return;
}
//****************************************************************************80

int SandiaSGMGA::sgmga_size 
( 
  int dim_num,
  double level_weight[],
  int level_max, 
  void ( SandiaRules2::*gw_compute_points[] ) ( int order, int dim, double x[] ),
  double tol,
  int growth,
  int ( SandiaRules::*gw_compute_order[] ) ( int level, int growth ) 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_SIZE sizes an SGMGA grid, discounting duplicates.
//
//  Discussion:
//
//    The sparse grid is the logical sum of product grids that satisfy
//    a particular constraint.
//
//    Depending on the 1D rules involved, there may be many duplicate points
//    in the sparse grid.
//
//    This function counts the unique points in the sparse grid.  It does this
//    in a straightforward way, by actually generating all the points, and
//    comparing them, with a tolerance for equality.
//
//    This function has been modified to automatically omit points for which
//    the "combinatorial coefficient" is zero, since such points would have
//    a weight of zero in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int dim, double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
//    Output, int SGMGA_SIZE, the number of unique points.
//
{
  double coef;
  int dim;
  int level;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  bool more_grids;
  bool more_points;
  int order;
  int *order_1d;
  int point;
  int *point_index;
  int point_num;
  int point_total_num;
  int point_total_num2;
  double *points;
  double q_max;
  double q_min;
  int seed;
  int *sparse_total_index;
  int *sparse_total_order;
  double *sparse_total_point;
//
//  Special cases.
//
  if ( level_max < 0 )
  {
    point_num = -1;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Get total number of points, including duplicates.
//
  point_total_num = sgmga_size_total ( dim_num, level_weight, 
    level_max, growth, gw_compute_order );
//
//  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
//  for the TOTAL set of points.
//
  sparse_total_order = new int[dim_num*point_total_num];
  sparse_total_index = new int[dim_num*point_total_num];

  point_total_num2 = 0;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  order_1d = new int[dim_num];
  point_index = new int[dim_num];
//
//  Initialization for SGMGA_VCN_ORDERED.
//
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - webbur->r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = false;
//
//  Seek all vectors LEVEL_1D which satisfy the constraint:
//
//    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
//      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
//      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
//
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, 
      level_1d, q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
//
//  Compute the combinatorial coefficient.
//
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d, 
      q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
//
//  Transform each 1D level to a corresponding 1D order.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      order_1d[dim] = ((*webbur).*gw_compute_order[dim])(level_1d[dim],growth);
    }
//
//  The inner loop generates a POINT of the GRID of the LEVEL.
//
    more_points = false;

    for ( ; ; )
    {
      webbur->vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

      if ( !more_points )
      {
        break;
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_order[dim+point_total_num2*dim_num] = order_1d[dim];
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_index[dim+point_total_num2*dim_num] = point_index[dim];
      }
      point_total_num2 = point_total_num2 + 1;
    }
  }
  delete [] level_1d;
  delete [] order_1d;
  delete [] point_index;
//
//  Now compute the coordinates of the TOTAL set of points.
//
  sparse_total_point = new double[dim_num*point_total_num];

  for ( point = 0; point < point_total_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_total_point[dim+point*dim_num] = webbur->r8_huge ( );
    }
  }
//
//  Compute the point coordinates.
//
  level_1d_max = new int[dim_num];
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }

    for ( level = 0; level <= level_1d_max[dim]; level++ )
    {
      order = ((*webbur).*gw_compute_order[dim]) ( level, growth );

      points = new double[order];

      ((*webbur2).*gw_compute_points[dim]) ( order, dim, points );

      for ( point = 0; point < point_total_num; point++ )
      {
        if ( sparse_total_order[dim+point*dim_num] == order )
        {
          sparse_total_point[dim+point*dim_num] = 
            points[sparse_total_index[dim+point*dim_num]-1];
        }
      }
      delete [] points;
    }
  }
//
//  Count the tolerably unique columns. 
//
  seed = 123456789;

  point_num = webbur->point_radial_tol_unique_count ( dim_num, point_total_num, 
    sparse_total_point, tol, &seed );

  delete [] level_1d_max;
  delete [] sparse_total_index;
  delete [] sparse_total_order;
  delete [] sparse_total_point;

  return point_num;
}
//****************************************************************************80

int SandiaSGMGA::sgmga_size_total 
( 
  int dim_num, 
  double level_weight[],
  int level_max, 
  int growth,
  int ( SandiaRules::*gw_compute_order[] ) ( int level, int growth )
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_SIZE_TOTAL sizes an SGMGA grid, counting duplicates.
//
//  Discussion:
//
//    This routine returns the total point count for an SGMGA
//    ( Sparse Grid of Mixed type with Growth rule and Anisotropic weights).
//
//    The sparse grid is the logical sum of product grids.
//
//    The sparse grid has an associated integer index LEVEL_MAX, whose lowest 
//    value is 0.  LEVEL_MAX = 0 indicates the sparse grid made up of one product 
//    grid, which in turn is the product of 1D factor grids of the lowest level.
//    This usually means the sparse grid with LEVEL_MAX equal to 0 is a
//    one point grid.
//
//    We can assign a level to each factor grid, and hence a LEVEL vector
//    to the corresponding product grid, and a weighted index
//    LEVEL_GRID (which will in general be a real number):
//
//      LEVEL_GRID = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL(I)
//
//    The product grid will participate in the formation of the sparse grid
//    if it satisfies the following weighted constraint:
//
//      LEVEL_MAX - DIM_NUM < LEVEL_GRID <= LEVEL_MAX
//
//    This routine determines the total number of abscissas in all the 
//    product rules used to form the SGMGA associated with the index LEVEL_MAX.
//    The count disregards duplication.  If the same multidimensional abcsissa
//    occurs in two different product rules that are part of the SGMGA, then
//    that single abcissa is counted twice. 
//
//    This computation is useful in cases where the entire set of abscissas
//    is going to be generated, preparatory to compression to finding, indexing
//    and merging the duplicate abcissass.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int GROWTH, the growth rule.  
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
//    Output, int SGMGA_SIZE_TOTAL, the number of points
//    including repetitions.
//
{
  double coef;
  int dim;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  bool more_grids;
  int *order_1d;
  int point_total_num;
  double q_max;
  double q_min;
//
//  Special case.
//
  if ( level_max == 0 )
  {
    point_total_num = 1;
    return point_total_num;
  }

  point_total_num = 0;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  order_1d = new int[dim_num];
//
//  Initialization for SGMGA_VCN_ORDERED.
//
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - webbur->r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = false;
//
//  Seek all vectors LEVEL_1D which satisfy the constraint:
//
//    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
//      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
//      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
//
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, 
      level_1d, q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
//
//  Compute the combinatorial coefficient.
//
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d, 
      q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
//
//  Transform each 1D level to a corresponding 1D order.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      order_1d[dim] = ((*webbur).*gw_compute_order[dim])(level_1d[dim],growth);
    }
    point_total_num = point_total_num + webbur->i4vec_product ( dim_num, 
      order_1d );
  }
  delete [] level_1d;
  delete [] level_1d_max;
  delete [] order_1d;

  return point_total_num;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_unique_index 
( 
  int dim_num,
  double level_weight[], 
  int level_max,
  void ( SandiaRules2::*gw_compute_points[] ) ( int order, int dim, double x[] ),
  double tol,
  int point_num,
  int point_total_num,
  int growth,
  int ( SandiaRules::*gw_compute_order[] ) ( int level, int growth ),
  int sparse_unique_index[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_UNIQUE_INDEX maps nonunique to unique points.
//
//  Discussion:
//
//    The sparse grid usually contains many points that occur in more
//    than one product grid.
//
//    When generating the point locations, it is easy to realize that a point
//    has already been generated.
//
//    But when it's time to compute the weights of the sparse grids, it is
//    necessary to handle situations in which weights corresponding to 
//    the same point generated in multiple grids must be collected together.
//
//    This routine generates ALL the points, including their multiplicities,
//    and figures out a mapping from them to the collapsed set of unique points.
//
//    This mapping can then be used during the weight calculation so that
//    a contribution to the weight gets to the right place.
//
//    The user must preallocate space for the output array SPARSE_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int dim, double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, int POINT_TOTAL_NUM, the total number of points 
//    in the grid. 
//
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//    Output, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists, 
//    for each (nonunique) point, the corresponding index of the same point in 
//    the unique listing.
//
{
  double coef;
  int dim;
  int level;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  bool more_grids;
  bool more_points;
  int order;
  int *order_1d;
  int point;
  int *point_index;
  int point_num2;
  int point_total_num2;
  double *points;
  double q_max;
  double q_min;
  int rep;
  int seed;
  int *sparse_total_index;
  int *sparse_total_order;
  double *sparse_total_point;
  int *undx;
//
//  Special cases.
//
  if ( level_max < 0 )
  {
    return;
  }

  if ( level_max == 0 )
  {
    sparse_unique_index[0] = 0;
    return;
  }
//
//  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
//  for the TOTAL set of points.
//
  sparse_total_order = new int[dim_num*point_total_num];
  sparse_total_index = new int[dim_num*point_total_num];

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  order_1d = new int[dim_num];
  point_index = new int[dim_num];

  point_total_num2 = 0;
//
//  Initialization for SGMGA_VCN_ORDERED.
//
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - webbur->r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = false;
//
//  Seek all vectors LEVEL_1D which satisfy the constraint:
//
//    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
//      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
//      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
//
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, 
      level_1d, q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
//
//  Compute the combinatorial coefficient.
//
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d, 
      q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
//
//  Transform each 1D level to a corresponding 1D order.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      order_1d[dim] = ((*webbur).*gw_compute_order[dim])(level_1d[dim],growth);
    }
//
//  The inner loop generates a POINT of the GRID of the LEVEL.
//
    more_points = false;

    for ( ; ; )
    {
      webbur->vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

      if ( !more_points )
      {
        break;
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_order[dim+point_total_num2*dim_num] = order_1d[dim];
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_index[dim+point_total_num2*dim_num] = point_index[dim];
      }
      point_total_num2 = point_total_num2 + 1;
    }
  }
  delete [] level_1d;
  delete [] level_1d_max;
  delete [] order_1d;
  delete [] point_index;
//
//  Now compute the coordinates of the TOTAL set of points.
//
  sparse_total_point = new double[dim_num*point_total_num];

  for ( point = 0; point < point_total_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_total_point[dim+point*dim_num] = webbur->r8_huge ( );
    }
  }
//
//  Compute the point coordinates.
//
  level_1d_max = new int[dim_num];
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }

    for ( level = 0; level <= level_1d_max[dim]; level++ )
    {
      order = ((*webbur).*gw_compute_order[dim]) ( level, growth );

      points = new double[order];

      ((*webbur2).*gw_compute_points[dim]) ( order, dim, points );

      for ( point = 0; point < point_total_num; point++ )
      {
        if ( sparse_total_order[dim+point*dim_num] == order )
        {
          sparse_total_point[dim+point*dim_num] = 
            points[sparse_total_index[dim+point*dim_num]-1];
        }
      }
      delete [] points;
    }
  }
//
//  Merge points that are too close.
//
  seed = 123456789;

  undx = new int[point_num];

  point_num2 = webbur->point_radial_tol_unique_index ( dim_num, point_total_num, 
    sparse_total_point, tol, &seed, undx, sparse_unique_index );

  for ( point = 0; point < point_total_num; point++ )
  {
    rep = undx[sparse_unique_index[point]];
    if ( point != rep )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_point[dim+point*dim_num] = sparse_total_point[dim+rep*dim_num];
      }
    }
  }
//
//  Construct an index that indicates the "rank" of the unique points.
//
  webbur->point_unique_index ( dim_num, point_total_num, sparse_total_point,
    point_num, undx, sparse_unique_index );

  delete [] sparse_total_index;
  delete [] sparse_total_order;
  delete [] sparse_total_point;
  delete [] undx;

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_vcn 
( 
  int n, 
  double w[],
  int x[],
  double q_min,
  double q_max, 
  bool *more 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN returns the next constrained vector.
//
//  Discussion:
//
//    This function is intended to replace the "naive" version, now called
//    SGMGA_VCN_NAIVE, which is too slow for high dimensional problems.
//
//    For nonnegative vectors X of dimension N, and nonnegative
//    weights W, we define:
//
//      Q = sum ( 1 <= I <= N ) W(I) * X(I)
//
//    and seek X satisfying the constraint:
//
//      Q_MIN < Q <= Q_MAX
//
//    This routine returns, one at a time exactly those X which satisfy
//    the constraint.  No attempt is made to return the X values in 
//    any particular order as far as Q goes.  
//
//  Example:
// 
//        W               4.0 3.0 5.0       
//      MIN     16.0       0   0   0
//      ---     ----      -----------
//        1     20.0       5   0   0
//        2     19.0       4   1   0
//        3     18.0       3   2   0
//        4     17.0       2   3   0
//        5     20.0       2   4   0
//        6     19.0       1   5   0
//        7     18.0       0   6   0
//        8     17.0       3   0   1
//        9     20.0       3   1   1
//       10     19.0       2   2   1
//       11     18.0       1   3   1
//       12     17.0       0   4   1
//       13     20.0       0   5   1
//       14     18.0       2   0   2
//       15     17.0       1   1   2
//       16     20.0       1   2   2
//       17     19.0       0   3   2
//       18     19.0       1   0   3
//       19     18.0       0   1   3
//       20     20.0       0   0   4
//      ---     ----      ----------
//      MAX     20.0       6   7   5         
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, double W[N], the weights, which should be nonnegative.
//    At least one weight must be positive.
//
//    Input/output, int X[N].  On first call, with 
//    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
//    the input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper limits on the sum.
//
//    Input/output, bool *MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  static int dir;
  int i;
  static int n2;
  static int nstart;
  double t;
  static int *xmax;
  static int xmin;
//
//  Initialization for first call.
//
//  Allocate XMAX to remember the currently maximum possible value for each X.
//
//  Locate NSTART, the index of the first nonzero weight.
//  The algorithm is easier to program if the last index we look at
//  has a nonzero weight, so that it can always make up the remainder.
//
  if ( !(*more) )
  {
    xmax = new int[n];

    nstart = - 1;

    for ( i = 0; i < n; i++ )
    {
      if ( 0.0 < w[i] )
      {
        nstart = i;
        break;
      }
    }
//
//  Theoretically, we could even handle the case where all weights are zero.
//  That case is ruled out elsewhere in this software, so I will not try
//  to deal with it here for now.
//
    if ( nstart == - 1 )
    {
      std::cerr << "\n";
      std::cerr << " SGMGA_VCN - Fatal error!\n";
      std::cerr << "  No weight is positive.\n";
      std::exit ( 1 );
    }
//
//  Initialize X to zero, even the indices we ignore.
//
    for ( i = 0; i < n; i++ )
    {
      x[i] = 0;
    }
//
//  N2 points to our current index of interest.
//
    n2 = n;
    dir = - 1;

    *more = true;
  }
//
//  Look for the next solution vector X.
//
  for ( ; ; )
  {
//
//  If no more, the search is terminated.
//
    if ( !(*more) )
    {
      break;
    }
//
//  DIR = -1, decrement N2, and, if possible, set X[N2] to XMIN.
//  DIR =  0, hold N2 at current value, and see if we can increment X[N2].
//
    else if ( dir == - 1 || dir == 0 )
    {
      if ( dir == - 1 )
      {
        n2 = n2 - 1;
      }

      if ( w[n2] == 0.0 )
      {
        xmin = 0;
        xmax[n2] = 0;
      }
      else if ( nstart < n2 )
      {
        xmin = 0;
        t = q_max;
        for ( i = n2 + 1; i < n; i++ )
        {
          t = t - w[i] * ( double ) x[i];
        }
        xmax[n2] = ( int ) ( webbur->r8_floor ( t / w[n2] ) );
      }
      else if ( n2 == nstart && dir == - 1 )
      {
        t = q_min;
        for ( i = n2 + 1; i < n; i++ )
        {
          t = t - w[i] * ( double ) x[i];
        }
        xmin = ( int ) ( webbur->r8_ceiling ( t / w[n2] ) );
        if ( xmin < 0 )
        {
          xmin = 0;
        }
        t = 0.0;
        for ( i = 0; i < n2; i++ )
        {
          t = t + w[i] * ( double ) x[i];
        }
        t = t + w[n2] * xmin;
        for ( i = n2 + 1; i < n; i++ )
        {
          t = t + w[i] * ( double ) x[i];
        }
        if ( t <= q_min )
        {
          xmin = xmin + 1;
        }
        x[n2] = xmin;
        t = q_max;
        for ( i = n2 + 1; i < n; i++ )
        {
          t = t - w[i] * ( double ) x[i];
        }
        xmax[n2] = ( int ) ( webbur->r8_floor ( t / w[n2] ) );
      }

      if ( xmax[n2] < xmin )
      {
        dir = + 1;
      }
      else
      {
        if ( n2 == nstart )
        {
          if ( dir == - 1 )
          {
            dir = 0;
            break;
          }
          else if ( dir == 0 )
          {
            x[n2] = x[n2] + 1;
            if ( x[n2] <= xmax[n2] )
            {
              break;
            }
            else
            {
              dir = + 1;
            }
          }
        }
        else
        {
          x[n2] = xmin;
        }
      }
    }
//
//  DIR = + 1:
//  Try moving backwards to find an index N2 whose X we can increment.
//
    else if ( dir == + 1 )
    {
      for ( ; ; )
      {
        if ( n2 == n - 1 )
        {
          dir = 0;
          *more = false;
          delete [] xmax;
          break;
        }

        n2 = n2 + 1;

        if ( 0.0 < w[n2] )
        {
          if ( x[n2] < xmax[n2] )
          {
            x[n2] = x[n2] + 1;
            dir = - 1;
            break;
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

double SandiaSGMGA::sgmga_vcn_coef 
( 
  int dim_num, 
  double level_weight[],
  int x[], 
  double q_max 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_COEF returns the "next" constrained vector's coefficient.
//
//  Discussion:
//
//    The related code "SGMGA_VCN_COEF_NAIVE" represents a "naive" 
//    approach to this calculation.  This code carries out the same calculation, but tries
//    to avoid the potential explosion in work that is exponential in the
//    spatial dimension for the naive approach.
//
//    We are considering integer vectors X of dimension DIM_NUM for which
//    the functional
//
//      Q(X) = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
//
//   satisfies the "Q" constraint:
//
//      Q_MIN < Q(X) <= Q_MAX
//
//    where LEVEL_WEIGHT is a vector of (essentially) positive weights.
//    Some, but not all of the entries of LEVEL_WEIGHT might be zero;
//    in that case, the corresponding values of X never vary, and do not
//    play a part in the following computation.
//
//    Supposing we have a suitable vector X, we now wish to count the 
//    number of distinct vectors Y which also satisfy the Q constraint
//    as well as the following "binary" constraint:
//
//      Y(I) = X(I) + B(I)
//
//    where every entry of B is 0 or 1.
//
//    Clearly, there are 2^DIM_NUM vectors Y which satisfy the binary
//    constraint, and a naive calculation would simply generate each 
//    possible Y, evaluate Q(Y), and if Y satisfies the Q constraint,
//    add it to the count.
//
//    But if we are considering even moderately large values of DIM_NUM, 
//    say 20 <= DIM_NUM, then the mere task of generating all possible 
//    Y vectors is burdensome.  If there are in fact likely to be only 
//    a few satisfactory Y vectors, (which depends on the values of 
//    Q_MAX and LEVEL_WEIGHT, of course) then it may be possible to
//    find and count them more rapidly.
//
//    This function attempts a more rapid computation by carrying out the
//    search in a particular order, and realizing that, in certain cases,
//    if a particular value Y* does not satisfy the Q constraint, then
//    a consecutive sequence of Y's following Y* also cannot satisfy the
//    constraint, and hence the search can jump over them.
//
//  Example:
//
//    DIM_NUM = 3
//    LEVEL_WEIGHT    3.0  2.0  1.0
//    Q_MAX    6.0
//
//    U = unsigned count 
//    S =   signed count returned as COEF
//                 
//    #   U  S   X1 X2 X3
//
//    1   8  0    0  0  0
//    2   7  1    0  0  1
//    3   6  0    0  0  2
//    4   5 -1    0  0  3
//    5   3 -1    0  0  4
//    6   2  0    0  0  5
//    7   1  1    0  0  6
//    8   6  0    0  1  0
//    9   5 -1    0  1  1
//   10   3 -1    0  1  2
//   11   2  0    0  1  3
//   12   1  1    0  1  4
//   13   3 -1    0  2  0
//   14   2  0    0  2  1
//   15   1  1    0  2  2
//   16   1  1    0  3  0
//   17   5 -1    1  0  0
//   18   3 -1    1  0  1
//   19   2  0    1  0  2
//   20   1  1    1  0  3
//   21   2  0    1  1  0
//   22   1  1    1  1  1
//   23   1  1    2  0  0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the weights.
//
//    Input, int X[DIM_NUM], satisfies the Q constraint.
//
//    Input, double Q_MAX, the Q constraint maximum.
//
//    Output, double SGMGA_VCN_COEF, the combinatorial coefficient.
//
{
  int *b;
  int b_sum;
  int c;
  double coef;
  int i;
  int j;
  double q;

  c = 0;
  b = new int[dim_num];

  for ( i = 0; i < dim_num; i++ )
  {
    b[i] = 0;
  }

  for ( ; ; )
  {
//
//  Generate the next binary perturbation.
//
    i = - 1;

    while ( i < dim_num - 1 )
    {
      i = i + 1;
//
//  If LEVEL_WEIGHT(I) == 0, B(I) is fixed at 0.  Next I.
//
      if ( level_weight[i] == 0.0 )
      {
      }
//
//  If B(I) is 1, set it to 0.  Next I.
//
      else if ( b[i] == 1 )
      {
        b[i] = 0;
      }
//
//  B(I) is 0.  Convert it to 1.
//
      else 
      {
        b[i] = 1;

        for ( ; ; )
        {
// 
//  Does X + B satisfy the Q_MAX constraint?
//
          q = 0.0;
          for ( j = 0; j < dim_num; j++ )
          {
            q = q + level_weight[j] * ( double ) ( x[j] + b[j] );
          }
          if ( q <= q_max )
          {
            break;
          }
//
//  If Q(X+B) now exceeds QMAX, B is rejected.  But we can also skip
//  all perturbations which agree with B through the I-th position.
//  To skip to the next "interesting" candidate, we essentially carry
//  out binary addition between B and a vector B' which has a single 1
//  in the I-th position.
//
          b[i] = 0;

          while ( i < dim_num - 1 )
          {
            i = i + 1;

            if ( level_weight[i] == 0.0 )
            {
            }
            else if ( b[i] == 1 )
            {
              b[i] = 0;
            }
            else
            {
              b[i] = 1;
              break;
            }
          }
        }
        break;
      }
    }
    b_sum = 0;
    for ( j = 0; j < dim_num; j++ )
    {
      b_sum = b_sum + b[j];
    }
//
//  X+B is another solution to be counted.
//
    c = c + 1 - 2 * ( b_sum % 2 );
//
//  We're done if we've got back to 0.
//
    if ( b_sum == 0 )
    {
      break;
    }
  }
  coef = ( double ) ( c );
  delete [] b;

  return coef;
}
//****************************************************************************80

double SandiaSGMGA::sgmga_vcn_coef_naive 
(
  int dim_num,
  double level_weight[], 
  int x_max[],
  int x[],
  double q_min,
  double q_max 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_COEF_NAIVE: "next" constrained vector's coefficient.
//
//  Discussion:
//
//    This function uses a naive approach to the computation, resulting in
//    a set of 2^DIM_NUM tasks.  Hence it is not suitable for cases where
//    DIM_NUM is moderately large.  The function SGMGA_VCN_COEF carries out
//    a more complicated but more efficient algorithm for the same computation.
//
//    We are given a vector X of dimension DIM_NUM which satisfies:
//    the following constraint:
//
//      sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
//
//    This routine computes the appropriate coefficient for X in the
//    anisotropic sparse grid scheme.
//
//    The coefficient is calculated as follows:
//
//      Let B be a binary vector of length DIM_NUM, and let ||B|| represent
//      the sum of the entries of B.
//
//      Coef = sum ( all B such that X+B satisfies constraints ) (-1)^||B||
//
//    Since X+0 satisfies the constraint, there is always at least one 
//    summand.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of components in the vector.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
//
//    Input, int X[DIM_NUM], a point which satisifies the constraints.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper limits on the sum.
//
//    Output, double SGMGA_VCN_COEF_NAIVE, the combinatorial coefficient.
//
{
  int *b;
  int b_sum;
  double coef;
  int i;
  double q;
  bool too_big;

  b = new int[dim_num];

  for ( i = 0; i < dim_num; i++ )
  {
    b[i] = 0;
  }
  coef = 1.0;

  for ( ; ; )
  {
//
//  Generate the next binary perturbation.
//
    webbur->binary_vector_next ( dim_num, b );
    b_sum = webbur->i4vec_sum ( dim_num, b );
//
//  We're done if we've got back to 0.
//
    if ( b_sum == 0 )
    {
      break;
    }
//
//  Does it satisfy the XMAX constraint?
//  (THIS CHECK IS SURPRISINGLY NECESSARY, PARTLY BECAUSE OF ZERO WEIGHT).
//
    too_big = false;
    for ( i = 0; i < dim_num; i++ )
    {
      if ( x_max[i] < x[i] + b[i] )
      {
        too_big = true;
        break;
      }
    }
    if ( too_big ) 
    {
      continue;
    }
//
//  Does it satisfy the Q_MAX constraint?
//
    q = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      q = q + level_weight[i] * ( double ) ( x[i] + b[i] );
    }

    if ( q <= q_max )
    {
      coef = coef + webbur->r8_mop ( b_sum );
    }
  }

  delete [] b;

  return coef;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_vcn_naive 
(
  int dim_num,
  double level_weight[],
  int x_max[], 
  int x[],
  double q_min,
  double q_max, 
  bool *more 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_NAIVE returns the next constrained vector.
//
//  Discussion:
//
//    This function uses a naive algorithm, which quickly becomes unsuitable
//    for higher dimensions.  The function SGMGA_VCN is an attempt at 
//    a more efficient calculation of the same quantities.
//
//    We consider vectors X of dimension DIM_NUM satisfying:
//
//      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
//
//    and define
//
//      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
//
//    and seek X satisfying the constraint:
//
//      Q_MIN < Q <= Q_MAX
//
//    For sparse grid applications, we compute
//
//      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
//
//    and assume there is an underlying LEVEL used to index the sets of 
//    constrained vectors, and that 
//
//      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
//      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
//      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
//
//    This routine returns, one at a time exactly those X which satisfy
//    the constraint.  No attempt is made to return the X values in 
//    any particular order as far as Q goes.  
//
//  Example:
//
//    LEVEL_WEIGHT:          1.000000        1.000000
//
//    Q_MIN:        0.000000
//    Q_MAX:        2.000000
//    X_MAX:                         2         2
//
//         1        1.000000         1         0
//         2        2.000000         2         0
//         3        1.000000         0         1
//         4        2.000000         1         1
//         5        2.000000         0         2
//
//    LEVEL_WEIGHT:          1.000000        2.000000
//
//    Q_MIN:       -1.000000
//    Q_MAX:        2.000000
//    X_MAX:                         2         1
//
//         1        0.000000         0         0
//         2        1.000000         1         0
//         3        2.000000         2         0
//         4        2.000000         0         1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of components in the vector.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
//
//    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper
//    limits on the sum.
//
//    Input/output, bool *MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;
  double q;

  if ( ! ( *more ) )
  {
    *more = true;
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = 0;
    }

    q = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      q = q + level_weight[i] * ( double ) ( x[i] );
    }

    if ( q_min < q && q <= q_max )
    {
      return;
    }
  }

  for ( ; ; )
  {
    j = 0;

    for ( ; ; )
    {
      if ( x[j] < x_max[j] )
      {
        break;
      }

      if ( dim_num - 1 <= j )
      {
        *more = false;
        return;
      }
      j = j + 1;
    }

    x[j] = x[j] + 1;
    for ( i = 0; i < j; i++ )
    {
      x[i] = 0;
    }

    q = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      q = q + level_weight[i] * ( double ) ( x[i] );
    }

    if ( q_min < q && q <= q_max )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_vcn_ordered 
(
  int dim_num,
  double level_weight[], 
  int x_max[],
  int x[],
  double q_min,
  double q_max,
  bool *more 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_ORDERED returns the "next" constrained vector, with ordering.
//
//  Discussion:
//
//    We consider vectors X of dimension DIM_NUM satisfying:
//
//      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
//
//    and define
//
//      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
//
//    and seek X's satisfying the constraint:
//
//      Q_MIN < Q <= Q_MAX
//
//    For sparse grid applications, we compute
//
//      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
//
//    and assume there is an underlying LEVEL used to index the sets of 
//    constrained vectors, and that 
//
//      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
//      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
//      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
//
//    This function returns, one at a time exactly those X which satisfy
//    the constraint.
//
//    A weak ordering is imposed on the solution vectors.  This function 
//    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
//    that the X vectors returned are roughly sorted (or at least binned) 
//    by Q value.
//
//  Example:
//
//    If the weights are also integral, then the X vectors are in fact SORTED 
//    by Q value:
//
//    LEVEL_WEIGHT:          1.000000        1.000000
//    Q_MIN:        0.000000
//    Q_MAX:        2.000000
//    X_MAX:                         2         2
//
//         1        1.000000         1         0
//         2        1.000000         0         1
//         3        2.000000         2         0
//         4        2.000000         1         1
//         5        2.000000         0         2
//
//    When the weights are not integral, then the X values are only BINNED
//    by Q value, that is, we first get all X's with Q values between Q_MIN
//    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
//
//    LEVEL_WEIGHT:             1.5               1
//    Q_MIN:  0.5
//    Q_MAX:  3
//    X_MAX:                           2         3
//
//           1             1.5         1         0
//           2               1         0         1
//           3             2.5         1         1
//           4               2         0         2
//           5               3         2         0
//           6               3         0         3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of components in the vector.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
//
//    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper
//    limits on the sum.
//
//    Input/output, bool *MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  double q;
  static double q_max2;
  static double q_min2;
//
//  On first call, initialize the subrange.
//
  if ( !(*more) )
  {
    q_min2 = q_min;
    q_max2 = webbur->r8_min ( q_min + 1.0, q_max );
  }
//
//  Call a lower level function to search the subrange.
//
  for ( ; ; )
  {
    sgmga_vcn ( dim_num, level_weight, x, q_min2, 
      q_max2, more );
//
//  If another solution was found, return it.
//
    if ( *more )
    {
      return;
    }
//
//  If the current subrange is exhausted, try to move to the next one.
//
    if ( q_max2 < q_max )
    {
      q_min2 = q_max2;
      q_max2 = webbur->r8_min ( q_max2 + 1.0, q_max );
    }
//
//  If there are no more subranges, we're done.
//
    else
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_vcn_ordered_naive 
(
  int dim_num,
  double level_weight[], 
  int x_max[],
  int x[],
  double q_min,
  double q_max,
  bool *more 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_ORDERED_NAIVE returns the "next" constrained vector, with ordering.
//
//  Discussion:
//
//    We consider vectors X of dimension DIM_NUM satisfying:
//
//      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
//
//    and define
//
//      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
//
//    and seek X's satisfying the constraint:
//
//      Q_MIN < Q <= Q_MAX
//
//    For sparse grid applications, we compute
//
//      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
//
//    and assume there is an underlying LEVEL used to index the sets of 
//    constrained vectors, and that 
//
//      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
//      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
//      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
//
//    This function returns, one at a time exactly those X which satisfy
//    the constraint.
//
//    A weak ordering is imposed on the solution vectors.  This function 
//    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
//    that the X vectors returned are roughly sorted (or at least binned) 
//    by Q value.
//
//  Example:
//
//    If the weights are also integral, then the X vectors are in fact SORTED 
//    by Q value:
//
//    LEVEL_WEIGHT:          1.000000        1.000000
//    Q_MIN:        0.000000
//    Q_MAX:        2.000000
//    X_MAX:                         2         2
//
//         1        1.000000         1         0
//         2        1.000000         0         1
//         3        2.000000         2         0
//         4        2.000000         1         1
//         5        2.000000         0         2
//
//    When the weights are not integral, then the X values are only BINNED
//    by Q value, that is, we first get all X's with Q values between Q_MIN
//    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
//
//    LEVEL_WEIGHT:             1.5               1
//    Q_MIN:  0.5
//    Q_MAX:  3
//    X_MAX:                           2         3
//
//           1             1.5         1         0
//           2               1         0         1
//           3             2.5         1         1
//           4               2         0         2
//           5               3         2         0
//           6               3         0         3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of components in the vector.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
//
//    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper
//    limits on the sum.
//
//    Input/output, bool *MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  double q;
  static double q_max2;
  static double q_min2;
//
//  On first call, initialize the subrange.
//
  if ( !(*more) )
  {
    q_min2 = q_min;
    q_max2 = webbur->r8_min ( q_min + 1.0, q_max );
  }
//
//  Call a lower level function to search the subrange.
//
  for ( ; ; )
  {
    sgmga_vcn_naive ( dim_num, level_weight, x_max, x, q_min2, 
      q_max2, more );
//
//  If another solution was found, return it.
//
    if ( *more )
    {
      return;
    }
//
//  If the current subrange is exhausted, try to move to the next one.
//
    if ( q_max2 < q_max )
    {
      q_min2 = q_max2;
      q_max2 = webbur->r8_min ( q_max2 + 1.0, q_max );
    }
//
//  If there are no more subranges, we're done.
//
    else
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void SandiaSGMGA::sgmga_weight 
(
  int dim_num,
  double level_weight[],
  int level_max, 
  void ( SandiaRules2::*gw_compute_weights[] ) ( int order, int dim, double w[] ),
  int point_num,
  int point_total_num,
  int sparse_unique_index[], 
  int growth, 
  int ( SandiaRules::*gw_compute_order[] ) ( int level, int growth ),
  double sparse_weight[] 
)

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_WEIGHT computes weights for an SGMGA grid.
//
//  Discussion:
//
//    The user must preallocate space for the output array SPARSE_WEIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int dim, double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, int POINT_NUM, the number of unique points 
//    in the grid. 
//
//    Input, int POINT_TOTAL_NUM, the total number of points 
//    in the grid. 
//
//    Input, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists, 
//    for each (nonunique) point, the corresponding index of the same point in 
//    the unique listing.
//
//    Input, int GROWTH, the growth rule.  
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
//    Output, double SPARSE_WEIGHT[POINT_NUM], the weights
//    associated with the sparse grid points.
//
{
  double coef;
  int dim;
  double *grid_weight;
  int level;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  bool more_grids;
  int order;
  int *order_1d;
  int order_nd;
  int point;
  int point_total;
  int point_unique;
  double q_max;
  double q_min;

  for ( point = 0; point < point_num; point++ )
  {
    sparse_weight[point] = 0.0;
  }

  point_total = 0;

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
//
//  Initialization for SGMGA_VCN_ORDERED.
//
  level_weight_min_pos = webbur->r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - webbur->r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = ( int ) ( webbur->r8_floor ( q_max / level_weight[dim] ) ) + 1;
      if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
      {
        level_1d_max[dim] = level_1d_max[dim] - 1;
      }
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = false;
//
//  Seek all vectors LEVEL_1D which satisfy the constraint:
//
//    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
//      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
//      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
//
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, 
      level_1d, q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
//
//  Compute the combinatorial coefficient.
//
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d, 
      q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
//
//  Transform each 1D level to a corresponding 1D order.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      order_1d[dim] = ((*webbur).*gw_compute_order[dim])(level_1d[dim],growth);
    }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
    order_nd = webbur->i4vec_product ( dim_num, order_1d );
//
//  Compute the weights for this grid.
//
//  The correct transfer of data from the product grid to the sparse grid
//  depends on the fact that the product rule weights are stored under colex
//  order of the points, and this is the same ordering implicitly used in
//  generating the SPARSE_UNIQUE_INDEX array.
//
    grid_weight = new double[order_nd];

    sgmga_product_weight ( dim_num, order_1d, order_nd, 
      gw_compute_weights, grid_weight );
//
//  Add these weights to the rule.
//
    for ( order = 0; order < order_nd; order++ )
    {
      point_unique = sparse_unique_index[point_total];

      point_total = point_total + 1;

      sparse_weight[point_unique] = sparse_weight[point_unique] 
        + coef * grid_weight[order];
    }

    delete [] grid_weight;
  }

  delete [] level_1d;
  delete [] level_1d_max;
  delete [] order_1d;

  return;
}

} // namespace ROL
