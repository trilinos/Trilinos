/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <stk_transfer/TransferP2P.hpp>

namespace {

template <unsigned DIM, typename Scalar>
Scalar distance_squared (const Scalar x[DIM], const Scalar y[DIM]) {
  Scalar d = 0;  
  for (unsigned i=0; i<DIM; ++i) d += (x[i]-y[i])*(x[i]-y[i]);
  return d;
}

}

// Template Specializations based on spatial dimension:
template void STK_TransferP2P::point_to_point_coarse_search<3>(IdentProcRelation &,
                                                               const PointMap    &,
                                                               const MDArray     &,
                                                               const BoundingBox::Data ,
                                                               const stk::ParallelMachine );
template void STK_TransferP2P::linear_interpolation<3> (MDArray &,
                                                        const MDArray &,
                                                        const IdentProcRelation &,
                                                        const MDArray &,
                                                        const MDArray &,
                                                        const stk::ParallelMachine  );

template void STK_TransferP2P::filter_with_fine_search<3>(IdentProcRelation &,
                                                          const PointMap &,
                                                          const MDArray &,
                                                          const stk::ParallelMachine  ) ;

template void STK_TransferP2P::delete_range_points_found<3>(PointMap &,
                                                         const IdentProcRelation &,
                                                         const stk::ParallelMachine  ) ;

template  void STK_TransferP2P::convert_to_map<3>(PointMap &, const MDArray &);


namespace STK_TransferP2P {


/*********************   LU_decomp()   *******************************
 *
 *     Decompose matrix A into the product of a unit lower triangular matrix, L,
 *     and an upper triangular matrix, U, with column pivoting.
 *             
 *     The matrix A is assumed to be stored in an array like this:
 *                     
 *                                    A[0]  A[1]  A[2]
 *                                    A[3]  A[4]  A[5]
 *                                    A[6]  A[7]  A[8]
 *     
 *     Upon completion, the entries A[3], A[6], A[7] of A are overwritten by
 *     the strictly lower portion of L.  And the rest of A is overwritten by
 *     the upper portion of U.  The pivots are stored in piv.  A plus one or
 *     minus one is returned in "sign" to allow computation of the determinant
 *     .. just multiply the trace times this value.
 *    
 *     A value of 0 (zero) is returned if the matrix is singular (assuming NDEBUG
 *     is not defined.)
 *                                                                                                          
 *     There are 3 divisions, 5 multiplications, and 5 add/sub
 *     plus pivoting compares and swaps.
 *                                                                                                                      
 *
 *                                             Drake 9-98
 */
int LU_decomp(double A[9], int piv[3], int* sign)
{
  piv[0] = 0; piv[1] = 1; piv[2] = 2;
  
  register double m;
  
#ifndef NDEBUG
  if ( A[0] == 0.0 && A[3] == 0.0 && A[6] == 0.0 ) return 0;
#endif
  
  (*sign) = 1;
  if ( (m = fabs(A[0])) < fabs(A[3]) ) { 
    if (fabs(A[3]) < fabs(A[6]))
    {   
      piv[0] = 2; piv[2] = 0;             // Switch rows 0 and 2.
      m = A[0]; A[0] = A[6]; A[6] = m;
      m = A[1]; A[1] = A[7]; A[7] = m;
      m = A[2]; A[2] = A[8]; A[8] = m;
    }   
    else
    {   
      piv[0] = 1; piv[1] = 0;             // Switch rows 0 and 1.
      m = A[0]; A[0] = A[3]; A[3] = m;
      m = A[1]; A[1] = A[4]; A[4] = m;
      m = A[2]; A[2] = A[5]; A[5] = m;
    }   
    (*sign) = -1; 
  }
  else if (m < fabs(A[6]))
    {   
      piv[0] = 2; piv[2] = 0;             // Switch rows 0 and 2.
      m = A[0]; A[0] = A[6]; A[6] = m;
      m = A[1]; A[1] = A[7]; A[7] = m;
      m = A[2]; A[2] = A[8]; A[8] = m;
      (*sign) = -1;
    }

  A[3] = m = A[3]/A[0];           /* multiplier m21 which is L21 */
  A[4] = A[4] - A[1] * m;         /* U22 */
  A[5] = A[5] - A[2] * m;         /* U23 */

  A[6] = m = A[6]/A[0];           /* multiplier m31 which is L31 */
  A[7] = A[7] - A[1] * m;
  A[8] = A[8] - A[2] * m;         /* intermediate U33 */

#ifndef NDEBUG
  if ( A[4] == 0.0 && A[7] == 0.0) return 0;
#endif

  if (fabs(A[4]) < fabs(A[7])) {
    int tmpPivot = piv[1];
    piv[1] = piv[2]; piv[2] = tmpPivot;
    m = A[3]; A[3] = A[6]; A[6] = m;     // Switch rows 1 and 2.
    m = A[4]; A[4] = A[7]; A[7] = m;
    m = A[5]; A[5] = A[8]; A[8] = m;
    (*sign) = -(*sign);
  }

  A[7] = m = A[7]/A[4];           /* multiplier m32 which is L32 */
  A[8] = A[8] - A[5] * m;         /* U33 */

#ifndef NDEBUG
  if ( A[8] == 0.0 ) return 0;
#endif

  return 1;
}

/*********************   LU_solve()   *******************************
 *     Uses a precomputed LU decomposition to find the solution to Ax = b for
 *         
 *     a length 3 vector b.
 *                 
 *     It solves Ly = b for vector y, then Ux = y for x.
 *                         
 *     The matrix A is assumed to be of the same form as the output from
 *     LU_decomp().
 *                                     
 *     The vector b is overwritten with the solution.
 *                                             
 *     The diagonals are checked for zeros if NDEBUG is not defined.
 *                                                     
 *     There are 3 div, 6 mult, and 6 add/sub plus compares and swaps to undo
 *     the pivoting.
 *                                                                 
 *     Drake 9/98
 *     
 */
int LU_solve(const double A[9], const int piv[3], double b[3])
{
#ifndef NDEBUG
  if (A[0] == 0.0 || A[4] == 0.0 || A[8] == 0.0) return 0;
#endif

  // Apply inverse of permutation matrix to RHS first.
    int i = piv[0];
  if (i != 0)
  {
    double t = b[0]; b[0] = b[i]; b[i] = t;
  }
  i = piv[2];
  if (i == 1 || (i == 0 && piv[1] == 2) )
  {
    double t = b[1]; b[1] = b[2]; b[2] = t;
  }
  b[1] = b[1] - A[3] * b[0];                 /* y2 = b2-m21*b1   (y1 = b[0]) */
  b[2] = b[2] - A[6] * b[0] - A[7] * b[1];   /* y3 = b3-m31*y1-m32*y2        */

  b[2] = b[2] / A[8];                        /* x3 = y3/a33                  */
  b[1] = (b[1] - A[5]*b[2]) / A[4];          /* x2 = (b2-a23*x3)/a22         */
  b[0] = (b[0] - A[1]*b[1] - A[2]*b[2]) / A[0];
                                             /* x1 = (b1-a12*x2-a13*x3)/a11  */

  return 1;
}


std::vector<double> solve_3_by_3_with_LU(const MDArray             M, 
                                         const std::vector<double> x) {
  double A[9];
  double b[3];
  int piv[3];
  int sign;
  for (unsigned i=0,k=0; i<3; ++i) 
    for (unsigned j=0; j<3; ++j) 
      A[k++] = M(i,j);
  for (unsigned i=0; i<3; ++i) b[i] = x[i];
  LU_decomp(A, piv, &sign);
  LU_solve (A, piv, b);
  std::vector<double> v(b,b+3);
  return v;
}

template <unsigned DIM>
void point_to_point_coarse_search(IdentProcRelation &RangeToDomain,
                                  const PointMap    &ToPoints,
                                  const MDArray     &FromPoints,
                                  const BoundingBox::Data radius,
                                  const stk::ParallelMachine comm) {
  
  const unsigned NumDomainPts = FromPoints.dimension(0);

  std::vector<BoundingBox> range_vector;
  std::vector<BoundingBox> domain_vector;

  BoundingBox::Data center[DIM];

  for (unsigned i=0; i<NumDomainPts; ++i) {
    for (unsigned j=0; j<DIM; ++j) {
      center[j] = FromPoints(i,j);
    }
    const BoundingBox::Key key(i, stk::parallel_machine_rank(comm));
    BoundingBox B(center, radius, key);
    domain_vector.push_back(B); 
  }

  for (PointMap::const_iterator i=ToPoints.begin(); i!=ToPoints.end(); ++i) {
    for (unsigned j=0; j<DIM; ++j) {
      center[j] = i->second[j];
    }
    const BoundingBox::Key key(i->first, stk::parallel_machine_rank(comm));
    BoundingBox B(center, radius, key);
    range_vector.push_back(B); 
  }

  stk::search::FactoryOrder order;
  order.m_communicator = comm;

  // Slightly confusing: coarse_search documentation has domain->range
  // relations sorted by domain key.  We want range->domain type relations
  // sorted on range key. It might appear we have the arguments revered
  // in coarse_search call, but really, this is what we want.
  stk::search::coarse_search(RangeToDomain, domain_vector, range_vector, order);
}



template <unsigned DIM>
void linear_interpolation (MDArray &ToValues,
                    const MDArray &FromValues,
                    const IdentProcRelation &RangeToDomain,
                    const MDArray &ToPoints,
                    const MDArray &FromPoints,
                    const stk::ParallelMachine  comm) {
  const unsigned numValues = ToValues.size();
  for (unsigned i=0; i<numValues; ++i)  {
    MDArray Corners(DIM+1,DIM);
    for (unsigned j=0; j<DIM+1; ++j) {
      const unsigned c = RangeToDomain[4*i+j].second.ident;
      for (unsigned k=0; k<DIM; ++k) Corners(j,k) = FromPoints(c,k);
    } 
    MDArray SpanVectors(DIM,DIM);
    for (unsigned j=1; j<DIM+1; ++j) 
      for (unsigned k=0; k<DIM; ++k) 
        SpanVectors(k,j-1) = Corners(j,k) - Corners(0,k);
    std::vector<double> point(DIM);
    const unsigned p = RangeToDomain[4*i].first.ident;
    for (unsigned k=0; k<DIM; ++k) point[k] = ToPoints(p,k)-Corners(0,k);

    std::vector<double> S;
    if (3==DIM) S = solve_3_by_3_with_LU(SpanVectors, point);
    else std::cerr<<__FILE__<<":"<<__LINE__<<" Only 3D supported."<<std::endl;

    std::vector<double> Values(DIM+1);
    for (unsigned j=0; j<DIM+1; ++j) {
      const unsigned c = RangeToDomain[4*i+j].second.ident;
      Values[j] = FromValues(c,0);
    }
    // So, we have choosen corner 0 as the base of the span
    // vectors and determined local (or parametric)
    // coordinates of the target point in the span vector
    // coordinate system.  These are stored in S. S is
    // dimension DIM and there are DIM+1 corner values.
    // The scalar used to scale the value at corner 0 
    // is 1-sum(S).
    double T=1;
    for (unsigned j=0; j<DIM; ++j) T -= S[j];
    double interpolated_value = T * Values[0];
    for (unsigned j=0; j<DIM; ++j)
      interpolated_value += S[j] * Values[j+1]; 
    ToValues[i] = interpolated_value; 
  }
}

  
// Want to find the N best elements to keep.  
template <unsigned DIM>
void filter_with_fine_search(IdentProcRelation &RangeToDomain,
                                const PointMap &ToPoints,
                                const MDArray &FromPoints,
                                const stk::ParallelMachine  comm) {

  typedef  std::map<double,unsigned,std::greater<double> > DIST_MAP;
  const unsigned NumRelations = RangeToDomain.size();

  DIST_MAP  smallest_distances;
  IdentProcRelation new_relations;

  MDArray::scalar_type from_center[DIM];
  MDArray::scalar_type   to_center[DIM];
  // One of those duel sorted indexing problems:
  for (IdentProc::Key i=0, j=0; j<NumRelations;) {
    const IdentProc::Key  range_index = RangeToDomain[j].first .ident;
    const IdentProc::Key domain_index = RangeToDomain[j].second.ident;
    if (i < range_index) {
      if (DIM+1 == smallest_distances.size()) {
        // Found enough points for linear interpolation.
        for (DIST_MAP::const_iterator d=smallest_distances.begin(); 
             d != smallest_distances.end(); ++d) {
          new_relations.push_back(RangeToDomain[d->second]);
        }
      } 
      smallest_distances.clear();
      ++i;
    } else {
      for (unsigned k=0; k<DIM; ++k) from_center[k] = FromPoints(domain_index,k);
      for (unsigned k=0; k<DIM; ++k)   to_center[k] =   ToPoints.find(range_index)->second[k];
      const double dist = distance_squared<DIM>(from_center, to_center);
      const DIST_MAP::value_type val(dist,j);

      // Is using a map too much memory allocation/deallocation? Maybe std::vector is better.
      if (smallest_distances.insert(val).second && DIM+1 < smallest_distances.size()) 
        smallest_distances.erase(smallest_distances.begin()); // DIST_MAP: Largest at begin()
      ++j;
    }
  }
  if (DIM+1 == smallest_distances.size()) {
    for (DIST_MAP::const_iterator d=smallest_distances.begin(); 
         d != smallest_distances.end(); ++d) {
      new_relations.push_back(RangeToDomain[d->second]);
    }
  } 
  RangeToDomain.swap(new_relations);
}

template <unsigned DIM>
void delete_range_points_found(PointMap &ToPoints,
                               const IdentProcRelation &RangeToDomain,
                               const stk::ParallelMachine  comm) {
  const unsigned NumRelations = RangeToDomain.size();

  for (unsigned i=0; i<NumRelations; i+=DIM+1)  {
    const PointMap::key_type ident = RangeToDomain[i].first.ident;
    ToPoints.erase(ident);
  }
}


template <unsigned DIM>
void convert_to_map(PointMap      &map_points, 
                    const MDArray &Points) {
  const unsigned NumPts = Points.dimension(0);
  PointMap::iterator at=map_points.begin();
  std::vector<double> v(DIM);
  for (unsigned i=0; i<NumPts; ++i) {
    for (unsigned j=0; j<DIM; ++j) v[j] = Points(i,j);
    const PointMap::value_type val(i,v);
    at = map_points.insert(at,val);
  }
}

}
