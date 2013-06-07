/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_TransferP2P_hpp
#define  STK_TransferP2P_hpp

#include <boost/smart_ptr/shared_ptr.hpp>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_mesh/base/Comm.hpp>


namespace stk {
namespace transfer {

template <class MESHA, class MESHB> class LinearInterpoate {

public :
typedef MESHA MeshA;
typedef MESHB MeshB;
typedef MeshA::EntityKey EntityKeyA;
typedef MeshB::EntityKey EntityKeyB;

typedef std::multimap<EntityKeyB, EntityKeyA> EntityKeyMap;

enum { Dim = MESHA::DIM };

static void filter_with_fine_search(EntityKeyMap  &rel,
                                    const MeshB   &ToPoints,
                                    const MeshA &FromPoints);

static void apply (MeshB        &meshb,
                    const MeshA         &mesha,
                    const EntityKeyMap &RangeToDomain);

private :

static int LU_decomp(double A[9], int piv[3], int* sign);
static int LU_solve(const double A[9], const int piv[3], double b[3])
static std::vector<double> solve_3_by_3_with_LU(const MDArray M, const std::vector<double> x);

enum { dim_eq = stk::StaticAssert<MeshB::Dimension==MeshA::Dimension::OK };
enum { dim_3  = stk::StaticAssert<               3==MeshA::Dimension::OK };

static double distance_squared(const double *x, const double *y) ;
static EntityKeyMap determine_best_fit(EntityKeyMap::const_iterator begin,
                               EntityKeyMap::const_iterator end,
                               const EntityKeyB          current,
                                 const MeshB             &ToPoints,
                                 const MeshA             &FromPoints) const;

};

template <class MESHA, class MESHB> double LinearInterpoate<MESHA,MESHB>::distance_squared(const double *x, const double *y) {
  double d = 0;
  for (unsigned i=0; i<Dim; ++i) d += (x[i]-y[i])*(x[i]-y[i]);
  return d;
}

template <class MESHA, class MESHB> EntityKeyMap LinearInterpoate<MESHA,MESHB>::determine_best_fit(
                               EntityKeyMap::const_iterator begin,
                               EntityKeyMap::const_iterator end,
                               const EntityKeyB          current,
                                 const MeshB             &ToPoints,
                                 const MeshA             &FromPoints) const {

  typedef  std::multimap<double,EntityKeyMap::value_type,std::greater<double> > DIST_MAP;
  DIST_MAP  smallest;

  const double *  to_coords =   ToPoints.coord(current);

  for (EntityKeyMap::const_iterator d=begin; d != end; ++d) {
    const EntityKeyA domain_index = d->second;
    const double *from_coords = FromPoints.coord(domain_index);

    const double dist = distance_squared(from_coords, to_coords);

    const DIST_MAP::value_type val(dist,*d);
    smallest.insert(val);
    // Is using a map too much memory allocation/deallocation? Maybe std::vector is better.
    if (Dim+1 < smallest.size()) smallest.erase(smallest.begin()); // DIST_MAP: Largest at begin()

  }
  EntityKeyMap ret;
  for (DIST_MAP::const_iterator d=smallest.begin(); d != smallest.end(); ++d) ret.insert(d->second);
  return ret;
}
template <class MESHA, class MESHB> void LinearInterpoate<MESHA,MESHB>::filter_with_fine_search(
                                    EntityKeyMap  &rel,
                                    const MeshB   &ToPoints,
                                    const MeshA &FromPoints) {

  for (EntityKeyMap::const_iterator j=rel.begin(); j!=rel.end(); ++j) {
    std::pair<EntityKeyMap::const_iterator,EntityKeyMap::const_iterator> keys=rel.equal_range(j->first);
    const unsigned num = distance(keys.first, keys.second);
    if (num <= Dim) {
      rel.erase(keys.first, keys.second);
    } else if (Dim+1 < num) {
      EntityKeyMap n = determine_best_fit(keys.first, keys.second, j->first, ToPoints, FromPoints);
      rel.erase(keys.first, keys.second);
      rel.insert(n.begin(), n.end()); 
    }   
    j = rel.second;
  }
}


template <class MESHA, class MESHB>  void LinearInterpoate<MESHA,MESHB>::apply (
                          MeshB        &meshb,
                    const MeshA         &mesha,
                    const EntityKeyMap &RangeToDomain){

  const unsigned dim  = MeshB::Dimension;
  const unsigned span = dim+1;

  MeshB::EntityKeyVec keysb;
  meshb.Keys(keysb);
  const unsigned numKeysb = keysb.size();
  const unsigned numValsb = meshb.num_values();

  for (unsigned i=0; i<numKeysb; ++i)  {
    const EntityKeyB to_key = keysb[i];
    const pair<EntityKeyMap::const_iterator, EntityKeyMap::const_iterator> from_keys = RangeToDomain.equal_range(to_key);
    const unsigned num_relations = distance(from_keys.first, from_keys.second);
    ThrowRequireMsg (span == num_relations,  
      __FILE__<<":"<<__LINE__<<" Expected "<<span<<" relations."<<" Found:"<<num_relations);

    MDArray Corners(span,dim);
    {
      unsigned n=0;
      for (EntityKeyMap const_iterator j=from_keys.first; j!= from_keys.second; ++j,++n) {
        const EntityKeyA from_key = j->second;
        const double *c = mesha.coord(from_key); 
        for (unsigned k=0; k<dim; ++k) Corners(n,k) = c[k];
      } 
    }
    MDArray SpanVectors(dim,dim);
    for (unsigned j=1; j<span; ++j) 
      for (unsigned k=0; k<dim; ++k) 
        SpanVectors(k,j-1) = Corners(j,k) - Corners(0,k);
    std::vector<double> point(dim);
    {
      const double  *c = meshb.coord(to_key);
      for (unsigned k=0; k<DIM; ++k) point[k] = c[k]-Corners(0,k);
    }

    std::vector<double> S;
    S = solve_3_by_3_with_LU(SpanVectors, point);

    std::vector<double> Values();
    for (unsigned f=0; f<numValsb; ++f)  {
      const unsigned   to_field_size =   meshb.value_size(to_key, f);
      const unsigned from_field_size =   mesha.value_size(*from_keys.first, f);
      const unsigned field_size = std::min(to_field_size, from_field_size);
      for (unsigned n=0; n<field_size; ++n) {
        unsigned m=0;
        for (EntityKeyMap const_iterator j=from_keys.first; j!= from_keys.second; ++j,++m) {
          const EntityKeyA from_key = j->second;
          const double *c = mesha.value(from_key,f);
          Values[m] = c[n];
        }
        // So, we have choosen corner 0 as the base of the span
        // vectors and determined local (or parametric)
        // coordinates of the target point in the span vector
        // coordinate system.  These are stored in S. S is
        // dimension DIM and there are DIM+1 corner values.
        // The scalar used to scale the value at corner 0 
        // is 1-sum(S).
        double T=1;
        for (unsigned j=0; j<dim; ++j) T -= S[j];
        double interpolated_value = T * Values[0];
        for (unsigned j=0; j<dim; ++j)
          interpolated_value += S[j] * Values[j+1]; 
        double  *c = meshb.value(to_key, f);
        c[n] = interpolated_value; 
      }
    }
  }
}


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
template <class MESHA, class MESHB> int LinearInterpoate<MESHA,MESHB>::LU_decomp(double A[9], int piv[3], int* sign)
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
template <class MESHA, class MESHB> int LinearInterpoate<MESHA,MESHB>::LU_solve(const double A[9], const int piv[3], double b[3])
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


template <class MESHA, class MESHB> std::vector<double> LinearInterpoate<MESHA,MESHB>::solve_3_by_3_with_LU(
         const MDArray             M,
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


} }
#endif

