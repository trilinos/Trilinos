// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef  STK_LINEARINTERPOLATE_HPP
#define  STK_LINEARINTERPOLATE_HPP

#include <iterator>

#include <Intrepid_FieldContainer.hpp>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/environment/ReportHandler.hpp>


namespace stk {
namespace transfer {

template <class MESHA, class MESHB> class LinearInterpolate {

public :
  typedef Intrepid::FieldContainer<double> MDArray;

  typedef MESHA                            MeshA;
  typedef MESHB                            MeshB;
  typedef typename MeshA::EntityKey        EntityKeyA;
  typedef typename MeshB::EntityKey        EntityKeyB;
  typedef typename MeshA::EntityProc       EntityProcA;
  typedef typename MeshB::EntityProc       EntityProcB;

  typedef std::pair    <EntityProcB, EntityProcA>         EntityProcRelation;
  typedef std::vector  <EntityProcRelation>               EntityProcRelationVec;

  typedef std::multimap<EntityKeyB, EntityKeyA> EntityKeyMap;

  enum { Dimension = 3 };

  static void post_coarse_search_filter(EntityProcRelationVec &range_to_domain,
                                        const MeshA     &mesha,
                                        const MeshB     &meshb);

  static void filter_to_nearest(EntityKeyMap    &BtoA,
                                const MeshA     &FromPoints,
                                const MeshB     &ToPoints);

  static void apply (MeshB               &meshb,
                     const MeshA         &mesha,
                     const EntityKeyMap &RangeToDomain);

  static std::vector<double> solve_3_by_3_with_LU(const MDArray &M, const std::vector<double> &x);
private :

  static int LU_decomp(double A[9], int piv[3], int* sign);
  static int LU_solve(const double A[9], const int piv[3], double b[3]);

  static double distance_squared(const double *x, const double *y) ;
  static EntityKeyMap determine_best_fit(const typename EntityKeyMap::const_iterator begin,
                                         const typename EntityKeyMap::const_iterator end,
                                         const EntityKeyB          current,
                                         const MeshB             &ToPoints,
                                         const MeshA             &FromPoints);
  static bool co_linear(const typename EntityKeyMap::value_type &e,
                        const EntityKeyMap &first_three,
                        const MeshA        &FromPoints);
  static bool co_planer(const typename EntityKeyMap::value_type &e,
                        const EntityKeyMap &first_three,
                        const MeshA        &FromPoints);

};

template <class MESHA, class MESHB> double LinearInterpolate<MESHA,MESHB>::distance_squared(const double *x, const double *y) {
  double d = 0;
  for (unsigned i=0; i<Dimension; ++i) d += (x[i]-y[i])*(x[i]-y[i]);
  return d;
}

template <class MESHA, class MESHB> bool LinearInterpolate<MESHA,MESHB>::co_linear(
                               const typename EntityKeyMap::value_type &e,
                               const EntityKeyMap &first_two,
                               const MeshA        &FromPoints) {

  const unsigned span = Dimension-1;
  const double *Corners[Dimension]={0};
  {
    unsigned n=0;
    for (typename EntityKeyMap::const_iterator j=first_two.begin(); j != first_two.end(); ++j,++n) {
      const EntityKeyA from_key = j->second;
      Corners[n] = FromPoints.coord(from_key);
    }
    Corners[n] = FromPoints.coord(e.second);
  }
  double S[span][Dimension];
  for (unsigned j=1; j<Dimension; ++j)
    for (unsigned k=0; k<Dimension; ++k)
      S[j-1][k] = Corners[j][k] - Corners[0][k];

  const double a_cross_b[Dimension] = { S[0][1]*S[1][2] - S[0][2]*S[1][1] ,
                                     -( S[0][0]*S[1][2] - S[0][2]*S[1][0]),
                                        S[0][0]*S[1][1] - S[0][1]*S[1][0]};
  const double norm_a_cross_b = a_cross_b[0]*a_cross_b[0] + a_cross_b[1]*a_cross_b[1] + a_cross_b[2]*a_cross_b[2];
  const bool ret = norm_a_cross_b < .00001;

  return ret;
}

template <class MESHA, class MESHB> bool LinearInterpolate<MESHA,MESHB>::co_planer(
                               const typename EntityKeyMap::value_type &e,
                               const EntityKeyMap &first_three,
                               const MeshA        &FromPoints) {

  const unsigned span = Dimension+1;
  const double *Corners[span]={0};
  {
    unsigned n=0;
    for (typename EntityKeyMap::const_iterator j=first_three.begin(); j != first_three.end(); ++j,++n) {
      const EntityKeyA from_key = j->second;
      Corners[n] = FromPoints.coord(from_key);
    }
    Corners[n] = FromPoints.coord(e.second);
  }
  double S[Dimension][Dimension];
  for (unsigned j=1; j<span; ++j)
    for (unsigned k=0; k<Dimension; ++k)
      S[j-1][k] = Corners[j][k] - Corners[0][k];

  const double a_cross_b[Dimension] = { S[0][1]*S[1][2] - S[0][2]*S[1][1] ,
                                     -( S[0][0]*S[1][2] - S[0][2]*S[1][0]),
                                        S[0][0]*S[1][1] - S[0][1]*S[1][0]};
  const double a_cross_b_dot_c = a_cross_b[0]*S[2][0] + a_cross_b[1]*S[2][1] + a_cross_b[2]*S[2][2];
  const bool ret = std::abs(a_cross_b_dot_c) < .00001;

  return ret;
}

template <class MESHA, class MESHB> typename LinearInterpolate<MESHA,MESHB>::EntityKeyMap LinearInterpolate<MESHA,MESHB>::determine_best_fit(
                               const typename EntityKeyMap::const_iterator begin,
                               const typename EntityKeyMap::const_iterator end,
                               const EntityKeyB          current,
                               const MeshB             &ToPoints,
                               const MeshA             &FromPoints) {

  typedef  std::multimap<double, typename EntityKeyMap::value_type, std::less<double> > DIST_MAP;
  DIST_MAP  sorted;

  const double *  to_coords =   ToPoints.coord(current);

  for (typename EntityKeyMap::const_iterator d=begin; d != end; ++d) {
    const EntityKeyA domain_index = d->second;
    const double *from_coords = FromPoints.coord(domain_index);

    const double dist = distance_squared(from_coords, to_coords);
    const typename DIST_MAP::value_type val(dist,*d);
    sorted.insert(val);
  }

  EntityKeyMap ret;
  for (typename DIST_MAP::const_iterator d=sorted.begin(); d != sorted.end() && ret.size() <= Dimension; ++d) {
    if      (ret.size()  < Dimension-1)            ret.insert(d->second);
    else if (ret.size() == Dimension-1) {
      if    (!co_linear(d->second,ret,FromPoints)) ret.insert(d->second);
    }else if(!co_planer(d->second,ret,FromPoints)) ret.insert(d->second);
  }
  if (ret.size() < Dimension+1) ret.clear();
  return ret;
}

template <class MESHA, class MESHB> void LinearInterpolate<MESHA,MESHB>::post_coarse_search_filter(
            EntityProcRelationVec &BtoA,
            const MeshA           &mesha,
            const MeshB           &meshb){
  const unsigned p_rank = parallel_machine_rank(mesha.comm());
  std::sort(BtoA.begin(), BtoA.end());
  typedef typename EntityProcRelationVec::iterator iterator;
  iterator k=BtoA.begin();
  for (iterator i=BtoA.begin(),j=BtoA.begin(); j!=BtoA.end();) {
    while (j!=BtoA.end() && i->first == j->first) {
      ++j;
    }
    const unsigned num = j-i;
    if (Dimension+2 < num || p_rank != i->first.proc()) {
      while (i!=j && i != BtoA.end()) {
        *k++ = *i++;
      }
    }
    else {
      i=j;
    }
  }
  BtoA.resize(k-BtoA.begin());
}

template <class MESHA, class MESHB> void LinearInterpolate<MESHA,MESHB>::filter_to_nearest(
                                    EntityKeyMap  &BtoA,
                                    const MeshA   &mesha,
                                    const MeshB   &meshb) {
  typedef typename EntityKeyMap::iterator iterator;
  for (iterator j=BtoA.begin(); j!=BtoA.end(); ) {
    std::pair<iterator, iterator> keys=BtoA.equal_range(j->first);
    const unsigned num = distance(keys.first, keys.second);
    ThrowRequireMsg (Dimension <  num,
      __FILE__<<":"<<__LINE__<<" Expected "<<Dimension+1<<" relations."<<" Found:"<<num<<" for Key:"<<j->first);
    EntityKeyMap n = determine_best_fit(keys.first, keys.second, j->first, meshb, mesha);
    BtoA.erase(keys.first, keys.second);
    BtoA.insert(n.begin(), n.end());
    j = keys.second;
  }
}


template <class MESHA, class MESHB>  void LinearInterpolate<MESHA,MESHB>::apply (
                          MeshB        &meshb,
                    const MeshA        &mesha,
                    const EntityKeyMap &RangeToDomain){

  typedef typename EntityKeyMap::const_iterator map_const_iterator;
  const unsigned span = Dimension+1;

  const unsigned numValsa = mesha.num_values();
  const unsigned numValsb = meshb.num_values();
    ThrowRequireMsg (numValsb == numValsa,
      __FILE__<<":"<<__LINE__<<" Found "<<numValsa<<" values for mesh a and "<<numValsb<<" for mesh b."
      <<" These should be the same.");


  for (map_const_iterator j,i=RangeToDomain.begin(); i!=RangeToDomain.end(); i=j)  {
    j = i;
    while (j!=RangeToDomain.end() && i->first == j->first) ++j;
    const unsigned num_relations = distance(i, j);

    const EntityKeyB to_key = i->first;
    ThrowRequireMsg (span == num_relations,
      __FILE__<<":"<<__LINE__<<" Expected "<<span<<" relations."<<" Found:"<<num_relations<<" for Key:"<<to_key);

    MDArray Corners(span,Dimension);
    {
      unsigned n=0;
      for (map_const_iterator k=i; k!=j; ++k,++n) {
        const EntityKeyA from_key = k->second;
        const double *c = mesha.coord(from_key);
        for (unsigned kk=0; kk<Dimension; ++kk) Corners(n,kk) = c[kk];
      }
    }

    MDArray SpanVectors(Dimension,Dimension);
    for (unsigned l=1; l<span; ++l)
      for (unsigned k=0; k<Dimension; ++k)
        SpanVectors(k,l-1) = Corners(l,k) - Corners(0,k);
    std::vector<double> point(Dimension);
    {
      const double  *c = meshb.coord(to_key);
      for (unsigned k=0; k<Dimension; ++k) point[k] = c[k]-Corners(0,k);
    }

    std::vector<double> S;
    S = solve_3_by_3_with_LU(SpanVectors, point);

    std::vector<double> Values(span);
    for (unsigned f=0; f<numValsb; ++f)  {
      const EntityKeyA from_key = i->second;
      const unsigned   to_field_size =   meshb.value_size(to_key,   f);
      const unsigned from_field_size =   mesha.value_size(from_key, f);
      const unsigned field_size = std::min(to_field_size, from_field_size);
      for (unsigned n=0; n<field_size; ++n) {
        unsigned m=0;
        for (map_const_iterator k=i; k!= j; ++k,++m) {
          const EntityKeyA this_from_key = k->second;
          const double *c = mesha.value(this_from_key,f);
          Values[m] = c[n];
        }
        // So, we have choosen corner 0 as the base of the span
        // vectors and determined local (or parametric)
        // coordinates of the target point in the span vector
        // coordinate system.  These are stored in S. S is
        // dimension Dimension and there are Dimension+1 corner values.
        // The scalar used to scale the value at corner 0
        // is 1-sum(S).
        double T=1;
        for (unsigned k=0; k<Dimension; ++k) T -= S[k];
        double interpolated_value = T * Values[0];
        for (unsigned k=0; k<Dimension; ++k)
          interpolated_value += S[k] * Values[k+1];
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
template <class MESHA, class MESHB> int LinearInterpolate<MESHA,MESHB>::LU_decomp(double A[9], int piv[3], int* sign)
{
  piv[0] = 0; piv[1] = 1; piv[2] = 2;

  double m;

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
template <class MESHA, class MESHB> int LinearInterpolate<MESHA,MESHB>::LU_solve(const double A[9], const int piv[3], double b[3])
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


template <class MESHA, class MESHB> std::vector<double> LinearInterpolate<MESHA,MESHB>::solve_3_by_3_with_LU(
         const MDArray             &M,
         const std::vector<double> &x) {
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
//check
double a[3]={0};
for (unsigned i=0; i<3; ++i)
for (unsigned j=0; j<3; ++j) a[i] += M(i,j)*v[j];
for (unsigned i=0; i<3; ++i) a[i] -= x[i];
const double norm = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
if (.00001 < norm) std::cout <<__FILE__<<":"<<__LINE__<<" Error in solve_3_by_3_with_LU"<<std::endl;


  return v;
}


} }
#endif

