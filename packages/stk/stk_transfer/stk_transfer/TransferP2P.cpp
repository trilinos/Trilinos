/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>

#include <stk_util/parallel/Parallel.hpp>

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
namespace STK_TransferP2P {


// Template Specializations based on spatial dimension:
template void point_to_point_coarse_search<3, MDArray>(IdentProcRelation &,
                                              const PointMap    &,
                                              const MDArray     &,
                                              const BoundingBox::Data ,
                                              const stk::ParallelMachine );
template void point_to_point_coarse_search<3, STKMesh>(IdentProcRelation &,
                                              const PointMap    &,
                                              const STKMesh     &,
                                              const BoundingBox::Data ,
                                              const stk::ParallelMachine );
template void linear_interpolation<3,MDArray> (MDArray &,
                                       const MDArray &,
                                       const IdentProcRelation &,
                                       const MDArray &,
                                       const MDArray &,
                                       const stk::ParallelMachine  );

template void linear_interpolation<3,STKMesh> (MDArray &,
                                       const STKMesh &,
                                       const IdentProcRelation &,
                                       const MDArray &,
                                       const STKMesh &,
                                       const stk::ParallelMachine  );

template void filter_with_fine_search<3,MDArray> (IdentProcRelation &,
                                                  const PointMap &,
                                                  const MDArray  &,
                                                  const stk::ParallelMachine) ;

template void filter_with_fine_search<3,STKMesh> (IdentProcRelation &,
                                                  const PointMap &,
                                                  const STKMesh  &,
                                                  const stk::ParallelMachine) ;
template void delete_range_points_found<3>(PointMap &,
                                           const IdentProcRelation &,
                                           const stk::ParallelMachine  ) ;

template void convert_to_map<3>(PointMap &, const MDArray &);



STKMesh::STKMesh(MetaDataPtr meta,
          BulkDataPtr bulk,
          EntityVec   &ent,
          VectorFieldType &coord,
          VectorFieldType &val) :
    meta_data        (meta),
    bulk_data        (bulk),
    entities         (ent),
    coordinates_field(coord),
    values_field     (val) {}
STKMesh::STKMesh(const STKMesh &M) :
    meta_data        (M.meta_data),
    bulk_data        (M.bulk_data),
    entities         (M.entities),
    coordinates_field(M.coordinates_field),
    values_field     (M.values_field) {}
STKMesh::~STKMesh(){}
BulkDataPtr &STKMesh::BulkData() {return bulk_data;}
STKMesh::iterator STKMesh::begin() const {return entities.begin();}
STKMesh::iterator   STKMesh::end() const {return entities.end();}
const double *STKMesh::Coord(const STKMesh::iterator i) const {
  const double *coord = bulk_data->field_data(coordinates_field, *i);
  return coord;
}
const double *STKMesh::Coord(const stk::mesh::EntityKey i) const {
  const stk::mesh::Entity  entity = bulk_data->get_entity(i);
  const double *coord = bulk_data->field_data(coordinates_field, entity);
  return  coord;
}
const double *STKMesh::Value(const stk::mesh::EntityKey i) const {
  const stk::mesh::Entity  entity = bulk_data->get_entity(i);
  const double *value = bulk_data->field_data(values_field, entity);
  return  value;
}
stk::mesh::EntityKey STKMesh::Key(STKMesh::iterator i) const {
  return bulk_data->entity_key(*i);
}
VectorFieldType &STKMesh::Coord(){return coordinates_field;}
VectorFieldType &STKMesh::Value(){return values_field;}



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

namespace {

template <unsigned DIM>
void boundingbox_vector(std::vector<BoundingBox> &vector,
                        const PointMap           &mesh,
                        const BoundingBox::Data   radius,
                        const stk::ParallelMachine comm) {
  BoundingBox::Data center[DIM];

  for (PointMap::const_iterator i=mesh.begin(); i!=mesh.end(); ++i) {
    for (unsigned j=0; j<DIM; ++j) center[j] = i->second[j];
    const IdentProc::Key Id(0,i->first);
    const BoundingBox::Key key(Id, stk::parallel_machine_rank(comm));
    BoundingBox B(center, radius, key);
    vector.push_back(B);
  }
}

template <unsigned DIM>
void boundingbox_vector(std::vector<BoundingBox> &vector,
                        const MDArray            &mesh,
                        const BoundingBox::Data   radius,
                        const stk::ParallelMachine comm) {

  const unsigned NumDomainPts = mesh.dimension(0);
  BoundingBox::Data center[DIM];

  for (unsigned i=0; i!=NumDomainPts; ++i) {
    for (unsigned j=0; j<DIM; ++j) center[j] = mesh(i,j);
    const IdentProc::Key Id(0,i);
    const BoundingBox::Key key(Id, stk::parallel_machine_rank(comm));
    BoundingBox B(center, radius, key);
    vector.push_back(B);
  }
}

template <unsigned DIM>
void boundingbox_vector(std::vector<BoundingBox> &vector,
                        const STKMesh            &mesh,
                        const BoundingBox::Data   radius,
                        const stk::ParallelMachine comm) {
  BoundingBox::Data center[DIM];

  for (typename STKMesh::iterator i=mesh.begin(); i!=mesh.end(); ++i) {
    const double *c = mesh.Coord(i);
    for (unsigned j=0; j<DIM; ++j) center[j] = c[j];
    const stk::mesh::EntityKey Id=mesh.Key(i);
    const BoundingBox::Key key(Id, stk::parallel_machine_rank(comm));
    BoundingBox B(center, radius, key);
    vector.push_back(B);
  }
}

}
template <unsigned DIM, class PointData>
void point_to_point_coarse_search(IdentProcRelation &RangeToDomain,
                                  const PointMap    &ToPoints,
                                  const PointData   &FromPoints,
                                  const BoundingBox::Data radius,
                                  const stk::ParallelMachine comm) {

  std::vector<BoundingBox> range_vector;
  std::vector<BoundingBox> domain_vector;

  boundingbox_vector<DIM>(domain_vector, FromPoints, radius, comm);
  boundingbox_vector<DIM>(range_vector,    ToPoints, radius, comm);

  stk::search::FactoryOrder order;
  order.m_communicator = comm;

  // Slightly confusing: coarse_search documentation has domain->range
  // relations sorted by domain key.  We want range->domain type relations
  // sorted on range key. It might appear we have the arguments revered
  // in coarse_search call, but really, this is what we want.
  stk::search::coarse_search(RangeToDomain, domain_vector, range_vector, order);
}

namespace {
template <class ArrayOrMesh>
const double *entity_coord(const ArrayOrMesh &FromPoints,
                           const stk::mesh::EntityKey &k);

template <class ArrayOrMesh>
const double *entity_value(const ArrayOrMesh &FromPoints,
                           const stk::mesh::EntityKey &k);

template<>
  const double *entity_coord<MDArray>(const MDArray &FromPoints,
                                      const stk::mesh::EntityKey &k) {
  const int j = k.id();
  const int i = FromPoints.getEnumeration(j,0);
  const double * d = &FromPoints[i];
  return d;
}

template<>
const double *entity_coord<STKMesh>(const STKMesh &FromMesh,
                                    const stk::mesh::EntityKey &k) {
  const double * d = FromMesh.Coord(k);
  return d;
}

template<>
  const double *entity_value<MDArray>(const MDArray &FromPoints,
                                      const stk::mesh::EntityKey &k) {
  const int j = k.id();
  const int i = FromPoints.getEnumeration(j,0);
  const double * d = &FromPoints[i];
  return d;
}
template<>
  const double *entity_value<STKMesh>(const STKMesh &FromMesh,
                                      const stk::mesh::EntityKey &k) {
  const double * d = FromMesh.Value(k);
  return d;
}
}


template <unsigned DIM, class MeshClass>
void linear_interpolation (MDArray          &ToValues,
                    const MeshClass         &FromValues,
                    const IdentProcRelation &RangeToDomain,
                    const MDArray           &ToPoints,
                    const MeshClass         &FromPoints,
                    const stk::ParallelMachine  comm) {
  const unsigned numValues = ToValues.dimension(0);
  const unsigned numFields = ToValues.dimension(1);
  for (unsigned i=0; i<numValues; ++i)  {
    MDArray Corners(DIM+1,DIM);
    for (unsigned j=0; j<DIM+1; ++j) {
      const IdentProc::Key key = RangeToDomain[4*i+j].second.ident;
      const double *c = entity_coord<MeshClass>(FromPoints, key);
      for (unsigned k=0; k<DIM; ++k) Corners(j,k) = c[k];
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
    for (unsigned f=0; f<numFields; ++f)  {
      for (unsigned j=0; j<DIM+1; ++j) {
        const IdentProc::Key k = RangeToDomain[4*i+j].second.ident;
        const double *c = entity_value(FromValues, k);
        Values[j] = c[f];
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
      ToValues(i,f) = interpolated_value;
    }
  }
}


// Want to find the N best elements to keep.
template <unsigned DIM, class MeshClass>
void filter_with_fine_search(IdentProcRelation  &RangeToDomain,
                                const PointMap  &ToPoints,
                                const MeshClass &FromPoints,
                                const stk::ParallelMachine  comm) {

  typedef  std::map<double,unsigned,std::greater<double> > DIST_MAP;
  const unsigned NumRelations = RangeToDomain.size();
  const unsigned my_rank = stk::parallel_machine_rank(comm);

  DIST_MAP  smallest_distances;
  IdentProcRelation new_relations;

  MDArray::scalar_type from_center[DIM];
  MDArray::scalar_type   to_center[DIM];
  // One of those duel sorted indexing problems:
  for (unsigned i=0, j=0; j<NumRelations;) {
    const unsigned        range_index = RangeToDomain[j].first .ident.id();
    const unsigned       processor_id = RangeToDomain[j].first .proc;
    const IdentProc::Key domain_index = RangeToDomain[j].second.ident;
    if (processor_id == my_rank) {
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
        const double *from_coords = entity_coord<MeshClass>(FromPoints, domain_index);
        for (unsigned k=0; k<DIM; ++k) from_center[k] = from_coords[k];
        for (unsigned k=0; k<DIM; ++k)   to_center[k] =   ToPoints.find(range_index)->second[k];
        const double dist = distance_squared<DIM>(from_center, to_center);
        const DIST_MAP::value_type val(dist,j);

        // Is using a map too much memory allocation/deallocation? Maybe std::vector is better.
        if (smallest_distances.insert(val).second && DIM+1 < smallest_distances.size())
          smallest_distances.erase(smallest_distances.begin()); // DIST_MAP: Largest at begin()
        ++j;
      }
    } else ++j;
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
    const PointMap::key_type id = RangeToDomain[i].first.ident.id();
    ToPoints.erase(id);
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

STKMesh convert_points_to_mesh(const MDArray &Coords,
                               const MDArray &Values,
                               const stk::ParallelMachine  comm) {
  const unsigned num_nodes   = Coords.dimension(0);
  const unsigned spatial_dim = Coords.dimension(1);
  const unsigned num_values  = Values.dimension(1);
  MetaDataPtr meta_data(new stk::mesh::MetaData(spatial_dim));
  BulkDataPtr bulk_data(new stk::mesh::BulkData(*meta_data.get(),comm));

  VectorFieldType &coordinates_field =
    meta_data->declare_field<VectorFieldType>("coordinates", spatial_dim);
  VectorFieldType &values_field =
    meta_data->declare_field<VectorFieldType>("values", num_values);


  stk::mesh::Part & allParts = meta_data->universal_part();
  stk::mesh::put_field(coordinates_field , stk::mesh::MetaData::NODE_RANK, allParts);
  stk::mesh::put_field(values_field ,      stk::mesh::MetaData::NODE_RANK, allParts);

  bulk_data->modification_begin();

  const unsigned entity_rank_count = meta_data->entity_rank_count();
  std::vector<size_t> requests(entity_rank_count,0);
  requests[0] = num_nodes;

  EntityVec entities;
  bulk_data->generate_new_entities(requests, entities);

  bulk_data->modification_end();

  for (unsigned i=0; i<num_nodes; ++i) {
    stk::mesh::Entity e = entities[i];
    double *coord = bulk_data->field_data(coordinates_field, e);
    double *value = bulk_data->field_data(values_field     , e);
    for (unsigned j=0; j<spatial_dim; ++j) coord[j] = Coords(i,j);
    for (unsigned j=0; j<num_values ; ++j) value[j] = Values(i,j);
  }

  STKMesh mesh(meta_data, bulk_data, entities, coordinates_field, values_field);

  return mesh;
}


void copy_domain_to_range_processors(STKMesh                   &Mesh,
                                     const IdentProcRelation   &RangeToDomain,
                                     const std::string         &transfer_name,
                                     const stk::ParallelMachine comm) {

  const unsigned my_rank = stk::parallel_machine_rank(comm);

  std::vector<stk::mesh::EntityProc> domain_to_ghost ;

  const IdentProcRelation::const_iterator end=RangeToDomain.end();
  for (IdentProcRelation::const_iterator i=RangeToDomain.begin(); i!=end; ++i) {

    const stk::mesh::EntityKey domain_entity_key = i->second.ident;
    const unsigned            domain_owning_rank = i->second.proc;
    const unsigned             range_owning_rank = i->first.proc;
    if (domain_owning_rank == my_rank && range_owning_rank != my_rank) {
      const stk::mesh::Entity entity = Mesh.BulkData()->get_entity(domain_entity_key);
      if (Mesh.BulkData()->parallel_owner_rank(entity) == (int)my_rank) { // should always be true by construction.
        const stk::mesh::EntityProc ep(entity, range_owning_rank);
        domain_to_ghost.push_back(ep);
      }
    }
  }

  Mesh.BulkData()->modification_begin();

  stk::mesh::Ghosting & transfer_domain_ghosting =
    Mesh.BulkData()->create_ghosting("Transfer "+transfer_name+" Ghosting");
  {
    std::vector<stk::mesh::EntityKey> receive;
    transfer_domain_ghosting.receive_list( receive );
    Mesh.BulkData()->change_ghosting( transfer_domain_ghosting ,
                                        domain_to_ghost ,
                                        receive );
  }
  Mesh.BulkData()->modification_end();

  // Copy coordinates to the newly ghosted nodes
  std::vector<const stk::mesh::FieldBase *> fields ;
  fields.push_back(&Mesh.Coord());
  fields.push_back(&Mesh.Value());

  stk::mesh::communicate_field_data( transfer_domain_ghosting , fields);
}



}
