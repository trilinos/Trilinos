// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#ifndef STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP
#define STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_precision.hpp>
#include <glm/gtx/transform.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
namespace stk { namespace mesh {

struct GetCoordiantes;

template <typename CoordinateFunctor>
class PeriodicBoundarySearch
{
public:
  enum CoordinatesTransform {
    TRANSLATION,
    ROTATIONAL,
    PROPER_RIGID
  };

  typedef double Scalar;
  typedef std::vector<std::pair<stk::mesh::Selector, stk::mesh::Selector> > SelectorPairVector;
  typedef stk::search::IdentProc<stk::mesh::EntityKey> SearchId;
  typedef stk::search::Point<Scalar> Point;
  typedef stk::search::Sphere<Scalar> Sphere;
  typedef std::vector< std::pair<Sphere,SearchId> > SphereIdVector;
  typedef std::vector<std::pair<SearchId,SearchId> > SearchPairVector;
  typedef std::vector<std::pair<stk::mesh::EntityKey,stk::mesh::EntityKey> > SearchPairSet;

  struct TransformHelper {

    CoordinatesTransform m_transform_type;
    glm::f64mat3x3 m_rotation;
    glm::f64vec3 m_translation;

    TransformHelper()
      : m_transform_type(TRANSLATION)
      , m_translation(0)
    {
      // Default is identity transform.
    }

    TransformHelper(const boost::array<double, 3> & trans_arg)
      : m_transform_type(TRANSLATION)
      , m_translation(trans_arg[0], trans_arg[1], trans_arg[2] )
    { }

    TransformHelper(double angle, const double axis[3])
      : m_transform_type(ROTATIONAL)
      , m_translation(0)
    {
      m_rotation = glm::f64mat3x3(glm::rotate(angle, axis[0], axis[1], axis[2]));
    }

    TransformHelper(double angle, const double axis[3], const double point[3])
      : m_transform_type(PROPER_RIGID)
      , m_translation(point[0], point[1], point[2])
    {
      m_rotation = glm::f64mat3x3(glm::rotate(angle, axis[0], axis[1], axis[2]));
      m_translation = m_translation - m_rotation*m_translation;  // Can't safely use -=.
    }

    template<typename RealType>
    bool getRowMajorRotation(std::vector<RealType> &buff) const
    {
      if (m_transform_type == TRANSLATION) {
        return false;
      }
      size_t dim = ((buff.size() == 4) ? 2 : 3);
      for (size_t col = 0; col < dim; ++col) {
        for (size_t row = 0; row < dim; ++row) {
          buff[dim * col + row] = m_rotation[row][col];
        }
      }
      return true;
    }
  };

public:

  PeriodicBoundarySearch( stk::mesh::BulkData & bulk_data,
                          const CoordinateFunctor & getCoordinates )
    : m_bulk_data(bulk_data),
      m_get_coordinates(getCoordinates),
      m_periodic_pairs(),
      m_search_results(),
      m_periodic_ghosts(NULL),
      m_firstCallToFindPeriodicNodes(true),
      m_hasRotationalPeriodicity(false)
  {}

  const SearchPairVector & get_pairs() const { return m_search_results; }

  stk::mesh::Selector get_domain_selector() const
  {
    stk::mesh::Selector result;
    for (SelectorPairVector::const_iterator it = m_periodic_pairs.begin(); it != m_periodic_pairs.end(); ++it)
    {
      result |= it->first;
    }
    return result;
  }

  stk::mesh::Selector get_range_selector() const
  {
    stk::mesh::Selector result;
    for (SelectorPairVector::const_iterator it = m_periodic_pairs.begin(); it != m_periodic_pairs.end(); ++it)
    {
      result |= it->second;
    }
    return result;
  }

  void find_periodic_nodes(stk::ParallelMachine parallel)
  {
    m_search_results.clear();
    m_unique_search_results.clear();

    //resolve multiple periodicity once
    if (m_firstCallToFindPeriodicNodes)
    {
      m_firstCallToFindPeriodicNodes = false;
      resolve_multi_periodicity();

      for (size_t i = 0; i < m_periodic_pairs.size(); ++i)
      {
        populate_transform(parallel, m_periodic_pairs[i].first, m_periodic_pairs[i].second, m_transforms[i]);
      }

      //here we need to alter the selectors as to remove the redundant entries
      remove_redundant_nodes();
    }

    for (size_t i = 0; i < m_periodic_pairs.size(); ++i)
    {
      find_periodic_nodes_for_given_pair(m_periodic_pairs[i].first, m_periodic_pairs[i].second, parallel,
                                         m_transforms[i]);
    }

    for (size_t i = 0; i < m_search_results.size(); ++i)
    {
      m_unique_search_results.push_back(std::make_pair(m_search_results[i].first.id(), m_search_results[i].second.id()));
    }

    std::sort(m_unique_search_results.begin(), m_unique_search_results.end());
    SearchPairSet::iterator itor = std::unique(m_unique_search_results.begin(), m_unique_search_results.end());
    m_unique_search_results.erase(itor, m_unique_search_results.end());
  }

  size_t size() const { return m_unique_search_results.size();}

  std::pair<stk::mesh::Entity, stk::mesh::Entity> get_node_pair(size_t i) const
  {
    return std::make_pair(m_bulk_data.get_entity(m_unique_search_results[i].first),
                          m_bulk_data.get_entity(m_unique_search_results[i].second));
  }

  template<typename RealType>
  bool get_search_row_major_rotation(std::vector<RealType> & buff) const
  {
    return m_transforms[0].getRowMajorRotation(buff);
  }

  void add_linear_periodic_pair(const stk::mesh::Selector & domain,
      const stk::mesh::Selector & range)
  {
    ThrowErrorIf(m_hasRotationalPeriodicity);
    m_periodic_pairs.push_back(std::make_pair(domain, range));
    m_transforms.push_back( TransformHelper() );
  }

  void add_rotational_periodic_pair(const stk::mesh::Selector & domain,
      const stk::mesh::Selector & range,
      const double theta,
      const double axis[],
      const double point[])
  {
    m_periodic_pairs.push_back(std::make_pair(domain, range));
    //only one periodic BC can exist with rotational periodicity
    ThrowRequire(m_periodic_pairs.size() == 1);

    m_hasRotationalPeriodicity = true;

    if (point[0]*point[0] + point[1]*point[1] + point[2]*point[2] == 0) {
      m_transforms.push_back( TransformHelper(theta, axis) );
    }
    else
      m_transforms.push_back(TransformHelper(theta, axis, point));
  }

  const stk::mesh::Ghosting & get_ghosting() const
  {
    ThrowRequire(m_periodic_ghosts);
    return *m_periodic_ghosts;
  }

  /* Ghost the search results as appropriate. The search results include results on all procs including
   * locally_owned OR globally_shared. We need to ghost the locally owned node (domain or range)
   * to the processor pointed to by the other node in the periodic pair
   */
  void create_ghosting(const std::string & name)
  {
    ThrowRequire(m_bulk_data.synchronized_state() == m_bulk_data.MODIFIABLE);
    const int parallel_rank = m_bulk_data.parallel_rank();
    std::vector<stk::mesh::EntityProc> send_nodes;
    for (size_t i=0, size=m_search_results.size(); i<size; ++i) {
        stk::mesh::Entity domain_node = m_bulk_data.get_entity(m_search_results[i].first.id());
        stk::mesh::Entity range_node = m_bulk_data.get_entity(m_search_results[i].second.id());

        bool isOwnedDomain = m_bulk_data.is_valid(domain_node) ? m_bulk_data.bucket(domain_node).owned() : false;
        bool isOwnedRange = m_bulk_data.is_valid(range_node) ? m_bulk_data.bucket(range_node).owned() : false;
        int domain_proc = m_search_results[i].first.proc();
        int range_proc = m_search_results[i].second.proc();

        //Shouldn't happen since one needs to be shared or owned (it was used in the search)

        if (isOwnedDomain && domain_proc == parallel_rank)
        {
          if (range_proc == parallel_rank) continue;

          ThrowRequire(m_bulk_data.parallel_owner_rank(domain_node) == domain_proc);

          send_nodes.push_back(stk::mesh::EntityProc(domain_node, range_proc));
//          std::cout << "On proc " << m_bulk_data.parallel_rank() << " we are sending domain node to range node "
//              << m_bulk_data.identifier(domain_node) << ":" << m_bulk_data.identifier(range_node)
//              << " since we own the domain and the range resides on proc " << range_proc << std::endl;
        }
        else if (isOwnedRange && range_proc == parallel_rank)
        {
          if (domain_proc == parallel_rank) continue;

          ThrowRequire(m_bulk_data.parallel_owner_rank(range_node) == range_proc);

          send_nodes.push_back(stk::mesh::EntityProc(range_node, domain_proc));
//          std::cout << "On proc " << m_bulk_data.parallel_rank() << " we are sending range node to domain node "
//              << m_bulk_data.identifier(domain_node) << ":" << m_bulk_data.identifier(range_node)
//              << " since we own the range and the domain resides on proc " << domain_proc << std::endl;
        }
//        else
//        {
//          std::cout << "On proc " << m_bulk_data.parallel_rank() << " we have both nodes unowned between procs: " << domain_proc << "," << range_proc << std::endl;
//        }
    }

    m_periodic_ghosts = &m_bulk_data.create_ghosting(name);
    m_bulk_data.change_ghosting(*m_periodic_ghosts, send_nodes);
  }

private:
  stk::mesh::BulkData & m_bulk_data;
  CoordinateFunctor m_get_coordinates;
  SelectorPairVector m_periodic_pairs;
  SearchPairVector m_search_results;
  SearchPairSet m_unique_search_results;

  stk::mesh::Ghosting * m_periodic_ghosts;
  typedef std::vector<TransformHelper > TransformVector;
  TransformVector m_transforms;
  bool m_firstCallToFindPeriodicNodes;
  bool m_hasRotationalPeriodicity;

  void resolve_multi_periodicity()
  {
    switch (m_periodic_pairs.size())
    {
      case 0:
      case 1:
        //nothing to do here
        break;
      case 2:
      {
        if ((m_transforms[0].m_transform_type != TRANSLATION)
            || (m_transforms[1].m_transform_type != TRANSLATION)) {
          ThrowErrorMsg("Rotation not supported when there are 2 periodic conditions.");
        }
        const stk::mesh::Selector & domainA = m_periodic_pairs[0].first;
        const stk::mesh::Selector & domainB = m_periodic_pairs[1].first;
        const stk::mesh::Selector domainIntersection = domainA & domainB;

        const stk::mesh::Selector & rangeA = m_periodic_pairs[0].second;
        const stk::mesh::Selector & rangeB = m_periodic_pairs[1].second;
        const stk::mesh::Selector rangeIntersection = rangeA & rangeB;

        //now add new pair with this
        m_periodic_pairs.push_back(std::make_pair(domainIntersection, rangeIntersection));
        m_transforms.push_back(TransformHelper());
        break;
      }
      case 3:
      {
        if ((m_transforms[0].m_transform_type != TRANSLATION)
            || (m_transforms[1].m_transform_type != TRANSLATION)
            || (m_transforms[2].m_transform_type != TRANSLATION)
            ) {
          ThrowErrorMsg("Rotation not supported when there are 3 periodic conditions.");
        }
        const stk::mesh::Selector domainA = m_periodic_pairs[0].first;
        const stk::mesh::Selector domainB = m_periodic_pairs[1].first;
        const stk::mesh::Selector domainC = m_periodic_pairs[2].first;

        const stk::mesh::Selector rangeA = m_periodic_pairs[0].second;
        const stk::mesh::Selector rangeB = m_periodic_pairs[1].second;
        const stk::mesh::Selector rangeC = m_periodic_pairs[2].second;

        //edges
        m_periodic_pairs.push_back(std::make_pair(domainA & domainB, rangeA & rangeB));
        m_transforms.push_back(TransformHelper());
        m_periodic_pairs.push_back(std::make_pair(domainB & domainC, rangeB & rangeC));
        m_transforms.push_back(TransformHelper());
        m_periodic_pairs.push_back(std::make_pair(domainA & domainC, rangeA & rangeC));
        m_transforms.push_back(TransformHelper());
        m_periodic_pairs.push_back(std::make_pair(domainA & domainB & domainC, rangeA & rangeB & rangeC));
        m_transforms.push_back(TransformHelper());
        break;
      }
      default:
        ThrowErrorMsg("Cannot handle more than 3 periodic pairs");
        break;
    }
  }

  void find_periodic_nodes_for_given_pair(stk::mesh::Selector side1,
                                          stk::mesh::Selector side2,
                                          stk::ParallelMachine parallel,
                                          TransformHelper & transform
                                          )
  {
    SearchPairVector search_results;
    SphereIdVector side_1_vector, side_2_vector;

    populate_search_vector(side1, side_1_vector);

    populate_search_vector(side2, side_2_vector);

    switch (transform.m_transform_type)
    {
      case TRANSLATION:
        translate_coordinates(side_1_vector,
                              side_2_vector,
                              transform.m_translation);
        break;
      case ROTATIONAL:
        //something
        rotate_coordinates(side_1_vector, side_2_vector, transform.m_rotation);
        break;
      case PROPER_RIGID:
        apply_affine_to_coordinates(side_1_vector, side_2_vector,
                                    transform.m_rotation, transform.m_translation);
        break;
      default:
        ThrowErrorMsg("Periodic transform method doesn't exist");
        break;
    }

    stk::search::coarse_search(side_1_vector, side_2_vector, stk::search::BOOST_RTREE, parallel, search_results);

    m_search_results.insert(m_search_results.end(), search_results.begin(), search_results.end());
  }

  void remove_redundant_nodes()
  {
    switch (m_periodic_pairs.size())
    {
      case 1:
        //do nothing
        break;
      case 3: //2 periodic pairs plus intersections
      {
        const stk::mesh::Selector domainA = m_periodic_pairs[0].first;
        const stk::mesh::Selector domainB = m_periodic_pairs[1].first;

        const stk::mesh::Selector rangeA = m_periodic_pairs[0].second;
        const stk::mesh::Selector rangeB = m_periodic_pairs[1].second;
        const stk::mesh::Selector rangeIntersection = rangeA & rangeB;

        m_periodic_pairs[0].second = rangeA - rangeIntersection;
        m_periodic_pairs[1].second = rangeB - rangeIntersection;

        break;
      }
      case 7:
      {
        const stk::mesh::Selector domainA = m_periodic_pairs[0].first;
        const stk::mesh::Selector domainB = m_periodic_pairs[1].first;
        const stk::mesh::Selector domainC = m_periodic_pairs[2].first;

        const stk::mesh::Selector rangeA = m_periodic_pairs[0].second;
        const stk::mesh::Selector rangeB = m_periodic_pairs[1].second;
        const stk::mesh::Selector rangeC = m_periodic_pairs[2].second;

        //point selector
        const stk::mesh::Selector rangeABC = rangeA & rangeB & rangeC;
        //edge selectors
        const stk::mesh::Selector rangeAB = rangeA & rangeB;
        const stk::mesh::Selector rangeAC = rangeA & rangeC;
        const stk::mesh::Selector rangeBC = rangeB & rangeC;

        //now we redefine the periodic pairs to remove redundancies
        m_periodic_pairs[0].second = rangeA - (rangeAB | rangeAC);
        m_periodic_pairs[1].second = rangeB - (rangeAB | rangeBC);
        m_periodic_pairs[2].second = rangeC - (rangeAC | rangeBC);

        //edges
        m_periodic_pairs[3].second = rangeAB - rangeABC;
        m_periodic_pairs[4].second = rangeBC - rangeABC;
        m_periodic_pairs[5].second = rangeAC - rangeABC;

        break;
      }
      default:
        ThrowErrorMsg("Cannot handle this number of periodic pairs");
        break;
    }

  }

  void populate_search_vector(stk::mesh::Selector side_selector
                              , SphereIdVector & aabb_vector
                             )
  {
    const double radius = 1e-10;
    const unsigned parallel_rank = m_bulk_data.parallel_rank();

    stk::mesh::BucketVector const& buckets = m_bulk_data.get_buckets( stk::topology::NODE_RANK, side_selector );

    size_t num_nodes = 0;

    for (size_t bindex = 0, num_buckets = buckets.size(); bindex < num_buckets; ++bindex) {
      stk::mesh::Bucket & b = *buckets[bindex];
      Point center;
      for (size_t ord =0, num_entities = b.size(); ord < num_entities; ++ord) {
        ++num_nodes;
        m_get_coordinates(b[ord], &center[0]);
        SearchId search_id( m_bulk_data.entity_key(b[ord]), parallel_rank);
        aabb_vector.push_back( std::make_pair( Sphere(center, radius), search_id));
      }
    }
  }

  void populate_transform(stk::ParallelMachine parallel,
      stk::mesh::Selector side_1,
      stk::mesh::Selector side_2,
      TransformHelper & transform)
  {

    switch(transform.m_transform_type )
    {
      case TRANSLATION:
      {
        //now get the centroids and populate the displacement
        double centroid_1[3];
        double centroid_2[3];

        calc_centroid(parallel, side_1 & m_bulk_data.mesh_meta_data().locally_owned_part(), centroid_1);
        calc_centroid(parallel, side_2 & m_bulk_data.mesh_meta_data().locally_owned_part(), centroid_2);
        for (int i = 0; i < 3; ++i)
          transform.m_translation[i] = centroid_2[i] - centroid_1[i];

        break;
      }
      case ROTATIONAL:
      case PROPER_RIGID:
      {
        //do nothing
        break;
      }
      default:
        ThrowErrorMsg("Unknown transform type.");
        break;
    }

  }

  void calc_centroid(stk::ParallelMachine parallel, stk::mesh::Selector side_selector, double * global_centroid)
  {
    const int spatial_dimension = m_bulk_data.mesh_meta_data().spatial_dimension();

    stk::mesh::BucketVector const& buckets =
      m_bulk_data.get_buckets( stk::topology::NODE_RANK, side_selector & m_bulk_data.mesh_meta_data().locally_owned_part());

    size_t num_nodes = 0;

    double local_centroid[3] = { };
    for (size_t bindex = 0, num_buckets = buckets.size(); bindex < num_buckets; ++bindex) {
      stk::mesh::Bucket & b = *buckets[bindex];
      double coords[3];
      for (size_t ord =0, num_entities = b.size(); ord < num_entities; ++ord) {
        ++num_nodes;
        m_get_coordinates(b[ord], coords);
        for (int i = 0; i < spatial_dimension; ++i) {
          local_centroid[i] += coords[i];
        }
      }
    }
    //do coordinate tranformation here
    size_t global_node_count = 0;

    stk::all_reduce_sum(parallel, &num_nodes, &global_node_count, 1);
    stk::all_reduce_sum(parallel, local_centroid, global_centroid, 3);
    for (int i = 0; i < 3; ++i)
    {
      global_centroid[i] = global_centroid[i] / global_node_count;
    }
  }

  void translate_coordinates(
      SphereIdVector & side_1_vector,
      SphereIdVector & side_2_vector,
      const glm::f64vec3 &translate) const
  {
    // translate domain to range, i.e. master to slave
    for (size_t i = 0, size = side_1_vector.size(); i < size; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        side_1_vector[i].first.center()[j] += translate[j];
      }
    }
  }

  void rotate_coordinates(
      SphereIdVector & side_1_vector,
      SphereIdVector & side_2_vector,
      const glm::f64mat3x3 & rotation) const
  {
    for (size_t iPoint = 0, size = side_1_vector.size(); iPoint < size; ++iPoint)
    {
      double *center = &side_1_vector[iPoint].first.center()[0];
      glm::f64vec3 ctr(center[0], center[1], center[2]);
      ctr = rotation * ctr;
      for (int i = 0; i < 3; ++i) {
        center[i] = ctr[i];
      }
    }
  }

  void apply_affine_to_coordinates(
      SphereIdVector & side_1_vector,
      SphereIdVector & side_2_vector,
      const glm::f64mat3x3 & rotation,
      const glm::f64vec3 & translation) const
  {
    for (size_t iPoint = 0, size = side_1_vector.size(); iPoint < size; ++iPoint)
    {
      double *center = &side_1_vector[iPoint].first.center()[0];
      glm::f64vec3 ctr(center[0], center[1], center[2]);
      ctr = rotation * ctr + translation; // Can't safely use LHS in RHS.
      for (int i = 0; i < 3; ++i) {
        center[i] = ctr[i];
      }
    }
  }

};

template <class CoordFieldType, typename Scalar = double>
struct GetCoordinates
{
  typedef void result_type;
  GetCoordinates(stk::mesh::BulkData & bulk_data, CoordFieldType & coords_field)
    : m_bulk_data(bulk_data),
      m_coords_field(coords_field)
  {}

  void operator()(stk::mesh::Entity e, Scalar * coords) const
  {
    const unsigned nDim = m_bulk_data.mesh_meta_data().spatial_dimension();
    const double * const temp_coords = stk::mesh::field_data(m_coords_field, e);
    for (unsigned i = 0; i < nDim; ++i) {
      coords[i] = temp_coords[i];
    }
  }

  stk::mesh::BulkData & m_bulk_data;
  CoordFieldType & m_coords_field;
};

template <class ModelCoordFieldType, class DispCoordFieldType, typename Scalar = double>
struct GetDisplacedCoordinates
{
  typedef void result_type;
  GetDisplacedCoordinates(stk::mesh::BulkData & bulk_data, ModelCoordFieldType & model_coord_field, DispCoordFieldType & disp_coord_field)
    : m_bulk_data(bulk_data),
      m_model_coord_field(model_coord_field),
      m_disp_coord_field(disp_coord_field)
  {}

  void operator()(stk::mesh::Entity e, Scalar * coords) const
  {
    const unsigned nDim = m_bulk_data.mesh_meta_data().spatial_dimension();
    const double * const temp_model_coords = stk::mesh::field_data(m_model_coord_field, e);
    const double * const temp_disp_coords = stk::mesh::field_data(m_disp_coord_field, e);
    for (unsigned i = 0; i < nDim; ++i) {
      coords[i] = temp_model_coords[i] + temp_disp_coords[i];
    }
  }

  stk::mesh::BulkData & m_bulk_data;
  ModelCoordFieldType & m_model_coord_field;
  DispCoordFieldType & m_disp_coord_field;
};

}} //namespace stk::mesh

namespace impl_hack {

void really_dumb_func();

} //namespace impl_hack

#endif /*STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP*/




