// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <stk_util/stk_config.h>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
namespace stk { namespace mesh {

class matrix3x3
{
    public:
        matrix3x3() : m_data(9,0) { }

        void setData(const int index, const double data)
        {
            STK_ThrowErrorIf(index<0 or index>8);
            m_data[index] = data;
        }

        double getData(const int index) const
        {
            STK_ThrowErrorIf(index<0 or index>8);
            return m_data[index];
        }

        double getData(const int row, const int col) const
        {
            STK_ThrowErrorIf(row<0 or row >2);
            STK_ThrowErrorIf(col<0 or col>2);
            int index = row*3 + col;
            return m_data[index];
        }

        int numEntries() const { return 9; }

        void transformVec(const double *in, double *result) const
        {
            result[0] = m_data[0] * in[0] + m_data[1] * in[1] + m_data[2] * in[2];
            result[1] = m_data[3] * in[0] + m_data[4] * in[1] + m_data[5] * in[2];
            result[2] = m_data[6] * in[0] + m_data[7] * in[1] + m_data[8] * in[2];
        }

        void transformVec(const std::vector<double>& in, std::vector<double> &result) const
        {
            result.resize(3);
            this->transformVec(in.data(), result.data());
        }

        ~matrix3x3() {}

        // matrix3x3& operator=(const matrix3x3& rhs); using default assignment operator
        // matrix3x3( const matrix3x3& ); Using default copy constructor

    private:
        std::vector<double> m_data;
};

inline void fillRotationMatrix(const double angleInRadians,  double axisX,  double axisY,  double axisZ, matrix3x3 &rotationMatrix)
{
    double magnitude = sqrt(axisX*axisX + axisY*axisY + axisZ*axisZ);
    STK_ThrowErrorIf(magnitude == 0);

    axisX /= magnitude;
    axisY /= magnitude;
    axisZ /= magnitude;

    double cosAngle = cos(angleInRadians);
    double oneMinusCosAngle = 1 - cosAngle;
    double sinAngle = sin(angleInRadians);

    rotationMatrix.setData(0,cosAngle + axisX*axisX*oneMinusCosAngle);
    rotationMatrix.setData(1,axisX*axisY*oneMinusCosAngle - axisZ*sinAngle);
    rotationMatrix.setData(2,axisX*axisZ*oneMinusCosAngle + axisY*sinAngle);

    rotationMatrix.setData(3,axisX*axisY*oneMinusCosAngle + axisZ*sinAngle);
    rotationMatrix.setData(4,cosAngle + axisY*axisY*oneMinusCosAngle);
    rotationMatrix.setData(5,axisY*axisZ*oneMinusCosAngle - axisX*sinAngle);

    rotationMatrix.setData(6,axisX*axisZ*oneMinusCosAngle - axisY*sinAngle);
    rotationMatrix.setData(7,axisY*axisZ*oneMinusCosAngle + axisX*sinAngle);
    rotationMatrix.setData(8,cosAngle + axisZ*axisZ*oneMinusCosAngle);
}

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

  struct SelectorMapping {
    SelectorMapping(const stk::mesh::Selector & domain, const stk::mesh::Selector & range) : m_domain(domain), m_range(range) {}
    stk::mesh::Selector m_domain;
    stk::mesh::Selector m_range;
    stk::mesh::Selector m_unique_range;
  };

  typedef double Scalar;
  typedef std::vector<SelectorMapping> SelectorMappingVector;
  typedef stk::search::IdentProc<stk::mesh::EntityKey> SearchId;
  typedef stk::search::Point<Scalar> Point;
  typedef stk::search::Sphere<Scalar> Sphere;
  typedef std::vector< std::pair<Sphere,SearchId> > SphereIdVector;
  typedef std::vector<std::pair<SearchId,SearchId> > SearchPairVector;
  typedef std::vector<std::pair<stk::mesh::EntityKey,stk::mesh::EntityKey> > SearchPairSet;

  struct TransformHelper {

    CoordinatesTransform m_transform_type;
    matrix3x3 m_rotation;
    std::vector<double> m_translation;

    TransformHelper()
      : m_transform_type(TRANSLATION)
      , m_translation(3,0)
    {
      // Default is identity transform.
    }

    TransformHelper(double angle, const double axis[3])
      : m_transform_type(ROTATIONAL)
      , m_translation(3,0)
    {
        fillRotationMatrix(angle, axis[0], axis[1], axis[2], m_rotation);
    }

    TransformHelper(double angle, const double axis[3], const double point[3])
      : m_transform_type(PROPER_RIGID)
      , m_translation(point, point+3)
    {
      fillRotationMatrix(angle, axis[0], axis[1], axis[2], m_rotation);

      std::vector<double> result;
      m_rotation.transformVec(m_translation, result);
      m_translation[0] = m_translation[0] - result[0];
      m_translation[1] = m_translation[1] - result[1];
      m_translation[2] = m_translation[2] - result[2];
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
          buff[dim * col + row] = m_rotation.getData(col, row);
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
      m_periodic_mappings(),
      m_search_results(),
      m_periodic_ghosts(NULL),
      m_firstCallToFindPeriodicNodes(true),
      m_hasRotationalPeriodicity(false)
  {}

  const SearchPairVector & get_pairs() const { return m_search_results; }

  const std::vector<TransformHelper>& get_transforms() const { return m_transforms; }

  stk::mesh::Selector get_domain_selector() const
  {
    stk::mesh::Selector result;
    for (typename SelectorMappingVector::const_iterator it = m_periodic_mappings.begin(); it != m_periodic_mappings.end(); ++it)
    {
      result |= it->m_domain;
    }
    return result;
  }

  stk::mesh::Selector get_range_selector() const
  {
    stk::mesh::Selector result;
    for (typename SelectorMappingVector::const_iterator it = m_periodic_mappings.begin(); it != m_periodic_mappings.end(); ++it)
    {
      result |= it->m_range;
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
    }

    for (size_t i = 0; i < m_periodic_mappings.size(); ++i)
    {
      populate_transform(parallel, m_periodic_mappings[i].m_domain, m_periodic_mappings[i].m_range, m_transforms[i]);
    }

    for (size_t i = 0; i < m_periodic_mappings.size(); ++i)
    {
      find_periodic_nodes_for_given_pair(m_periodic_mappings[i].m_domain, m_periodic_mappings[i].m_unique_range, parallel,
                                         m_transforms[i], m_search_tolerances[i]);
    }

    for (size_t i = 0; i < m_search_results.size(); ++i)
    {
      m_unique_search_results.emplace_back(m_search_results[i].first.id(), m_search_results[i].second.id());
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
                                const stk::mesh::Selector & range,
                                const double search_tol = 1.e-10)
  {
    STK_ThrowErrorIf(m_hasRotationalPeriodicity);
    m_periodic_mappings.emplace_back(domain, range);
    m_transforms.push_back( TransformHelper() );
    m_search_tolerances.push_back(search_tol);
  }

  void add_rotational_periodic_pair(const stk::mesh::Selector & domain,
      const stk::mesh::Selector & range,
      const double theta,
      const double axis[],
      const double point[],
      const double search_tol = 1.e-10)
  {
    m_periodic_mappings.emplace_back(domain, range);
    m_search_tolerances.push_back(search_tol);

    //only one periodic BC can exist with rotational periodicity
    STK_ThrowRequire(m_periodic_mappings.size() == 1);

    m_hasRotationalPeriodicity = true;

    if (point[0]*point[0] + point[1]*point[1] + point[2]*point[2] == 0) {
      m_transforms.push_back( TransformHelper(theta, axis) );
    }
    else
      m_transforms.push_back(TransformHelper(theta, axis, point));
  }

  const stk::mesh::Ghosting & get_ghosting() const
  {
    STK_ThrowRequire(m_periodic_ghosts);
    return *m_periodic_ghosts;
  }

  /* Ghost the search results as appropriate. The search results include results on all procs including
   * locally_owned OR globally_shared. We need to ghost the locally owned node (domain or range)
   * to the processor pointed to by the other node in the periodic pair
   */
  void create_ghosting(const std::string & name)
  {
    STK_ThrowRequire(m_bulk_data.in_modifiable_state());
    const int parallel_rank = m_bulk_data.parallel_rank();
    std::vector<stk::mesh::EntityProc> send_nodes;
    for (size_t i=0; i<m_search_results.size(); ++i) {
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

          STK_ThrowRequire(m_bulk_data.parallel_owner_rank(domain_node) == domain_proc);

          send_nodes.emplace_back(domain_node, range_proc);
//          std::cout << "On proc " << m_bulk_data.parallel_rank() << " we are sending domain node to range node "
//              << m_bulk_data.identifier(domain_node) << ":" << m_bulk_data.identifier(range_node)
//              << " since we own the domain and the range resides on proc " << range_proc << std::endl;
        }
        else if (isOwnedRange && range_proc == parallel_rank)
        {
          if (domain_proc == parallel_rank) continue;

          STK_ThrowRequire(m_bulk_data.parallel_owner_rank(range_node) == range_proc);

          send_nodes.emplace_back(range_node, domain_proc);
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
  SelectorMappingVector m_periodic_mappings;
  SearchPairVector m_search_results;
  SearchPairSet m_unique_search_results;

  stk::mesh::Ghosting * m_periodic_ghosts;
  typedef std::vector<TransformHelper > TransformVector;
  TransformVector m_transforms;
  std::vector<double> m_search_tolerances;
  bool m_firstCallToFindPeriodicNodes;
  bool m_hasRotationalPeriodicity;

  void resolve_multi_periodicity()
  {
    switch (m_periodic_mappings.size())
    {
      case 0:
        break;
      case 1:
        m_periodic_mappings[0].m_unique_range = m_periodic_mappings[0].m_range;
        break;
      case 2:
      {
        if ((m_transforms[0].m_transform_type != TRANSLATION)
            || (m_transforms[1].m_transform_type != TRANSLATION)) {
          STK_ThrowErrorMsg("Rotation not supported when there are 2 periodic conditions.");
        }
        const stk::mesh::Selector & domainA = m_periodic_mappings[0].m_domain;
        const stk::mesh::Selector & domainB = m_periodic_mappings[1].m_domain;
        const stk::mesh::Selector domainIntersection = domainA & domainB;

        const stk::mesh::Selector & rangeA = m_periodic_mappings[0].m_range;
        const stk::mesh::Selector & rangeB = m_periodic_mappings[1].m_range;
        const stk::mesh::Selector rangeIntersection = rangeA & rangeB;

        //now add new pair with this
        m_periodic_mappings.emplace_back(domainIntersection, rangeIntersection);
        m_search_tolerances.push_back(1.e-10);
        m_transforms.push_back(TransformHelper());

        //remove redundant entries from unique range
        m_periodic_mappings[0].m_unique_range = m_periodic_mappings[0].m_range - rangeIntersection;
        m_periodic_mappings[1].m_unique_range = m_periodic_mappings[1].m_range - rangeIntersection;
        m_periodic_mappings[2].m_unique_range = m_periodic_mappings[2].m_range;
        break;
      }
      case 3:
      {
        if ((m_transforms[0].m_transform_type != TRANSLATION)
            || (m_transforms[1].m_transform_type != TRANSLATION)
            || (m_transforms[2].m_transform_type != TRANSLATION)
            ) {
          STK_ThrowErrorMsg("Rotation not supported when there are 3 periodic conditions.");
        }
        const stk::mesh::Selector & domainA = m_periodic_mappings[0].m_domain;
        const stk::mesh::Selector & domainB = m_periodic_mappings[1].m_domain;
        const stk::mesh::Selector & domainC = m_periodic_mappings[2].m_domain;

        const stk::mesh::Selector & rangeA = m_periodic_mappings[0].m_range;
        const stk::mesh::Selector & rangeB = m_periodic_mappings[1].m_range;
        const stk::mesh::Selector & rangeC = m_periodic_mappings[2].m_range;

        //point selector
        const stk::mesh::Selector domainABC = domainA & domainB & domainC;
        const stk::mesh::Selector rangeABC = rangeA & rangeB & rangeC;
        //edge selectors
        const stk::mesh::Selector domainAB = domainA & domainB;
        const stk::mesh::Selector domainAC = domainA & domainC;
        const stk::mesh::Selector domainBC = domainB & domainC;
        const stk::mesh::Selector rangeAB = rangeA & rangeB;
        const stk::mesh::Selector rangeAC = rangeA & rangeC;
        const stk::mesh::Selector rangeBC = rangeB & rangeC;

        //edges and point
        m_periodic_mappings.push_back(SelectorMapping(domainAB, rangeAB));
        m_search_tolerances.push_back(1.e-10);
        m_transforms.push_back(TransformHelper());
        m_periodic_mappings.push_back(SelectorMapping(domainBC, rangeBC));
        m_search_tolerances.push_back(1.e-10);
        m_transforms.push_back(TransformHelper());
        m_periodic_mappings.push_back(SelectorMapping(domainAC, rangeAC));
        m_search_tolerances.push_back(1.e-10);
        m_transforms.push_back(TransformHelper());
        m_periodic_mappings.push_back(SelectorMapping(domainABC, rangeABC));
        m_search_tolerances.push_back(1.e-10);
        m_transforms.push_back(TransformHelper());

        //remove redundant entries from unique range
        m_periodic_mappings[0].m_unique_range = m_periodic_mappings[0].m_range - (rangeAB | rangeAC);
        m_periodic_mappings[1].m_unique_range = m_periodic_mappings[1].m_range - (rangeAB | rangeBC);
        m_periodic_mappings[2].m_unique_range = m_periodic_mappings[2].m_range - (rangeAC | rangeBC);
        m_periodic_mappings[3].m_unique_range = m_periodic_mappings[3].m_range - rangeABC;
        m_periodic_mappings[4].m_unique_range = m_periodic_mappings[4].m_range - rangeABC;
        m_periodic_mappings[5].m_unique_range = m_periodic_mappings[5].m_range - rangeABC;
        m_periodic_mappings[6].m_unique_range = m_periodic_mappings[6].m_range;
        break;
      }
      default:
        STK_ThrowErrorMsg("Cannot handle more than 3 periodic pairs");
        break;
    }
  }

  void find_periodic_nodes_for_given_pair(stk::mesh::Selector side1,
                                          stk::mesh::Selector side2,
                                          stk::ParallelMachine parallel,
                                          TransformHelper & transform,
                                          const double search_tolerance
                                          )
  {
    SearchPairVector search_results;
    SphereIdVector side_1_vector, side_2_vector;

    populate_search_vector(side1, side_1_vector, search_tolerance);

    populate_search_vector(side2, side_2_vector, search_tolerance);

    switch (transform.m_transform_type)
    {
      case TRANSLATION:
        translate_coordinates(side_1_vector,
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
        STK_ThrowErrorMsg("Periodic transform method doesn't exist");
        break;
    }

    stk::search::coarse_search(side_1_vector, side_2_vector, stk::search::KDTREE, parallel, search_results);

    m_search_results.insert(m_search_results.end(), search_results.begin(), search_results.end());
  }

  void populate_search_vector(stk::mesh::Selector side_selector,
                              SphereIdVector & aabb_vector,
                              const double search_tolerance
                             )
  {
    const unsigned parallel_rank = m_bulk_data.parallel_rank();

    stk::mesh::BucketVector const& buckets = m_bulk_data.get_buckets( stk::topology::NODE_RANK, side_selector );

    for (size_t bindex = 0, num_buckets = buckets.size(); bindex < num_buckets; ++bindex) {
      stk::mesh::Bucket & b = *buckets[bindex];
      Point center;
      for (size_t ord =0, num_entities = b.size(); ord < num_entities; ++ord) {
        m_get_coordinates(b[ord], &center[0]);
        SearchId search_id( m_bulk_data.entity_key(b[ord]), parallel_rank);
        aabb_vector.emplace_back( Sphere(center, search_tolerance), search_id);
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
        STK_ThrowErrorMsg("Unknown transform type.");
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
      const std::vector<double> &translate) const
  {
    for (auto && side_1 : side_1_vector)
    {
      for (int j = 0; j < 3; ++j)
      {
        side_1.first.center()[j] += translate[j];
      }
    }
  }

  void rotate_coordinates(
      SphereIdVector & side_1_vector,
      SphereIdVector & side_2_vector,
      const matrix3x3 & rotation) const
  {
    for (auto && side_1 : side_1_vector)
    {
      double *center = &(side_1.first.center()[0]);
      std::vector<double> ctr(center, center+3);
      rotation.transformVec(ctr.data(), center);
    }
  }

  void apply_affine_to_coordinates(
      SphereIdVector & side_1_vector,
      SphereIdVector & side_2_vector,
      const matrix3x3 & rotation,
      const std::vector<double> & translation) const
  {
    for (auto && side_1 : side_1_vector)
    {
      double *center = &side_1.first.center()[0];
      std::vector<double> ctr(center, center+3);
      rotation.transformVec(ctr.data(), center);
      for (int i = 0; i < 3; ++i) {
        center[i] += translation[i];
      }
    }
  }

};

template <class CoordFieldType, typename Scalar = double>
struct GetCoordinates
{
  typedef void result_type;
  GetCoordinates(stk::mesh::BulkData & bulk_data, const CoordFieldType & coords_field)
    : m_bulk_data(bulk_data),
      m_coords_field(coords_field)
  {}

  void operator()(stk::mesh::Entity e, Scalar * coords) const
  {
    const unsigned nDim = m_bulk_data.mesh_meta_data().spatial_dimension();
    const double * const temp_coords = reinterpret_cast<Scalar *>(stk::mesh::field_data(m_coords_field, e));
    for (unsigned i = 0; i < nDim; ++i) {
      coords[i] = temp_coords[i];
    }
  }

  stk::mesh::BulkData & m_bulk_data;
  const CoordFieldType & m_coords_field;
};

}} //namespace stk::mesh

namespace impl_hack {

void really_dumb_func();

} //namespace impl_hack

#endif /*STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP*/




