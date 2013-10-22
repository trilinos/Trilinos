#ifndef STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP
#define STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk { namespace mesh {

struct GetCoordiantes;

template <typename CoordinateFunctor>
class PeriodicBoundarySearch
{
public:
  enum CoordinatesTransform {
    TRANSLATION,
    ROTATIONAL
  };

  struct TransformHelper {
    TransformHelper()
      : translation(),
        transformType(TRANSLATION)
    {
      rotation[0] = 1.0;
      rotation[1] = 0.0;
      rotation[2] = 0.0;
      rotation[3] = 0.0;
      rotation[4] = 1.0;
      rotation[5] = 0.0;
      rotation[6] = 0.0;
      rotation[7] = 0.0;
      rotation[8] = 1.0;
    }

    TransformHelper(const boost::array<double, 3> & translation_arg)
      : translation(translation_arg),
        transformType(TRANSLATION)
    {
      rotation[0] = 1.0;
      rotation[1] = 0.0;
      rotation[2] = 0.0;
      rotation[3] = 0.0;
      rotation[4] = 1.0;
      rotation[5] = 0.0;
      rotation[6] = 0.0;
      rotation[7] = 0.0;
      rotation[8] = 1.0;
    }

    TransformHelper(const boost::array<double, 9> & rotation_arg)
      : rotation(rotation_arg),
        translation(),
        transformType(ROTATIONAL)
        {}

    boost::array<double, 9> rotation;
    boost::array<double, 3> translation;
    CoordinatesTransform transformType;

  };

public:
  typedef double Scalar;
  typedef std::vector<std::pair<stk::mesh::Selector, stk::mesh::Selector> > SelectorPairVector;
  typedef stk::search::ident::IdentProc<stk::mesh::EntityKey,unsigned> SearchId;
  typedef stk::search::box::SphereBoundingBox<SearchId,Scalar,3> AABB;
  typedef std::vector<AABB> AABBVector;
  typedef std::vector<std::pair<SearchId,SearchId> > SearchPairVector;


  PeriodicBoundarySearch( stk::mesh::BulkData & bulk_data,
                          const CoordinateFunctor & getCoordinates )
    : m_bulk_data(bulk_data),
      m_get_coordinates(getCoordinates),
      m_periodic_pairs(),
      m_search_results(),
      m_periodic_ghosts(NULL),
      m_initMultiPeriodic(true)
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

    if (m_initMultiPeriodic)
    {
      m_initMultiPeriodic = false;

      resolve_multi_periodicity();
    }

    for (size_t i = 0; i < m_periodic_pairs.size(); ++i)
    {
      stk::mesh::Selector & side1 = m_periodic_pairs[i].first;
      stk::mesh::Selector & side2 = m_periodic_pairs[i].second;
      const CoordinatesTransform transform_method = m_transforms[i].transformType;

      find_periodic_nodes_for_given_pair(side1, side2, transform_method, parallel, m_transforms[i]);
    }
  }

  size_t size() const { return m_search_results.size();}

  std::pair<stk::mesh::Entity, stk::mesh::Entity> get_node_pair(size_t i) const
  {
    return std::make_pair(m_bulk_data.get_entity(m_search_results[i].first.ident),
        m_bulk_data.get_entity(m_search_results[i].second.ident));
  }


  void add_linear_periodic_pair(const stk::mesh::Selector & domain,
      const stk::mesh::Selector & range)
  {
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

    //create rotation matrix from domain to range
    boost::array<double, 9> rotationMatrix;
    rotationMatrix[0] = cos(theta);   rotationMatrix[3] = -sin(theta);    rotationMatrix[6] = 0.0;
    rotationMatrix[1] = sin(theta);   rotationMatrix[4] = cos(theta);     rotationMatrix[7] = 0.0;
    rotationMatrix[2] = 0.0;          rotationMatrix[5] = 0.0;            rotationMatrix[8] = 1.0;
    m_transforms.push_back( TransformHelper(rotationMatrix) );

  }

  const stk::mesh::Ghosting & get_ghosting() const
  {
    ThrowRequire(m_periodic_ghosts);
    return *m_periodic_ghosts;
  }

  void create_ghosting(const std::string & name)
  {
    ThrowRequire(m_bulk_data.synchronized_state() == m_bulk_data.MODIFIABLE);
    const int parallel_rank = m_bulk_data.parallel_rank();
    size_t num_constraints = 0;
    std::vector<stk::mesh::EntityProc> send_nodes;
    for (size_t i=0, size=m_search_results.size(); i<size; ++i) {
        stk::mesh::Entity domain_node = m_bulk_data.get_entity(m_search_results[i].first.ident);
        int domain_proc = m_search_results[i].first.proc;
        stk::mesh::Entity range_node = m_bulk_data.get_entity(m_search_results[i].second.ident);
        int range_proc = m_search_results[i].second.proc;

        if (parallel_rank == domain_proc) ++num_constraints;

        if ((parallel_rank != domain_proc) && (parallel_rank == range_proc)) {
          send_nodes.push_back(stk::mesh::EntityProc(range_node, domain_proc));
        }
        else if ((parallel_rank == domain_proc) && (parallel_rank != range_proc)) {
          send_nodes.push_back(stk::mesh::EntityProc(domain_node, range_proc));
        }
        ThrowAssert((parallel_rank == domain_proc) || (parallel_rank == range_proc));
    }

    m_periodic_ghosts = &m_bulk_data.create_ghosting(name);
    m_bulk_data.change_ghosting(*m_periodic_ghosts, send_nodes);
  }

private:
  stk::mesh::BulkData & m_bulk_data;
  CoordinateFunctor m_get_coordinates;
  SelectorPairVector m_periodic_pairs;
  SearchPairVector m_search_results;
  stk::mesh::Ghosting * m_periodic_ghosts;
  typedef std::vector<TransformHelper > TransformVector;
  TransformVector m_transforms;
  bool m_initMultiPeriodic;


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
        const stk::mesh::Selector & domainA = m_periodic_pairs[0].first;
        const stk::mesh::Selector & domainB = m_periodic_pairs[1].first;
        const stk::mesh::Selector domainIntersection = domainA & domainB;

        const stk::mesh::Selector & rangeA = m_periodic_pairs[0].second;
        const stk::mesh::Selector & rangeB = m_periodic_pairs[1].second;
        const stk::mesh::Selector rangeIntersection = rangeA & rangeB;

        //now add new pair with this
        m_periodic_pairs.push_back(std::make_pair(domainIntersection, rangeIntersection));
        //TODO: need logic to handle various combinations of transforms
        m_transforms.push_back(TransformHelper());
        break;
      }
      case 3:
      {
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
        ThrowRequireMsg(false, "Cannot handle this number of periodic pairs");
        break;
    }
  }

  void find_periodic_nodes_for_given_pair(stk::mesh::Selector side1,
      stk::mesh::Selector side2,
      CoordinatesTransform transform_method,
      stk::ParallelMachine parallel,
      TransformHelper & transform)
  {
    SearchPairVector search_results;
    AABBVector side_1_vector, side_2_vector;
    double local_centroid[6] =
    { 0 };
    size_t local_node_count[2];

    local_node_count[0] = populate_search_vector(side1, side_1_vector, local_centroid);

    local_node_count[1] = populate_search_vector(side2, side_2_vector, local_centroid + 3);

    switch (transform_method)
    {
      case TRANSLATION:
        translate_coordinates(parallel,
                              local_node_count,
                              local_centroid,
                              side_1_vector,
                              side_2_vector,
                              transform.translation);
        break;
      case ROTATIONAL:
        //something
        rotate_coordinates(parallel, local_node_count, local_centroid, side_1_vector, side_2_vector, transform.rotation);
        break;
      default:
        ThrowRequireMsg(false, "Periodic transform method doesn't exist");
        break;
    }
    stk::search::FactoryOrder order;
    order.m_communicator = parallel;
    stk::search::coarse_search(search_results, side_2_vector, side_1_vector, order);

    m_search_results.insert(m_search_results.end(), search_results.begin(), search_results.end());
  }

  size_t populate_search_vector(stk::mesh::Selector side_selector
                              , AABBVector & aabb_vector
                              , double * centroid
                             )
  {
    const int spatial_dimension = m_bulk_data.mesh_meta_data().spatial_dimension();
    const double radius = 1e-10;
    const unsigned parallel_rank = m_bulk_data.parallel_rank();

    stk::mesh::BucketVector buckets;

    stk::mesh::get_buckets( side_selector
                            ,m_bulk_data.buckets(stk::topology::NODE_RANK)
                            ,buckets
                            );

    size_t num_nodes = 0;

    for (size_t bindex = 0, num_buckets = buckets.size(); bindex < num_buckets; ++bindex) {
      stk::mesh::Bucket & b = *buckets[bindex];
      double coords[3];
      for (size_t ord =0, num_entities = b.size(); ord < num_entities; ++ord) {
        ++num_nodes;
        m_get_coordinates(b[ord], coords);
        SearchId search_id( m_bulk_data.entity_key(b[ord]), parallel_rank);
        aabb_vector.push_back(AABB( coords, radius, search_id));
        for (int i=0; i<spatial_dimension; ++i) {
          centroid[i] += coords[i];
        }
      }
    }
    return num_nodes;
  }

  void translate_coordinates(stk::ParallelMachine parallel,
      size_t local_node_count[2],
      double local_centroid[6],
      AABBVector & side_1_vector,
      AABBVector & side_2_vector,
      boost::array<double, 3> & translate) const
  {
    //do coordinate tranformation here
    size_t global_node_count[2] =
    { };
    double global_centroid[6] =
    { };
    stk::all_reduce_sum(parallel, local_node_count, global_node_count, 2);
    stk::all_reduce_sum(parallel, local_centroid, global_centroid, 6);
    for (int i = 0; i < 3; ++i)
    {
      translate[i] = (global_centroid[i + 3] / global_node_count[1]) - (global_centroid[i] / global_node_count[0]);
    }
    // translate domain to range, i.e. master to slave
    for (size_t i = 0, size = side_1_vector.size(); i < size; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        side_1_vector[i].center[j] += translate[j];
      }
    }
  }


  void rotate_coordinates(stk::ParallelMachine parallel,
      size_t local_node_count[2],
      double local_centroid[6],
      AABBVector & side_1_vector,
      AABBVector & side_2_vector,
      boost::array<double, 9> & rotation) const
  {
    //TODO: now assuming that periodicity is around z axis
    for (size_t iPoint = 0, size = side_1_vector.size(); iPoint < size; ++iPoint)
    {
      double r[3];
      for (int i = 0; i < 3; ++i)
      {
        r[i] = side_1_vector[iPoint].center[i];
        side_1_vector[iPoint].center[i] = 0.0;
      }

      for (int j = 0; j < 3; ++j)
      {
        for (int i = 0; i < 3; ++i)
        {
          side_1_vector[iPoint].center[i] += rotation[i + 3*j]*r[j];
        }
      }
    }


  }

};

template <class CoordFieldType, typename Scalar = double, unsigned SpatialDimension = 3>
struct GetCoordinates
{
  typedef void result_type;
  GetCoordinates(stk::mesh::BulkData & bulk_data, CoordFieldType & coords_field)
    : m_bulk_data(bulk_data),
      m_coords_field(coords_field)
  {}

  void operator()(stk::mesh::Entity e, Scalar * coords) const
  {
    const double * const temp_coords = m_bulk_data.field_data(m_coords_field, e);
    for (unsigned i = 0; i < SpatialDimension; ++i) {
      coords[i] = temp_coords[i];
    }
  }

  stk::mesh::BulkData & m_bulk_data;
  CoordFieldType & m_coords_field;
};

template <class ModelCoordFieldType, class DispCoordFieldType, typename Scalar = double, unsigned SpatialDimension = 3>
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
    const double * const temp_model_coords = m_bulk_data.field_data(m_model_coord_field, e);
    const double * const temp_disp_coords = m_bulk_data.field_data(m_disp_coord_field, e);
    for (unsigned i = 0; i < SpatialDimension; ++i) {
      coords[i] = temp_model_coords[i] + temp_disp_coords[i];
    }
  }

  stk::mesh::BulkData & m_bulk_data;
  ModelCoordFieldType & m_model_coord_field;
  DispCoordFieldType & m_disp_coord_field;
};




}}//namespace stk::mesh

#endif /*STK_SEARCH_UTIL_STK_MESH_PERIODIC_BOUNDARY_SEARCH_HPP*/




