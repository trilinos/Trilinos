#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/Sphere.hpp>

namespace stk {
namespace tools {

void fix_node_sharing_via_search(stk::mesh::BulkData& bulkData, double tolerance)
{
  const int psize = bulkData.parallel_size();
  const int prank = bulkData.parallel_rank();

  const stk::mesh::FieldBase *coordsField = bulkData.mesh_meta_data().coordinate_field();
  const unsigned dim = bulkData.mesh_meta_data().spatial_dimension();

  if(psize > 1)
  {
    typedef stk::search::Sphere<double> Sphere;
    typedef stk::search::IdentProc<stk::mesh::EntityId, int> IdentProc;
    typedef stk::search::Point<double> Point;
    typedef std::vector<std::pair<Sphere, IdentProc> > SphereVec;
    typedef std::vector<std::pair<IdentProc, IdentProc> > Results;

    SphereVec sourceBboxVector;

    const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::NODE_RANK);
    for(const stk::mesh::Bucket *bucket : buckets)
    {
      if (bucket->owned())
      {
        stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(*coordsField,
          [&](auto& coordsData) {
            for (size_t k = 0; k < bucket->size(); ++k) {
              stk::mesh::Entity node = (*bucket)[k];
              stk::mesh::EntityId id = bulkData.identifier(node);
              auto coords = coordsData.entity_values(node);
              IdentProc searchId = IdentProc(id, prank);
              Point point(coords(0_comp), dim >= 2 ? coords(1_comp) : 0.0, dim == 3 ? coords(2_comp) : 0.0);
              Sphere sphere(point, tolerance);
              sourceBboxVector.push_back(std::make_pair(sphere, searchId));
            }
          }
        );
      }
    }

    Results searchResults;
    stk::search::SearchMethod searchMethod = stk::search::KDTREE;
    stk::search::coarse_search(sourceBboxVector, sourceBboxVector, searchMethod, bulkData.parallel(), searchResults);

    for(const auto &resultPair : searchResults)
    {
      stk::mesh::EntityId id = resultPair.second.id();
      int proc = resultPair.second.proc();
      if(proc != prank && id == resultPair.first.id())
      {
        stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, id);
        bulkData.add_node_sharing(node, proc);
      }
    }
  }
}

}}

