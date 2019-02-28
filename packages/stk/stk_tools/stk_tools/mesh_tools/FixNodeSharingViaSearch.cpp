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

void fix_node_sharing_via_search(stk::mesh::BulkData& bulkData)
{
    const int psize = bulkData.parallel_size();
    const int prank = bulkData.parallel_rank();

    const stk::mesh::FieldBase *coordsField = bulkData.mesh_meta_data().coordinate_field();

    if(psize > 1)
    {
        typedef stk::search::Sphere<double> Sphere;
        typedef stk::search::IdentProc<stk::mesh::EntityId, int> IdentProc;
        typedef stk::search::Point<double> Point;
        typedef std::vector<std::pair<Sphere, IdentProc> > SphereVec;
        typedef std::vector<std::pair<IdentProc, IdentProc> > Results;

        SphereVec sourceBboxVector;

        const double radius = 1.0e-6;
        const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::NODE_RANK);
        for(const stk::mesh::Bucket *bucket : buckets)
        {
            if(bucket->owned())
            {
                for(size_t k = 0; k < bucket->size(); ++k)
                {
                    stk::mesh::Entity node = (*bucket)[k];
                    stk::mesh::EntityId id = bulkData.identifier(node);
                    const double* coords = static_cast<double*>(stk::mesh::field_data(*coordsField, node));
                    IdentProc searchId = IdentProc(id, prank);
                    Point point(coords[0], coords[1], bulkData.mesh_meta_data().spatial_dimension() == 3 ? coords[2] : 0.0);
                    Sphere sphere(point, radius);
                    sourceBboxVector.push_back(std::make_pair(sphere, searchId));
                }
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

