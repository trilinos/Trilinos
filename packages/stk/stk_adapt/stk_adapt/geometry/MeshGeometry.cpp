#include "MeshGeometry.hpp"

MeshGeometry::MeshGeometry(GeometryKernel* geom)
{
    geomKernel = geom;
}

MeshGeometry::~MeshGeometry()
{

}

void MeshGeometry::add_evaluator(GeometryEvaluator* evaluator)
{
    geomEvaluators.push_back(evaluator);
}

void MeshGeometry::add_evaluators(std::vector<GeometryEvaluator*> evaluators)
{
    geomEvaluators.insert(geomEvaluators.end(), evaluators.begin(), evaluators.end());
}

const std::vector<GeometryEvaluator*>& MeshGeometry::getGeomEvaluators()
{
  return geomEvaluators;
}

void MeshGeometry::snap_points_to_geometry(PerceptMesh* mesh_data)
{
    BulkData& bulkData = *mesh_data->getBulkData();
    VectorFieldType* coordField = mesh_data->getCoordinatesField();

    const std::vector<Bucket*> & buckets = bulkData.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

    for ( std::vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        size_t s;
        for (s=0; s<geomEvaluators.size(); s++)
        {
          if (geomEvaluators[s]->mMesh(**k))
            {
               Bucket & bucket = **k ;
               const unsigned num_nodes_in_bucket = bucket.size();

               for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                 {
                   Entity& node = bucket[iNode];

                   double * coord = stk::mesh::field_data( *coordField , node );
                   geomKernel->snap_to(coord, geomEvaluators[s]->mGeometry);
                 }
            }
        }
      }
}
