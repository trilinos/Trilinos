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

void MeshGeometry::snap_points_to_geometry()
{
    size_t i;
    for (i=0; i<geomEvaluators.size(); i++)
    {
        size_t j;
        for (j=0; j<geomEvaluators[i]->mesh.size()/3; j++)
        {
            double *point = &(geomEvaluators[i]->mesh[j*3]);
            double *new_point = geomKernel->snap_to(point, geomEvaluators[i]->geometry);
        }
    }
}
