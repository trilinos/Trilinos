#include "GeometryFactory.hpp"
#include "PerceptMesh.hpp"

GeometryFactory::GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry)
{
    geomKernel = kernel;
    geomDatabase = geometry;
}

GeometryFactory::~GeometryFactory()
{

}

bool GeometryFactory::read_file(const std::string& filename)
{
    return geomKernel->read_file(filename, geometryHandles);
}

void GeometryFactory::create_evaluators(PerceptMesh* mesh_data)
{
    size_t i;
    for (i=0; i<geometryHandles.size(); i++)
    {
        int mesh = geomKernel->get_attribute(std::string("mesh_sideset"),
                                             geometryHandles[i]);
        GeometryEvaluator* eval = new GeometryEvaluator;
        geomDatabase->add_evaluator(eval);
        eval->geometry = geometryHandles[i];
        eval->mesh = get_nodes_in_sideset(mesh_data, mesh);
    }
}

PointSet& GeometryFactory::get_nodes_in_sideset(PerceptMesh* mesh_data, int sideset)
{
    PointSet points;

    return points;
}
