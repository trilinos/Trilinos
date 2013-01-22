#include "GeometryFactory.hpp"
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

using namespace stk;
using namespace mesh;
using namespace percept;

GeometryFactory::GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry)
{
    geomKernel = kernel;
    geomDatabase = geometry;
}

GeometryFactory::~GeometryFactory()
{
#if 1
  const std::vector<GeometryEvaluator*>& evaluators = geomDatabase->getGeomEvaluators();
  for (unsigned i = 0; i < evaluators.size(); i++)
    {
      delete evaluators[i];
    }
#endif
}

bool GeometryFactory::read_file(const std::string& filename, PerceptMesh* mesh_data)
{
    std::vector<GeometryHandle> geometry_entities;
    if (!geomKernel->read_file(filename, geometry_entities))
        return false;
    for (size_t i=0; i<geometry_entities.size(); i++)
    {
        std::string str = geomKernel->get_attribute(geometry_entities[i]);
        bool partial_string_match_ok = true;
        Part* part = mesh_data->get_non_const_part(str, partial_string_match_ok);
#if DEBUG_GEOM_SNAP
        std::cout << "tmp geom 3dm attribute = " << str << " part= " << part 
                  << " part name= " << (part?part->name():"null") << std::endl;
#endif
        if (part)
        {
            GeometryEvaluator* eval = new GeometryEvaluator(part);
            eval->mGeometry = geometry_entities[i];
            //eval->mGeometry = i;
            geomDatabase->add_evaluator(eval);
        }
    }
    return true;
}
