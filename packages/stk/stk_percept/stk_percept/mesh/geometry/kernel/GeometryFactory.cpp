#include "GeometryFactory.hpp"
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <boost/algorithm/string.hpp>    


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

static stk::mesh::Part* 
getPart(PerceptMesh *eMesh, std::string part_name, bool partial_string_match_ok)
{
  stk::mesh::Part* found_part =0;
  boost::algorithm::to_lower(part_name);
  if (!partial_string_match_ok)
    {
      found_part = eMesh->get_fem_meta_data()->get_part(part_name);
    }
  else
    {
      const stk::mesh::PartVector & parts = eMesh->get_fem_meta_data()->get_parts();
      unsigned nparts = parts.size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          size_t found = part.name().find(part_name);
          if (found != std::string::npos)
            {
              // skip "old" parts
              if (part.name().find(PerceptMesh::s_omit_part) == std::string::npos)
                {
                  found_part = &part;
                  break;
                }
            }
        }
    }
  const bool error_check = true;
  if (error_check && !found_part)
    {
      std::ostringstream msg;
      msg << "GeometryFactor::getPart() couldn't find part with name = " << part_name;
      std::cout << msg.str() << std::endl;
      const stk::mesh::PartVector & parts = eMesh->get_fem_meta_data()->get_parts();
      unsigned nparts = parts.size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          if (!stk::mesh::is_auto_declared_part(part))
            std::cout << "part= " << part.name() << std::endl;
        }
      throw std::runtime_error(msg.str());
    }
  return found_part;
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
        //Part* part = mesh_data->get_non_const_part(str, partial_string_match_ok);
        Part* part = getPart(mesh_data, str, partial_string_match_ok);
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
