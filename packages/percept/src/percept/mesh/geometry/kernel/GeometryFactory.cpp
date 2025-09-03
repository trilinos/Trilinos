// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "GeometryFactory.hpp"
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <stk_util/diag/StringUtil.hpp>

namespace percept {

GeometryFactory::GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry)
{
    geomKernel = kernel;
    geomDatabase = geometry;
}

GeometryFactory::~GeometryFactory()
{
#if 0
  const std::vector<GeometryEvaluator*>& evaluators = geomDatabase->getGeomEvaluators();
  for (unsigned i = 0; i < evaluators.size(); i++)
    {
      delete evaluators[i];
    }
#endif
}


static stk::mesh::Part* 
getPart(stk::mesh::MetaData *meta_data, std::string part_name, bool partial_string_match_ok)
{
  stk::mesh::Part* found_part =0;
  sierra::make_lower(part_name);
  if (!partial_string_match_ok)
    {
      found_part = meta_data->get_part(part_name);
    }
  else
    {
      const stk::mesh::PartVector & parts = meta_data->get_mesh_parts();
      unsigned nparts = parts.size();

      unsigned pname_len_min = std::numeric_limits<unsigned>::max();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          std::string stk_part_name = sierra::make_lower(part.name());
          size_t found = stk_part_name.find(part_name);
          if (found != std::string::npos)
            {
              // skip "old" parts
              if (stk_part_name.find(PerceptMesh::s_omit_part) == std::string::npos)
                {
                  if (stk_part_name.length() < pname_len_min)
                    {
                      found_part = &part;
                      pname_len_min = stk_part_name.length();
                    }
                }
            }
        }
    }

  const bool error_check = true;
  if (error_check && !found_part && part_name != "edgeseams")
    {
      std::ostringstream msg;
      msg << "GeometryFactor::getPart() couldn't find part with name = " << part_name;
      std::cout << msg.str() << std::endl;
      const stk::mesh::PartVector & parts = meta_data->get_parts();
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

bool GeometryFactory::read_file(const std::string& filename, stk::mesh::MetaData *meta_data)
{
  std::vector<GeometryHandle> geometry_entities;
  if (!geomKernel->read_file(filename, geometry_entities))
    return false;
  for (size_t i=0; i<geometry_entities.size(); i++)
    {
      std::string str = geometry_entities[i].attribute;
      bool partial_string_match_ok = true;

      if( (str.substr(0,1))=="#"){

        std::string character = "";
        std::string newstr1 = "";
        int pos = 1;
        while ( character != " "){ //parse out first partname
          character = str.substr(pos,1);
          if (character == " ")
            break;
          newstr1 = newstr1 + character;
          pos++;
        }
        pos++;

        std::string newstr2 = str.substr(pos, str.length()-pos);
        stk::mesh::Part* part = getPart(meta_data, newstr1, partial_string_match_ok);
        if (part)
          {
            GeometryEvaluator* eval = new GeometryEvaluator(part);
            eval->mGeometry = geometry_entities[i];
            geomDatabase->add_evaluator(eval);
          }
        stk::mesh::Part* part2 = getPart(meta_data, newstr2, partial_string_match_ok);
        if (part2)
          {
            GeometryEvaluator* eval = new GeometryEvaluator(part2);
            eval->mGeometry = geometry_entities[i];
            geomDatabase->add_evaluator(eval);
          }

      }
      else{
        stk::mesh::Part* part = getPart(meta_data, str, partial_string_match_ok);
        if (part)
          {
            GeometryEvaluator* eval = new GeometryEvaluator(part);
            eval->mGeometry = geometry_entities[i];
            geomDatabase->add_evaluator(eval);
          }
      }
    }
  return true;
}

} // namespace percept
