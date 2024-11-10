// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Surface_Parser.hpp>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_Parser.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

namespace {

Sphere *
parse_sphere(const Parser::Node & ic_node)
{
  std::vector<double> center;
  if (ic_node.get_if_present("center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC sphere.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for IC sphere.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for sphere must be -1 or 1.\n";
    }
  }

  double radius = 0.0;
  if (!ic_node.get_if_present("radius", radius))
  {
    stk::RuntimeDoomedAdHoc() << "Missing radius for IC sphere.\n";
  }

  return new Sphere(stk::math::Vector3d(center.data()), radius, sign);
}

Ellipsoid *
parse_ellipsoid(const Parser::Node & ic_node)
{
  std::vector<double> center;
  if (ic_node.get_if_present("center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC sphere.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for IC sphere.\n";
  }

  std::vector<double> semiaxes;
  if (ic_node.get_if_present("semiaxes", semiaxes))
  {
    if (semiaxes.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for semiaxes for IC ellipsoid.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing semiaxes for IC ellipsoid.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for ellipsoid must be -1 or 1.\n";
    }
  }

  std::vector<double> rotationVec;
  if (ic_node.get_if_present("rotation", rotationVec))
  {
    if (semiaxes.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for rotation for IC ellipsoid.\n";
    }
  }

  return new Ellipsoid(center, semiaxes, rotationVec, sign);
}

Plane *
parse_plane(const Parser::Node & ic_node)
{
  std::vector<double> normal;
  if (ic_node.get_if_present("normal", normal))
  {
    if (normal.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC plane.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing normal for IC plane.\n";
  }

  double multiplier = 1.0;
  ic_node.get_if_present("multiplier", multiplier);

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC plane must be -1 or 1.\n";
    }
  }
  multiplier *= sign;

  double offset = 0.0;
  if (!ic_node.get_if_present("offset", offset))
  {
    stk::RuntimeDoomedAdHoc() << "Missing offset for IC plane.\n";
  }

  return new Plane(normal.data(), offset, multiplier);
}

Cylinder *
parse_cylinder(const Parser::Node & ic_node)
{
  std::vector<double> p1;
  if (ic_node.get_if_present("p1", p1))
  {
    if (p1.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for p1 for IC cylinder.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing normal for IC cylinder.\n";
  }

  std::vector<double> p2;
  if (ic_node.get_if_present("p2", p2))
  {
    if (p2.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC cylinder.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing p2 for IC cylinder.\n";
  }

  double radius = 0.0;
  if (!ic_node.get_if_present("radius", radius))
  {
    stk::RuntimeDoomedAdHoc() << "Missing radius for IC cylinder.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC cylinder must be -1 or 1.\n";
    }
  }

  return new Cylinder(p1.data(), p2.data(), radius, sign);
}

Surface *
parse_facets(const Parser::Node & ic_node, const stk::diag::Timer &parent_timer)
{
  std::string surface_name;
  ic_node.get_if_present("name", surface_name);

  std::string facet_filename;
  if (!ic_node.get_if_present("filename", facet_filename))
  {
    stk::RuntimeDoomedAdHoc() << "Missing filename for IC facets.\n";
  }

  std::string facet_format;
  if (!ic_node.get_if_present("format", facet_format))
  {
    stk::RuntimeDoomedAdHoc() << "Missing format for IC facets.\n";
  }

  bool scaleSpecified = false;
  double scale = 1.0;
  if (ic_node.get_if_present("scale", scale))
  {
    scaleSpecified = true;
    if (scale <= 0.0)
    {
      stk::RuntimeDoomedAdHoc() << "Scale for IC facets must be >= 0.\n";
    }
  }

  stk::math::Vector3d scaleVec{scale, scale, scale};
  std::vector<std::string> scaleComponents = {"scaleX","scaleY","scaleZ"};
  for (int i=0; i<3; i++)
  {
    if (ic_node.get_if_present(scaleComponents[i], scaleVec[i]))
    {
      if (scaleSpecified)
        stk::RuntimeDoomedAdHoc() << "Cannot specify both scale and " << scaleComponents[i] << "\n";
      if (scaleVec[i] <= 0.0)
        stk::RuntimeDoomedAdHoc() << scaleComponents[i] << " for IC facets must be >= 0.\n";
    }
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC facets must be -1 or 1.\n";
    }
  }

  std::transform(facet_format.begin(), facet_format.end(), facet_format.begin(), ::toupper);
  Surface * surface = nullptr;

  if (facet_format == "STL")
    surface = new STLSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "FAC")
    surface = new FACSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "PLY")
    surface = new PLYSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "EXO")
    surface = new EXOSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else
    stk::RuntimeDoomedAdHoc() << "Unrecognized facet format: " << facet_format;

  return surface;
}

Surface *
parse_mesh_surface(const Parser::Node & ic_node, const stk::mesh::MetaData & meta)
{
  std::string surface_name;
  STK_ThrowRequire(ic_node.get_if_present("mesh", surface_name));
  const AuxMetaData & auxMeta = AuxMetaData::get(meta);

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC facets must be -1 or 1.\n";
    }
  }

  STK_ThrowErrorMsgIf(!auxMeta.has_part(surface_name), "Could not locate a surface named " << surface_name);
  const stk::mesh::Part & io_part = auxMeta.get_part(surface_name);

  const stk::mesh::Field<double>* coords = reinterpret_cast<const stk::mesh::Field<double>*>(&auxMeta.get_current_coordinates().field());
  STK_ThrowRequire(nullptr != coords);

  const stk::mesh::Selector surfaceSelector = stk::mesh::Selector(io_part);

  if (meta.spatial_dimension() == 2)
    return new MeshSurface<Facet2d>(meta, *coords, surfaceSelector, sign);
  return new MeshSurface<Facet3d>(meta, *coords, surfaceSelector, sign);
}

LevelSet_String_Function *
parse_string_function(const Parser::Node & ic_node)
{
  std::string expression;
  if (!ic_node.get_if_present("expression", expression))
  {
    stk::RuntimeDoomedAdHoc() << "Missing expression for string_function.\n";
  }

  LevelSet_String_Function * surf = new LevelSet_String_Function(expression);

  std::vector<double> bounds;
  if (ic_node.get_if_present("bounding_box", bounds))
  {
    if (bounds.size() == 6)
    {
      const BoundingBox surfBbox( stk::math::Vector3d(bounds[0],bounds[1],bounds[2]), stk::math::Vector3d(bounds[3],bounds[4],bounds[5]) );
      surf->set_bounding_box(surfBbox);
    }
    else
    {
      stk::RuntimeDoomedAdHoc() << "bounding_box for string_function must be a vector of length 6 (for both 2D or 3D) (xmin,ymin,zmin, xmax,ymax,zmax).\n";
    }
  }

  return surf;
}

}

Surface *
Surface_Parser::parse(const Parser::Node & parserNode, const stk::mesh::MetaData & meta, const stk::diag::Timer &parentTimer)
{
  if (parserNode.is_scalar())
  {
    // support random specified as Scalar (no :)
    std::string ic_type = parserNode.as<std::string>();
    if (ic_type == "random")
    {
      return new Random(0);
    }
    return nullptr;
  }

  if ( parserNode.get_null_if_present("sphere") )
  {
    return parse_sphere(parserNode);
  }
  else if ( parserNode.get_null_if_present("ellipsoid") )
  {
    return parse_ellipsoid(parserNode);
  }
  else if ( parserNode.get_null_if_present("plane") )
  {
    return parse_plane(parserNode);
  }
  else if ( parserNode.get_null_if_present("cylinder") )
  {
    return parse_cylinder(parserNode);
  }
  else if ( parserNode.get_null_if_present("string_function") )
  {
    return parse_string_function(parserNode);
  }
  else if ( parserNode.get_null_if_present("facets") )
  {
    return parse_facets(parserNode, parentTimer);
  }
  else if ( parserNode.get_scalar_if_present("mesh") )
  {
    return parse_mesh_surface(parserNode, meta);
  }
  else if ( parserNode.get_null_if_present("random") )
  {
    int seed = 0;
    parserNode.get_if_present("seed", seed);
    return new Random(seed);
  }
  return nullptr;
}

} // namespace krino
