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
#include <Akri_MeshSurface.hpp>
#include <Akri_Parser.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <cmath>

namespace krino {

namespace {

void validate_sign(double sign, const std::string & surfaceName)
{
  if (sign != -1.0 && sign != 1.0)
  {
    stk::RuntimeDoomedAdHoc()
        << "Sign for " << surfaceName << " must be -1 or 1.\n";
  }
}

bool is_finite_vector(const std::vector<double> & values)
{
  for (const double value : values)
    if (!std::isfinite(value))
      return false;
  return true;
}

bool had_new_parse_error(const unsigned oldDoomedCount)
{
  return stk::get_doomed_count() > oldDoomedCount;
}

std::shared_ptr<Sphere>
parse_sphere(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::vector<double> center;
  if (ic_node.get_if_present("center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for sphere.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for sphere.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    validate_sign(sign, "sphere");
  }

  double radius = 0.0;
  if (!ic_node.get_if_present("radius", radius))
  {
    stk::RuntimeDoomedAdHoc() << "Missing radius for sphere.\n";
  }
  else if (radius <= 0.0)
  {
    stk::RuntimeDoomedAdHoc() << "Radius for sphere must be > 0.\n";
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;
  return std::make_shared<Sphere>(stk::math::Vector3d(center.data()), radius, sign);
}

std::shared_ptr<Ellipsoid>
parse_ellipsoid(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::vector<double> center;
  if (ic_node.get_if_present("center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for ellipsoid.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for ellipsoid.\n";
  }

  std::vector<double> semiaxes;
  if (ic_node.get_if_present("semiaxes", semiaxes))
  {
    if (semiaxes.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for semiaxes for ellipsoid.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing semiaxes for ellipsoid.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    validate_sign(sign, "ellipsoid");
  }

  std::vector<double> rotationVec;
  if (ic_node.get_if_present("rotation", rotationVec))
  {
    if (rotationVec.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for rotation for ellipsoid.\n";
    }
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;
  return std::make_shared<Ellipsoid>(center, semiaxes, rotationVec, sign);
}

std::shared_ptr<Cuboid>
parse_cuboid(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::vector<double> center;
  if (ic_node.get_if_present("center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for cuboid.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for cuboid.\n";
  }

  std::vector<double> dimensions;
  if (ic_node.get_if_present("dimensions", dimensions))
  {
    if (dimensions.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for dimensions for cuboid.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing dimensions for cuboid.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    validate_sign(sign, "cuboid");
  }

  std::vector<double> rotation;
  if (ic_node.get_if_present("rotation", rotation))
  {
    if (rotation.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for rotation for cuboid.\n";
    }

    if (had_new_parse_error(oldDoomedCount))
      return nullptr;
    return std::make_shared<Cuboid>(stk::math::Vector3d(center.data()), stk::math::Vector3d(dimensions.data()), stk::math::Vector3d(rotation.data()), sign);
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;
  return std::make_shared<Cuboid>(stk::math::Vector3d(center.data()), stk::math::Vector3d(dimensions.data()), sign);
}

std::shared_ptr<Plane>
parse_plane(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::vector<double> normal;
  if (ic_node.get_if_present("normal", normal))
  {
    if (normal.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for normal for plane.\n";
    }
    else
    {
      if (!is_finite_vector(normal))
      {
        stk::RuntimeDoomedAdHoc() << "Plane normal values must be finite.\n";
      }

      const double normalMagnitudeSquared =
          normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];

      if (normalMagnitudeSquared == 0.0)
      {
        stk::RuntimeDoomedAdHoc() << "Plane normal must have nonzero magnitude.\n";
      }
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing normal for plane.\n";
  }

  double multiplier = 1.0;
  if (ic_node.get_if_present("multiplier", multiplier))
  {
    if (!std::isfinite(multiplier))
    {
      stk::RuntimeDoomedAdHoc() << "Plane multiplier must be finite.\n";
    }
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    validate_sign(sign, "plane");
  }
  multiplier *= sign;

  double offset = 0.0;
  const bool hasOffset = ic_node.get_if_present("offset", offset);
  if (hasOffset && !std::isfinite(offset))
  {
    stk::RuntimeDoomedAdHoc() << "Plane offset must be finite.\n";
  }

  std::vector<double> point;
  const bool hasPoint = ic_node.get_if_present("point", point);
  if (hasPoint)
  {
    if (point.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for point for plane.\n";
    }
    else if (!is_finite_vector(point))
    {
      stk::RuntimeDoomedAdHoc() << "Plane point values must be finite.\n";
    }
  }

  double fractionOfMeshExtent = 0.0;
  const bool hasFractionOfMeshExtent = ic_node.get_if_present("fraction_of_mesh_extent", fractionOfMeshExtent);
  if (hasFractionOfMeshExtent && !std::isfinite(fractionOfMeshExtent))
  {
    stk::RuntimeDoomedAdHoc() << "Plane fraction_of_mesh_extent must be finite.\n";
  }

  const int numPlacements = int(hasOffset) + int(hasPoint) + int(hasFractionOfMeshExtent);
  if (numPlacements > 1)
  {
    stk::RuntimeDoomedAdHoc() << "Specify only one of offset, point or fraction_of_mesh_extent for plane.\n";
  }
  else if (numPlacements == 0)
  {
    stk::RuntimeDoomedAdHoc() << "Missing offset, point or fraction_of_mesh_extent for plane.\n";
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;

  if (hasPoint)
    return std::make_shared<Plane>(normal.data(), point.data(), multiplier);

  if (hasFractionOfMeshExtent)
    return Plane::at_fraction_of_extent(normal.data(), fractionOfMeshExtent, multiplier);

  return std::make_shared<Plane>(normal.data(), offset, multiplier);
}

std::shared_ptr<Cylinder>
parse_cylinder(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::vector<double> p1;
  if (ic_node.get_if_present("p1", p1))
  {
    if (p1.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for p1 for cylinder.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing p1 for cylinder.\n";
  }

  std::vector<double> p2;
  if (ic_node.get_if_present("p2", p2))
  {
    if (p2.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for p2 for cylinder.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing p2 for cylinder.\n";
  }

  double radius = 0.0;
  if (!ic_node.get_if_present("radius", radius))
  {
    stk::RuntimeDoomedAdHoc() << "Missing radius for cylinder.\n";

  }
  else if (radius <= 0.0)
  {
    stk::RuntimeDoomedAdHoc() << "Radius for cylinder must be > 0.\n";
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    validate_sign(sign, "cylinder");
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;
  return std::make_shared<Cylinder>(stk::math::Vector3d(p1.data()), stk::math::Vector3d(p2.data()), radius, sign);
}

std::shared_ptr<Faceted_Surface_From_File>
parse_facets(const Parser::Node & ic_node, const stk::diag::Timer &parent_timer)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::string facet_filename;
  if (!ic_node.get_if_present("filename", facet_filename))
  {
    stk::RuntimeDoomedAdHoc() << "Missing filename for facets.\n";
  }

  std::string facet_format;
  if (!ic_node.get_if_present("format", facet_format))
  {
    stk::RuntimeDoomedAdHoc() << "Missing format for facets.\n";
  }

  bool scaleSpecified = false;
  double scale = 1.0;
  if (ic_node.get_if_present("scale", scale))
  {
    scaleSpecified = true;
    if (scale <= 0.0)
    {
      stk::RuntimeDoomedAdHoc() << "Scale for facets must be > 0.\n";
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
        stk::RuntimeDoomedAdHoc() << scaleComponents[i] << " for facets must be >= 0.\n";
    }
  }

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    validate_sign(sign, "facets");
  }

  bool removeCoincidentFacets = false;
  ic_node.get_if_present("remove_coincident_facets", removeCoincidentFacets);

  bool checkWatertight = false;
  ic_node.get_if_present("check_watertight", checkWatertight);

  std::vector<double> watertightBounds;
  if (ic_node.get_if_present("watertight_bounding_box", watertightBounds))
  {
    if (watertightBounds.size() != 6)
    {
      stk::RuntimeDoomedAdHoc() << "watertight_bounding_box for facets must be a vector of length 6 (for both 2D or 3D) (xmin,ymin,zmin, xmax,ymax,zmax).\n";
    }
  }

  std::transform(facet_format.begin(), facet_format.end(), facet_format.begin(), ::toupper);

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;

  std::shared_ptr<Faceted_Surface_From_File> surface;
  if (facet_format == "STL")
    surface = std::make_shared<STLSurface>(parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "FAC")
    surface = std::make_shared<FACSurface>(parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "PLY")
    surface = std::make_shared<PLYSurface>(parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "EXO")
    surface = std::make_shared<EXOSurface>(parent_timer, facet_filename, sign, scaleVec);
  else
  {
    stk::RuntimeDoomedAdHoc() << "Unrecognized facet format: " << facet_format;
    return nullptr;
  }

  surface->set_remove_coincident_facets(removeCoincidentFacets);
  surface->set_check_watertight(checkWatertight);
  if (watertightBounds.size() == 6)
    surface->set_watertight_bounding_box(BoundingBox(stk::math::Vector3d(watertightBounds[0],watertightBounds[1],watertightBounds[2]), stk::math::Vector3d(watertightBounds[3],watertightBounds[4],watertightBounds[5])));
  return surface;
}

std::shared_ptr<FacetedSurfaceBase>
parse_mesh_surface(const Parser::Node & ic_node, const stk::mesh::MetaData & meta)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::string surface_name;

  if (!ic_node.get_if_present("mesh", surface_name) &&
      !ic_node.get_if_present("name", surface_name))
  {
    stk::RuntimeDoomedAdHoc()
        << "Missing mesh surface name. Expected 'mesh' or 'name'.\n";
  }

  const AuxMetaData & auxMeta = AuxMetaData::get(meta);

  double sign = 1.0;
  if (ic_node.get_if_present("sign", sign))
  {
    if (sign != -1.0 && sign != 1.0)
    {
      stk::RuntimeDoomedAdHoc()
          << "Sign for mesh surface must be -1 or 1.\n";
    }
  }

  bool removeCoincidentFacets = false;
  ic_node.get_if_present("remove_coincident_facets", removeCoincidentFacets);

  bool checkWatertight = false;
  ic_node.get_if_present("check_watertight", checkWatertight);

  std::vector<double> watertightBounds;
  if (ic_node.get_if_present("watertight_bounding_box", watertightBounds))
  {
    if (watertightBounds.size() != 6)
    {
      stk::RuntimeDoomedAdHoc() << "watertight_bounding_box for mesh surface must be a vector of length 6 (for both 2D or 3D) (xmin,ymin,zmin, xmax,ymax,zmax).\n";
    }
  }

  if (!auxMeta.has_part(surface_name))
  {
    stk::RuntimeDoomedAdHoc() << "Could not locate a surface named " << surface_name;
  }

  const stk::mesh::Field<double>* coords =
      reinterpret_cast<const stk::mesh::Field<double>*>(
          &auxMeta.get_current_coordinates().field());

  if (coords == nullptr)
  {
    stk::RuntimeDoomedAdHoc()
        << "Could not locate current coordinates field for mesh surface.\n";
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;

  const stk::mesh::Part & io_part = auxMeta.get_part(surface_name);
  const stk::mesh::Selector surfaceSelector = stk::mesh::Selector(io_part);

  std::shared_ptr<FacetedSurfaceBase> surface;
  if (meta.spatial_dimension() == 2)
    surface = std::make_shared<MeshSurface<Facet2d>>(meta, *coords, surfaceSelector, sign);
  else
    surface = std::make_shared<MeshSurface<Facet3d>>(meta, *coords, surfaceSelector, sign);

  surface->set_remove_coincident_facets(removeCoincidentFacets);
  surface->set_check_watertight(checkWatertight);
  if (watertightBounds.size() == 6)
    surface->set_watertight_bounding_box(BoundingBox(stk::math::Vector3d(watertightBounds[0],watertightBounds[1],watertightBounds[2]), stk::math::Vector3d(watertightBounds[3],watertightBounds[4],watertightBounds[5])));
  return surface;
}

std::shared_ptr<LevelSet_String_Function>
parse_string_function(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::string expression;
  if (!ic_node.get_if_present("expression", expression))
  {
    stk::RuntimeDoomedAdHoc() << "Missing expression for string_function.\n";
  }

  std::vector<double> bounds;
  if (ic_node.get_if_present("bounding_box", bounds))
  {
    if (bounds.size() != 6)
    {
      stk::RuntimeDoomedAdHoc() << "bounding_box for string_function must be a vector of length 6 (for both 2D or 3D) (xmin,ymin,zmin, xmax,ymax,zmax).\n";
    }
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;

  auto surf = std::make_shared<LevelSet_String_Function>(expression);
  if (!bounds.empty())
  {
    const BoundingBox surfBbox( stk::math::Vector3d(bounds[0],bounds[1],bounds[2]), stk::math::Vector3d(bounds[3],bounds[4],bounds[5]) );
    surf->set_bounding_box(surfBbox);
  }
  return surf;
}

std::shared_ptr<Implicit_Neural_Representation>
parse_implicit_neural_representation(const Parser::Node & ic_node)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::string model_filename;
  if (!ic_node.get_if_present("filename", model_filename))
  {
    stk::RuntimeDoomedAdHoc() << "Missing model filename for implicit neural representation.\n";
  }

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;
  return std::make_shared<Implicit_Neural_Representation>(model_filename);
}

}

std::shared_ptr<Surface>
Surface_Parser::parse(
    const Parser::Node & parserNode,
    const stk::mesh::MetaData & meta,
    const stk::diag::Timer & parentTimer)
{
  if (parserNode.is_scalar())
  {
    std::string ic_type = parserNode.as<std::string>();
    if (ic_type == "random")
    {
      return std::make_shared<Random>(0);
    }
    stk::RuntimeDoomedAdHoc()
        << "Unrecognized scalar surface specification '" << ic_type << "'. "
        << "The only supported scalar surface specification is 'random'. "
        << "For other surfaces, use a map such as 'sphere:', 'cylinder:', "
        << "'cuboid:', 'plane:', etc. "
        << parserNode.info();

    return nullptr;
  }

  Parser::Node surfaceNode;

  surfaceNode = parserNode.get_options_if_present("sphere");
  if (surfaceNode)
  {
    return parse_sphere(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("ellipsoid");
  if (surfaceNode)
  {
    return parse_ellipsoid(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("cuboid");
  if (surfaceNode)
  {
    return parse_cuboid(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("plane");
  if (surfaceNode)
  {
    return parse_plane(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("cylinder");
  if (surfaceNode)
  {
    return parse_cylinder(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("string_function");
  if (surfaceNode)
  {
    return parse_string_function(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("facets");
  if (surfaceNode)
  {
    return parse_facets(surfaceNode, parentTimer);
  }

  surfaceNode = parserNode.get_options_if_present("implicit_neural_representation");
  if (surfaceNode)
  {
    return parse_implicit_neural_representation(surfaceNode);
  }

  surfaceNode = parserNode.get_options_if_present("random");
  if (surfaceNode)
  {
    int seed = 0;
    surfaceNode.get_if_present("seed", seed);
    return std::make_shared<Random>(seed);
  }

  /*
   * Mesh is a little different because the existing flat syntax is:
   *
   *   mesh: surface_part_name
   *   sign: 1.0
   *
   * so handle it separately.
   */
  if (parserNode.get_scalar_if_present("mesh"))
  {
    return parse_mesh_surface(parserNode, meta);
  }

  const Parser::Node meshNode = parserNode.get_map_if_present("mesh");
  if (meshNode)
  {
    return parse_mesh_surface(meshNode, meta);
  }

  return nullptr;
}

} // namespace krino
