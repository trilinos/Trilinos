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
#include <Akri_FacetFileReader.hpp>
#include <Akri_FacetFileSurface.hpp>
#include <Akri_MeshBlockSkinSurface.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_Parser.hpp>
#include <Akri_Transformation.hpp>
#include <Akri_TransformedSurface.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <algorithm>
#include <cmath>
#include <optional>

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

// The sharp feature options of a faceted surface: sharp_feature_angle (in degrees) turns on its
// sharp features, and open_boundary_edges_are_sharp_features also makes its open boundary edges
// features.
struct SharpFeatureOptions
{
  std::optional<double> featureAngleInDegrees;
  bool areOpenBoundaryEdgesFeatures{false};
};

static SharpFeatureOptions
parse_sharp_feature_options(const Parser::Node & ic_node, const std::string & surfaceType)
{
  SharpFeatureOptions options;
  double featureAngle = 0.;
  if (ic_node.get_if_present("sharp_feature_angle", featureAngle))
  {
    if (!(featureAngle > 0. && featureAngle < 180.))
      stk::RuntimeDoomedAdHoc() << "sharp_feature_angle for " << surfaceType
                                << " must be more than 0 and less than 180 degrees.\n";
    options.featureAngleInDegrees = featureAngle;
  }
  ic_node.get_if_present(
      "open_boundary_edges_are_sharp_features", options.areOpenBoundaryEdgesFeatures);
  if (options.areOpenBoundaryEdgesFeatures && !options.featureAngleInDegrees)
    stk::RuntimeDoomedAdHoc() << "open_boundary_edges_are_sharp_features for " << surfaceType
                              << " needs a sharp_feature_angle.\n";
  return options;
}

static void
set_sharp_feature_options(MeshFacetedSurface & surface, const SharpFeatureOptions & options)
{
  if (options.featureAngleInDegrees)
    surface.set_sharp_feature_angle(
        *options.featureAngleInDegrees, options.areOpenBoundaryEdgesFeatures);
}

static std::shared_ptr<Surface> parse_facets(const Parser::Node & ic_node)
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
        stk::RuntimeDoomedAdHoc() << scaleComponents[i] << " for facets must be > 0.\n";
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

  // Facets are read into a connected surface (MeshFacetedSurface), by default with coincident
  // vertices merged so that facets meeting at a point share its vertex.
  bool mergeCoincidentVertices = true;
  ic_node.get_if_present("merge_coincident_vertices", mergeCoincidentVertices);

  const SharpFeatureOptions sharpFeatureOptions = parse_sharp_feature_options(ic_node, "facets");
  if (sharpFeatureOptions.featureAngleInDegrees && !mergeCoincidentVertices)
    stk::RuntimeDoomedAdHoc()
        << "sharp_feature_angle for facets needs merge_coincident_vertices.\n";

  std::transform(facet_format.begin(), facet_format.end(), facet_format.begin(), ::toupper);
  if (!is_facet_file_format(facet_format))
    stk::RuntimeDoomedAdHoc() << "Unrecognized facet format: " << facet_format;

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;

  auto surface = std::make_shared<FacetFileSurface>(stk::EnvData::parallel_comm(),
      facet_format,
      facet_filename,
      sign,
      scaleVec,
      mergeCoincidentVertices);
  surface->set_remove_coincident_facets(removeCoincidentFacets);
  surface->set_check_watertight(checkWatertight);
  if (watertightBounds.size() == 6)
    surface->set_watertight_bounding_box(BoundingBox(
        stk::math::Vector3d(watertightBounds[0], watertightBounds[1], watertightBounds[2]),
        stk::math::Vector3d(watertightBounds[3], watertightBounds[4], watertightBounds[5])));
  set_sharp_feature_options(*surface, sharpFeatureOptions);
  return surface;
}

// The names of the blocks of a mesh_block_skin: one, as block, or a list, as blocks.
static std::vector<std::string> parse_mesh_block_skin_block_names(const Parser::Node & skinNode)
{
  std::vector<std::string> blockNames;
  std::string blockName;
  const bool hasBlock = skinNode.get_if_present("block", blockName);
  if (hasBlock)
    blockNames.push_back(blockName);
  std::vector<std::string> listedBlockNames;
  if (skinNode.get_if_present("blocks", listedBlockNames))
  {
    if (hasBlock)
      stk::RuntimeDoomedAdHoc() << "mesh_block_skin takes either block or blocks, not both.\n";
    blockNames = listedBlockNames;
  }
  if (blockNames.empty())
    stk::RuntimeDoomedAdHoc() << "mesh_block_skin needs a block, or a list of blocks.\n";
  return blockNames;
}

// The skin of blocks of an Exodus mesh file (see MeshBlockSkinSurface): filename, block or blocks, sign
// and the sharp feature options.  The file is read when the surface is first prepared.
static std::shared_ptr<Surface> parse_mesh_block_skin(const Parser::Node & skinNode, const stk::mesh::MetaData & meta)
{
  const auto oldDoomedCount = stk::get_doomed_count();

  std::string filename;
  if (!skinNode.get_if_present("filename", filename))
    stk::RuntimeDoomedAdHoc() << "Missing filename for mesh_block_skin.\n";

  const std::vector<std::string> blockNames = parse_mesh_block_skin_block_names(skinNode);

  double sign = 1.0;
  if (skinNode.get_if_present("sign", sign))
    validate_sign(sign, "mesh_block_skin");

  const SharpFeatureOptions sharpFeatureOptions = parse_sharp_feature_options(skinNode, "mesh_block_skin");

  if (had_new_parse_error(oldDoomedCount))
    return nullptr;

  auto surface = std::make_shared<MeshBlockSkinSurface>(meta.spatial_dimension(),
      get_block_skin_mesh_file(stk::EnvData::parallel_comm(), filename),
      blockNames,
      int(sign));
  set_sharp_feature_options(*surface, sharpFeatureOptions);
  return surface;
}

std::shared_ptr<MeshSurface>
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

  const SharpFeatureOptions sharpFeatureOptions =
      parse_sharp_feature_options(ic_node, "mesh surface");

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

  auto surface = std::make_shared<MeshSurface>(meta, *coords, surfaceSelector, sign);

  surface->set_remove_coincident_facets(removeCoincidentFacets);
  surface->set_check_watertight(checkWatertight);
  if (watertightBounds.size() == 6)
    surface->set_watertight_bounding_box(BoundingBox(stk::math::Vector3d(watertightBounds[0],watertightBounds[1],watertightBounds[2]), stk::math::Vector3d(watertightBounds[3],watertightBounds[4],watertightBounds[5])));
  set_sharp_feature_options(*surface, sharpFeatureOptions);
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

// A vector of 3 reals for key in a motion, or of 2 with z = 0 for those that are not rotations.
static std::optional<stk::math::Vector3d> parse_motion_vector(const Parser::Node & motionNode,
    const std::string & key,
    const bool isRotation,
    const std::string & surfaceType)
{
  std::vector<double> values;
  if (!motionNode.get_if_present(key, values))
    return std::nullopt;
  if (values.size() == 3)
    return stk::math::Vector3d(values.data());
  if (!isRotation && values.size() == 2)
    return stk::math::Vector3d(values[0], values[1], 0.);
  stk::RuntimeDoomedAdHoc() << "Expecting " << (isRotation ? "3" : "2 or 3") << " real values for " << key
                            << " in the motion of " << surfaceType << ".\n";
  return std::nullopt;
}

// The motion of a surface, given by the map motion in its options, if present (see Transformation).
static std::optional<Transformation> parse_motion(const Parser::Node & surfaceNode, const std::string & surfaceType)
{
  const Parser::Node motionNode = surfaceNode.get_if_present("motion");
  if (!motionNode)
    return std::nullopt;
  if (!motionNode.is_map())
  {
    stk::RuntimeDoomedAdHoc() << "The motion of " << surfaceType << " must be a map of "
                              << "translational_velocity, rotational_velocity, reference_point, initial_displacement and initial_rotation.\n";
    return std::nullopt;
  }

  static const std::vector<std::string> motionKeys{"translational_velocity", "rotational_velocity", "reference_point", "initial_displacement", "initial_rotation"};
  for (auto & key : motionNode.keys())
    if (std::find(motionKeys.begin(), motionKeys.end(), key) == motionKeys.end())
      stk::RuntimeDoomedAdHoc() << "Unrecognized " << key << " in the motion of " << surfaceType << ".\n";

  Transformation motion;
  if (const auto v = parse_motion_vector(motionNode, "translational_velocity", false, surfaceType))
    motion.set_translational_velocity(*v);
  if (const auto v = parse_motion_vector(motionNode, "rotational_velocity", true, surfaceType))
    motion.set_rotational_velocity(*v);
  if (const auto v = parse_motion_vector(motionNode, "reference_point", false, surfaceType))
    motion.set_reference_point(*v);
  if (const auto v = parse_motion_vector(motionNode, "initial_displacement", false, surfaceType))
    motion.set_initial_displacement(*v);
  if (const auto v = parse_motion_vector(motionNode, "initial_rotation", true, surfaceType))
    motion.set_initial_rotation(*v);
  return motion;
}

// The surface, moved by the motion in its options if there is one.
static std::shared_ptr<Surface> with_motion_if_present(std::shared_ptr<Surface> surface,
    const Parser::Node & surfaceNode,
    const std::string & surfaceType)
{
  const auto oldDoomedCount = stk::get_doomed_count();
  const std::optional<Transformation> motion = parse_motion(surfaceNode, surfaceType);
  if (!surface || had_new_parse_error(oldDoomedCount))
    return nullptr;
  if (!motion)
    return surface;

  const Plane * plane = dynamic_cast<const Plane *>(surface.get());
  if (nullptr != plane && plane->needs_offset_from_extent())
  {
    stk::RuntimeDoomedAdHoc() << "Motion is not supported for a plane at a fraction_of_mesh_extent.\n";
    return nullptr;
  }
  return std::make_shared<TransformedSurface>(surface, *motion);
}

std::shared_ptr<Surface> Surface_Parser::parse(const Parser::Node & parserNode,
    const stk::mesh::MetaData & meta,
    const stk::diag::Timer & /*parentTimer*/)
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
    return with_motion_if_present(parse_sphere(surfaceNode), surfaceNode, "sphere");
  }

  surfaceNode = parserNode.get_options_if_present("ellipsoid");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_ellipsoid(surfaceNode), surfaceNode, "ellipsoid");
  }

  surfaceNode = parserNode.get_options_if_present("cuboid");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_cuboid(surfaceNode), surfaceNode, "cuboid");
  }

  surfaceNode = parserNode.get_options_if_present("plane");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_plane(surfaceNode), surfaceNode, "plane");
  }

  surfaceNode = parserNode.get_options_if_present("cylinder");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_cylinder(surfaceNode), surfaceNode, "cylinder");
  }

  surfaceNode = parserNode.get_options_if_present("string_function");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_string_function(surfaceNode), surfaceNode, "string_function");
  }

  surfaceNode = parserNode.get_options_if_present("facets");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_facets(surfaceNode), surfaceNode, "facets");
  }

  surfaceNode = parserNode.get_options_if_present("mesh_block_skin");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_mesh_block_skin(surfaceNode, meta), surfaceNode, "mesh_block_skin");
  }

  surfaceNode = parserNode.get_options_if_present("implicit_neural_representation");
  if (surfaceNode)
  {
    return with_motion_if_present(parse_implicit_neural_representation(surfaceNode), surfaceNode, "implicit_neural_representation");
  }

  surfaceNode = parserNode.get_options_if_present("random");
  if (surfaceNode)
  {
    int seed = 0;
    surfaceNode.get_if_present("seed", seed);
    return with_motion_if_present(std::make_shared<Random>(seed), surfaceNode, "random");
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
    return with_motion_if_present(parse_mesh_surface(parserNode, meta), parserNode, "mesh");
  }

  const Parser::Node meshNode = parserNode.get_map_if_present("mesh");
  if (meshNode)
  {
    return with_motion_if_present(parse_mesh_surface(meshNode, meta), meshNode, "mesh");
  }

  return nullptr;
}

} // namespace krino
