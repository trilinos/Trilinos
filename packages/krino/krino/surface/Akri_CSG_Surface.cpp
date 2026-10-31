// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CSG_Surface.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_CSG_Tree_Output.hpp>
#include <Akri_Filtered_CSG_Tree.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <limits>

namespace krino{

CSG_Surface::CSG_Surface(const CSGTree & tree)
{
  set_tree(tree);
}

void CSG_Surface::set_tree(const CSGTree & tree)
{
  myAllTree = tree;
  myAllSurfaces = myAllTree.gather_surfaces();

  int subsurfaces_might_produce_wrong_sign = false;
  for ( auto&& surface : myAllSurfaces )
  {
    if (surface->truncated_distance_may_have_wrong_sign())
    {
      subsurfaces_might_produce_wrong_sign = true;
      break;
    }
  }
  int local_subsurfaces_might_produce_wrong_sign = subsurfaces_might_produce_wrong_sign;
  stk::all_reduce_max(stk::EnvData::parallel_comm(), &local_subsurfaces_might_produce_wrong_sign, &subsurfaces_might_produce_wrong_sign, 1);
  my_subsurfaces_might_produce_wrong_sign = subsurfaces_might_produce_wrong_sign;
}

CSG_Surface::~CSG_Surface() {}

size_t CSG_Surface::storage_size() const
{
  size_t tot_size = 0;
  for (auto && surf : myAllSurfaces)
    tot_size += surf->storage_size();
  return tot_size;
}

void CSG_Surface::set_transformation(Transformation * trans)
{
  Surface::set_transformation(trans);
  for (auto && subsurf : myAllSurfaces)
    subsurf->set_transformation(trans);
}

bool CSG_Surface::empty() const
{
  return myAllTree.empty();
}

void CSG_Surface::clear()
{
  myAllTree.clear();
  myNearbyTree.clear();
}

void
CSG_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  for ( auto&& surface : myAllSurfaces )
  {
    surface->prepare_to_compute(time, point_bbox, truncation_length);
  }

  if (truncation_length > 0.)
  {
    BoundingBox padded_point_bbox = point_bbox;
    padded_point_bbox.pad(truncation_length);

    myNearbySurfaces.clear();
    for ( auto&& surface : myAllSurfaces )
    {
      // keep nearby surfaces
      if (FACETED_SURFACE_2D == surface->type() ||
          FACETED_SURFACE_3D == surface->type() ||
          surface->does_intersect(padded_point_bbox))
      {
        myNearbySurfaces.emplace_back(surface);
      }
    }
    myNearbyTree = filtered_CSG_tree(myAllTree, myNearbySurfaces);
  }
  else
  {
    myNearbySurfaces = myAllSurfaces;
    myNearbyTree = myAllTree;
  }
}

static double evaluate_levelset(
    const std::shared_ptr<CSGTreeNode>& node,
    const stk::math::Vector3d &x,
    const double narrowBand,
    const double farFieldValue,
    const bool isInsideNegative)
{
    STK_ThrowAssert(node != nullptr);

    if (node->node_operator() == CSGOperation::Leaf)
        return node->surface()->possibly_truncated_point_signed_distance(x, narrowBand, farFieldValue);

    switch (node->node_operator())
    {
        case CSGOperation::Union:
        {
            const auto & operands = node->vector_operands();
            STK_ThrowAssert(!operands.empty());
            double phi = evaluate_levelset(operands[0], x, narrowBand, farFieldValue, isInsideNegative);
            for (unsigned i=1; i<operands.size(); ++i)
            {
              double phiNode = evaluate_levelset(operands[i], x, narrowBand, farFieldValue, isInsideNegative);
              phi = isInsideNegative ? std::min(phi, phiNode) : std::max(phi, phiNode);
            }
            return phi;
        }

        case CSGOperation::Intersection:
        {
            const auto & operands = node->vector_operands();
            STK_ThrowAssert(!operands.empty());
            double phi = evaluate_levelset(operands[0], x, narrowBand, farFieldValue, isInsideNegative);
            for (unsigned i=1; i<operands.size(); ++i)
            {
              double phiNode = evaluate_levelset(operands[i], x, narrowBand, farFieldValue, isInsideNegative);
              phi = isInsideNegative ? std::max(phi, phiNode) : std::min(phi, phiNode);
            }
            return phi;
        }

        case CSGOperation::Difference:
        {
            const double phiLeft = evaluate_levelset(node->left_operand(), x, narrowBand, farFieldValue, isInsideNegative);
            const double phiRight = evaluate_levelset(node->right_operand(), x, narrowBand, farFieldValue, isInsideNegative);

            return isInsideNegative ? std::max(phiLeft, -phiRight) : std::min(phiLeft, -phiRight);
        }

        default:
        {
            throw std::runtime_error("Unknown CSG operation.");
        }
    }
}

double
CSG_Surface::possibly_truncated_point_signed_distance(const stk::math::Vector3d &x, const double narrowBandSize, const double farFieldValue) const
{
  if (myNearbyTree.empty())
    return farFieldValue;

  if (myIsInsideNegative)
  {
    STK_ThrowRequireMsg(farFieldValue >= narrowBandSize, "CSG_Surface has a specific requirement for farFieldValue due to min/max operations.");
    const double limitVal = (narrowBandSize == 0.) ? std::numeric_limits<double>::max() : farFieldValue;
    return std::min(limitVal, evaluate_levelset(myNearbyTree.requireRoot(), x, narrowBandSize, farFieldValue, myIsInsideNegative));
  }

  STK_ThrowRequireMsg(farFieldValue <= -narrowBandSize, "CSG_Surface has a specific requirement for farFieldValue due to min/max operations.");
  const double limitVal = (narrowBandSize == 0.) ? -std::numeric_limits<double>::max() : farFieldValue;
  return std::max(limitVal, evaluate_levelset(myNearbyTree.requireRoot(), x, narrowBandSize, farFieldValue, myIsInsideNegative));
}

void CSG_Surface::insert_into(BoundingBox & bbox) const
{
  for (auto && subsurf : myAllSurfaces)
    subsurf->insert_into(bbox);
}

bool CSG_Surface::does_intersect(const BoundingBox & bbox) const
{
  for (auto && subsurf : myAllSurfaces)
    if (subsurf->does_intersect(bbox))
      return true;
  return false;
}

std::string CSG_Surface::output() const
{
    return output_CSG_surface_tree(myAllTree);
}


} // namespace krino



