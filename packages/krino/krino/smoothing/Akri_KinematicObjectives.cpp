// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_KinematicObjectives.hpp>

#include <vector>

#include <Akri_KinematicElement.hpp>
#include <Akri_KinematicUtils.hpp>
#include <Akri_ModifiedNeoHookeanModel.hpp>
#include <Akri_SethHillModel.hpp>
#include <Akri_SizeShapeModel.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino
{

namespace
{

// The kinematic element objective of a model: the element energy and its sensitivity to the node
// coordinates, for tri3 elements in 2D and tet4 elements in 3D.
template <class MODEL> class KinematicModelObjective : public KinematicElementObjective
{
public:
  KinematicModelObjective() = default;
  explicit KinematicModelObjective(const MODEL & model) : myModel(model) {}

  virtual double compute_element_objective(const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords) const override
  {
    STK_ThrowAssert((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4));
    if (spatialDim == 2) return tri3_2d_element_energy(myModel, refSize, elemNodeCoords);
    return tet4_element_energy(myModel, refSize, elemNodeCoords);
  }

  virtual void fill_element_sensitivity(const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const override
  {
    STK_ThrowAssert((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4));
    if (spatialDim == 2)
      kinematicUtils::fill_element_sensitivity(
          tri3_2d_element_forces(myModel, refSize, elemNodeCoords), elemGradContrib);
    else
      kinematicUtils::fill_element_sensitivity(
          tet4_element_forces(myModel, refSize, elemNodeCoords), elemGradContrib);
  }

  virtual double compute_element_objective(const unsigned spatialDim,
      const unsigned npe,
      const double * metricMapping,
      const double * const * elemNodeCoords) const override
  {
    STK_ThrowRequireMsg((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4),
        "Anisotropic kinematic objectives support only tri3 elements in 2D and tet4 elements in 3D "
        "(tri6 and tet10 elements are given as their vertices).");
    if (spatialDim == 2) return tri3_2d_element_energy(myModel, metricMapping, elemNodeCoords);
    return tet4_element_energy(myModel, metricMapping, elemNodeCoords);
  }

  virtual void fill_element_sensitivity(const unsigned spatialDim,
      const unsigned npe,
      const double * metricMapping,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const override
  {
    STK_ThrowRequireMsg((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4),
        "Anisotropic kinematic objectives support only tri3 elements in 2D and tet4 elements in 3D "
        "(tri6 and tet10 elements are given as their vertices).");
    if (spatialDim == 2)
      kinematicUtils::fill_element_sensitivity(
          tri3_2d_element_forces(myModel, metricMapping, elemNodeCoords), elemGradContrib);
    else
      kinematicUtils::fill_element_sensitivity(
          tet4_element_forces(myModel, metricMapping, elemNodeCoords), elemGradContrib);
  }

private:
  MODEL myModel;
};

} // namespace

std::shared_ptr<KinematicElementObjective> make_size_shape_objective()
{
  return std::make_shared<KinematicModelObjective<SizeShapeModel>>();
}

std::shared_ptr<KinematicElementObjective> make_modified_neo_hookean_objective()
{
  return std::make_shared<KinematicModelObjective<ModifiedNeoHookeanModel>>();
}

std::shared_ptr<KinematicElementObjective> make_seth_hill_objective(const SethHillParams & params)
{
  return std::make_shared<KinematicModelObjective<SethHillModel>>(SethHillModel(params));
}

std::shared_ptr<KinematicElementObjective> make_simplified_seth_hill_objective()
{
  return std::make_shared<KinematicModelObjective<SimplifiedSethHillModel>>();
}

} // namespace krino
