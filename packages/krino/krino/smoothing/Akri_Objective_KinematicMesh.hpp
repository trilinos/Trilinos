// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_KINEMATICMESH_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_KINEMATICMESH_HPP_
#include <Akri_Objective_Mesh.hpp>
#include <Akri_SethHillConfig.hpp>

namespace krino {

class KinematicElementObjective
{
public:
  virtual ~KinematicElementObjective() {}
  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords) const = 0;

  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const = 0;
};

class KinematicMeshObjective : public MeshObjectiveBase
{
public:
  KinematicMeshObjective(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const KinematicElementObjective & elemSmoothingObj)
  : MeshObjectiveBase(mesh, coordsField, elemSelector), myElemObj(elemSmoothingObj) {}

  virtual double compute_element_objective(
          const unsigned spatialDim,
          const unsigned npe,
          const size_t elemIndex,
          const double * const * elemNodeCoords) const override;
  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const size_t elemIndex,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const override;

  void set_constant_element_ref_size(const double);
  void set_element_ref_size(const CoordinatesFieldRef coordsField, const SmoothReference & reference);
  void set_element_ref_sizes(const std::vector<double> & elemRefSizes);

private:
  const KinematicElementObjective & myElemObj;
  std::vector<double> myElemRefSizes;
};

}


#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_KINEMATICMESH_HPP_ */
