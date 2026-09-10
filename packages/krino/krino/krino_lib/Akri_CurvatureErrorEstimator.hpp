// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CURVATUREERRORESTIMATOR_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CURVATUREERRORESTIMATOR_HPP_
#include <Akri_FieldRef.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

class InterfaceGeometry;
class RefinementSupport;

class CurvatureErrorEstimator
{
public:
  virtual void precompute(const stk::mesh::BulkData& mesh,
      const FieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const FieldRef nodeMarkerField) = 0;
  virtual double estimate_element_error(const stk::mesh::BulkData& mesh, const stk::mesh::Entity elem) const = 0;
  virtual double refine_tolerance() const = 0;
  virtual double unrefine_tolerance() const = 0;
  virtual ~CurvatureErrorEstimator() {}
};

class CurvatureTimesElementSizeError : public CurvatureErrorEstimator
{
public:
  CurvatureTimesElementSizeError(const InterfaceGeometry & interfaceGeometry, const double refinementTol, const stk::diag::Timer & parentTimer)
  : myInterfaceGeometry(interfaceGeometry), myRefinementTol(refinementTol), myTimer("Curvature Est", parentTimer)  {}
  void precompute(const stk::mesh::BulkData& mesh,
      const FieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const FieldRef nodeMarkerField) override {}
  double estimate_element_error(const stk::mesh::BulkData& mesh, const stk::mesh::Entity elem) const override;
  double refine_tolerance() const override { return myRefinementTol; }
  double unrefine_tolerance() const override { return 0.25*myRefinementTol; } // appropriate for error that converges linearly
private:
  const InterfaceGeometry & myInterfaceGeometry;
  double myRefinementTol;
  mutable stk::diag::Timer myTimer;
};

class FacetAngleCurvatureError : public CurvatureErrorEstimator
{
public:
  FacetAngleCurvatureError(const InterfaceGeometry & interfaceGeometry, const double refinementTol, const stk::diag::Timer & parentTimer)
  : myInterfaceGeometry(interfaceGeometry), myRefinementTol(refinementTol), myTimer("Curvature Angle Est", parentTimer)  {}
  void precompute(const stk::mesh::BulkData& mesh,
      const FieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const FieldRef nodeMarkerField) override;
  double estimate_element_error(const stk::mesh::BulkData& mesh, const stk::mesh::Entity elem) const override;
  double refine_tolerance() const override { return myRefinementTol; }
  double unrefine_tolerance() const override { return 0.25*myRefinementTol; } // appropriate for error that converges linearly
private:
  const InterfaceGeometry & myInterfaceGeometry;
  double myRefinementTol;
  std::vector<double> myErrorByElemOffset;
  mutable stk::diag::Timer myTimer;
};

class FacetChordalErrorFromGradientJump : public CurvatureErrorEstimator
{
// This class estimates the chordal error (distance from faceted surface to underlying curved surface).
// It uses the jump in the gradient across the sides of the element.  This jump can
// be used to closely approximate the angle between the facets of the surface.
public:
  FacetChordalErrorFromGradientJump(const InterfaceGeometry & interfaceGeometry, const double refinementTol, const stk::diag::Timer & parentTimer)
  : myInterfaceGeometry(interfaceGeometry), myRefinementTol(refinementTol), myTimer("Chordal Err Est", parentTimer)  {}
  void precompute(const stk::mesh::BulkData& mesh,
      const FieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const FieldRef nodeMarkerField) override;
  double estimate_element_error(const stk::mesh::BulkData& mesh, const stk::mesh::Entity elem) const override;
  double refine_tolerance() const override { return myRefinementTol; }
  double unrefine_tolerance() const override { return 0.0625*myRefinementTol; } // appropriate for error that converges quadratically
private:
  const InterfaceGeometry & myInterfaceGeometry;
  double myRefinementTol;
  std::vector<double> myErrorByElemOffset;
  mutable stk::diag::Timer myTimer;
};

std::unique_ptr<CurvatureErrorEstimator> build_curvature_error_estimator(const RefinementSupport & refinementSupport,
    const InterfaceGeometry & interfaceGeometry,
    const stk::diag::Timer & parentTimer);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CURVATUREERRORESTIMATOR_HPP_ */
