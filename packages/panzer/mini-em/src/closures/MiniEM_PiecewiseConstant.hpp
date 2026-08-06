// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_PIECEWISECONSTANT_DECL_HPP
#define MINIEM_PIECEWISECONSTANT_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace mini_em {

  using panzer::Cell;
  using panzer::Point;

  /** \brief Computes a piecewise-constant scalar field: value0 inside
    * the axis-aligned box [xl,xr]x[yl,yr][x[zl,zr]], and value1
    * outside it.
    */
  template<typename EvalT, typename Traits>
  class PiecewiseConstant : public panzer::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    /** \brief Constructor.
      *
      * \param[in] name name of the field to output.
      * \param[in] ir integration rule to evaluate at.
      * \param[in] fl field layout library used to look up the field's data layout.
      * \param[in] value0 value inside the box [xl,xr]x[yl,yr][x[zl,zr]].
      * \param[in] value1 value outside the box.
      * \param[in] xl lower x bound of the box.
      * \param[in] xr upper x bound of the box.
      * \param[in] yl lower y bound of the box.
      * \param[in] yr upper y bound of the box.
      * \param[in] zl lower z bound of the box (3D only).
      * \param[in] zr upper z bound of the box (3D only).
      * \param[in] DoF_ name of the DOF whose basis supplies the coordinates to evaluate at.
      */
    PiecewiseConstant(const std::string & name,
                      const panzer::IntegrationRule & ir,
                      const panzer::FieldLayoutLibrary & fl,
                      const double value0,
                      const double value1,
                      const double xl,
                      const double xr,
                      const double yl,
                      const double yr,
                      const double zl,
                      const double zr,
                      const std::string& DoF_);

    /// \brief Computes the piecewise-constant field for this workset per the class description.
    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT,Cell,Point> values;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree, ir_dim;
    double value0_, value1_;
    double xl_, xr_, yl_, yr_, zl_, zr_;
  };

}

#include "MiniEM_PiecewiseConstant_impl.hpp"

#endif
