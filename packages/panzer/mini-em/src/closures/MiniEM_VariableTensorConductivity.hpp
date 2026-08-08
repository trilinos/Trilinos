// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_VARIABLETENSORCONDUCTIVITY_DECL_HPP
#define MINIEM_VARIABLETENSORCONDUCTIVITY_DECL_HPP

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

  /** \brief Computes a piecewise-constant magnetized (Hall-effect)
    * conductivity tensor field, using one of three (sigma,beta)
    * parameter sets depending on which of two nested, hardcoded
    * axis-aligned boxes centered at the origin the point falls in:
    * (sigma0,beta0) inside [0.4,0.6]^dim, (sigma1,beta1) inside
    * [0.2,0.8]^dim but outside the first box, and (sigma2,beta2)
    * outside both boxes. Each tensor has the form
    * sigma*(I + beta⊗beta - [beta]_x) (see
    * mini_em::TensorConductivity for the constant single-region
    * version).
    */
  template<typename EvalT, typename Traits>
  class VariableTensorConductivity : public panzer::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    /** \brief Constructor.
      *
      * \param[in] name name of the field to output.
      * \param[in] ir integration rule to evaluate at.
      * \param[in] fl field layout library used to look up the field's data layout.
      * \param[in] sigma0_ conductivity of the inner box region [0.4,0.6]^dim.
      * \param[in] sigma1_ conductivity of the middle region, inside [0.2,0.8]^dim but outside the inner box.
      * \param[in] sigma2_ conductivity outside both boxes.
      * \param[in] betax0_ x-component of the inner region's Hall parameter vector.
      * \param[in] betay0_ y-component of the inner region's Hall parameter vector.
      * \param[in] betaz0_ z-component of the inner region's Hall parameter vector.
      * \param[in] betax1_ x-component of the middle region's Hall parameter vector.
      * \param[in] betay1_ y-component of the middle region's Hall parameter vector.
      * \param[in] betaz1_ z-component of the middle region's Hall parameter vector.
      * \param[in] betax2_ x-component of the outer region's Hall parameter vector.
      * \param[in] betay2_ y-component of the outer region's Hall parameter vector.
      * \param[in] betaz2_ z-component of the outer region's Hall parameter vector.
      * \param[in] DoF_ name of the DOF whose basis supplies the coordinates to evaluate at.
      */
    VariableTensorConductivity(const std::string & name,
                               const panzer::IntegrationRule & ir,
                               const panzer::FieldLayoutLibrary & fl,
                               const double & sigma0_,
                               const double & sigma1_,
                               const double & sigma2_,
                               const double & betax0_,
                               const double & betay0_,
                               const double & betaz0_,
                               const double & betax1_,
                               const double & betay1_,
                               const double & betaz1_,
                               const double & betax2_,
                               const double & betay2_,
                               const double & betaz2_,
                               const std::string& DoF_);

    /// \brief Looks up the integration rule index needed by evaluateFields(). Called once before the first evaluation.
    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    /// \brief Computes the conductivity tensor field for this workset per the class description.
    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT,Cell,Point,Dim,Dim> conductivity;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree, ir_dim, ir_index;
    double sigma0, betax0, betay0, betaz0;
    double sigma1, betax1, betay1, betaz1;
    double sigma2, betax2, betay2, betaz2;
  };

}

#include "MiniEM_VariableTensorConductivity_impl.hpp"

#endif
