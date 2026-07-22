// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __myEquationSet_hpp__
#define   __myEquationSet_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>
#include <vector>

// Panzer
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"

// Phalanx
#include "Phalanx_FieldManager.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

/**
 * 	\brief Our equation set.
 *
 * 	This class represents our equation set, which contains only a single
 * 	equation:
 * 	\f[
 * 	  -\Delta u(x, y) + (1 - 8 \pi^2) u(x, y) = \sin(2 \pi x) \sin(2 \pi y).
 * 	\f]
 */
template <typename EvalT>
class MyEquationSet
  :
  public panzer::EquationSet_DefaultImpl<EvalT>
{
  public:
    
    /**
     * 	\brief Default Constructor.
     *
     * 	This first validates the parameters that are passed in (params), and
     * 	then creates the degree of freedom ("U") and its gradient, and finally
     * 	adds the ID of the closure model that will be used to supply the source
     * 	term.
     *
     * 	\param[in] params                  An unnamed `ParameterList` that is
     * 	                                   an entry of our physics block
     * 	                                   `ParameterList` in the input XML
     * 	                                   file.
     * 	\param[in] defaultIntegrationOrder The integration order to use if none
     * 	                                   is supplied in the input XML file.
     * 	\param[in] cellData                This is passed to the
     * 	                                   `EquationSet_DefaultImpl`
     * 	                                   constructor.
     * 	\param[in] gd                      This is passed to the
     * 	                                   `EquationSet_DefaultImpl`
     * 	                                   constructor.
     * 	\param[in] buildTransientSupport   This is passed to the
     * 	                                   `EquationSet_DefaultImpl`
     * 	                                   constructor.
     */
    MyEquationSet(
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      const int&                                  defaultIntegrationOrder,
      const panzer::CellData&                     cellData,
      const Teuchos::RCP<panzer::GlobalData>&     gd,
      const bool                                  buildTransientSupport);

    /**
     * 	\brief Build and register the equation set evaluators.
     *
     * 	This routine builds up the directed acyclic graph bit by bit.  We'll
     * 	add an `Integrator_GradBasisDotVector` for the Laplacian term, and then
     * 	two `Integrator_BasisTimesScalar`s for the Helmholtz and source terms.
     * 	These three integrals will then be summed up to give us the residual
     * 	representation of our equation.
     *
     * 	\param[in/out] fm           The object that stores all the fields we
     * 	                            have available to build up the directed
     * 	                            acylcic graph, to which we'll add all the
     * 	                            evaluators that make up this equation set.
     * 	\param[in]     fieldLibrary This is unused in this routine, though it
     * 	                            is part of the `EquationSet_DefaultImpl`
     * 	                            interface.
     * 	\param[in]     userData     This is unused in this routine, though it
     * 	                            is part of the `EquationSet_DefaultImpl`
     * 	                            interface.
     */
    void
    buildAndRegisterEquationSetEvaluators(
      PHX::FieldManager<panzer::Traits>& fm,
      const panzer::FieldLibrary&        fieldLibrary,
      const Teuchos::ParameterList&      userData) const;

  private:
    
    /**
     * 	\brief The name of the degree of freedom.
     */
    std::string dofName_;

}; // end of class MyEquationSet

#include "myEquationSetImpl.hpp"

#endif // __myEquationSet_hpp__
