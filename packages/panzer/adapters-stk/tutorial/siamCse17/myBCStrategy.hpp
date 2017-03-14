// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef   __myBCStrategy_hpp__
#define   __myBCStrategy_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>
#include <vector>

// Panzer
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_Traits.hpp"

// Phalanx
#include "Phalanx_FieldManager.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

/**
 *  \brief Our zero Dirichlet boundary condition.
 *
 *  This class represents our zero Dirichlet boundary condition.  The name
 *  (MyBCStrategy) could have been anything we liked, and perhaps something
 *  more informative like ZeroDirichletBC would have been better.
 */
template <typename EvalT>
class MyBCStrategy
  :
  public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>
{
  public:

    /**
     *  \brief Default Constructor
     *
     *  Passes `bc` and `globalData` to the `BCStrategy_Dirichlet_DefaultImpl`
     *  constructor and ensures that the boundary condition "Strategy" from the
     *  input XML file matches this BC.
     *
     *  \param[in] bc         A `panzer::BC` constructed from an entry in the
     *                        "Boundary Conditions" `ParameterList` in the
     *                        input XML file.
     *  \param[in] globalData The global data for the problem, which holds the
     *                        parameter library and the default output stream.
     *
     *  \throws std::logic_error If "Strategy" isn't "MyBCStrategy".
     */
    MyBCStrategy(
      const panzer::BC&                       bc,
      const Teuchos::RCP<panzer::GlobalData>& globalData);

    /**
     *  \brief Set up this boundary condition.
     *
     *  Alert the Panzer library of a degree of freedom that is required by
     *  this boundary condition, and then find the basis corresponding to that
     *  degree of freedom.
     *
     *  \param[in] sidePB   The physics block to which this boundary condition
     *                      will be applied.
     *  \param[in] userData This is unused in this routine, though it is part
     *                      of the `BCStrategy_Dirichlet_DefaultImpl`
     *                      interface.
     *
     *  \throws std::runtime_error If we cannot find a basis for the degree of
     *                             freedom.
     */
    void
    setup(
      const panzer::PhysicsBlock&   sidePB,
      const Teuchos::ParameterList& userData);

    /**
     *  \brief Build and register evaluators.
     *
     *  Build any `Evaluator`s needed for this boundary condition (in this
     *  particular case, just a `panzer::Constant`), and register them with the
     *  `FieldManager`.
     *
     *  \param[in/out] fm       The object that holds all the fields that we
     *                          can use to build up the directed acyclic graph
     *                          for our problem.
     *  \param[in]     pb       This is unused in this routine, though it is
     *                          part of the `BCStrategy_Dirichlet_DefaultImpl`
     *                          interface.
     *  \param[in]     factory  This is unused in this routine, though it is
     *                          part of the `BCStrategy_Dirichlet_DefaultImpl`
     *                          interface.
     *  \param[in]     models   This is unused in this routine, though it is
     *                          part of the `BCStrategy_Dirichlet_DefaultImpl`
     *                          interface.
     *  \param[in]     userData This is unused in this routine, though it is
     *                          part of the `BCStrategy_Dirichlet_DefaultImpl`
     *                          interface.
     */
    void
    buildAndRegisterEvaluators(
      PHX::FieldManager<panzer::Traits>&                       fm,
      const panzer::PhysicsBlock&                              pb,
      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>&
                                                               factory,
      const Teuchos::ParameterList&                            models,
      const Teuchos::ParameterList&                            userData) const;

    /**
     *  \brief The basis for this boundary condition.
     *
     *  This is the basis corresponding to the degree of freedom to which this
     *  boundary condition applies.
     */
    Teuchos::RCP<panzer::PureBasis> basis;
}; // end of class MyBCStrategy

#include "myBCStrategyImpl.hpp"

#endif // __myBCStrategy_hpp__
