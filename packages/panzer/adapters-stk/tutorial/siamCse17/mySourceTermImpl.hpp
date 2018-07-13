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

#ifndef   __mySourceTermImpl_hpp__
#define   __mySourceTermImpl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <cmath>

// Kokkos
#include "Kokkos_Parallel.hpp"

// Panzer
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
MySourceTerm<EvalT, Traits>::
MySourceTerm(
  const std::string&             name,
  const panzer::IntegrationRule& ir)
  :
  irDegree_(ir.cubature_degree)
{
  using panzer::Cell;
  using panzer::Point;
  using PHX::MDField;

  // This is the field that will be storing the values of our source function
  // at every cell and integration point.
  result_ = MDField<ScalarT, Cell, Point>(name, ir.dl_scalar);

  // This specifies the field that this Evaluator will be evaluating.  Since
  // our source term is an analytic function of x and y, it has no
  // dependencies, so we have no need to call addDependentField().
  this->addEvaluatedField(result_);

  // Multiple instances of the same Evaluator can be created with different 
  // "name"s passed into the constructor.  This is admittedly not useful, in
  // our particular case, but could easily be useful if we were to generalize
  // our source term to take multiple arguments, for instance, k_x and k_y,
  // instead of having those both hard-coded to 2\pi.
  this->setName("MySourceTerm(" + name + ")");
} // end of Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
MySourceTerm<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  using panzer::getIntegrationRuleIndex;
  irIndex_ = getIntegrationRuleIndex(irDegree_, (*sd.worksets_)[0]);
} // end of postRegistrationSetup()

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
MySourceTerm<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using Kokkos::parallel_for;
  using panzer::index_t;

  const auto& coords = workset.int_rules[irIndex_]->ip_coordinates;

  // Loop over the cells in the workset.
  parallel_for(workset.num_cells, [=] (const index_t cell)
  {
    // Loop over the integration points in the cell.
    for (int point(0); point < result_.extent_int(1); ++point)
    {
      // Get the (x,y) location corresponding to this integration point in this
      // cell.
      const double& x(coords(cell, point, 0));
      const double& y(coords(cell, point, 1));

      // Compute the value of the source term at this (x,y) location.
      result_(cell, point) = sin(2 * M_PI * x) * sin(2 * M_PI * y);
    } // end loop over the integration points
  }); // end loop over the cells in the workset
} // end of evaluateFields()

#endif // __mySourceTermImpl_hpp__
