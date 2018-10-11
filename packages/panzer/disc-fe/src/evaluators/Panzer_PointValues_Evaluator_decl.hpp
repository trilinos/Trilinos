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

#ifndef PANZER_POINT_VALUES_EVALUATOR_DECL_HPP
#define PANZER_POINT_VALUES_EVALUATOR_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
template<typename EvalT, typename Traits>
class PointValues_Evaluator
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    PointValues_Evaluator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // is anything other than ScalarT really needed here?
  PointValues2<double> pointValues;
 
  PHX::MDField<double,NODE,Dim> refPointArray;

  bool useBasisValuesRefArray; // if true then basis is non-null
  Teuchos::RCP<const panzer::PureBasis> basis;
  std::size_t basis_index;

  //! Initialization method to unify the constructors.
  template <typename ArrayT>
  void initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                  const Teuchos::Ptr<const ArrayT> & userArray,
                  // const Teuchos::Ptr<const Kokkos::DynRankView<double,PHX::Device> > & userArray,
                  const Teuchos::RCP<const panzer::PureBasis> & pureBasis);

public:
  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Kokkos::DynRankView<double,PHX::Device> & userArray);

  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const PHX::MDField<double, panzer::IP, panzer::Dim> & userArray);

  //! This builds a point rule from the basis function reference points in the workset
  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Teuchos::RCP<const panzer::PureBasis> & pureBasis);

}; // end of class PointValues_Evaluator


}

#endif
