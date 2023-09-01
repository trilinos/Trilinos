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

#ifndef PANZER_STK_PROJECT_FIELD_DECL_HPP
#define PANZER_STK_PROJECT_FIELD_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer_stk {

/** \brief Given a field, perform a local L2 projection onto the desired basis.
 * 
 * \note Requires that orientations be given in the \c postRegistrationSetup phase.
*/
template<typename EvalT, typename Traits> 
class ProjectField
  : public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
    {
   
public:

  /**
   * Constructor for the ProjectField evaluator. 
   * 
   * \param[in] inName Name of the source MDField
   * \param[in] src Basis of the source field
   * \param[in] dst Target basis
   * \param[in] outname (Optional) Name for the projected MDField
   */
  
  ProjectField(const std::string & inName, Teuchos::RCP<panzer::PureBasis> src,
               Teuchos::RCP<panzer::PureBasis> dst, std::string outname = "");
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& fm);
  
  void evaluateFields(typename Traits::EvalData d);

private:

  using ScalarT = typename EvalT::ScalarT;

  const std::string field_name_;
  Teuchos::RCP<const panzer::PureBasis> srcBasis_;
  Teuchos::RCP<const panzer::PureBasis> dstBasis_;
  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations_;

  PHX::MDField<double,panzer::Cell,panzer::BASIS> result_;
  PHX::MDField<double,panzer::Cell,panzer::BASIS> source_;

  Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> local_orts_;

};

}

// **************************************************************
#endif
