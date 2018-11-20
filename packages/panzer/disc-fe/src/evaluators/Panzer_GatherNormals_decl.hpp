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

#ifndef PANZER_EVALUATOR_GATHER_NORMALS_DECL_HPP
#define PANZER_EVALUATOR_GATHER_NORMALS_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_PointValues2.hpp"

namespace panzer {

/** \brief Gathers tangent vectors per field from the global indexer and
    stores them in the field manager.
*/
template<typename EvalT, typename Traits> 
class GatherNormals
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>,
    public CloneableEvaluator  {
   
public:
  
  GatherNormals(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  {return Teuchos::rcp(new GatherNormals<EvalT,Traits>(pl));}
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  std::string dof_name;

  PHX::MDField<ScalarT,Cell,NODE,Dim> gatherFieldNormals;

  GatherNormals();

  Teuchos::RCP<const PureBasis> basis;
  Teuchos::RCP<const PointRule> pointRule;
  Kokkos::DynRankView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> faceNormal; // face normals
  Kokkos::DynRankView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> refFaceNormal; // reference face normals

  PointValues2<double> pointValues;
  PHX::MDField<const double, Cell, IP, Dim, Dim, void, void, void, void>
    constJac_;

  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations;
};

}

// **************************************************************
#endif
