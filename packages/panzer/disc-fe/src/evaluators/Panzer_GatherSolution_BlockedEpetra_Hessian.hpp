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

#ifndef __Panzer_GatherSolution_BlockedEpetra_Hessian_hpp__
#define __Panzer_GatherSolution_BlockedEpetra_Hessian_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// Epetra gather solution file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class GatherSolution_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers)
     : indexers_(indexers) {}

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                               const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>(indexers_,pl)); }

private:
  typedef typename panzer::Traits::Hessian EvalT;
  typedef typename panzer::Traits::Hessian::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  std::vector<std::string> indexerNames_;

  std::string globalDataKey_; // what global data does this fill?

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > indexers_;
  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  bool useTimeDerivativeSolutionVector_;

  std::string sensitivitiesName_; // This sets which gather operations have sensitivities

  // first derivative fields
  int gatherSeedIndex_;              // what gather seed in the workset to use
                                     // if less than zero then use alpha or beta
                                     // as appropriate
  bool firstSensitivitiesAvailable_; // Turn on the first derivative sensitivities 
                                     // to turn on/off a certain set of sensitivities
  bool firstApplySensitivities_;     // This is a local variable that is used by evaluateFields
                        
  // handle second derivatives                         
  std::string sensitivities2ndPrefix_; // Prefix for field containing the sensitivities
  bool secondSensitivitiesAvailable_;  // Turn on the second derivative sensitivities 
  bool secondApplySensitivities_;      // This is a local variable that is used by evaluateFields

  Teuchos::RCP<Thyra::ProductVectorBase<double> > x_;
  Teuchos::RCP<Thyra::ProductVectorBase<double> > dx_;

  GatherSolution_BlockedEpetra();
};

}

#endif // end hessian support

#endif
