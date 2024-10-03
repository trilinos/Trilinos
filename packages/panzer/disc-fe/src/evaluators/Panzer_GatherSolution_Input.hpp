// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_GatherSolution_Input_hpp__
#define __Panzer_GatherSolution_Input_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

namespace panzer {

// forward declaration
class PureBasis;

/** This class parses the input parameters for the gather solution 
  * evaluators. Its goal is to unify the input for all of those 
  * evaluators. 
  */
class GatherSolution_Input : Teuchos::ParameterListAcceptorDefaultBase {
public:
  GatherSolution_Input();

  /** Set the parameter list, this is the complete state. This will modify 
    * the list.
    */
  void setParameterList(const Teuchos::ParameterList & pl);

  /** Set the parameter list, this is the complete state. This will modify 
    * the list.
    */
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
  
  //! Get valid parameters
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  // Accessor functions

  //! The names of the DOFs to be gathered (all types)
  const std::vector<std::string> & getDofNames() const { return  dofNames_; }

  /** The indexer names used to gather the DOFs. Will be same length as getDofNames
    * (all types).
    */
  const std::vector<std::string> & getIndexerNames() const { return indexerNames_; }

  //! Basis definiting layout of dof names (all types)
  Teuchos::RCP<const PureBasis> getBasis() const { return basis_; }

  //! Gather a time derivative vector?  (all types)
  bool useTimeDerivativeSolutionVector() const { return useTimeDerivSolnVec_; }
  
  //! Name of the global evaluation data container to use for the source vector (all types)
  std::string getGlobalDataKey() const { return globalDataKey_; }
  
  // tangent
  
  //! Get the name of the tangent fields (tangent only)
  const std::vector<std::vector<std::string> > & getTangentNames() const { return tangentNames_; }

  // jacobian
  
  //! The name of the sensitivities. Enables sensitivities at "preEvaluate" time (Jacobian and Hessian)
  std::string getSensitivitiesName() const { return sensName_; }

  //! What index to use for initializing the seed (Jacobian and Hessian)
  int getGatherSeedIndex() const { return gatherSeedIndex_; }

  //! Are first derivative sensitivities enabled or disabled? (Jacobian and Hessian)
  bool firstSensitivitiesAvailable() { return firstSensAvail_; }

  // hessian
  
  //! Are second derivative sensitivies enabled or disabled (Hessian only)
  bool secondSensitivitiesAvailable() { return secondSensAvail_; }

  //! What prefix to use for the GEDC for second derivative sensitivity direction (Hessian only)
  std::string getSecondSensitivityDataKeyPrefix() { return secondSensDataKeyPrefix_; }

private:
  GatherSolution_Input(const GatherSolution_Input &); // hide me

  // residual
  std::vector<std::string> dofNames_;   
  std::vector<std::string> indexerNames_;   
  Teuchos::RCP<const PureBasis> basis_;
  bool useTimeDerivSolnVec_;
  std::string globalDataKey_;
  
  // tangent
  std::vector<std::vector<std::string> > tangentNames_;   

  // jacobian
  std::string sensName_;
  int gatherSeedIndex_;
  bool firstSensAvail_;

  // hessian
  bool secondSensAvail_;
  std::string secondSensDataKeyPrefix_;
};

}

#endif
