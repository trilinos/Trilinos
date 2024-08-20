// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_HELPERS_HPP
#define MINIEM_HELPERS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_String_Utilities.hpp"

#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_STK_IOClosureModel_Factory_TemplateBuilder.hpp"

#include "MiniEM_AddFieldsToMesh.hpp"

namespace mini_em {

  enum linearAlgebraType {
    linAlgTpetra,
    linAlgEpetra
  };

  enum physicsType {
    MAXWELL,
    DARCY
  };

  enum solverType {
    AUGMENTATION,
    MUELU,
    ML,
    CG,
    GMRES
  };

  void getMesh(Teuchos::ParameterList &mesh_pl,
               std::string &meshFile,
               int &x_elements,
               int &y_elements,
               int &z_elements,
               int &basis_order,
               Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
               Teuchos::RCP<panzer_stk::STK_Interface> &mesh,
               Teuchos::RCP<panzer_stk::STK_MeshFactory> &mesh_factory,
               double &mesh_size);

  Teuchos::RCP<Teuchos::ParameterList> getSolverParameters(linearAlgebraType linAlgebra,
                                                           physicsType physics,
                                                           solverType solver,
                                                           int dim,
                                                           Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
                                                           Teuchos::RCP<Teuchos::FancyOStream> &out,
                                                           std::string &xml,
                                                           int basis_order,
                                                           const bool preferTPLs = false,
                                                           const bool truncateMueLuHierarchy = false);

  void setClosureParameters(physicsType physics,
                            Teuchos::ParameterList &physicsEqSet,
                            Teuchos::ParameterList &closure_models,
                            double dt,
                            std::string &auxModelID);

  void setAuxiliaryOperatorParameters(physicsType physics,
                                      solverType solver,
                                      int basis_order,
                                      std::string pCoarsenScheduleStr,
                                      bool matrixFree,
                                      Teuchos::ParameterList &input_params,
                                      Teuchos::ParameterList &lin_solver_pl,
                                      std::string &auxFieldOrder);

  void createExodusFile(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                        Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory,
                        Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                        const bool & exodus_out,
                        Teuchos::RCP<const Teuchos::MpiComm<int> > comm,
                        physicsType physics);

  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
  buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
                            const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
                            const Teuchos::RCP<panzer::GlobalIndexer> & globalIndexer,
                            const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                            const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                            const Teuchos::ParameterList & closure_model_pl,
                            physicsType physics);
}


#endif
