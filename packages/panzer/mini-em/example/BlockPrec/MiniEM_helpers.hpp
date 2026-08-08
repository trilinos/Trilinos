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

  //! Which linear algebra backend to build the linear system with.
  enum linearAlgebraType {
    linAlgTpetra,
    linAlgEpetra
  };

  //! Which physics (equation set) the mini-em example is solving.
  enum physicsType {
    MAXWELL,
    DARCY
  };

  //! Which preconditioning/solver strategy to use.
  enum solverType {
    AUGMENTATION,
    MUELU,
    ML,
    CG,
    GMRES,
    MAXWELL1_RS,
    MAXWELL1_SA_RS,
    MAXWELL1_EMIN,
    DIRECT
  };

  /** \brief Builds the mesh and mesh factory from a mesh ParameterList, either by reading an Exodus file or building an inline mesh (Quad/Tri for 2D, Hex/Tet for 3D, per meshType).
    *
    * \param[in] mesh_pl mesh ParameterList; its "Source" entry selects Exodus/Pamgen/inline mesh generation.
    * \param[in] meshFile if non-empty, an Exodus file to read the mesh from, overriding mesh_pl's "Source".
    * \param[in,out] x_elements number of elements in x for an inline mesh; if positive, overrides mesh_pl, otherwise set on return to the actual value used.
    * \param[in,out] y_elements number of elements in y for an inline mesh; if positive, overrides mesh_pl, otherwise set on return to the actual value used.
    * \param[in,out] z_elements number of elements in z for an inline mesh (3D only); if positive, overrides mesh_pl, otherwise set on return to the actual value used.
    * \param[in] meshType if non-empty, the inline mesh type ("tet"/"hex"/"tri"/"quad"), overriding mesh_pl's "Mesh Type".
    * \param[in] basis_order unused by this function.
    * \param[in] comm MPI communicator to build the mesh on.
    * \param[out] mesh set to the built (uncommitted) mesh.
    * \param[out] mesh_factory set to the mesh factory used to build mesh.
    * \param[out] mesh_size set to the inline mesh's characteristic element size (1/max element count per axis); unset for Exodus/Pamgen meshes.
    */
  void getMesh(Teuchos::ParameterList &mesh_pl,
               std::string &meshFile,
               int &x_elements,
               int &y_elements,
               int &z_elements,
               std::string meshType,
               int &basis_order,
               Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
               Teuchos::RCP<panzer_stk::STK_Interface> &mesh,
               Teuchos::RCP<panzer_stk::STK_MeshFactory> &mesh_factory,
               double &mesh_size);

  /** \brief Builds the Stratimikos/Teko solver ParameterList for the given linear algebra backend, physics, and solver strategy, optionally reading overrides from an xml file.
    *
    * \param[in] linAlgebra which linear algebra backend (Tpetra/Epetra) the solver is configured for.
    * \param[in] physics which physics (Maxwell/Darcy) the solver is configured for.
    * \param[in] solver which solver/preconditioner strategy to load the base configuration for.
    * \param[in] dim spatial dimension of the problem; selects the 2D solver configuration override when 2.
    * \param[in] comm MPI communicator, passed through to updateParams() when loading each xml file.
    * \param[in] out output stream, passed through to updateParams() when loading each xml file.
    * \param[in] xml if non-empty, an xml file to load instead of the solver's default configuration.
    * \param[in] basis_order FE basis order; if greater than 1, layers on the higher-order MueLu configuration.
    * \param[in] preferTPLs if true, layers on the "prefer third-party libraries" MueLu configuration override.
    * \param[in] useBarriers if true, layers on the MueLu configuration override that adds MPI barriers for profiling.
    * \param[in] truncateMueLuHierarchy if true, layers on the MueLu configuration override that truncates the multigrid hierarchy.
    */
  Teuchos::RCP<Teuchos::ParameterList> getSolverParameters(linearAlgebraType linAlgebra,
                                                           physicsType physics,
                                                           solverType solver,
                                                           int dim,
                                                           Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
                                                           Teuchos::RCP<Teuchos::FancyOStream> &out,
                                                           std::string &xml,
                                                           int basis_order,
                                                           const bool preferTPLs = false,
                                                           const bool useBarriers = false,
                                                           const bool truncateMueLuHierarchy = false);

  /** \brief Derives and fills in additional closure model parameters (e.g. inverse permeability/diffusivity, 1/dt) needed by the physics equation set and copies the relevant material parameters into a separate closure model list (auxModelID) used by the auxiliary equation sets.
    *
    * \param[in] physics which physics (Maxwell/Darcy) is being solved; selects which parameters are derived/copied.
    * \param[in] physicsEqSet the physics equation set's ParameterList; may be updated with a derived "Inverse Permeability" entry (Maxwell only).
    * \param[in,out] closure_models the closure model ParameterList; filled in with derived entries (e.g. "1/dt", "1/mu") and the auxiliary model sublist.
    * \param[in] dt the time step, used to derive the "dt"/"1/dt" closure model entries.
    * \param[out] auxModelID set to the name of the auxiliary closure model sublist ("electromagnetics_aux" or "darcy_aux") that was filled in.
    */
  void setClosureParameters(physicsType physics,
                            Teuchos::ParameterList &physicsEqSet,
                            Teuchos::ParameterList &closure_models,
                            double dt,
                            std::string &auxModelID);

  /** \brief Fills in the auxiliary equation set and preconditioner ParameterLists needed to build the requested solver's auxiliary operators (mass matrices, Schur complements, weak gradient, etc.), given the basis order, p-multigrid coarsening schedule, and matrix-free option.
    *
    * \param[in] physics which physics (Maxwell/Darcy) is being solved; selects which auxiliary model ID and operators are set up.
    * \param[in] solver which solver/preconditioner strategy is being configured; selects which auxiliary operators are needed.
    * \param[in] basis_order FE basis order; prepended to the p-multigrid coarsening schedule if not already present.
    * \param[in] pCoarsenScheduleStr comma-separated p-multigrid coarsening schedule, used as a fallback if "Assembly"/"p coarsen schedule" is not set in input_params.
    * \param[in] matrixFree if true, configures the auxiliary operators to be built matrix-free instead of assembled.
    * \param[in] input_params the top-level input ParameterList; its "Physics Blocks", "Assembly", "Auxiliary Operators", and "Closure Models" sublists are read and updated with the auxiliary equation set configuration.
    * \param[in,out] lin_solver_pl the solver ParameterList; updated with the p-multigrid coarsening schedule where required.
    * \param[out] auxFieldOrder set to the blocked field order string for the auxiliary DOF manager.
    */
  void setAuxiliaryOperatorParameters(physicsType physics,
                                      solverType solver,
                                      int basis_order,
                                      std::string pCoarsenScheduleStr,
                                      bool matrixFree,
                                      Teuchos::ParameterList &input_params,
                                      Teuchos::ParameterList &lin_solver_pl,
                                      std::string &auxFieldOrder);

  /** \brief Registers the physics blocks' DOF fields with the mesh for output, completes mesh construction, and (if exodus_out is true) sets up the mesh for writing to "mesh_output.exo".
    *
    * \param[in] physicsBlocks the physics blocks whose DOF fields should be registered with the mesh for output.
    * \param[in] mesh_factory mesh factory used to complete mesh construction.
    * \param[in] mesh the mesh to register fields with and complete construction of.
    * \param[in] exodus_out if true, sets up the mesh for writing to "mesh_output.exo".
    * \param[in] comm MPI communicator to complete mesh construction on.
    * \param[in] physics which physics (Maxwell/Darcy) is being solved; Maxwell additionally registers a "J" cell-average vector output.
    */
  void createExodusFile(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                        Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory,
                        Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                        const bool & exodus_out,
                        Teuchos::RCP<const Teuchos::MpiComm<int> > comm,
                        physicsType physics);

  /** \brief Builds and registers the evaluators for a response library that writes the solution and requested closure-model fields out to the STK mesh.
    *
    * \param[in] physicsBlocks the physics blocks to build the response library over.
    * \param[in] linObjFactory linear object factory used to construct the response library.
    * \param[in] wkstContainer workset container the response library evaluates over.
    * \param[in] globalIndexer global indexer describing the DOF layout.
    * \param[in] cm_factory closure model factory used to build the closure-model fields to be written out.
    * \param[in] mesh the STK mesh to write the solution/closure-model fields to.
    * \param[in] closure_model_pl the closure model ParameterList to build fields from.
    * \param[in] physics which physics (Maxwell/Darcy) is being solved.
    */
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
