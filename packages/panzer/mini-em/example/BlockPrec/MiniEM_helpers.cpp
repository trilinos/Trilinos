#include "MiniEM_helpers.hpp"

namespace mini_em {

  void getMesh(Teuchos::ParameterList &mesh_pl,
               std::string &meshFile,
               int &x_elements,
               int &y_elements,
               int &z_elements,
               int &basis_order,
               Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
               Teuchos::RCP<panzer_stk::STK_Interface> &mesh,
               Teuchos::RCP<panzer_stk::STK_MeshFactory> &mesh_factory,
               double &mesh_size) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Teuchos::TimeMonitor tMmesh(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: build mesh")));

    int dim = 3;

    if (mesh_pl.get<std::string>("Source") == "Exodus File" || meshFile != "") { // Exodus file reader...
      RCP<Teuchos::ParameterList> pl;
      RCP<Teuchos::ParameterList> input_pl = rcp(new Teuchos::ParameterList(mesh_pl.sublist("Exodus File")));
      if (meshFile == "") {
        pl = rcp(new Teuchos::ParameterList(input_pl->sublist("Exodus Parameters")));
      } else {
        pl = rcp(new Teuchos::ParameterList());
        pl->set("File Name",meshFile);
      }
      mesh_factory = Teuchos::RCP<panzer_stk::STK_MeshFactory>(new panzer_stk::STK_ExodusReaderFactory());
      mesh_factory->setParameterList(pl);
      // build mesh
      mesh = mesh_factory->buildUncommitedMesh((*comm->getRawMpiComm())());
    } else if (mesh_pl.get<std::string>("Source") ==  "Pamgen Mesh") { // Pamgen mesh generator
      Teuchos::ParameterList & pamgen_pl = mesh_pl.sublist("Pamgen Mesh");
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList(pamgen_pl.sublist("Pamgen Parameters")));
      pl->set("File Type","Pamgen");
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory());
      mesh_factory->setParameterList(pl);
      // build mesh
      mesh = mesh_factory->buildUncommitedMesh((*comm->getRawMpiComm())());
    } else if (mesh_pl.get<std::string>("Source") == "Inline Mesh") { // Inline mesh generator
      // set mesh factory parameters
      Teuchos::ParameterList & inline_gen_pl = mesh_pl.sublist("Inline Mesh");
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList(inline_gen_pl.sublist("Mesh Factory Parameter List")));
      dim = inline_gen_pl.get<int>("Mesh Dimension");

      // overrides from command line
      if (x_elements > 0)
        pl->set<int>("X Elements",x_elements);
      if (y_elements > 0)
        pl->set<int>("Y Elements",y_elements);
      if (dim == 3 && z_elements > 0)
        pl->set<int>("Z Elements",z_elements);

      // build mesh
      if (dim == 3) {
        if (inline_gen_pl.get<std::string>("Mesh Type") == "tet")
          mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory());
        else if (inline_gen_pl.get<std::string>("Mesh Type") == "quad")
          mesh_factory = rcp(new panzer_stk::CubeHexMeshFactory());
        else
          throw;
      } else if (dim == 2) {
        if (inline_gen_pl.get<std::string>("Mesh Type") == "tet")
          mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory());
        else if (inline_gen_pl.get<std::string>("Mesh Type") == "quad")
          mesh_factory = rcp(new panzer_stk::SquareQuadMeshFactory());
        else
          throw;
      }
      mesh_factory->setParameterList(pl);
      mesh = mesh_factory->buildUncommitedMesh((*comm->getRawMpiComm())());

      x_elements = pl->get<int>("X Elements");
      y_elements = pl->get<int>("Y Elements");
      if (dim == 3) {
        z_elements = pl->get<int>("Z Elements");
        mesh_size = 1.0/std::max(x_elements,std::max(y_elements,z_elements));
      } else
        mesh_size = 1.0/std::max(x_elements,y_elements);
    } else
      throw;
  }


  void updateParams(const std::string & xml,
                    Teuchos::RCP<Teuchos::ParameterList> pl,
                    const Teuchos::RCP<const Teuchos::MpiComm<int> > comm,
                    const Teuchos::RCP<Teuchos::FancyOStream> out) {
    *out << "Loading solver config from " << xml << std::endl;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xml,pl.ptr(),*comm);
  }


  Teuchos::RCP<Teuchos::ParameterList> getSolverParameters(linearAlgebraType linAlgebra,
                                                           physicsType physics,
                                                           solverType solver,
                                                           int dim,
                                                           Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
                                                           Teuchos::RCP<Teuchos::FancyOStream> &out,
                                                           std::string &xml) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Not yet implemented
    if (solver == MUELU_DARCY)
      throw;

    if ((solver == AUGMENTATION) ||
        (solver == ML_REFMAXWELL) ||
        (solver == MUELU_REFMAXWELL) ||
        (solver == MUELU_MAXWELL_HO)) {
      TEUCHOS_ASSERT(physics == MAXWELL);
    } else if (solver == MUELU_DARCY) {
      TEUCHOS_ASSERT(physics == DARCY);
    }

    RCP<Teuchos::ParameterList> lin_solver_pl = Teuchos::rcp(new Teuchos::ParameterList("Linear Solver"));
    {
      if (xml == "") {
        // Load a solver configuration
        // This input deck choice depends on
        // * chosen solver
        // * linear algebra library
        // * spatial dimension
        // * node type
        if (solver == AUGMENTATION)
          if (linAlgebra == linAlgTpetra)
            updateParams("solverAugmentation.xml", lin_solver_pl, comm, out);
          else
            updateParams("solverAugmentationEpetra.xml", lin_solver_pl, comm, out);
        else if (solver == CG)
          if (linAlgebra == linAlgTpetra)
            updateParams("solverCG.xml", lin_solver_pl, comm, out);
          else
            throw;
        else if (solver == GMRES)
          if (linAlgebra == linAlgTpetra)
            updateParams("solverGMRES.xml", lin_solver_pl, comm, out);
          else
            throw;
        else if (solver == ML_REFMAXWELL) {
          updateParams("solverMLRefMaxwell.xml", lin_solver_pl, comm, out);
        } else if (solver == MUELU_REFMAXWELL || solver == MUELU_MAXWELL_HO) {
          if (linAlgebra == linAlgTpetra) {
            updateParams("solverMueLuRefMaxwell.xml", lin_solver_pl, comm, out);

            if (dim == 2)
              updateParams("solverMueLuRefMaxwell2D.xml", lin_solver_pl, comm, out);

#ifdef KOKKOS_ENABLE_OPENMP
            if (typeid(panzer::TpetraNodeType).name() == typeid(Tpetra::KokkosCompat::KokkosOpenMPWrapperNode).name()) {
              if (linAlgebra == linAlgTpetra)
                updateParams("solverMueLuRefMaxwellOpenMP.xml", lin_solver_pl, comm, out);
              else {
                std::cout << std::endl
                          << "WARNING" << std::endl
                          << "MueLu RefMaxwell + Epetra + OpenMP does currently not work." << std::endl
                          << "The Xpetra-Epetra interface is missing \"setAllValues\" with kokkos views." << std::endl << std::endl;
                throw;
              }
            }
#endif
#ifdef KOKKOS_ENABLE_CUDA
            if (typeid(panzer::TpetraNodeType).name() == typeid(Tpetra::KokkosCompat::KokkosCudaWrapperNode).name())
              updateParams("solverMueLuRefMaxwellCuda.xml", lin_solver_pl, comm, out);
#endif
          } else {
            updateParams("solverMueLuRefMaxwellEpetra.xml", lin_solver_pl, comm, out);

            if (dim == 2)
              updateParams("solverMueLuRefMaxwell2D.xml", lin_solver_pl, comm, out);
          }
          if (solver == MUELU_MAXWELL_HO) {
            RCP<Teuchos::ParameterList> lin_solver_pl_lo = lin_solver_pl;
            lin_solver_pl = rcp(new Teuchos::ParameterList("Linear Solver"));
            updateParams("solverMueLuMaxwellHO.xml", lin_solver_pl, comm, out);
#ifdef KOKKOS_ENABLE_CUDA
            if (typeid(panzer::TpetraNodeType).name() == typeid(Tpetra::KokkosCompat::KokkosCudaWrapperNode).name()) {
              updateParams("solverMueLuMaxwellHOCuda.xml", lin_solver_pl, comm, out);
            }
#endif
            Teuchos::ParameterList& mueluList = lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("Maxwell").sublist("S_E Preconditioner").sublist("Preconditioner Types").sublist("MueLu");
            if (mueluList.isParameter("coarse: type") && mueluList.get<std::string>("coarse: type") == "RefMaxwell")
              mueluList.set("coarse: params", lin_solver_pl_lo->sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("Maxwell").sublist("S_E Preconditioner").sublist("Preconditioner Types").sublist("MueLuRefMaxwell"));
          }
        } else if (solver == MUELU_DARCY) {
          updateParams("solverMueLuRefDarcy.xml", lin_solver_pl, comm, out);
        }
      } else
        updateParams(xml, lin_solver_pl, comm, out);
    }

    return lin_solver_pl;
  }


  void setClosureParameters(physicsType physics,
                            Teuchos::ParameterList &physicsEqSet,
                            Teuchos::ParameterList &closure_models,
                            double dt,
                            std::string &auxModelID) {
    std::string modelID = physicsEqSet.get<std::string>("Model ID");
    if (physics == MAXWELL) {
      auxModelID = "electromagnetics_aux";
      { // Set up closure model lists
        if (!physicsEqSet.isType<std::string>("Inverse Permeability") && !closure_models.sublist(modelID).isSublist("1/mu")) {
          std::string paramLabel = physicsEqSet.get<std::string>("Permeability");
          if (closure_models.sublist(modelID).sublist(paramLabel).isType<double>("Value")) {
            double mu = closure_models.sublist(modelID).sublist(paramLabel).get<double>("Value");
            closure_models.sublist(modelID).sublist("1/mu").set("Value",1.0/mu);
          } else
            TEUCHOS_ASSERT(false);
          physicsEqSet.set("Inverse Permeability", "1/mu");
        }

        closure_models.sublist(modelID).sublist("1/dt").set<double>("Value",1.0/dt);
        if (physicsEqSet.isType<std::string>("Current")) {
          std::string paramLabel = physicsEqSet.get<std::string>("Current");
          if (closure_models.sublist(modelID).sublist(paramLabel).get<std::string>("Type") == "GAUSSIAN PULSE")
            closure_models.sublist(modelID).sublist(paramLabel).set<double>("dt",dt); // set pulse width such that dt resolves it
        }

        // copy over entries to closure model for solver
        std::vector<std::string> parameters = {"Permittivity", "Permeability", "Conductivity", "Inverse Permeability"};
        for (auto it = parameters.begin(); it != parameters.end(); ++it) {
          std::string paramLabel = physicsEqSet.get<std::string>(*it);
          closure_models.sublist(auxModelID).sublist(paramLabel) = closure_models.sublist(modelID).sublist(paramLabel);
          if (closure_models.sublist(auxModelID).sublist(paramLabel).isType<std::string>("DoF Name"))
            closure_models.sublist(auxModelID).sublist(paramLabel).set("DoF Name", "AUXILIARY_EDGE");
        }
      } // Set up closure model lists
    }
    else if (physics == DARCY) {
      auxModelID = "darcy_aux";

      closure_models.sublist(modelID).sublist("1/dt").set<double>("Value",1.0/dt);

      // copy over entries to closure model for solver
      std::vector<std::string> parameters = {"Inverse Diffusivity"};
      for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        std::string paramLabel = physicsEqSet.get<std::string>(*it);
        closure_models.sublist(auxModelID).sublist(paramLabel) = closure_models.sublist(modelID).sublist(paramLabel);
        if (closure_models.sublist(auxModelID).sublist(paramLabel).isType<std::string>("DoF Name"))
          closure_models.sublist(auxModelID).sublist(paramLabel).set("DoF Name", "AUXILIARY_FACE");
      }
    }
    closure_models.sublist(auxModelID).sublist("dt").set<double>("Value",dt);
    closure_models.sublist(auxModelID).sublist("1/dt").set<double>("Value",1.0/dt);
  }


  void setAuxiliaryOperatorParameters(physicsType physics,
                                      solverType solver,
                                      int basis_order,
                                      std::string pCoarsenScheduleStr,
                                      bool matrixFree,
                                      Teuchos::ParameterList &input_params,
                                      Teuchos::ParameterList &lin_solver_pl,
                                      std::string &auxFieldOrder) {

    Teuchos::ParameterList & physicsBlock_pl = input_params.sublist("Physics Blocks");
    Teuchos::ParameterList & assembly_pl     = input_params.sublist("Assembly");
    Teuchos::ParameterList & aux_ops_pl      = input_params.sublist("Auxiliary Operators");

    auxFieldOrder = "blocked:";

    std::string auxModelID;
    if (physics == MAXWELL)
      auxModelID = "electromagnetics_aux";
    else if (physics == DARCY)
      auxModelID = "darcy_aux";

    pCoarsenScheduleStr = assembly_pl.get<std::string>("p coarsen schedule", pCoarsenScheduleStr);
    std::vector<std::string> pCoarsenScheduleVecStr;
    std::vector<int> pCoarsenSchedule;
    panzer::StringTokenizer(pCoarsenScheduleVecStr, pCoarsenScheduleStr, ",");
    panzer::TokensToInts(pCoarsenSchedule, pCoarsenScheduleVecStr);
    if (basis_order > 1)
      pCoarsenSchedule.insert(pCoarsenSchedule.begin(), basis_order);

    { // Check that this is a valid schedule.
      auto it = pCoarsenSchedule.begin();
      int p = *it;
      TEUCHOS_ASSERT_EQUALITY(p, basis_order);
      ++it;
      while (it != pCoarsenSchedule.end()) {
        int q = *it;
        TEUCHOS_ASSERT(q < p);
        ++it;
        p = q;
      }
      TEUCHOS_ASSERT_EQUALITY(pCoarsenSchedule.back(), 1);
    }

    if (lin_solver_pl.sublist("Preconditioner Types").isSublist("Teko") &&
        lin_solver_pl.sublist("Preconditioner Types").sublist("Teko").isSublist("Inverse Factory Library")) {
      if (lin_solver_pl.sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").isSublist("Maxwell"))
        lin_solver_pl.sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("Maxwell").set("p coarsen schedule",std::to_string(basis_order)+","+pCoarsenScheduleStr);
      if (lin_solver_pl.sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").isSublist("Darcy"))
        lin_solver_pl.sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("Darcy").set("p coarsen schedule",std::to_string(basis_order)+","+pCoarsenScheduleStr);
    }

    Teuchos::ParameterList& auxPhysicsBlocksPL = physicsBlock_pl.sublist("Auxiliary Physics Block");

    for (auto it = pCoarsenSchedule.begin(); it != pCoarsenSchedule.end(); ++it) {

      if (physics == MAXWELL) {
        std::string auxNodalField, auxEdgeField, opPostfix;
        int polynomialOrder = *it;
        // Are we setting up lower order operators?
        if (polynomialOrder != basis_order) {
          auxNodalField = "AUXILIARY_NODE_" + std::to_string(polynomialOrder);
          auxEdgeField = "AUXILIARY_EDGE_" + std::to_string(polynomialOrder);
          opPostfix = " "+std::to_string(polynomialOrder);
        } else {
          auxNodalField = "AUXILIARY_NODE";
          auxEdgeField = "AUXILIARY_EDGE";
          opPostfix = "";
        }

        if (solver == MUELU_REFMAXWELL || solver == ML_REFMAXWELL || solver == MUELU_MAXWELL_HO)
          auxFieldOrder += " "+auxNodalField+" "+auxEdgeField;
        else
          auxFieldOrder += " "+auxEdgeField;

        if (solver == MUELU_REFMAXWELL || solver == ML_REFMAXWELL || solver == MUELU_MAXWELL_HO) {
          // discrete gradient
          auto gradPL = Teuchos::ParameterList();
          gradPL.set("Source", auxNodalField);
          gradPL.set("Target", auxEdgeField);
          gradPL.set("Op", "grad");
          gradPL.set("matrix-free", polynomialOrder != 1 ? matrixFree : false);
          aux_ops_pl.sublist("Discrete Gradient"+opPostfix) = gradPL;
        }

        // Schur complement
        auto schurComplementPL = Teuchos::ParameterList();
        schurComplementPL.set("Type", "Auxiliary SchurComplement");
        schurComplementPL.set("DOF Name", auxEdgeField);
        schurComplementPL.set("Basis Type", "HCurl");
        schurComplementPL.set("Model ID", auxModelID);
        schurComplementPL.set("Permittivity", "epsilon");
        schurComplementPL.set("Conductivity", "sigma");
        schurComplementPL.set("Inverse Permeability", "1/mu");
        schurComplementPL.set("Basis Order", polynomialOrder);
        schurComplementPL.set("Integration Order", 2*polynomialOrder);
        auxPhysicsBlocksPL.sublist("Auxiliary Edge SchurComplement Physics"+opPostfix) = schurComplementPL;

        if (solver == MUELU_MAXWELL_HO) {
          // Projected Schur complement
          auto projectedSchurComplementPL = Teuchos::ParameterList();
          projectedSchurComplementPL.set("Type", "Auxiliary ProjectedSchurComplement");
          projectedSchurComplementPL.set("DOF Name", auxNodalField);
          projectedSchurComplementPL.set("Basis Type", "HGrad");
          projectedSchurComplementPL.set("Model ID", auxModelID);
          projectedSchurComplementPL.set("Permittivity", "epsilon");
          projectedSchurComplementPL.set("Conductivity", "sigma");
          projectedSchurComplementPL.set("Basis Order", polynomialOrder);
          projectedSchurComplementPL.set("Integration Order", 2*polynomialOrder);
          auxPhysicsBlocksPL.sublist("Auxiliary Node ProjectedSchurComplement"+opPostfix) = projectedSchurComplementPL;
        }
      } else if (physics == DARCY) {
        std::string auxFaceField, auxEdgeField, opPostfix;
        int polynomialOrder = *it;
        // Are we setting up lower order operators?
        if (polynomialOrder != basis_order) {
          auxFaceField = "AUXILIARY_FACE_" + std::to_string(polynomialOrder);
          auxEdgeField = "AUXILIARY_EDGE_" + std::to_string(polynomialOrder);
          opPostfix = " "+std::to_string(polynomialOrder);
        } else {
          auxFaceField = "AUXILIARY_FACE";
          auxEdgeField = "AUXILIARY_EDGE";
          opPostfix = "";
        }

        if (solver == MUELU_DARCY)
          auxFieldOrder += " "+auxEdgeField + " "+auxFaceField;
        else
          auxFieldOrder += " "+auxFaceField;

        if (solver == MUELU_DARCY) {
          // discrete curl
          auto curlPL = Teuchos::ParameterList();
          curlPL.set("Source", auxEdgeField);
          curlPL.set("Target", auxFaceField);
          curlPL.set("Op", "curl");
          curlPL.set("matrix-free", polynomialOrder != 1 ? matrixFree : false);
          aux_ops_pl.sublist("Discrete Curl"+opPostfix) = curlPL;
        }

        // Schur complement
        auto schurComplementPL = Teuchos::ParameterList();
        schurComplementPL.set("Type", "Auxiliary DarcySchurComplement");
        schurComplementPL.set("DOF Name", auxFaceField);
        schurComplementPL.set("Basis Type", "HDiv");
        schurComplementPL.set("Model ID", auxModelID);
        schurComplementPL.set("Inverse Diffusivity", "1/kappa");
        schurComplementPL.set("Basis Order", polynomialOrder);
        schurComplementPL.set("Integration Order", 2*polynomialOrder);
        auxPhysicsBlocksPL.sublist("Auxiliary Face DarcySchurComplement Physics"+opPostfix) = schurComplementPL;
      }

    }

    // Set up additional mass matrices for RefMaxwell
    if (solver == MUELU_REFMAXWELL || solver == ML_REFMAXWELL || solver == MUELU_MAXWELL_HO) {
      std::string auxNodalField, auxEdgeField, opPostfix;
      if (basis_order != 1) {
        auxNodalField = "AUXILIARY_NODE_" + std::to_string(1);
        auxEdgeField = "AUXILIARY_EDGE_" + std::to_string(1);
        opPostfix = " "+std::to_string(1);
      } else {
        auxNodalField = "AUXILIARY_NODE";
        auxEdgeField = "AUXILIARY_EDGE";
        opPostfix = "";
      }

      // Edge mass matrix with unit weight
      auto massEdgePL = Teuchos::ParameterList();
      massEdgePL.set("Type", "Auxiliary Mass Matrix");
      massEdgePL.set("DOF Name", auxEdgeField);
      massEdgePL.set("Basis Type", "HCurl");
      massEdgePL.set("Model ID", auxModelID);
      massEdgePL.set("Basis Order", 1);
      massEdgePL.set("Integration Order", 2);
      auxPhysicsBlocksPL.sublist("Auxiliary Edge Mass Physics"+opPostfix) = massEdgePL;

      // Edge mass matrix with 1/mu weight
      auto massEdgeWeightedPL = Teuchos::ParameterList();
      massEdgeWeightedPL.set("Type", "Auxiliary Mass Matrix");
      massEdgeWeightedPL.set("DOF Name", auxEdgeField);
      massEdgeWeightedPL.set("Basis Type", "HCurl");
      massEdgeWeightedPL.set("Model ID", auxModelID);
      massEdgeWeightedPL.set("Field Multipliers", "1/mu");
      massEdgeWeightedPL.set("Basis Order", 1);
      massEdgeWeightedPL.set("Integration Order", 2);
      massEdgeWeightedPL.set("Operator Label", "weighted ");
      auxPhysicsBlocksPL.sublist("Auxiliary Edge Mass Physics weighted"+opPostfix) = massEdgeWeightedPL;

      // Nodal mass matrix
      auto massNodePL = Teuchos::ParameterList();
      massNodePL.set("Type", "Auxiliary Mass Matrix");
      massNodePL.set("DOF Name", auxNodalField);
      massNodePL.set("Basis Type", "HGrad");
      massNodePL.set("Model ID", auxModelID);
      massNodePL.set("Field Multipliers", "mu,1/dt");
      massNodePL.set("Basis Order", 1);
      massNodePL.set("Integration Order", 2);
      auxPhysicsBlocksPL.sublist("Auxiliary Node Mass Physics"+opPostfix) = massNodePL;

    } else if (solver == MUELU_DARCY) {

      std::string auxEdgeField, auxFaceField, opPostfix;
      if (basis_order != 1) {
        auxEdgeField = "AUXILIARY_EDGE_" + std::to_string(1);
        auxFaceField = "AUXILIARY_FACE_" + std::to_string(1);
        opPostfix = " "+std::to_string(1);
      } else {
        auxEdgeField = "AUXILIARY_EDGE";
        auxFaceField = "AUXILIARY_FACE";
        opPostfix = "";
      }

      // Face mass matrix with unit weight
      auto massFacePL = Teuchos::ParameterList();
      massFacePL.set("Type", "Auxiliary Mass Matrix");
      massFacePL.set("DOF Name", auxFaceField);
      massFacePL.set("Basis Type", "HDiv");
      massFacePL.set("Model ID", auxModelID);
      massFacePL.set("Basis Order", 1);
      massFacePL.set("Integration Order", 2);
      auxPhysicsBlocksPL.sublist("Auxiliary Face Mass Physics"+opPostfix) = massFacePL;

      // Face mass matrix with 1/mu weight
      // auto massFaceWeightedPL = Teuchos::ParameterList();
      // massFaceWeightedPL.set("Type", "Auxiliary Mass Matrix");
      // massFaceWeightedPL.set("DOF Name", auxFaceField);
      // massFaceWeightedPL.set("Basis Type", "HDiv");
      // massFaceWeightedPL.set("Model ID", auxModelID);
      // massFaceWeightedPL.set("Field Multipliers", "1/mu");
      // massFaceWeightedPL.set("Basis Order", 1);
      // massFaceWeightedPL.set("Integration Order", 2);
      // massFaceWeightedPL.set("Operator Label", "weighted ");
      // auxPhysicsBlocksPL.sublist("Auxiliary Face Mass Physics weighted"+opPostfix) = massFaceWeightedPL;

      // Edge mass matrix
      auto massEdgePL = Teuchos::ParameterList();
      massEdgePL.set("Type", "Auxiliary Mass Matrix");
      massEdgePL.set("DOF Name", auxEdgeField);
      massEdgePL.set("Basis Type", "HCurl");
      massEdgePL.set("Model ID", auxModelID);
      massEdgePL.set("Field Multipliers", "1/dt");
      massEdgePL.set("Basis Order", 1);
      massEdgePL.set("Integration Order", 2);
      auxPhysicsBlocksPL.sublist("Auxiliary Edge Mass Physics"+opPostfix) = massEdgePL;

    }

    // Set up interpolations between levels
    auto it = pCoarsenSchedule.begin();
    int p = *it;
    ++it;
    while (it != pCoarsenSchedule.end()) {
      int q = *it;

      std::string space, space2;
      if (physics == MAXWELL) {
        space = "AUXILIARY_EDGE";
        space2 = "Hcurl";
      } else if (physics == DARCY) {
        space = "AUXILIARY_FACE";
        space2 = "Hdiv";
      }

      auto interpPL = Teuchos::ParameterList();
      interpPL.set("Source", space+"_"+std::to_string(q));
      interpPL.set("Target", p != basis_order ? space+"_"+std::to_string(p) : space);
      interpPL.set("Op", "value");
      interpPL.set("matrix-free", matrixFree);
      aux_ops_pl.sublist("Interpolation " + space2 + " " + std::to_string(q) + "->" + std::to_string(p)) = interpPL;

      p = q;
      ++it;
    }
  }


  void createExodusFile(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                        Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory,
                        Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                        const bool & exodus_out,
                        Teuchos::RCP<const Teuchos::MpiComm<int> > comm,
                        physicsType physics) {
    for(std::size_t i=0;i<physicsBlocks.size();i++) {
      Teuchos::RCP<panzer::PhysicsBlock> pb = physicsBlocks[i]; // we are assuming only one physics block

      const std::vector<panzer::StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end());

      // build string for modifiying vectors
      std::vector<std::string> dimenStr(3);
      dimenStr[0] = "X"; dimenStr[1] = "Y"; dimenStr[2] = "Z";

      // add basis to DOF manager: block specific
      std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator fieldItr;
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
        Teuchos::RCP<const panzer::PureBasis> basis = fieldItr->second;
        if(basis->getElementSpace()==panzer::PureBasis::HGRAD)
          mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
        else if(basis->getElementSpace()==panzer::PureBasis::CONST )
          mesh->addCellField(fieldItr->first,pb->elementBlockID());
        else if(basis->getElementSpace()==panzer::PureBasis::HCURL ||
                basis->getElementSpace()==panzer::PureBasis::HDIV    ) {
          for(int dim=0;dim<basis->dimension();++dim)
            mesh->addCellField(fieldItr->first+dimenStr[dim],pb->elementBlockID());
        } else if(basis->getElementSpace()==panzer::PureBasis::HVOL)
          mesh->addCellField(fieldItr->first,pb->elementBlockID());
      }

      std::vector<std::string> block_names;
      mesh->getElementBlockNames(block_names);

      Teuchos::ParameterList output_pl("Output");
      output_pl.sublist("Cell Average Quantities");
      Teuchos::ParameterList& cell_avg_v = output_pl.sublist("Cell Average Vectors");
      if (physics == MAXWELL)
        cell_avg_v.set(block_names[0],"J");
      output_pl.sublist("Cell Quantities");
      output_pl.sublist("Nodal Quantities");
      output_pl.sublist("Allocate Nodal Quantities");
      mini_em::addFieldsToMesh(*mesh,output_pl);
    }
    mesh_factory->completeMeshConstruction(*mesh,(*comm->getRawMpiComm())());

    if (exodus_out)
      mesh->setupExodusFile("mesh_output.exo");
  }

  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
  buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
                            const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
                            const Teuchos::RCP<panzer::GlobalIndexer> & globalIndexer,
                            const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                            const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                            const Teuchos::ParameterList & closure_model_pl,
                            physicsType physics) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
      = rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,globalIndexer,linObjFactory));

    // get a vector of all the element blocks
    std::vector<std::string> eBlocks;
    mesh->getElementBlockNames(eBlocks);

    panzer_stk::RespFactorySolnWriter_Builder builder;
    builder.mesh = mesh;

    stkIOResponseLibrary->addResponse("Main Field Output",eBlocks,builder);

    std::vector<std::string> block_names;
    mesh->getElementBlockNames(block_names);

    // this automatically adds in the nodal fields
    Teuchos::ParameterList output_pl("Output");
    output_pl.sublist("Cell Average Quantities");
    Teuchos::ParameterList& cell_avg_v = output_pl.sublist("Cell Average Vectors");
    if (physics == MAXWELL)
      cell_avg_v.set(block_names[0],"J");

    output_pl.sublist("Cell Quantities");
    output_pl.sublist("Nodal Quantities");
    output_pl.sublist("Allocate Nodal Quantities");
    panzer_stk::IOClosureModelFactory_TemplateBuilder<panzer::Traits> io_cm_builder(cm_factory,mesh,
                                                                                    output_pl);
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> io_cm_factory;
    io_cm_factory.buildObjects(io_cm_builder);

    stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                                  io_cm_factory,
                                                  closure_model_pl,
                                                  Teuchos::ParameterList());

    return stkIOResponseLibrary;
  }

}
