//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_MAXWELL1_DEF_HPP
#define MUELU_MAXWELL1_DEF_HPP

#include <sstream>

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "MueLu_Maxwell1_decl.hpp"
#include "MueLu_Maxwell_Utils.hpp"

#include "MueLu_ReitzingerPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_HierarchyManager.hpp"
#include <MueLu_HierarchyUtils.hpp>
# include "MueLu_Utilities_kokkos.hpp"
#include "MueLu_VerbosityLevel.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_ML2MueLuParameterTranslator.hpp>
#include <MueLu_RefMaxwellSmoother.hpp>

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif


namespace MueLu {


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {
    return SM_Matrix_->getDomainMap();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {
    return SM_Matrix_->getRangeMap();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameters(Teuchos::ParameterList& list) {

    if (list.isType<std::string>("parameterlist: syntax") && list.get<std::string>("parameterlist: syntax") == "ml") {
      list.remove("parameterlist: syntax");
      Teuchos::ParameterList newList;

      // interpret ML list
      newList.sublist("maxwell1: 22list") = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list,"Maxwell"));


      // Hardwiring options to ensure ML compatibility
      newList.sublist("maxwell1: 22list").set("use kokkos refactor", false);
      newList.sublist("maxwell1: 22list").set("tentative: constant column sums", false);
      newList.sublist("maxwell1: 22list").set("tentative: calculate qr", false);

      newList.sublist("maxwell1: 11list").set("use kokkos refactor", false);
      newList.sublist("maxwell1: 11list").set("tentative: constant column sums", false);
      newList.sublist("maxwell1: 11list").set("tentative: calculate qr", false);

      if(list.isParameter("aggregation: damping factor") && list.get<double>("aggregation: damping factor") == 0.0)
        newList.sublist("maxwell1: 11list").set("multigrid algorithm", "unsmoothed reitzinger");
      else
        newList.sublist("maxwell1: 11list").set("multigrid algorithm", "smoothed reitzinger");
      newList.sublist("maxwell1: 11list").set("aggregation: type", "uncoupled");

      newList.sublist("maxwell1: 22list").set("multigrid algorithm", "unsmoothed");
      newList.sublist("maxwell1: 22list").set("aggregation: type", "uncoupled");

      if (newList.sublist("maxwell1: 22list").isType<std::string>("verbosity"))
        newList.set("verbosity", newList.sublist("maxwell1: 22list").get<std::string>("verbosity"));

      // Move coarse solver and smoother stuff to 11list
      std::vector<std::string> convert = {"coarse:", "smoother:", "smoother: pre", "smoother: post"};
      for (auto it = convert.begin(); it != convert.end(); ++it) {
        if (newList.sublist("maxwell1: 22list").isType<std::string>(*it + " type")) {
          newList.sublist("maxwell1: 11list").set(*it+" type", newList.sublist("maxwell1: 22list").get<std::string>(*it+" type"));
          newList.sublist("maxwell1: 22list").remove(*it+" type");
        }
        if (newList.sublist("maxwell1: 22list").isSublist(*it+" params")) {
          newList.sublist("maxwell1: 11list").set(*it+" params", newList.sublist("maxwell1: 22list").sublist(*it+" params"));
          newList.sublist("maxwell1: 22list").remove(*it+" params");
        }
      }

      newList.sublist("maxwell1: 22list").set("smoother: type", "none");
      newList.sublist("maxwell1: 22list").set("coarse: type", "none");

      list = newList;
    }
    std::string  mode_string   = list.get("maxwell1: mode",                  MasterList::getDefault<std::string>("maxwell1: mode"));
    applyBCsTo22_              = list.get("maxwell1: apply BCs to 22",       false);
    dump_matrices_             = list.get("maxwell1: dump matrices",         MasterList::getDefault<bool>("maxwell1: dump matrices"));

    // Default smoother.  We'll copy this around.
    Teuchos::ParameterList defaultSmootherList;
    defaultSmootherList.set("smoother: type", "CHEBYSHEV");
    defaultSmootherList.sublist("smoother: params").set("chebyshev: degree",2);
    defaultSmootherList.sublist("smoother: params").set("chebyshev: ratio eigenvalue",7.0);
    defaultSmootherList.sublist("smoother: params").set("chebyshev: eigenvalue max iterations",30);

    // Make sure verbosity gets passed to the sublists
    std::string verbosity = list.get("verbosity","high");
    VerboseObject::SetDefaultVerbLevel(toVerbLevel(verbosity));

    // Check the validity of the run mode
    if(mode_ != MODE_GMHD_STANDARD) {
      if(mode_string == "standard")          mode_ = MODE_STANDARD;
      else if(mode_string == "refmaxwell")   mode_ = MODE_REFMAXWELL;
      else if(mode_string == "edge only")    mode_ = MODE_EDGE_ONLY;
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Must use mode 'standard', 'refmaxwell', 'edge only', or use the GMHD constructor.");
      }
    }

    // If we're in edge only or standard modes, then the (2,2) hierarchy gets built without smoothers.
    // Otherwise we use the user's smoothers (defaulting to Chebyshev if unspecified)
    if(list.isSublist("maxwell1: 22list"))
      precList22_     =  list.sublist("maxwell1: 22list");
    else if(list.isSublist("refmaxwell: 22list"))
      precList22_     =  list.sublist("refmaxwell: 22list");
    if(mode_ == MODE_EDGE_ONLY || mode_ == MODE_STANDARD || mode_ == MODE_GMHD_STANDARD)
      precList22_.set("smoother: pre or post","none");
    else if(!precList22_.isType<std::string>("Preconditioner Type") &&
       !precList22_.isType<std::string>("smoother: type") &&
       !precList22_.isType<std::string>("smoother: pre type") &&
       !precList22_.isType<std::string>("smoother: post type")) {
      precList22_ = defaultSmootherList;
    }
    precList22_.set("verbosity",precList22_.get("verbosity",verbosity));

    

    // For the (1,1) hierarchy we'll use Hiptmair (STANDARD) or Chebyshev (EDGE_ONLY / REFMAXWELL) if
    // the user doesn't specify things
    if(list.isSublist("maxwell1: 11list"))
      precList11_     =  list.sublist("maxwell1: 11list");
    else if(list.isSublist("refmaxwell: 11list"))
      precList11_     =  list.sublist("refmaxwell: 11list");

    if(mode_ == MODE_GMHD_STANDARD) {
      precList11_.set("smoother: pre or post","none");
      precList11_.set("smoother: type", "none");
    }
    if(!precList11_.isType<std::string>("Preconditioner Type") &&
       !precList11_.isType<std::string>("smoother: type") &&
       !precList11_.isType<std::string>("smoother: pre type") &&
       !precList11_.isType<std::string>("smoother: post type")) {
      if(mode_ == MODE_EDGE_ONLY || mode_ == MODE_REFMAXWELL) {
        precList11_ = defaultSmootherList;
      }
      if (mode_ == MODE_STANDARD)  {
        precList11_.set("smoother: type", "HIPTMAIR");
        precList11_.sublist("hiptmair: smoother type 1","CHEBYSHEV");
        precList11_.sublist("hiptmair: smoother type 2","CHEBYSHEV");
        precList11_.sublist("hiptmair: smoother list 1") = defaultSmootherList;
        precList11_.sublist("hiptmair: smoother list 2") = defaultSmootherList;
      }
    }
    precList11_.set("verbosity",precList11_.get("verbosity",verbosity));

    // Reuse support
    if (enable_reuse_ &&
        !precList11_.isType<std::string>("Preconditioner Type") &&
        !precList11_.isParameter("reuse: type"))
      precList11_.set("reuse: type", "full");
    if (enable_reuse_ &&
        !precList22_.isType<std::string>("Preconditioner Type") &&
        !precList22_.isParameter("reuse: type"))
      precList22_.set("reuse: type", "full");


    // Are we using Kokkos?
# ifdef HAVE_MUELU_SERIAL
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosSerialWrapperNode).name())
      useKokkos_ = false;
# endif
# ifdef HAVE_MUELU_OPENMP
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosOpenMPWrapperNode).name())
      useKokkos_ = true;
# endif
# ifdef HAVE_MUELU_CUDA
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name())
      useKokkos_ = true;
# endif
    useKokkos_ = list.get("use kokkos refactor",useKokkos_);

    parameterList_ = list;

  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GMHDSetupHierarchy(Teuchos::ParameterList& List) const {
    Teuchos::ParameterList precListGmhd;

    MueLu::HierarchyUtils<SC,LO,GO,NO>::CopyBetweenHierarchies(*Hierarchy11_,*HierarchyGmhd_, "P",  "Psubblock", "RCP<Matrix>");

    HierarchyGmhd_->GetLevel(0)->Set("A", GmhdA_Matrix_);
    GmhdA_Matrix_->setObjectLabel("GmhdA");

    TEUCHOS_TEST_FOR_EXCEPTION( !List.isSublist("maxwell1: Gmhdlist"), Exceptions::RuntimeError, "Must provide maxwell1: Gmhdlist for GMHD setup");    
    precListGmhd     =  List.sublist("maxwell1: Gmhdlist");
    precListGmhd.set("coarse: max size",1);
    precListGmhd.set("max levels",HierarchyGmhd_->GetNumLevels());
    RCP<MueLu::HierarchyManager<SC,LO,GO,NO> > mueLuFactory = rcp(new MueLu::ParameterListInterpreter<SC,LO,GO,NO>(precListGmhd,GmhdA_Matrix_->getDomainMap()->getComm()));
    HierarchyGmhd_->setlib(GmhdA_Matrix_->getDomainMap()->lib());
    HierarchyGmhd_->SetProcRankVerbose(GmhdA_Matrix_->getDomainMap()->getComm()->getRank());
    mueLuFactory->SetupHierarchy(*HierarchyGmhd_);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute(bool reuse) {
    /* Algorithm overview for Maxwell1 construction:

       1) Create a nodal auxillary hierarchy based on (a) the user's nodal matrix or (b) a matrix constructed
       by D0^T A D0 if the user doesn't provide a nodal matrix.  We call this matrix "NodeAggMatrix."

       2)  If the user provided a node matrix, we use the prolongators from the auxillary nodal hierarchy
       to generate matrices for smoothers on all levels.  We call this "NodeMatrix."  Otherwise NodeMatrix = NodeAggMatrix

       3) We stick all of the nodal P matrices and NodeMatrix objects on the final (1,1) hierarchy and then use the
       ReitzingerPFactory to generate the edge P and A matrices.
     */


#ifdef HAVE_MUELU_CUDA
    if (parameterList_.get<bool>("maxwell1: cuda profile setup", false)) cudaProfilerStart();
#endif

    std::string timerLabel;
    if (reuse)
      timerLabel = "MueLu Maxwell1: compute (reuse)";
    else
      timerLabel = "MueLu Maxwell1: compute";
    RCP<Teuchos::TimeMonitor> tmCompute = getTimer(timerLabel);

    ////////////////////////////////////////////////////////////////////////////////
    // Generate Kn and apply BCs (if needed)
    bool have_generated_Kn = false;
    if(Kn_Matrix_.is_null()) {
      GetOStream(Runtime0) << "Maxwell1::compute(): Kn not provided.  Generating." << std::endl;
      Kn_Matrix_ = generate_kn();
      have_generated_Kn = true;
    }

    if (parameterList_.get<bool>("rap: fix zero diagonals", true)) {
      magnitudeType threshold;
      if (parameterList_.isType<magnitudeType>("rap: fix zero diagonals threshold"))
        threshold = parameterList_.get<magnitudeType>("rap: fix zero diagonals threshold",
                                                      Teuchos::ScalarTraits<double>::eps());
      else
        threshold = Teuchos::as<magnitudeType>(parameterList_.get<double>("rap: fix zero diagonals threshold",
                                                                          Teuchos::ScalarTraits<double>::eps()));
      Scalar replacement = Teuchos::as<Scalar>(parameterList_.get<double>("rap: fix zero diagonals replacement",
                                                                          MasterList::getDefault<double>("rap: fix zero diagonals replacement")));
      Xpetra::MatrixUtils<SC,LO,GO,NO>::CheckRepairMainDiagonal(Kn_Matrix_, true, GetOStream(Warnings1), threshold, replacement);
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Generate the (2,2) Hierarchy
    Kn_Matrix_->setObjectLabel("Maxwell1 (2,2)");
    
    /* Critical ParameterList changes */
    if (!Coords_.is_null())
      precList22_.sublist("user data").set("Coordinates",Coords_);

    /* Repartitioning *must* be in sync between hierarchies, but the
     only thing we need to watch here is the subcomms, since ReitzingerPFactory
     won't look at all the other stuff */
    if(precList22_.isParameter("repartition: enable")) {
      bool repartition = precList22_.get<bool>("repartition: enable");
      precList11_.set("repartition: enable",repartition);

      // If we're repartitioning (2,2), we need to rebalance for (1,1) to do the right thing
      if(repartition)
        precList22_.set("repartition: rebalance P and R",true);

      if (precList22_.isParameter("repartition: use subcommunicators")) {
        precList11_.set("repartition: use subcommunicators", precList22_.get<bool>("repartition: use subcommunicators"));
        
        // We do not want (1,1) and (2,2) blocks being repartitioned seperately, so we specify the map that
        // is going to be used (this is generated in ReitzingerPFactory)
        if(precList11_.get<bool>("repartition: use subcommunicators")==true)
          precList11_.set("repartition: use subcommunicators in place",true);
      }
      else {
        // We'll have Maxwell1 default to using subcommunicators if you don't specify
        precList11_.set("repartition: use subcommunicators", true);
        precList22_.set("repartition: use subcommunicators", true);
        
        // We do not want (1,1) and (2,2) blocks being repartitioned seperately, so we specify the map that
        // is going to be used (this is generated in ReitzingerPFactory)
        precList11_.set("repartition: use subcommunicators in place",true);
      }        
        
    }
    else
      precList11_.remove("repartition: enable", false);
   

    ////////////////////////////////////////////////////////////////////////////////
    // Remove explicit zeros from matrices
    /*
    Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(parameterList_,D0_Matrix_,SM_Matrix_);


    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());;
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo",          true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*SM_Matrix_, "SM_Matrix", params);
    }
    */


    ////////////////////////////////////////////////////////////////////////////////
    // Detect Dirichlet boundary conditions
    if (!reuse) {
      magnitudeType rowSumTol = precList11_.get("aggregation: row sum drop tol",-1.0);
      Maxwell_Utils<SC,LO,GO,NO>::detectBoundaryConditionsSM(SM_Matrix_,D0_Matrix_,rowSumTol,
                                                             useKokkos_,BCrowsKokkos_,BCcolsKokkos_,BCdomainKokkos_,
                                                             BCedges_,BCnodes_,BCrows_,BCcols_,BCdomain_,
                                                             allEdgesBoundary_,allNodesBoundary_);
      if (IsPrint(Statistics2)) {
        GetOStream(Statistics2) << "MueLu::Maxwell1::compute(): Detected " << BCedges_ << " BC rows and " << BCnodes_ << " BC columns." << std::endl;
      }
    }

    if (allEdgesBoundary_) {
      // All edges have been detected as boundary edges.
      // Do not attempt to construct sub-hierarchies, but just set up a single level preconditioner.
      GetOStream(Warnings0) << "All edges are detected as boundary edges!" << std::endl;
      mode_ = MODE_EDGE_ONLY;

      // Generate single level hierarchy for the edge
      precList22_.set("max levels", 1);
    }
     
 
    if (allNodesBoundary_) {
      // All Nodes have been detected as boundary nodes.
      // Do not attempt to construct sub-hierarchies, but just set up a single level preconditioner.
      GetOStream(Warnings0) << "All nodes are detected as boundary nodes!" << std::endl;
      mode_ = MODE_EDGE_ONLY;

      // Generate single level hierarchy for the edge
      precList22_.set("max levels", 1);
    }
                                              
    ////////////////////////////////////////////////////////////////////////////////
    // Build (2,2) hierarchy
    Hierarchy22_ = MueLu::CreateXpetraPreconditioner(Kn_Matrix_, precList22_);


    ////////////////////////////////////////////////////////////////////////////////
    // Apply boundary conditions to D0 (if needed)
    if(!reuse) {
      D0_Matrix_->resumeFill();
      Scalar replaceWith;
      if (D0_Matrix_->getRowMap()->lib() == Xpetra::UseEpetra)
        replaceWith= Teuchos::ScalarTraits<SC>::eps();
      else
        replaceWith = Teuchos::ScalarTraits<SC>::zero();
      
      if(applyBCsTo22_) {
        GetOStream(Runtime0) << "Maxwell1::compute(): nuking BC rows/cols of D0" << std::endl;        
        if (useKokkos_) {
          Utilities_kokkos::ZeroDirichletCols(D0_Matrix_,BCcolsKokkos_,replaceWith);
        } else {
          Utilities::ZeroDirichletCols(D0_Matrix_,BCcols_,replaceWith);
        }
      }
      else {
        GetOStream(Runtime0) << "Maxwell1::compute(): nuking BC rows of D0" << std::endl;        
        if (useKokkos_) {
          Utilities_kokkos::ZeroDirichletRows(D0_Matrix_,BCrowsKokkos_,replaceWith);
        } else {
          Utilities::ZeroDirichletRows(D0_Matrix_,BCrows_,replaceWith);
        }
      }

      D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
    }


    ////////////////////////////////////////////////////////////////////////////////
    // What ML does is generate nodal prolongators with an auxillary hierarchy based on the 
    // user's (2,2) matrix.  The actual nodal matrices for smoothing are generated by the
    // Hiptmair smoother construction.  We're not going to do that --- we'll 
    // do as we insert them into the final (1,1) hierarchy.

    // Level 0
    RCP<Matrix> Kn_Smoother_0;
    if(have_generated_Kn) {
      Kn_Smoother_0 = Kn_Matrix_;
    }
    else {
      Kn_Smoother_0 = generate_kn();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Copy the relevant (2,2) data to the (1,1) hierarchy
    Hierarchy11_ = rcp(new Hierarchy("Maxwell1 (1,1)"));
    RCP<Matrix> OldSmootherMatrix;
    for(int i=0; i<Hierarchy22_->GetNumLevels(); i++) {
      Hierarchy11_->AddNewLevel();
      RCP<Level> NodeL = Hierarchy22_->GetLevel(i);
      RCP<Level> EdgeL = Hierarchy11_->GetLevel(i);
      RCP<Operator> NodeOp      = NodeL->Get<RCP<Operator> >("A");
      RCP<Matrix> NodeAggMatrix  = rcp_dynamic_cast<Matrix>(NodeOp);
      std::string labelstr = FormattingHelper::getColonLabel(EdgeL->getObjectLabel());      

      if(i==0) {
        EdgeL->Set("A", SM_Matrix_);
        EdgeL->Set("D0", D0_Matrix_);

        EdgeL->Set("NodeAggMatrix",NodeAggMatrix);
        EdgeL->Set("NodeMatrix",Kn_Smoother_0);
        OldSmootherMatrix = Kn_Smoother_0;
      }
      else {
        // Set the Nodal P
        auto NodalP = NodeL->Get<RCP<Matrix> >("P");
        EdgeL->Set("Pnodal",NodalP);

        // If we repartition a processor away, a RCP<Operator> is stuck
        // on the level instead of an RCP<Matrix>
        if(!NodeAggMatrix.is_null()) {
          EdgeL->Set("NodeAggMatrix",NodeAggMatrix);
          if(!have_generated_Kn) {
            // The user gave us a Kn, so we'll need to create the smoother matrix via RAP 
            RCP<Matrix> NewKn = Maxwell_Utils<SC,LO,GO,NO>::PtAPWrapper(OldSmootherMatrix,NodalP,parameterList_,labelstr);
            EdgeL->Set("NodeMatrix",NewKn);
            OldSmootherMatrix = NewKn;
          }
          else {          
            // The user didn't give us a Kn, so we aggregate and smooth with the same matrix
            EdgeL->Set("NodeMatrix",NodeAggMatrix);
          }
        }
        else {
          // We've partitioned things away.
          EdgeL->Set("NodeMatrix",NodeOp);
          EdgeL->Set("NodeAggMatrix",NodeOp);
        }
      }

      // Get the importer if we have one (for repartitioning)
      // This will get used in ReitzingerPFactory
      if(NodeL->IsAvailable("Importer")) {
        auto importer = NodeL->Get<RCP<const Import> >("Importer");
        EdgeL->Set("NodeImporter",importer);
      }
    }// end Hierarchy22 loop


    ////////////////////////////////////////////////////////////////////////////////
    // Generating the (1,1) Hierarchy
    std::string fine_smoother = precList11_.get<std::string>("smoother: type");
    {
      SM_Matrix_->setObjectLabel("A(1,1)");
      precList11_.set("coarse: max size",1);
      precList11_.set("max levels",Hierarchy22_->GetNumLevels());

      const bool refmaxwellCoarseSolve = (precList11_.get<std::string>("coarse: type",
                                                                       MasterList::getDefault<std::string>("coarse: type")) == "RefMaxwell");
      if (refmaxwellCoarseSolve) {
        GetOStream(Runtime0) << "Maxwell1::compute(): Will set up RefMaxwell coarse solver" << std::endl;
        precList11_.set("coarse: type", "none");
      }

      // Rip off non-serializable data before validation
      Teuchos::ParameterList nonSerialList11, processedPrecList11;
      MueLu::ExtractNonSerializableData(precList11_, processedPrecList11, nonSerialList11);
      RCP<HierarchyManager<SC,LO,GO,NO> > mueLuFactory = rcp(new ParameterListInterpreter<SC,LO,GO,NO>(processedPrecList11,SM_Matrix_->getDomainMap()->getComm()));
      Hierarchy11_->setlib(SM_Matrix_->getDomainMap()->lib());
      Hierarchy11_->SetProcRankVerbose(SM_Matrix_->getDomainMap()->getComm()->getRank());
      // Stick the non-serializible data on the hierarchy.
      HierarchyUtils<SC,LO,GO,NO>::AddNonSerializableDataToHierarchy(*mueLuFactory,*Hierarchy11_, nonSerialList11);
      mueLuFactory->SetupHierarchy(*Hierarchy11_);

      if (refmaxwellCoarseSolve) {
        GetOStream(Runtime0) << "Maxwell1::compute(): Setting up RefMaxwell coarse solver" << std::endl;
        auto coarseLvl = Hierarchy11_->GetLevel(Hierarchy11_->GetNumLevels()-1);
        auto coarseSolver = rcp(new MueLu::RefMaxwellSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node>("RefMaxwell",
                                                                                                      precList11_.sublist("coarse: params")));
        coarseSolver->Setup(*coarseLvl);
        coarseLvl->Set("PreSmoother",
                       rcp_dynamic_cast<SmootherBase>(coarseSolver, true));
      }
      
      if(mode_ == MODE_REFMAXWELL) {
        if(Hierarchy11_->GetNumLevels() > 1) {
          RCP<Level> EdgeL = Hierarchy11_->GetLevel(1);
          P11_ = EdgeL->Get<RCP<Matrix> >("P");
        }
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Allocate temporary MultiVectors for solve (only needed for RefMaxwell style)
    allocateMemory(1);

    describe(GetOStream(Runtime0));


#ifdef MUELU_MAXWELL1_DEBUG
    for(int i=0; i<Hierarchy11_->GetNumLevels(); i++) {
      RCP<Level> L = Hierarchy11_->GetLevel(i);
      RCP<Matrix> EdgeMatrix = rcp_dynamic_cast<Matrix>(L->Get<RCP<Operator> >("A"));
      RCP<Matrix> NodeMatrix = rcp_dynamic_cast<Matrix>(L->Get<RCP<Operator> >("NodeMatrix"));
      RCP<Matrix> NodeAggMatrix = rcp_dynamic_cast<Matrix>(L->Get<RCP<Operator> >("NodeAggMatrix"));
      RCP<Matrix> D0         =rcp_dynamic_cast<Matrix>( L->Get<RCP<Operator> >("D0"));
      
      auto nrmE = EdgeMatrix->getFrobeniusNorm();
      auto nrmN = NodeMatrix->getFrobeniusNorm();
      auto nrmNa = NodeAggMatrix->getFrobeniusNorm();
      auto nrmD0= D0->getFrobeniusNorm();

      std::cout<<"DEBUG: Norms on Level "<<i<<" E/N/NA/D0 = "<<nrmE<<" / "<<nrmN <<" / "<<nrmNa<<" / "<< nrmD0 <<std::endl;
      std::cout<<"DEBUG: NNZ on Level    "<<i<<" E/N/NA/D0 = "<< 
        EdgeMatrix->getGlobalNumEntries()<<" / "<<
        NodeMatrix->getGlobalNumEntries()<<" / "<<
        NodeAggMatrix->getGlobalNumEntries()<<" / "<<
        D0->getGlobalNumEntries()<<std::endl;      
    }
#endif

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::generate_kn() const {
    // NOTE: This does not nicely support reuse, but the relevant code can be copied from
    // RefMaxwell when we decide we want to do this.

    // NOTE: Boundary conditions OAZ are handled via the "rap: fix zero diagonals threshold"
    RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu Maxwell1: Build Kn");
    
    Level fineLevel, coarseLevel;
    fineLevel.SetFactoryManager(null);
    coarseLevel.SetFactoryManager(null);
    coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
    fineLevel.SetLevelID(0);
    coarseLevel.SetLevelID(1);
    fineLevel.Set("A",SM_Matrix_);
    coarseLevel.Set("P",D0_Matrix_);
    //coarseLevel.Set("Coordinates",Coords_);
    
    coarseLevel.setlib(SM_Matrix_->getDomainMap()->lib());
    fineLevel.setlib(SM_Matrix_->getDomainMap()->lib());
    coarseLevel.setObjectLabel("Maxwell1 (2,2)");
    fineLevel.setObjectLabel("Maxwell1 (2,2)");
    
    RCP<RAPFactory> rapFact = rcp(new RAPFactory());
    ParameterList rapList = *(rapFact->GetValidParameterList());
    rapList.set("transpose: use implicit", true);
    rapList.set("rap: triple product", parameterList_.get<bool>("rap: triple product", false));
    rapFact->SetParameterList(rapList);
    coarseLevel.Request("A", rapFact.get());
    if (enable_reuse_) {
      coarseLevel.Request("AP reuse data", rapFact.get());
      coarseLevel.Request("RAP reuse data", rapFact.get());
    }
    
    RCP<Matrix> Kn_Matrix = coarseLevel.Get< RCP<Matrix> >("A", rapFact.get());
    Kn_Matrix->setObjectLabel("A(2,2)");

    return Kn_Matrix;
  }



  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::allocateMemory(int numVectors) const {
    if(mode_ == MODE_REFMAXWELL) {
      RCP<Teuchos::TimeMonitor> tmAlloc = getTimer("MueLu Maxwell1: Allocate MVs");

      residualFine_ = MultiVectorFactory::Build(SM_Matrix_->getRangeMap(), numVectors);
     
      if(!Hierarchy11_.is_null() && Hierarchy11_->GetNumLevels() > 1) {
        RCP<Level> EdgeL = Hierarchy11_->GetLevel(1);
        RCP<Matrix> A = EdgeL->Get<RCP<Matrix> >("A");
        residual11c_ = MultiVectorFactory::Build(A->getRangeMap(), numVectors);
        update11c_   = MultiVectorFactory::Build(A->getDomainMap(), numVectors);
      }
      
      if(!Hierarchy22_.is_null()) {
        residual22_ = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(), numVectors);
        update22_   = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(), numVectors);
      }

    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const Matrix& A, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      Xpetra::IO<SC, LO, GO, NO>::Write(name, A);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const MultiVector& X, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      Xpetra::IO<SC, LO, GO, NO>::Write(name, X);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dumpCoords(const RealValuedMultiVector& X, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      Xpetra::IO<coordinateType, LO, GO, NO>::Write(name, X);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const Teuchos::ArrayRCP<bool>& v, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      std::ofstream out(name);
      for (size_t i = 0; i < Teuchos::as<size_t>(v.size()); i++)
        out << v[i] << "\n";
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const Kokkos::View<bool*, typename Node::device_type>& v, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      std::ofstream out(name);
      auto vH = Kokkos::create_mirror_view (v);
          Kokkos::deep_copy(vH , v);
          for (size_t i = 0; i < vH.size(); i++)
            out << vH[i] << "\n";
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Teuchos::TimeMonitor> Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTimer(std::string name, RCP<const Teuchos::Comm<int> > comm) const {
    if (IsPrint(Timings)) {
      if (!syncTimers_)
        return Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));
      else {
        if (comm.is_null())
          return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name), SM_Matrix_->getRowMap()->getComm().ptr()));
        else
          return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name), comm.ptr()));
      }
    } else
      return Teuchos::null;
  }




  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetMatrix(RCP<Matrix> SM_Matrix_new, bool ComputePrec) {
    bool reuse = !SM_Matrix_.is_null();
    SM_Matrix_ = SM_Matrix_new;
    dump(*SM_Matrix_, "SM.m");
    if (ComputePrec)
      compute(reuse);
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseRefMaxwellAdditive(const MultiVector& RHS, MultiVector& X) const {
    // make sure that we have enough temporary memory
    const SC one = Teuchos::ScalarTraits<SC>::one();
    if (!allEdgesBoundary_ && X.getNumVectors() != residualFine_->getNumVectors())
      allocateMemory(X.getNumVectors());

    TEUCHOS_TEST_FOR_EXCEPTION(Hierarchy11_.is_null() || Hierarchy11_->GetNumLevels() == 0, Exceptions::RuntimeError, "(1,1) Hiearchy is null.");

    // 1) Run fine pre-smoother using Hierarchy11
    RCP<Level> Fine = Hierarchy11_->GetLevel(0);
    if (Fine->IsAvailable("PreSmoother")) {
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu Maxwell1: PreSmoother");
      RCP<SmootherBase> preSmoo = Fine->Get< RCP<SmootherBase> >("PreSmoother");
      preSmoo->Apply(X, RHS, true);
   }

    // 2) Compute residual
    {
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu Maxwell1: residual calculation");
      Utilities::Residual(*SM_Matrix_, X, RHS,*residualFine_);
    }

    // 3a) Restrict residual to (1,1) Hierarchy's level 1 and execute (1,1) hierarchy (use startLevel and InitialGuessIsZero)
    if(!P11_.is_null()) {
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu Maxwell1: (1,1) correction");
      P11_->apply(*residualFine_,*residual11c_,Teuchos::TRANS);
      Hierarchy11_->Iterate(*residual11c_,*update11c_,true,1);
    }

    // 3b) Restrict residual to (2,2) Hierarchy's level 0 and execute (2,2) hierarchy (use InitialGuessIsZero)
    if (!allNodesBoundary_) {
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu Maxwell1: (2,2) correction");
      D0_Matrix_->apply(*residualFine_,*residual22_,Teuchos::TRANS);
      Hierarchy22_->Iterate(*residual22_,*update22_,true,0);
    }

    // 4) Prolong both updates back into X-vector (Need to do both the P11 null and not null cases
    {
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu Maxwell1: Orolongation");
      if(!P11_.is_null())
        P11_->apply(*update11c_,X,Teuchos::NO_TRANS,one,one);
      if (!allNodesBoundary_)
        D0_Matrix_->apply(*update22_,X,Teuchos::NO_TRANS,one,one);
    }

    // 5) Run fine post-smoother using Hierarchy11
    if (Fine->IsAvailable("PostSmoother")) {
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu Maxwell1: PostSmoother");
      RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
      postSmoo->Apply(X, RHS, false);
    }


  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseStandard(const MultiVector& RHS, MultiVector& X) const {
    Hierarchy11_->Iterate(RHS,X,1,true);
  }
 
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply (const MultiVector& RHS, MultiVector& X,
                                                                Teuchos::ETransp /* mode */,
                                                                  Scalar /* alpha */,
                                                                  Scalar /* beta */) const {
    RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu Maxwell1: solve");
    if(mode_ == MODE_GMHD_STANDARD)
      HierarchyGmhd_->Iterate(RHS,X,1,true);
    else if(mode_ == MODE_STANDARD || mode_ == MODE_EDGE_ONLY)
      applyInverseStandard(RHS,X);
    else if(mode_ == MODE_REFMAXWELL)
      applyInverseRefMaxwellAdditive(RHS,X);
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Must use mode 'standard', 'refmaxwell' or 'edge only' when not doing GMHD.");
  }


 
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasTransposeApply() const {
    return false;
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  initialize(const Teuchos::RCP<Matrix> & D0_Matrix,
             const Teuchos::RCP<Matrix> & Kn_Matrix,
             const Teuchos::RCP<MultiVector>  & Nullspace,
             const Teuchos::RCP<RealValuedMultiVector>  & Coords,
             Teuchos::ParameterList& List)
  {
    // some pre-conditions
    TEUCHOS_ASSERT(D0_Matrix!=Teuchos::null);

#ifdef HAVE_MUELU_DEBUG
    if(!Kn_Matrix.is_null()) {
      TEUCHOS_ASSERT(Kn_Matrix->getDomainMap()->isSameAs(*D0_Matrix->getDomainMap()));
      TEUCHOS_ASSERT(Kn_Matrix->getRangeMap()->isSameAs(*D0_Matrix->getDomainMap()));
    }


    TEUCHOS_ASSERT(D0_Matrix->getRangeMap()->isSameAs(*D0_Matrix->getRowMap()));
#endif

    Hierarchy11_   = Teuchos::null;
    Hierarchy22_   = Teuchos::null;
    HierarchyGmhd_  = Teuchos::null;
    if (mode_ != MODE_GMHD_STANDARD) mode_ = MODE_STANDARD;

    // Default settings
    useKokkos_=false;
    allEdgesBoundary_=false;
    allNodesBoundary_=false;
    dump_matrices_ = false;
    enable_reuse_=false;
    syncTimers_=false;
    applyBCsTo22_ = false;

    
    // set parameters
    setParameters(List);

    if (D0_Matrix->getRowMap()->lib() == Xpetra::UseTpetra) {
      // We will remove boundary conditions from D0, and potentially change maps, so we copy the input.
      // Fortunately, D0 is quite sparse.
      // We cannot use the Tpetra copy constructor, since it does not copy the graph.

      RCP<Matrix> D0copy = MatrixFactory::Build(D0_Matrix->getRowMap(), D0_Matrix->getColMap(), 0);
      RCP<CrsMatrix> D0copyCrs = rcp_dynamic_cast<CrsMatrixWrap>(D0copy,true)->getCrsMatrix();
      ArrayRCP<const size_t> D0rowptr_RCP;
      ArrayRCP<const LO>     D0colind_RCP;
      ArrayRCP<const SC>     D0vals_RCP;
      rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix,true)->getCrsMatrix()->getAllValues(D0rowptr_RCP,
                                                                                    D0colind_RCP,
                                                                                    D0vals_RCP);

      ArrayRCP<size_t> D0copyrowptr_RCP;
      ArrayRCP<LO>     D0copycolind_RCP;
      ArrayRCP<SC>     D0copyvals_RCP;
      D0copyCrs->allocateAllValues(D0vals_RCP.size(),D0copyrowptr_RCP,D0copycolind_RCP,D0copyvals_RCP);
      D0copyrowptr_RCP.deepCopy(D0rowptr_RCP());
      D0copycolind_RCP.deepCopy(D0colind_RCP());
      D0copyvals_RCP.deepCopy(D0vals_RCP());
      D0copyCrs->setAllValues(D0copyrowptr_RCP,
                              D0copycolind_RCP,
                              D0copyvals_RCP);
      D0copyCrs->expertStaticFillComplete(D0_Matrix->getDomainMap(), D0_Matrix->getRangeMap(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix,true)->getCrsMatrix()->getCrsGraph()->getImporter(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix,true)->getCrsMatrix()->getCrsGraph()->getExporter());
      D0_Matrix_ = D0copy;
    } else
      D0_Matrix_ = MatrixFactory::BuildCopy(D0_Matrix);


    Kn_Matrix_    = Kn_Matrix;
    Coords_       = Coords;
    Nullspace_    = Nullspace;

    if(!Kn_Matrix_.is_null()) dump(*Kn_Matrix_, "Kn.m");
    if (!Nullspace_.is_null())    dump(*Nullspace_, "nullspace.m");
    if (!Coords_.is_null())       dumpCoords(*Coords_, "coords.m");

  }



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell1<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel /* verbLevel */) const {

    std::ostringstream oss;

    RCP<const Teuchos::Comm<int> > comm = SM_Matrix_->getDomainMap()->getComm();

#ifdef HAVE_MPI
    int root;
    if (!Kn_Matrix_.is_null())
      root = comm->getRank();
    else
      root = -1;

    int actualRoot;
    reduceAll(*comm, Teuchos::REDUCE_MAX, root, Teuchos::ptr(&actualRoot));
    root = actualRoot;
#endif


    oss << "\n--------------------------------------------------------------------------------\n" <<
      "---                            Maxwell1 Summary                            ---\n"
      "--------------------------------------------------------------------------------" << std::endl;
    oss << std::endl;

    GlobalOrdinal numRows;
    GlobalOrdinal nnz;

    SM_Matrix_->getRowMap()->getComm()->barrier();

    numRows = SM_Matrix_->getGlobalNumRows();
    nnz = SM_Matrix_->getGlobalNumEntries();

    Xpetra::global_size_t tt = numRows;
    int rowspacer = 3; while (tt != 0) { tt /= 10; rowspacer++; }
    tt = nnz;
    int nnzspacer = 2; while (tt != 0) { tt /= 10; nnzspacer++; }

    oss  << "block " << std::setw(rowspacer) << " rows " << std::setw(nnzspacer) << " nnz " << std::setw(9) << " nnz/row" << std::endl;
    oss << "(1, 1)" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;

    if (!Kn_Matrix_.is_null()) {
      // ToDo: make sure that this is printed correctly
      numRows = Kn_Matrix_->getGlobalNumRows();
      nnz = Kn_Matrix_->getGlobalNumEntries();

      oss << "(2, 2)" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;
    }

    oss << std::endl;


    std::string outstr = oss.str();

#ifdef HAVE_MPI
    RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    MPI_Comm rawComm = (*mpiComm->getRawMpiComm())();

    int strLength = outstr.size();
    MPI_Bcast(&strLength, 1, MPI_INT, root, rawComm);
    if (comm->getRank() != root)
      outstr.resize(strLength);
    MPI_Bcast(&outstr[0], strLength, MPI_CHAR, root, rawComm);
#endif

    out << outstr;

    if (!Hierarchy11_.is_null())
      Hierarchy11_->describe(out, GetVerbLevel());

    if (!Hierarchy22_.is_null())
      Hierarchy22_->describe(out, GetVerbLevel());

    if (!HierarchyGmhd_.is_null())
      HierarchyGmhd_->describe(out, GetVerbLevel());

    if (IsPrint(Statistics2)) {
      // Print the grid of processors
      std::ostringstream oss2;

      oss2 << "Sub-solver distribution over ranks" << std::endl;
      oss2 << "( (1,1) block only is indicated by '1', (2,2) block only by '2', and both blocks by 'B' and none by '.')" << std::endl;

      int numProcs = comm->getSize();
#ifdef HAVE_MPI
      RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);

      RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();
#endif

      char status = 0;
      if (!Kn_Matrix_.is_null())
        status += 1;
      std::vector<char> states(numProcs, 0);
#ifdef HAVE_MPI
      MPI_Gather(&status, 1, MPI_CHAR, &states[0], 1, MPI_CHAR, 0, *rawMpiComm);
#else
      states.push_back(status);
#endif

      int rowWidth = std::min(Teuchos::as<int>(ceil(sqrt(numProcs))), 100);
      for (int proc = 0; proc < numProcs; proc += rowWidth) {
        for (int j = 0; j < rowWidth; j++)
          if (proc + j < numProcs)
            if (states[proc+j] == 0)
              oss2 << ".";
            else if (states[proc+j] == 1)
              oss2 << "1";
            else if (states[proc+j] == 2)
              oss2 << "2";
            else
              oss2 << "B";
          else
            oss2 << " ";

        oss2 << "      " << proc << ":" << std::min(proc + rowWidth, numProcs) - 1 << std::endl;
      }
      oss2 << std::endl;
      GetOStream(Statistics2) << oss2.str();
    }


  }



} // namespace

#define MUELU_MAXWELL1_SHORT
#endif //ifdef MUELU_MAXWELL1_DEF_HPP
