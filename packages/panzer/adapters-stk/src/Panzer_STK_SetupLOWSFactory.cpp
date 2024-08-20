// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_STK_SetupLOWSFactory.hpp"
#include "Panzer_STK_ParameterListCallback.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"

#include "Teuchos_AbstractFactoryStd.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "EpetraExt_VectorOut.h"
#endif // PANZER_HAVE_EPETRA_STACK

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#ifdef PANZER_HAVE_TEKO
#include "Teko_StratimikosFactory.hpp"
#endif

#ifdef PANZER_HAVE_MUELU
#include "Stratimikos_MueLuHelpers.hpp"
//#include "MatrixMarket_Tpetra.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#endif

namespace panzer_stk {

namespace {

  bool
  determineCoordinateField(const panzer::GlobalIndexer & globalIndexer,std::string & fieldName)
  {
    std::vector<std::string> elementBlocks;
    globalIndexer.getElementBlockIds(elementBlocks);

    // grab fields for first block
    std::set<int> runningFields;
    {
      const std::vector<int> & fields = globalIndexer.getBlockFieldNumbers(elementBlocks[0]);
      runningFields.insert(fields.begin(),fields.end());
    }

    // grab fields for first block
    for(std::size_t i=1;i<elementBlocks.size();i++) {
      const std::vector<int> & fields = globalIndexer.getBlockFieldNumbers(elementBlocks[i]);

      std::set<int> currentFields(runningFields);
      runningFields.clear();
      std::set_intersection(fields.begin(),fields.end(),
                            currentFields.begin(),currentFields.end(),
                            std::inserter(runningFields,runningFields.begin()));
    }

    if(runningFields.size()<1)
      return false;

    fieldName = globalIndexer.getFieldString(*runningFields.begin());
    return true;
  }

  void
  fillFieldPatternMap(const panzer::DOFManager & globalIndexer,
                      const std::string & fieldName,
                      std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > & fieldPatterns)
  {
     std::vector<std::string> elementBlocks;
     globalIndexer.getElementBlockIds(elementBlocks);

     for(std::size_t e=0;e<elementBlocks.size();e++) {
        std::string blockId = elementBlocks[e];

        if(globalIndexer.fieldInBlock(fieldName,blockId))
           fieldPatterns[blockId] =
              Teuchos::rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(globalIndexer.getFieldPattern(blockId,fieldName),true);
     }
  }

  void
  fillFieldPatternMap(const panzer::GlobalIndexer & globalIndexer,
                      const std::string & fieldName,
                      std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > & fieldPatterns)
  {
    using Teuchos::Ptr;
    using Teuchos::ptrFromRef;
    using Teuchos::ptr_dynamic_cast;
    using panzer::DOFManager;

    // first standard dof manager
    {
      Ptr<const DOFManager> dofManager = ptr_dynamic_cast<const DOFManager>(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
  }
} // end anonymous namespace

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer,
                   const Teuchos::RCP<panzer_stk::STKConnManager> & stkConn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo,
                   const Teuchos::RCP<const panzer::GlobalIndexer> & auxGlobalIndexer,
                   bool useCoordinates
                   )
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    // Note if you want to use new solvers within Teko they have to be added to the solver builer
    // before teko is added. This is because Teko steals its defaults from the solver its being injected
    // into!

    #ifdef PANZER_HAVE_MUELU
    {
      Stratimikos::enableMueLu<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(linearSolverBuilder,"MueLu");
      Stratimikos::enableMueLuRefMaxwell<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(linearSolverBuilder,"MueLuRefMaxwell");
      Stratimikos::enableMueLuMaxwell1<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(linearSolverBuilder,"MueLuMaxwell1");
      #ifndef PANZER_HIDE_DEPRECATED_CODE
      // the next two are only for backwards compatibility
      Stratimikos::enableMueLu<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(linearSolverBuilder,"MueLu-Tpetra");
      Stratimikos::enableMueLuRefMaxwell<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(linearSolverBuilder,"MueLuRefMaxwell-Tpetra");
      #endif
    }
    #endif // MUELU

    #ifdef PANZER_HAVE_TEKO
    RCP<Teko::RequestHandler> reqHandler_local = reqHandler;

    if(!blockedAssembly) {

       std::string fieldName;

       // try to set request handler from member variable; This is a potential segfault
       // if its internally stored data (e.g. callback) gets released and all its data
       // required by ML or whatever gets hosed
       if(reqHandler_local==Teuchos::null)
          reqHandler_local = rcp(new Teko::RequestHandler);

       // add in the coordinate parameter list callback handler
       if(determineCoordinateField(*globalIndexer,fieldName)) {
          std::map<std::string,RCP<const panzer::Intrepid2FieldPattern> > fieldPatterns;
          fillFieldPatternMap(*globalIndexer,fieldName,fieldPatterns);

          RCP<panzer_stk::ParameterListCallback> callback = rcp(new
                panzer_stk::ParameterListCallback(fieldName,fieldPatterns,stkConn_manager,
                rcp_dynamic_cast<const panzer::GlobalIndexer>(globalIndexer)));
          reqHandler_local->addRequestCallback(callback);

          // determine if you want rigid body null space modes...currently an extremely specialized case!
          if(strat_params->sublist("Preconditioner Types").isSublist("ML")) {
/*           COMMENTING THIS OUT FOR NOW, this is causing problems with some of the preconditioners in optimization...not sure why

             Teuchos::ParameterList & ml_params = strat_params->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");

             {
                // force parameterlistcallback to build coordinates
                callback->preRequest(Teko::RequestMesg(rcp(new Teuchos::ParameterList())));

                // extract coordinate vectors
                std::vector<double> & xcoords = const_cast<std::vector<double> & >(callback->getXCoordsVector());
                std::vector<double> & ycoords = const_cast<std::vector<double> & >(callback->getYCoordsVector());
                std::vector<double> & zcoords = const_cast<std::vector<double> & >(callback->getZCoordsVector());

                ml_params.set<double*>("x-coordinates",&xcoords[0]);
                ml_params.set<double*>("y-coordinates",&ycoords[0]);
                ml_params.set<double*>("z-coordinates",&zcoords[0]);
             }
*/
/*
             bool useRigidBodyNullSpace = false;
             if(ml_params.isType<std::string>("null space: type"))
               useRigidBodyNullSpace = ml_params.get<std::string>("null space: type") == "pre-computed";

             if(useRigidBodyNullSpace) {
                // force parameterlistcallback to build coordinates
                callback->preRequest(Teko::RequestMesg(rcp(new Teuchos::ParameterList())));

                RCP<std::vector<double> > rbm = rcp(new std::vector<double>);
                std::vector<double> & rbm_ref = *rbm;

                // extract coordinate vectors
                std::vector<double> & xcoords = const_cast<std::vector<double> & >(callback->getXCoordsVector());
                std::vector<double> & ycoords = const_cast<std::vector<double> & >(callback->getYCoordsVector());
                std::vector<double> & zcoords = const_cast<std::vector<double> & >(callback->getZCoordsVector());

                // use ML to build the null space modes for ML
                int Nnodes     = Teuchos::as<int>(xcoords.size());
                int NscalarDof = 0;
                int Ndof       = spatialDim;
                int nRBM       = spatialDim==3 ? 6 : (spatialDim==2 ? 3 : 1);
                int rbmSize    = Nnodes*(nRBM+NscalarDof)*(Ndof+NscalarDof);
                rbm_ref.resize(rbmSize);

                ML_Coord2RBM(Nnodes,&xcoords[0],&ycoords[0],&zcoords[0],&rbm_ref[0],Ndof,NscalarDof);

                ml_params.set<double*>("null space: vectors",&rbm_ref[0]);
                ml_params.set<int>("null space: dimension",nRBM);

                callback->storeExtraVector(rbm);
             }
*/
          }

          if(writeCoordinates) {
#ifdef PANZER_HAVE_EPETRA_STACK
             // force parameterlistcallback to build coordinates
             callback->preRequest(Teko::RequestMesg(rcp(new Teuchos::ParameterList())));

             // extract coordinate vectors
             const std::vector<double> & xcoords = callback->getXCoordsVector();
             const std::vector<double> & ycoords = callback->getYCoordsVector();
             const std::vector<double> & zcoords = callback->getZCoordsVector();

             // use epetra to write coordinates to matrix market files
             Epetra_MpiComm ep_comm(*mpi_comm->getRawMpiComm()); // this is OK access to RawMpiComm becase its declared on the stack?
                                                                 // and all users of this object are on the stack (within scope of mpi_comm
             Epetra_Map map(-1,xcoords.size(),0,ep_comm);

             RCP<Epetra_Vector> vec;
             switch(spatialDim) {
             case 3:
                vec = rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("zcoords.mm",*vec);
                // Intentional fall-through.
             case 2:
                vec = rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
                EpetraExt::VectorToMatrixMarketFile("ycoords.mm",*vec);
                // Intentional fall-through.
             case 1:
                vec = rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("xcoords.mm",*vec);
                break;
             default:
                TEUCHOS_ASSERT(false);
             }
#else
             TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: Panzer_STK_SetupLOWSFactory.cpp - writeCoordinates not implemented for Tpetra yet!");
#endif
          }

#ifdef PANZER_HAVE_MUELU
          if(rcp_dynamic_cast<const panzer::GlobalIndexer>(globalIndexer)!=Teuchos::null
             && useCoordinates) {
             if(!writeCoordinates)
                callback->preRequest(Teko::RequestMesg(rcp(new Teuchos::ParameterList())));

             typedef Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> Map;
             typedef Tpetra::MultiVector<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> MV;

             // extract coordinate vectors and modify strat_params to include coordinate vectors
             unsigned dim = Teuchos::as<unsigned>(spatialDim);
             RCP<MV> coords;
             for(unsigned d=0;d<dim;d++) {
               const std::vector<double> & coord = callback->getCoordsVector(d);

               // no coords vector has been build yet, build one
               if(coords==Teuchos::null) {
                 if(globalIndexer->getNumFields()==1) {
                   RCP<const panzer::GlobalIndexer> ugi
                       = rcp_dynamic_cast<const panzer::GlobalIndexer>(globalIndexer);
                   std::vector<panzer::GlobalOrdinal> ownedIndices;
                   ugi->getOwnedIndices(ownedIndices);
                   RCP<const Map> coords_map = rcp(new Map(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),ownedIndices,0,mpi_comm));
                   coords = rcp(new MV(coords_map,dim));
                 }
                 else {
                   RCP<const Map> coords_map = rcp(new Map(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),coord.size(),0,mpi_comm));
                   coords = rcp(new MV(coords_map,dim));
                 }
               }

               // sanity check the size
               TEUCHOS_ASSERT(coords->getLocalLength()==coord.size());

               // fill appropriate coords vector
               Teuchos::ArrayRCP<double> dest = coords->getDataNonConst(d);
               for(std::size_t i=0;i<coord.size();i++)
                 dest[i] = coord[i];
             }

             // inject coordinates into parameter list
             Teuchos::ParameterList & muelu_params = strat_params->sublist("Preconditioner Types").sublist("MueLu");
             muelu_params.set<RCP<MV> >("Coordinates",coords);
          }
          #endif
       }
       // else write_out_the_mesg("Warning: No unique field determines the coordinates, coordinates unavailable!")

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler_local);
    }
    else {
       // try to set request handler from member variable
       if(reqHandler_local==Teuchos::null)
          reqHandler_local = rcp(new Teko::RequestHandler);

       std::string fieldName;
       if(determineCoordinateField(*globalIndexer,fieldName)) {
          RCP<const panzer::BlockedDOFManager> blkDofs =
             rcp_dynamic_cast<const panzer::BlockedDOFManager>(globalIndexer);
          RCP<const panzer::BlockedDOFManager> auxBlkDofs =
             rcp_dynamic_cast<const panzer::BlockedDOFManager>(auxGlobalIndexer);
          RCP<panzer_stk::ParameterListCallbackBlocked> callback =
                rcp(new panzer_stk::ParameterListCallbackBlocked(stkConn_manager,blkDofs,auxBlkDofs));
          reqHandler_local->addRequestCallback(callback);
       }

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler_local);

       if(writeCoordinates) {
#ifdef PANZER_HAVE_EPETRA_STACK
          RCP<const panzer::BlockedDOFManager> blkDofs =
             rcp_dynamic_cast<const panzer::BlockedDOFManager>(globalIndexer);

          // loop over blocks
          const std::vector<RCP<panzer::GlobalIndexer>> & dofVec
             = blkDofs->getFieldDOFManagers();
          for(std::size_t i=0;i<dofVec.size();i++) {

            // add in the coordinate parameter list callback handler
            TEUCHOS_ASSERT(determineCoordinateField(*dofVec[i],fieldName));

            std::map<std::string,RCP<const panzer::Intrepid2FieldPattern> > fieldPatterns;
            fillFieldPatternMap(*dofVec[i],fieldName,fieldPatterns);
            panzer_stk::ParameterListCallback plCall(fieldName,fieldPatterns,stkConn_manager,dofVec[i]);
            plCall.buildArrayToVector();
            plCall.buildCoordinates();

            // extract coordinate vectors
            const std::vector<double> & xcoords = plCall.getXCoordsVector();
            const std::vector<double> & ycoords = plCall.getYCoordsVector();
            const std::vector<double> & zcoords = plCall.getZCoordsVector();

            // use epetra to write coordinates to matrix market files
            Epetra_MpiComm ep_comm(*mpi_comm->getRawMpiComm()); // this is OK access to RawMpiComm becase its declared on the stack?
                                                                // and all users of this object are on the stack (within scope of mpi_comm
            Epetra_Map map(-1,xcoords.size(),0,ep_comm);

            RCP<Epetra_Vector> vec;
            switch(spatialDim) {
            case 3:
               vec = rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_zcoords.mm").c_str(),*vec);
               // Intentional fall-through.
            case 2:
               vec = rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_ycoords.mm").c_str(),*vec);
               // Intentional fall-through.
            case 1:
               vec = rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_xcoords.mm").c_str(),*vec);
               break;
            default:
               TEUCHOS_ASSERT(false);
            }

            // TODO add MueLu code...
            #ifdef PANZER_HAVE_MUELU
            if(useCoordinates) {

              typedef Xpetra::Map<int,panzer::GlobalOrdinal> Map;
              typedef Xpetra::MultiVector<double,int,panzer::GlobalOrdinal> MV;

              // TODO This is Epetra-specific
              RCP<const Map> coords_map = Xpetra::MapFactory<int,panzer::GlobalOrdinal>::Build(Xpetra::UseEpetra,
                                            Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),
                  //Teuchos::ArrayView<GO>(ownedIndices),
                  xcoords.size(),
                  0,
                  mpi_comm
              );

              unsigned dim = Teuchos::as<unsigned>(spatialDim);

              RCP<MV> coords = Xpetra::MultiVectorFactory<double,int,panzer::GlobalOrdinal>::Build(coords_map,spatialDim);

              for(unsigned d=0;d<dim;d++) {
                // sanity check the size
                TEUCHOS_ASSERT(coords->getLocalLength()==xcoords.size());

                // fill appropriate coords vector
                Teuchos::ArrayRCP<double> dest = coords->getDataNonConst(d);
                for(std::size_t j=0;j<coords->getLocalLength();++j) {
                  if (d == 0) dest[j] = xcoords[j];
                  if (d == 1) dest[j] = ycoords[j];
                  if (d == 2) dest[j] = zcoords[j];
                }
              }

              // TODO This is Epetra-specific
              // inject coordinates into parameter list
              Teuchos::ParameterList & muelu_params = strat_params->sublist("Preconditioner Types").sublist("MueLu");
              muelu_params.set<RCP<MV> >("Coordinates",coords);

            }
            #endif

          } /* end loop over all physical fields */
#else // PANZER_HAVE EPETRA
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: Panzer_STK_SetupLOWSFactory - writeCoordinates not implemented for Tpetra yet!")
#endif
       }

       if(writeTopo) {
          /*
          RCP<const panzer::BlockedDOFManager> blkDofs =
             rcp_dynamic_cast<const panzer::BlockedDOFManager>(globalIndexer);

          writeTopology(*blkDofs);
          */
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                     "Topology writing is no longer implemented. It needs to be reimplemented for the "
                                     "default DOFManager (estimate 2 days with testing)");
       }
    }
    #endif

    linearSolverBuilder.setParameterList(strat_params);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    return lowsFactory;
  }


  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer,
                   const Teuchos::RCP<panzer::ConnManager> & conn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo,
                   const Teuchos::RCP<const panzer::GlobalIndexer> & auxGlobalIndexer,
                   bool useCoordinates
                   )
  {
    #ifdef PANZER_HAVE_TEKO
    Teuchos::RCP<Teko::RequestHandler> reqHandler_local = reqHandler;
    if(reqHandler_local==Teuchos::null)
      reqHandler_local = Teuchos::rcp(new Teko::RequestHandler);
    #endif

    auto stk_conn_manager = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager>(conn_manager,true);

    return buildLOWSFactory(blockedAssembly,globalIndexer,stk_conn_manager,spatialDim,mpi_comm,strat_params,
#ifdef PANZER_HAVE_TEKO
                            reqHandler_local,
#endif
                            writeCoordinates,
                            writeTopo,
                            auxGlobalIndexer,
                            useCoordinates
                            );

    // should never reach this
    TEUCHOS_ASSERT(false);
    return Teuchos::null;
  }
}
