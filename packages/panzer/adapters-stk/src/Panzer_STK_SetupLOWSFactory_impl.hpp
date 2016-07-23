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

#ifndef __Panzer_STK_SetupLOWSFactory_impl_hpp__
#define __Panzer_STK_SetupLOWSFactory_impl_hpp__

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_ParameterListCallback.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"

#include "Teuchos_AbstractFactoryStd.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "EpetraExt_VectorOut.h"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#ifdef PANZER_HAVE_TEKO
#include "Teko_StratimikosFactory.hpp"
#endif

#ifdef PANZER_HAVE_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#include "MatrixMarket_Tpetra.hpp"
#endif

#ifdef PANZER_HAVE_IFPACK2
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif

namespace panzer_stk {

namespace {

  bool 
  determineCoordinateField(const panzer::UniqueGlobalIndexerBase & globalIndexer,std::string & fieldName)
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

#ifdef PANZER_HAVE_FEI
  template<typename GO>
  void 
  fillFieldPatternMap(const panzer::DOFManagerFEI<int,GO> & globalIndexer,
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
#endif

  template<typename GO>
  void 
  fillFieldPatternMap(const panzer::DOFManager<int,GO> & globalIndexer,
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
  fillFieldPatternMap(const panzer::UniqueGlobalIndexerBase & globalIndexer,
                      const std::string & fieldName,
                      std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > & fieldPatterns)
  {
    using Teuchos::Ptr;
    using Teuchos::ptrFromRef;
    using Teuchos::ptr_dynamic_cast;
    using panzer::DOFManager;
#ifdef PANZER_HAVE_FEI
    using panzer::DOFManagerFEI;
#endif

    // first standard dof manager
    {
      Ptr<const DOFManager<int,int> > dofManager = ptr_dynamic_cast<const DOFManager<int,int> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
    {
      Ptr<const DOFManager<int,panzer::Ordinal64> > dofManager = ptr_dynamic_cast<const DOFManager<int,panzer::Ordinal64> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }

#ifdef PANZER_HAVE_FEI
    // now FEI dof manager
    {
      Ptr<const DOFManagerFEI<int,int> > dofManager = ptr_dynamic_cast<const DOFManagerFEI<int,int> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
    {
      Ptr<const DOFManagerFEI<int,panzer::Ordinal64> > dofManager = ptr_dynamic_cast<const DOFManagerFEI<int,panzer::Ordinal64> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
#endif
  }
}

  template<typename GO>
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer_stk::STKConnManager<GO> > & stkConn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo
                   )
  {
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    // Note if you want to use new solvers within Teko they have to be added to the solver builer
    // before teko is added. This is because Teko steals its defaults from the solver its being injected
    // into!

    #ifdef PANZER_HAVE_MUELU
    {
      // TAW: the following is probably not optimal but it corresponds to what have been there before...
      Stratimikos::enableMueLu(linearSolverBuilder,"MueLu");
      Stratimikos::enableMueLu<int,panzer::Ordinal64,panzer::TpetraNodeType>(linearSolverBuilder,"MueLu-Tpetra");
    }
    #endif // MUELU
    #ifdef PANZER_HAVE_IFPACK2
    {
      typedef Thyra::PreconditionerFactoryBase<double> Base;
      typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double, int, panzer::Ordinal64,panzer::TpetraNodeType> > Impl;

      linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
    }
    #endif // MUELU


    #ifdef PANZER_HAVE_TEKO
    Teuchos::RCP<Teko::RequestHandler> reqHandler_local = reqHandler;

    if(!blockedAssembly) {

       std::string fieldName;

       // try to set request handler from member variable
       if(reqHandler_local==Teuchos::null)
          reqHandler_local = Teuchos::rcp(new Teko::RequestHandler);

       // add in the coordinate parameter list callback handler
       if(determineCoordinateField(*globalIndexer,fieldName)) {
          std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > fieldPatterns;
          fillFieldPatternMap(*globalIndexer,fieldName,fieldPatterns);

          Teuchos::RCP<panzer_stk::ParameterListCallback<int,GO> > callback = Teuchos::rcp(new
                panzer_stk::ParameterListCallback<int,GO>(fieldName,fieldPatterns,stkConn_manager,
                Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,GO> >(globalIndexer)));
          reqHandler_local->addRequestCallback(callback);

          if(writeCoordinates) {
             // force parameterlistcallback to build coordinates
             callback->preRequest(Teko::RequestMesg(Teuchos::rcp(new Teuchos::ParameterList())));

             // extract coordinate vectors
             const std::vector<double> & xcoords = callback->getXCoordsVector();
             const std::vector<double> & ycoords = callback->getYCoordsVector();
             const std::vector<double> & zcoords = callback->getZCoordsVector();

             // use epetra to write coordinates to matrix market files
             Epetra_MpiComm ep_comm(*mpi_comm->getRawMpiComm()); // this is OK access to RawMpiComm becase its declared on the stack?
                                                                 // and all users of this object are on the stack (within scope of mpi_comm
             Epetra_Map map(-1,xcoords.size(),0,ep_comm);

             Teuchos::RCP<Epetra_Vector> vec;
             switch(spatialDim) {
             case 3:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("zcoords.mm",*vec);
             case 2:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
                EpetraExt::VectorToMatrixMarketFile("ycoords.mm",*vec);
             case 1:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("xcoords.mm",*vec);
                break;
             default:
                TEUCHOS_ASSERT(false);
             }
          }

          #ifdef PANZER_HAVE_MUELU
          if(Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> >(globalIndexer)!=Teuchos::null) {
             if(!writeCoordinates)
                callback->preRequest(Teko::RequestMesg(Teuchos::rcp(new Teuchos::ParameterList())));

             typedef Tpetra::Map<int,panzer::Ordinal64,panzer::TpetraNodeType> Map;
             typedef Tpetra::MultiVector<double,int,panzer::Ordinal64,panzer::TpetraNodeType> MV;

             // extract coordinate vectors and modify strat_params to include coordinate vectors
             unsigned dim = Teuchos::as<unsigned>(spatialDim);
             Teuchos::RCP<MV> coords;
             for(unsigned d=0;d<dim;d++) {
               const std::vector<double> & coord = callback->getCoordsVector(d);

               // no coords vector has been build yet, build one
               if(coords==Teuchos::null) {
                 if(globalIndexer->getNumFields()==1) {
                   Teuchos::RCP<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> > ugi
                       = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> >(globalIndexer);
                   std::vector<panzer::Ordinal64> ownedIndices;
                   ugi->getOwnedIndices(ownedIndices);
                   Teuchos::RCP<const Map> coords_map = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::Ordinal64>::invalid(),ownedIndices,0,mpi_comm));
                   coords = Teuchos::rcp(new MV(coords_map,dim));
                 }
                 else {
                   Teuchos::RCP<const Map> coords_map = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::Ordinal64>::invalid(),coord.size(),0,mpi_comm));
                   coords = Teuchos::rcp(new MV(coords_map,dim));
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
             Teuchos::ParameterList & muelu_params = strat_params->sublist("Preconditioner Types").sublist("MueLu-Tpetra");
             muelu_params.set<Teuchos::RCP<MV> >("Coordinates",coords);
          }
          #endif
       }
       // else write_out_the_mesg("Warning: No unique field determines the coordinates, coordinates unavailable!")

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler_local);
    }
    else {
       // try to set request handler from member variable
       if(reqHandler_local==Teuchos::null)
          reqHandler_local = Teuchos::rcp(new Teko::RequestHandler);

       std::string fieldName;
       if(determineCoordinateField(*globalIndexer,fieldName)) {
          Teuchos::RCP<const panzer::BlockedDOFManager<int,GO> > blkDofs =
             Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager<int,GO> >(globalIndexer);
          Teuchos::RCP<panzer_stk::ParameterListCallbackBlocked<int,GO> > callback =
                Teuchos::rcp(new panzer_stk::ParameterListCallbackBlocked<int,GO>(stkConn_manager,blkDofs));
          reqHandler_local->addRequestCallback(callback);
       }

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler_local);

       if(writeCoordinates) {
          Teuchos::RCP<const panzer::BlockedDOFManager<int,GO> > blkDofs =
             Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager<int,GO> >(globalIndexer);

          // loop over blocks
          const std::vector<Teuchos::RCP<panzer::UniqueGlobalIndexer<int,GO> > > & dofVec
             = blkDofs->getFieldDOFManagers();
          for(std::size_t i=0;i<dofVec.size();i++) {
            std::string fieldName;

            // add in the coordinate parameter list callback handler
            TEUCHOS_ASSERT(determineCoordinateField(*dofVec[i],fieldName));

            std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > fieldPatterns;
            fillFieldPatternMap(*dofVec[i],fieldName,fieldPatterns);
            panzer_stk::ParameterListCallback<int,GO> plCall(fieldName,fieldPatterns,stkConn_manager,dofVec[i]);
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

            Teuchos::RCP<Epetra_Vector> vec;
            switch(spatialDim) {
            case 3:
               vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_zcoords.mm").c_str(),*vec);
            case 2:
               vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_ycoords.mm").c_str(),*vec);
            case 1:
               vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_xcoords.mm").c_str(),*vec);
               break;
            default:
               TEUCHOS_ASSERT(false);
            }
          }
       }

       if(writeTopo) {
          Teuchos::RCP<const panzer::BlockedDOFManager<int,GO> > blkDofs =
             Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager<int,GO> >(globalIndexer);

          writeTopology(*blkDofs);
       }
    }
    #endif

    linearSolverBuilder.setParameterList(strat_params);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    return lowsFactory;
  }

  template<typename GO>
  void 
  writeTopology(const panzer::BlockedDOFManager<int,GO> & blkDofs)
  {
    using Teuchos::RCP;

    // loop over each field block
    const std::vector<RCP<panzer::UniqueGlobalIndexer<int,GO> > > & blk_dofMngrs = blkDofs.getFieldDOFManagers();
    for(std::size_t b=0;b<blk_dofMngrs.size();b++) {
#ifdef PANZER_HAVE_FEI
      RCP<panzer::DOFManagerFEI<int,GO> > dofMngr = Teuchos::rcp_dynamic_cast<panzer::DOFManagerFEI<int,GO> >(blk_dofMngrs[b],true);

      std::vector<std::string> eBlocks;
      dofMngr->getElementBlockIds(eBlocks);

      // build file name
      std::stringstream fileName;
      fileName << "elements_" << b;
      std::ofstream file(fileName.str().c_str());

      // loop over each element block, write out topology
      for(std::size_t e=0;e<eBlocks.size();e++)
        writeTopology(*dofMngr,eBlocks[e],file);
#else
      TEUCHOS_ASSERT(false);
#endif
    }
  }

#ifdef PANZER_HAVE_FEI
  template <typename GO>
  void 
  writeTopology(const panzer::DOFManagerFEI<int,GO> & dofs,const std::string & block,std::ostream & os)
  {
    std::vector<std::string> fields(dofs.getElementBlockGIDCount(block));

    const std::set<int> & fieldIds = dofs.getFields(block);
    for(std::set<int>::const_iterator itr=fieldIds.begin();itr!=fieldIds.end();++itr) {
      std::string field = dofs.getFieldString(*itr);

      // get the layout of each field
      const std::vector<int> & fieldOffsets = dofs.getGIDFieldOffsets(block,*itr);
      for(std::size_t f=0;f<fieldOffsets.size();f++)
        fields[fieldOffsets[f]] = field;

    }

    // print the layout of the full pattern
    os << "#" << std::endl;
    os << "# Element Block \"" << block << "\"" << std::endl;
    os << "#   field pattern = [ " << fields[0];
    for(std::size_t f=1;f<fields.size();f++)
      os << ", " << fields[f];
    os << " ]" << std::endl;
    os << "#" << std::endl;

    const std::vector<int> & elements = dofs.getElementBlock(block);
    for(std::size_t e=0;e<elements.size();e++) {
      std::vector<GO> gids;
      dofs.getElementGIDs(elements[e],gids,block);

      // output gids belonging to this element
      os << "[ " << gids[0];
      for(std::size_t g=1;g<gids.size();g++)
        os << ", " << gids[g];
      os << " ]" << std::endl;
    }
  }
#endif

}

#endif
