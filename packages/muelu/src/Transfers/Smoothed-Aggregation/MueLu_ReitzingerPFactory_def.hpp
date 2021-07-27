// @HEADER
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
#ifndef MUELU_REITZINGERPFACTORY_DEF_HPP
#define MUELU_REITZINGERPFACTORY_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_ReitzingerPFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))


#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Pnodal",             Teuchos::null, "Generating factory of the matrix P");
    validParamList->set< RCP<const FactoryBase> >("D0",                 Teuchos::null, "Generating factory of the matrix D0");
    validParamList->set< RCP<const FactoryBase> >("NodeMatrix",         Teuchos::null, "Generating factory of the matrix NodeMatrix");
    //    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    //    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    //    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory of the coordinates");

    // Make sure we don't recursively validate options for the matrixmatrix kernels
    ParameterList norecurse;
    norecurse.disableRecursiveValidation();
    validParamList->set<ParameterList> ("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {

    const ParameterList& pL = GetParameterList();

    // NOTE: This guy can only either be 'Nullspace' or 'Scaled Nullspace' or else the validator above will cause issues
    std::string nspName = "Nullspace";
    if(pL.isParameter("Nullspace name")) nspName = pL.get<std::string>("Nullspace name");

    Input(fineLevel, "A");
    Input(coarseLevel, "Pnodal");
    Input(fineLevel, "D0");
    Input(fineLevel, "NodeMatrix");
    //    Input(fineLevel, "Aggregates");
    //    Input(fineLevel, "UnAmalgamationInfo");
    /*
    if( fineLevel.GetLevelID() == 0 &&
        fineLevel.IsAvailable("Coordinates", NoFactory::get()) &&     // we have coordinates (provided by user app)
        pL.get<bool>("tentative: build coarse coordinates") ) {       // and we want coordinates on other levels
      bTransferCoordinates_ = true;                                   // then set the transfer coordinates flag to true
      Input(fineLevel, "Coordinates");
    } else if (bTransferCoordinates_) {
      Input(fineLevel, "Coordinates");
    }
    */
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);
    //    typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinate_type;
    //    typedef Xpetra::MultiVector<coordinate_type,LO,GO,NO> RealValuedMultiVector;
    //    typedef Xpetra::MultiVectorFactory<coordinate_type,LO,GO,NO> RealValuedMultiVectorFactory;
    using MT  = typename Teuchos::ScalarTraits<SC>::magnitudeType;
    using XMM = Xpetra::MatrixMatrix<SC,LO,GO,NO>;
    Teuchos::FancyOStream& out0=GetBlackHole();
    const ParameterList& pL = GetParameterList();
    RCP<Matrix>                EdgeMatrix    = Get< RCP<Matrix> >               (fineLevel, "A");
    RCP<Matrix>                Pn            = Get< RCP<Matrix> >               (coarseLevel, "Pnodal");
    RCP<Matrix>                D0            = Get< RCP<Matrix> >               (fineLevel, "D0");
    RCP<Matrix>                NodeMatrix    = Get< RCP<Matrix> >               (fineLevel, "NodeMatrix");
    //    RCP<Aggregates>            aggregates    = Get< RCP<Aggregates> >           (fineLevel, "Aggregates");

    RCP<CrsMatrixWrap> Pn_crs = rcp_dynamic_cast<CrsMatrixWrap>(Pn);

    // Matrix matrix params
    RCP<ParameterList> mm_params = rcp(new ParameterList);;
    if(pL.isSublist("matrixmatrix: kernel params"))
      mm_params->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

    /* Get all of the aggregation information */
    //    LO numAggs      = aggregates->GetNumAggregates();
    //    ArrayRCP<LO> aggStart;
    //    ArrayRCP<LO> aggToRowMapLO;
    //    ArrayRCP<GO> aggToRowMapGO;
    //    if (goodMap) {
    //      amalgInfo->UnamalgamateAggregatesLO(*aggregates, aggStart, aggToRowMapLO);
    //      GetOStream(Runtime1) << "Column map is consistent with the row map, good." << std::endl;
    //    }
    //    else {
    //      TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError, "MueLu::ReitzingerPFactory: Needs a good map.");
    //    }

    // FIXME: We probably need to ghost the aggregates.  Do we?

    // FIXME: We need to make sure Pn isn't normalized

    // FIXME: We need to look through and see which of these really need importers and which ones don't

    /* Generate the Pn * D0 matrix and its transpose */
    RCP<Matrix> D0_Pn, PnT_D0T;
    {
      RCP<Matrix> dummy;
      SubFactoryMonitor m2(*this, "Generate D0*Pn", coarseLevel);
      D0_Pn = XMM::Multiply(*D0,false,*Pn,false,dummy,out0,true,true,"D0*Pn",mm_params);
    }
    {
      // FIXME: Do we need to optimize the transpose?
      SubFactoryMonitor m2(*this, "Transpose D0*Pn", coarseLevel);
      PnT_D0T = Utilities::Transpose(*D0_Pn, true);
    }

    // FIXME: This is using deprecated interfaces
    ArrayView<const LO>     colind_E, colind_N;
    ArrayView<const SC>     values_E, values_N;

    size_t Ne=EdgeMatrix->getNodeNumRows();
    size_t Nn=NodeMatrix->getNodeNumRows();

    // Upper bound on local number of coarse edges
    size_t max_edges = (NodeMatrix->getNodeNumEntries() + Nn +1) / 2;      
    ArrayRCP<size_t>  D0_rowptr(Ne+1);
    ArrayRCP<LO>      D0_colind(max_edges);
    ArrayRCP<SC>      D0_values(max_edges);
    D0_rowptr[0] = 0;


    LO current = 0;
    for(LO i=0; i<(LO)Nn; i++) {
     
      // FIXME: We don't really want an std::map here.  This is just a first cut implementation
      using value_type = std::pair<LO,SC>;
      std::map<LO, value_type> ce_map;
      
      // FIXME: This is using deprecated interfaces
      PnT_D0T->getLocalRowView(i,colind_E,values_E);
      
      for(LO j=0; j<(LO)colind_E.size(); j++) {
        LO j_row = colind_E[j];
        D0_Pn->getLocalRowView(j_row,colind_N,values_N);

        // Skip incomplete rows
        if(colind_N.size() != 2) continue;

        // FIXME:  Do we even need the values?????

        for(LO k=0; k<(LO)colind_N.size(); k++) {
          // FIXME: Should this be done in GID space?  Probably
          if(colind_N[k] > i) {
            ce_map.emplace(std::make_pair(colind_N[k],std::make_pair(colind_E[j],values_N[k])));
          }
        }//end for k < colind_N.size()
      }// end for j < colind_E.size()
      
      
      // std::map is sorted, so we'll just iterate through this
      for(auto iter=ce_map.begin(); iter != ce_map.end(); iter++) {
        LO col = iter->first;
        //        printf("[%d] Adding edge %d => %d\n",i,i,col);
        D0_colind[current]  = i;
        D0_values[current] = -1;
        current++;
        D0_colind[current]  = col;
        D0_values[current] = 1;
        current++;
        D0_rowptr[current / 2] = current;
      }

    }// end for i < Nn

    LO num_coarse_edges = current / 2;
    D0_rowptr.resize(num_coarse_edges+1);
    D0_colind.resize(current);
    D0_values.resize(current);



    // Count the total number of edges
    // FIXME: Since we're not using global ids above, this won't work in parallel.  Need to fix this.
    // I suspect we'll have to conditionally number edges based on owners.  Somehow.
    GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    RCP<const Map> ownedCoarseEdgeMap = Xpetra::MapFactory<LO,GO,NO>::Build(EdgeMatrix->getRowMap()->lib(), GO_INVALID, num_coarse_edges,EdgeMatrix->getRowMap()->getIndexBase(),EdgeMatrix->getRowMap()->getComm());
    RCP<const Map> ownedPlusSharedCoarseEdgeMap = ownedCoarseEdgeMap;

    // FIXME: Not sure this is right in parallel.  Does owning an edge guarantee that it is in the node's colmap
    RCP<const Map> ownedCoarseNodeMap = Pn->getDomainMap();
    RCP<const Map> ownedPlusSharedCoarseNodeMap  = Pn_crs->getColMap();

    // Create the coarse D0
    RCP<CrsMatrix> D0_coarse;
    {
      SubFactoryMonitor m2(*this, "Build D0", coarseLevel);
      // FIXME: We can be smarter with memory here
      // FIXME: Can we cheeseball the importer somehow?
      D0_coarse = CrsMatrixFactory::Build(ownedCoarseEdgeMap,ownedPlusSharedCoarseNodeMap,0);
      TEUCHOS_TEST_FOR_EXCEPTION(D0_coarse.is_null(), Exceptions::RuntimeError, "MueLu::ReitzingerPFactory: CrsMatrixFatory failed.");

      // FIXME: Deprecated code
      ArrayRCP<size_t>  ia;
      ArrayRCP<LO>      ja;
      ArrayRCP<SC>     val;
      D0_coarse->allocateAllValues(current, ia, ja, val);
      //      printf("D0_rowptr.size() = %d D0_colind.size() = %d D0_values.size() = %d\n",(int)D0_rowptr.size(),(int)D0_colind.size(),(int)D0_values.size());

      //      printf("Before: ia.size() = %d ja.size() = %d val.size() = %d\n",(int)ia.size(),(int)ja.size(),(int)val.size());
      std::copy(D0_rowptr.begin(),D0_rowptr.end(),ia.begin());
      std::copy(D0_colind.begin(),D0_colind.end(),ja.begin());
      std::copy(D0_values.begin(),D0_values.end(),val.begin());
      //      printf("After: ia.size() = %d ja.size() = %d val.size() = %d\n",(int)ia.size(),(int)ja.size(),(int)val.size());
      D0_coarse->setAllValues(ia, ja, val);

#if 0
      printf("D0: ia  :");
      for(int i=0; i<(int)ia.size(); i++)
        printf("%d ",(int)ia[i]);
      printf("\nD0: ja  :");
      for(int i=0; i<(int)ja.size(); i++)
        printf("%d ",ja[i]);
      printf("\n");
#endif

      D0_coarse->expertStaticFillComplete(ownedCoarseNodeMap,ownedCoarseEdgeMap);
    }
    RCP<Matrix> D0_coarse_m = rcp(new CrsMatrixWrap(D0_coarse));

    // Create the Pe matrix, but with the extra entries.  From ML's notes:
    /* The general idea is that the matrix                              */
    /*                        T_h P_n T_H^*                             */
    /* is almost Pe. If we make sure that P_n contains 1's and -1's, the*/
    /* matrix triple product will yield a matrix with +/- 1 and +/- 2's.*/
    /* If we remove all the 1's and divide the 2's by 2. we arrive at Pe*/
    RCP<Matrix> Pe;
    {
      SubFactoryMonitor m2(*this, "Generate Pe (pre-fix)", coarseLevel);
      RCP<Matrix> dummy;
      RCP<Matrix> Pn_D0cT = XMM::Multiply(*Pn,false,*D0_coarse_m,true,dummy,out0,true,true,"Pn*D0c'",mm_params);
      Pe = XMM::Multiply(*D0,false,*Pn_D0cT,false,dummy,out0,true,true,"D0*(Pn*D0c')",mm_params);
    }  

    /* Weed out the +/- entries */
    { 
      SubFactoryMonitor m2(*this, "Generate Pe (post-fix)", coarseLevel);
      Pe->resumeFill();
      SC one = Teuchos::ScalarTraits<SC>::one();
      MT two = 2*Teuchos::ScalarTraits<MT>::one();
      SC zero = Teuchos::ScalarTraits<SC>::zero();
      SC neg_one = - one;
      // FIXME: Deprecated code
      for(LO i=0; i<(LO) Ne; i++) {
        ArrayView<const LO>  columns;
        ArrayView<const SC>  values;
        Pe->getLocalRowView(i,columns,values);
        // FIXME: This won't work on fancy nodes
        ArrayView<SC> values_nc = Teuchos::av_const_cast<SC>(values);
        for (LO j=0; j<(LO)values.size(); j++)
          if ((values_nc[j] == one || values_nc[j] == neg_one))
            values_nc[j] = zero;       
          else
            values_nc[j] /= two;
      }//end for i < Ne
      Pe->fillComplete(Pe->getDomainMap(),Pe->getRangeMap());
    }


    /* Check commuting property */
    CheckCommutingProperty(*Pe,*D0_coarse_m,*D0,*Pn);
    
    /* Set output on the level */
    Set(coarseLevel,"P",Pe);
    Set(coarseLevel,"D0",D0_coarse_m);
    

  }// end Build

 template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
 void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
 CheckCommutingProperty(const Matrix & Pe, const Matrix & D0_c, const Matrix& D0_f, const Matrix & Pn) const {   
   if(IsPrint(Statistics0)) {
     using XMM = Xpetra::MatrixMatrix<SC,LO,GO,NO>;
     using MT   = typename Teuchos::ScalarTraits<SC>::magnitudeType;
     SC one = Teuchos::ScalarTraits<SC>::one();
     SC zero = Teuchos::ScalarTraits<SC>::zero();

     RCP<Matrix> dummy;
     Teuchos::FancyOStream &out0=GetBlackHole();     
     RCP<Matrix> left  = XMM::Multiply(Pe,false,D0_c,false,dummy,out0);
     RCP<Matrix> right = XMM::Multiply(D0_f,false,Pn,false,dummy,out0);
     
     // We need a non-FC matrix for the add, sadly
     RCP<CrsMatrix> sum_c = CrsMatrixFactory::Build(left->getRowMap(),left->getNodeMaxNumRowEntries()+right->getNodeMaxNumRowEntries());    
     RCP<Matrix> summation = rcp(new CrsMatrixWrap(sum_c));
     XMM::TwoMatrixAdd(*left,  false, one, *summation, zero);
     XMM::TwoMatrixAdd(*right, false, -one, *summation, one);
     
     MT norm = summation->getFrobeniusNorm();     
     GetOStream(Statistics0) << "CheckCommutingProperty: ||Pe D0_c - D0_f Pn || = "<<norm<<std::endl;
   }

 }//end CheckCommutingProperty



} //namespace MueLu



#define MUELU_REITZINGERPFACTORY_SHORT
#endif // MUELU_REITZINGERPFACTORY_DEF_HPP
