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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_TENTATIVEPFACTORY_DEF_HPP
#define MUELU_TENTATIVEPFACTORY_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_TentativePFactory_decl.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_NullspaceFactory.hpp" //FIXME
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // validParamList->set< bool >("QR",                                         true, "Use QR factorization"); Not implemented for QR=false

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",          Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",          Teuchos::null, "Generating factory of the coarse map");
    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Aggregates");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "UnAmalgamationInfo");
    Input(fineLevel, "CoarseMap");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level & fineLevel, Level & coarseLevel) const {

    FactoryMonitor m(*this, "Build", coarseLevel);

    RCP<Matrix>           A          = Get< RCP<Matrix> >          (fineLevel, "A");
    RCP<Aggregates>       aggregates = Get< RCP<Aggregates> >      (fineLevel, "Aggregates");
    RCP<AmalgamationInfo> amalgInfo  = Get< RCP<AmalgamationInfo> >(fineLevel, "UnAmalgamationInfo");
    RCP<MultiVector>      nullspace  = Get< RCP<MultiVector> >     (fineLevel, "Nullspace");
    RCP<const Map>        coarseMap  = Get< RCP<const Map> >       (fineLevel, "CoarseMap");

    // Build
    RCP<MultiVector> coarseNullspace; RCP<Matrix> Ptentative; // output of MakeTentative()

    MakeTentative(*A, *aggregates, *amalgInfo, *nullspace, coarseMap, coarseNullspace, Ptentative);

    // Level Set
    Set(coarseLevel, "Nullspace", coarseNullspace);
    Set(coarseLevel, "P",         Ptentative);

    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics0,0) << Utils::PrintMatrixInfo(*Ptentative, "Ptent", params);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MakeTentative(
                     const Matrix& fineA, const Aggregates& aggregates, const AmalgamationInfo& amalgInfo,
                     const MultiVector & fineNullspace, RCP<const Map> coarseMap,
                     RCP<MultiVector> & coarseNullspace, RCP<Matrix> & Ptentative) const
  {
    RCP<const Teuchos::Comm<int> > comm = fineA.getRowMap()->getComm();
    LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    // number of aggregates
    GO numAggs = aggregates.GetNumAggregates();

    // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
    // aggStart is a pointer into aggToRowMap
    // aggStart[i]..aggStart[i+1] are indices into aggToRowMap
    // aggToRowMap[aggStart[i]]..aggToRowMap[aggStart[i+1]-1] are the DOFs in aggregate i
    ArrayRCP<LO> aggStart;
    ArrayRCP< GO > aggToRowMap;
    AmalgamationFactory::UnamalgamateAggregates(aggregates, amalgInfo, aggStart, aggToRowMap);

    // find size of the largest aggregate
    LO maxAggSize=0;
    for (GO i=0; i<numAggs; ++i) {
      LO sizeOfThisAgg = aggStart[i+1] - aggStart[i];
      if (sizeOfThisAgg > maxAggSize) maxAggSize = sizeOfThisAgg;
    }

    // dimension of fine level nullspace
    const size_t NSDim = fineNullspace.getNumVectors();

    // index base for coarse Dof map (usually 0)
    GO indexBase=fineA.getRowMap()->getIndexBase();

    const RCP<const Map> nonUniqueMap = AmalgamationFactory::ComputeUnamalgamatedImportDofMap(aggregates, amalgInfo);
    const RCP<const Map> uniqueMap    = fineA.getDomainMap();
    RCP<const Import> importer = ImportFactory::Build(uniqueMap, nonUniqueMap);
    RCP<MultiVector> fineNullspaceWithOverlap = MultiVectorFactory::Build(nonUniqueMap,NSDim);
    fineNullspaceWithOverlap->doImport(fineNullspace,*importer,Xpetra::INSERT);

    // Pull out the nullspace vectors so that we can have random access.
    ArrayRCP< ArrayRCP<const SC> > fineNS(NSDim);
    for (size_t i=0; i<NSDim; ++i)
      fineNS[i] = fineNullspaceWithOverlap->getData(i);

    //Allocate storage for the coarse nullspace.
    coarseNullspace = MultiVectorFactory::Build(coarseMap,NSDim);

    ArrayRCP< ArrayRCP<SC> > coarseNS(NSDim);
    for (size_t i=0; i<NSDim; ++i)
      if (coarseMap->getNodeNumElements() > 0) coarseNS[i] = coarseNullspace->getDataNonConst(i);


    //This makes the rowmap of Ptent the same as that of fineA.
    //This requires moving some parts of some local Q's to other processors
    //because aggregates can span processors.
    RCP<const Map > rowMapForPtent = fineA.getRowMap();
    size_t numRowsForPtent = rowMapForPtent->getNodeNumElements();
    const Map& rowMapForPtentRef = *rowMapForPtent;

    // prerequisites: rowMapForPtent, NSDim

    // Set up storage for the rows of the local Qs that belong to other processors.
    // FIXME This is inefficient and could be done within the main loop below with std::vector's.
    RCP<const Map> colMap = fineA.getColMap();
    Array<GO> ghostGIDs;
    for (LO j=0; j<numAggs; ++j) {
      for (LO k=aggStart[j]; k<aggStart[j+1]; ++k) {
        if (rowMapForPtentRef.isNodeGlobalElement(aggToRowMap[k]) == false) {
          ghostGIDs.push_back(aggToRowMap[k]);
        }
      }
    }
    RCP<const Map > ghostQMap = MapFactory::Build(fineA.getRowMap()->lib(),
                                                  Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                  ghostGIDs,
                                                  indexBase, fineA.getRowMap()->getComm()); //JG:Xpetra::global_size_t>?
    //Vector to hold bits of Q that go to other processors.
    RCP<MultiVector> ghostQvalues = MultiVectorFactory::Build(ghostQMap,NSDim);
    //Note that Epetra does not support MultiVectors templated on Scalar != double.
    //So to work around this, we allocate an array of Vectors.  This shouldn't be too
    //expensive, as the number of Vectors is NSDim.
    Array<RCP<Xpetra::Vector<GO,LO,GO,Node> > > ghostQcolumns(NSDim);
    for (size_t i=0; i<NSDim; ++i)
      ghostQcolumns[i] = Xpetra::VectorFactory<GO,LO,GO,Node>::Build(ghostQMap);
    RCP<Xpetra::Vector<GO,LO,GO,Node> > ghostQrowNums = Xpetra::VectorFactory<GO,LO,GO,Node>::Build(ghostQMap);
    ArrayRCP< ArrayRCP<SC> > ghostQvals;
    ArrayRCP< ArrayRCP<GO> > ghostQcols;
    ArrayRCP< GO > ghostQrows;
    if (ghostQvalues->getLocalLength() > 0) {
      ghostQvals.resize(NSDim);
      ghostQcols.resize(NSDim);
      for (size_t i=0; i<NSDim; ++i) {
        ghostQvals[i] = ghostQvalues->getDataNonConst(i);
        ghostQcols[i] = ghostQcolumns[i]->getDataNonConst(0);
      }
      ghostQrows = ghostQrowNums->getDataNonConst(0);
    }

    //importer to handle moving Q
    importer = ImportFactory::Build(ghostQMap, fineA.getRowMap());

    // Dense QR solver
    Teuchos::SerialQRDenseSolver<LO,SC> qrSolver;

    //Allocate temporary storage for the tentative prolongator.
    Array<GO> globalColPtr(maxAggSize*NSDim,0);
    Array<LO> localColPtr(maxAggSize*NSDim,0);
    Array<SC> valPtr(maxAggSize*NSDim,0.);

    //Create column map for Ptent, estimate local #nonzeros in Ptent,  and create Ptent itself.
    const Map& coarseMapRef = *coarseMap;
    size_t nzEstimate = numRowsForPtent*NSDim;

    // For the 3-arrays constructor
    ArrayRCP<size_t>  ptent_rowptr;
    ArrayRCP<LO>      ptent_colind;
    ArrayRCP<Scalar>  ptent_values;

    // Because ArrayRCPs are slow...
    ArrayView<size_t> rowptr_v;
    ArrayView<LO>     colind_v;
    ArrayView<Scalar> values_v;

    // For temporary usage
    Array<size_t>    rowptr_temp;
    Array<LO>        colind_temp;
    Array<Scalar>    values_temp;

    RCP<CrsMatrix> PtentCrs;

    if (aggregates.AggregatesCrossProcessors()) {
      RCP<CrsMatrixWrap> PtentCrsWrap = rcp(new CrsMatrixWrap(rowMapForPtent, NSDim, Xpetra::StaticProfile));
      PtentCrs   = PtentCrsWrap->getCrsMatrix();
      Ptentative = PtentCrsWrap;
    }
    else {
      // Note: NNZ/row estimate is zero since we're using ESFC.
      RCP<CrsMatrixWrap> PtentCrsWrap = rcp(new CrsMatrixWrap(rowMapForPtent, coarseMap, 0, Xpetra::StaticProfile));
      PtentCrs   = PtentCrsWrap->getCrsMatrix();
      Ptentative = PtentCrsWrap;
      // Since the QR will almost certainly have zeros (NSDim>1) or we might have boundary conditions (NSDim==1),
      // we'll use our own temp storage for the colind/values
      // arrays, since we will need to shrink the arrays down later.
      // Perform initial seeding of the rowptr
      rowptr_temp.resize(numRowsForPtent+1,0);
      rowptr_temp[0]=0;
      for(size_t i=1; i < numRowsForPtent+1; ++i)
	rowptr_temp[i] = rowptr_temp[i-1] + NSDim;

      colind_temp.resize(nzEstimate,INVALID);
      values_temp.resize(nzEstimate,Teuchos::ScalarTraits<Scalar>::zero());

      // Alias the ArrayViews for these guys
      rowptr_v = rowptr_temp();
      colind_v = colind_temp();
      values_v = values_temp();
    }

    //*****************************************************************
    //Loop over all aggregates and calculate local QR decompositions.
    //*****************************************************************
    GO qctr=0; //for indexing into Ptent data vectors
    const Map& nonUniqueMapRef = *nonUniqueMap;

    size_t total_nnz_count=0;

    for (GO agg=0; agg<numAggs; ++agg)
    {
      LO myAggSize = aggStart[agg+1]-aggStart[agg];
      // For each aggregate, extract the corresponding piece of the nullspace and put it in the flat array,
      // "localQR" (in column major format) for the QR routine.
      Teuchos::SerialDenseMatrix<LO,SC> localQR(myAggSize, NSDim);
      for (size_t j=0; j<NSDim; ++j) {
        bool bIsZeroNSColumn = true;
        for (LO k=0; k<myAggSize; ++k)
        {
          // aggToRowMap[aggPtr[i]+k] is the kth DOF in the ith aggregate
          // fineNS[j][n] is the nth entry in the jth NS vector
          try{
            SC nsVal = fineNS[j][ nonUniqueMapRef.getLocalElement(aggToRowMap[aggStart[agg]+k]) ]; // extract information from fine level NS
            localQR(k,j) = nsVal;
            if (nsVal != 0.0) bIsZeroNSColumn = false;
          }
          catch(...) {
            std::cout << "length of fine level nsp: " << fineNullspace.getGlobalLength() << std::endl;
            std::cout << "length of fine level nsp w overlap: " << fineNullspaceWithOverlap->getGlobalLength() << std::endl;
            std::cout << "(local?) aggId=" << agg << std::endl;
            std::cout << "aggSize=" << myAggSize << std::endl;
            std::cout << "agg DOF=" << k << std::endl;
            std::cout << "NS vector j=" << j << std::endl;
            std::cout << "j*myAggSize + k = " << j*myAggSize + k << std::endl;
            std::cout << "aggToRowMap["<<agg<<"][" << k << "] = " << aggToRowMap[aggStart[agg]+k] << std::endl;
            std::cout << "id aggToRowMap[agg][k]=" << aggToRowMap[aggStart[agg]+k] << " is global element in nonUniqueMap = " <<
nonUniqueMapRef.isNodeGlobalElement(aggToRowMap[aggStart[agg]+k]) << std::endl;
            std::cout << "colMap local id aggToRowMap[agg][k]=" << nonUniqueMapRef.getLocalElement(aggToRowMap[aggStart[agg]+k]) << std::endl;
            std::cout << "fineNS...=" << fineNS[j][ nonUniqueMapRef.getLocalElement(aggToRowMap[aggStart[agg]+k]) ] << std::endl;
            std::cerr << "caught an error!" << std::endl;
          }
        } //for (LO k=0 ...
        TEUCHOS_TEST_FOR_EXCEPTION(bIsZeroNSColumn == true, Exceptions::RuntimeError, "MueLu::TentativePFactory::MakeTentative: fine level NS part has a zero column. Error.");
      } //for (LO j=0 ...

#if 0
      std::cout << "Input" << std::endl;
      std::cout << "myAggSize " << myAggSize << std::endl;
      // loop over rows
      for (size_t i=0; i<myAggSize; i++) {
        // loop over cols
        for (size_t j=0; j<NSDim; j++) {
          std::cout << localQR(i,j); std::cout << "\t";
        }
        std::cout << std::endl;
      }
#endif

      Xpetra::global_size_t offset=agg*NSDim;

      if(myAggSize >= Teuchos::as<LocalOrdinal>(NSDim)) {
        // calculate QR decomposition (standard)
        // R is stored in localQR (size: myAggSize x NSDim)

        // Householder multiplier
        SC tau = localQR(0,0);

        if (NSDim == 1) {
          // Only one nullspace vector, so normalize by hand
          Magnitude dtemp=0;
          for (size_t k=0; k<static_cast<size_t>(myAggSize); ++k) {
	    Magnitude tmag = Teuchos::ScalarTraits<SC>::magnitude(localQR(k,0));
            dtemp += tmag*tmag;
          }
          dtemp = Teuchos::ScalarTraits<Magnitude>::squareroot(dtemp);
          tau = localQR(0,0);
          localQR(0,0) = dtemp;
        } else {
          qrSolver.setMatrix( Teuchos::rcp(&localQR, false) );
          qrSolver.factor();
        }

         // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
         // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
         // This stores the (offset+k)th entry only if it is local according to the coarseMap.
         for (size_t j=0; j<NSDim; ++j) {
           for (size_t k=0; k<=j; ++k) {
             try {
               if (coarseMapRef.isNodeLocalElement(offset+k)) {
                 coarseNS[j][offset+k] = localQR(k, j); //TODO is offset+k the correct local ID?!
               }
             }
             catch(...) {
               std::cout << "caught error in coarseNS insert, j="<<j<<", offset+k = "<<offset+k<<std::endl;
             }
           }
         }

         // Calculate Q, the tentative prolongator.
         // The Lapack GEQRF call only works for myAggsize >= NSDim

         if (NSDim == 1) {
           // Only one nullspace vector, so calculate Q by hand
           Magnitude dtemp = Teuchos::ScalarTraits<SC>::magnitude(localQR(0,0));
           localQR(0,0) = tau;
           dtemp = 1 / dtemp;
           for (LocalOrdinal i=0; i<myAggSize; ++i) {
             localQR(i,0) *= dtemp ;
           }
         } else {
           qrSolver.formQ();
           Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SC> > qFactor = qrSolver.getQ();
           for (size_t j=0; j<NSDim; j++) {
             for (size_t i=0; i<static_cast<size_t>(myAggSize); i++) {
               localQR(i,j) = (*qFactor)(i,j);
             }
           }
         }

         // end default case (myAggSize >= NSDim)
      } else {  // special handling for myAggSize < NSDim (i.e. 1pt nodes)
        // construct R by hand, i.e. keep first myAggSize rows untouched
        //GetOStream(Warnings0,0) << "TentativePFactory (WARNING): aggregate with " << myAggSize << " DOFs and nullspace dim " << NSDim << ". special handling of QR decomposition." << std::endl;

        localQR.reshape(NSDim,NSDim);
        for (size_t i=myAggSize; i<NSDim; i++) {
          localQR(i,i) = Teuchos::ScalarTraits<SC>::one();
        }

        // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
        // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
        // This stores the (offset+k)th entry only if it is local according to the coarseMap.

        for (size_t j=0; j<NSDim; ++j) {
          for (size_t k=0; k<=j; ++k) {
            try {
              if (coarseMapRef.isNodeLocalElement(offset+k)) {
                coarseNS[j][offset+k] = localQR(k,j); // agg has only one node
              }
            }
            catch(...) {
              std::cout << "caught error in coarseNS insert, j="<<j<<", offset+k = "<<offset+k<<std::endl;
            }
          }
        }

        // Calculate Q, the tentative prolongator.
        // The Lapack GEQRF call only works for myAggsize >= NSDim
        // special handling for very small aggregates (with myAggsize < NSDim)

        // calculate identity matrix with zero columns for the last NSDim-myAggsize columns.
        for (size_t i=0; i<Teuchos::as<size_t>(myAggSize); i++) {
          // loop over cols
          for (size_t j=0; j<NSDim; j++) {
            if (j==i) localQR(i,j) = Teuchos::ScalarTraits<SC>::one();
            else localQR(i,j) = Teuchos::ScalarTraits<SC>::zero();
          }
        }
      } // end else (special handling for 1pt aggregates)

      //Process each row in the local Q factor.  If the row is local to the current processor
      //according to the rowmap, insert it into Ptentative.  Otherwise, save it in ghostQ
      //to be communicated later to the owning processor.
      //FIXME -- what happens if maps are blocked?
      for (GO j=0; j<myAggSize; ++j) {
        //This loop checks whether row associated with current DOF is local, according to rowMapForPtent.
        //If it is, the row is inserted.  If not, the row number, columns, and values are saved in
        //MultiVectors that will be sent to other processors.
        GO globalRow = aggToRowMap[aggStart[agg]+j];
	LO localRow  = rowMapForPtent->getLocalElement(globalRow); // CMS: There has to be an efficient way to do this...

        //TODO is the use of Xpetra::global_size_t below correct?
        if( rowMapForPtentRef.isNodeGlobalElement(globalRow) == false )
        {
          ghostQrows[qctr] = globalRow;
          for (size_t k=0; k<NSDim; ++k) {
            ghostQcols[k][qctr] = coarseMapRef.getGlobalElement(agg*NSDim+k);
            ghostQvals[k][qctr] = localQR(j,k);
          }
          ++qctr;
        } else {
          size_t nnz=0;
          for (size_t k=0; k<NSDim; ++k) {
            try{
              if (localQR(j,k) != Teuchos::ScalarTraits<SC>::zero()) {
		localColPtr[nnz]  = agg * NSDim + k;
		globalColPtr[nnz] = coarseMapRef.getGlobalElement(localColPtr[nnz]);
                valPtr[nnz] = localQR(j,k);
		++total_nnz_count;
                ++nnz;
              }
            }
            catch(...) {
              std::cout << "caught error in colPtr/valPtr insert, current index="<<nnz<<std::endl;
            }
          } //for (size_t k=0; k<NSDim; ++k)

          try{
	    if (aggregates.AggregatesCrossProcessors()) {
	      Ptentative->insertGlobalValues(globalRow,globalColPtr.view(0,nnz),valPtr.view(0,nnz));
	    }
	    else {
	      // Copy all of the *active* cols/vals into the colind_v/values_v arrays.
	      size_t start = rowptr_v[localRow];
	      for(size_t i=0; i<nnz; ++i) {
		colind_v[start+i] = localColPtr[i];
		values_v[start+i] = valPtr[i];
	      }
	    }
          }
          catch(...) {
            std::cout << "pid " << fineA.getRowMap()->getComm()->getRank()
                      << "caught error during Ptent row insertion, global row "
                      << globalRow << std::endl;
          }
        } //if (rowMapForPtent->getGlobalElement(localRow) == ...
      } //for (GO j=0; j<myAggSize; ++j)

    } // for (LO agg=0; agg<numAggs; ++agg)


    // ***********************************************************
    // ************* end of aggregate-wise QR ********************
    // ***********************************************************

    if (!aggregates.AggregatesCrossProcessors()) {
      GetOStream(Runtime1,0) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;
      // The QR will have plenty of zeros if we're NSDim > 1.  If NSDim==1, we still might have a handful of BCs.
      // We need to shrink down the CRS arrays and copy them over the the "real" CrsArrays.

      // Now allocate the final arrays
      PtentCrs->allocateAllValues(total_nnz_count,ptent_rowptr,ptent_colind,ptent_values);

      // Because ArrayRCPs are slow...
      ArrayView<size_t> rowptr_new = ptent_rowptr();
      ArrayView<LO>     colind_new = ptent_colind();
      ArrayView<Scalar> values_new = ptent_values();
      size_t count=0;
      // Collapse and copy
      for(size_t i=0; i<numRowsForPtent; i++) {
	rowptr_new[i]=count;
	for(size_t j=rowptr_v[i]; j<rowptr_v[i+1] && colind_v[j]!=INVALID; j++){
	  colind_new[count] = colind_v[j];
	  values_new[count] = values_v[j];
	    count++;
	}
      }
      rowptr_new[numRowsForPtent]=count;

      if(count!=total_nnz_count) throw std::runtime_error("MueLu error in data copy!");

      // Regardless, it is time to call setAllValues & ESFC
      PtentCrs->setAllValues(ptent_rowptr,ptent_colind,ptent_values);
      PtentCrs->expertStaticFillComplete(coarseMap,fineA.getDomainMap());
    } else {
      GetOStream(Runtime1,0) << "TentativePFactory : aggregates may cross process boundaries" << std::endl;
      // Import ghost parts of Q factors and insert into Ptentative.
      // First import just the global row numbers.
      RCP<Xpetra::Vector<GO,LO,GO,Node> > targetQrowNums = Xpetra::VectorFactory<GO,LO,GO,Node>::Build(rowMapForPtent);
      targetQrowNums->putScalar(-1);
      targetQrowNums->doImport(*ghostQrowNums,*importer,Xpetra::INSERT);
      ArrayRCP< GO > targetQrows = targetQrowNums->getDataNonConst(0);

      // Now create map based on just the row numbers imported.
      Array<GO> gidsToImport;
      gidsToImport.reserve(targetQrows.size());
      for (typename ArrayRCP<GO>::iterator r=targetQrows.begin(); r!=targetQrows.end(); ++r) {
        if (*r > -1) {
          gidsToImport.push_back(*r);
        }
      }
      RCP<const Map > reducedMap = MapFactory::Build( fineA.getRowMap()->lib(),
                                                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                      gidsToImport, indexBase, fineA.getRowMap()->getComm()    );

      // Import using the row numbers that this processor will receive.
      importer = ImportFactory::Build(ghostQMap, reducedMap);

      Array<RCP<Xpetra::Vector<GO,LO,GO,Node> > > targetQcolumns(NSDim);
      for (size_t i=0; i<NSDim; ++i) {
        targetQcolumns[i] = Xpetra::VectorFactory<GO,LO,GO,Node>::Build(reducedMap);
        targetQcolumns[i]->doImport(*(ghostQcolumns[i]),*importer,Xpetra::INSERT);
      }
      RCP<MultiVector> targetQvalues = MultiVectorFactory::Build(reducedMap,NSDim);
      targetQvalues->doImport(*ghostQvalues,*importer,Xpetra::INSERT);

      ArrayRCP< ArrayRCP<SC> > targetQvals;
      ArrayRCP<ArrayRCP<GO> > targetQcols;
      if (targetQvalues->getLocalLength() > 0) {
        targetQvals.resize(NSDim);
        targetQcols.resize(NSDim);
        for (size_t i=0; i<NSDim; ++i) {
          targetQvals[i] = targetQvalues->getDataNonConst(i);
          targetQcols[i] = targetQcolumns[i]->getDataNonConst(0);
        }
      }

      valPtr = Array<SC>(NSDim,0.);
      globalColPtr = Array<GO>(NSDim,0);
      for (typename Array<GO>::iterator r=gidsToImport.begin(); r!=gidsToImport.end(); ++r) {
        if (targetQvalues->getLocalLength() > 0) {
          for (size_t j=0; j<NSDim; ++j) {
            valPtr[j] = targetQvals[j][reducedMap->getLocalElement(*r)];
            globalColPtr[j] = targetQcols[j][reducedMap->getLocalElement(*r)];
          }
          Ptentative->insertGlobalValues(*r, globalColPtr.view(0,NSDim), valPtr.view(0,NSDim));
        } //if (targetQvalues->getLocalLength() > 0)
      }

      Ptentative->fillComplete(coarseMap,fineA.getDomainMap()); //(domain,range) of Ptentative
    } //if (!aggregatesAreLocal)



//    RCP<const Map> realColMap = Ptentative->getColMap();
//    sleep(1);
//    comm->barrier();
//    fos->setOutputToRootOnly(0);
//    *fos << "====================\nthe real col map\n======================" << std::endl;
//    fos->setOutputToRootOnly(-1);
//    realColMap->describe(*fos,Teuchos::VERB_EXTREME);

    // if available, use striding information of fine level matrix A for range map and coarseMap as domain map
    // otherwise use plain range map of Ptent = plain range map of A for range map and coarseMap as domain map.
    // Note: the latter is not really safe, since there is no striding information for the range map. This is not
    // really a problem, since striding information is always available on the intermedium levels and the coarsest levels.
    if(fineA.IsView("stridedMaps") == true) {
      Ptentative->CreateView("stridedMaps", fineA.getRowMap("stridedMaps"), coarseMap);
    } else Ptentative->CreateView("stridedMaps", Ptentative->getRangeMap(), coarseMap);

  } //MakeTentative()

} //namespace MueLu

//TODO: noQR_

// TODO ReUse: If only P or Nullspace is missing, TentativePFactory can be smart and skip part of the computation.

#define MUELU_TENTATIVEPFACTORY_SHORT
#endif // MUELU_TENTATIVEPFACTORY_DEF_HPP
