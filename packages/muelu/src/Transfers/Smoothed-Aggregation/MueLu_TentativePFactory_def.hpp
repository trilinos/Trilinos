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

#include "MueLu_TentativePFactory_decl.hpp"
#include "MueLu_QR_Interface.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_NullspaceFactory.hpp" //FIXME
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TentativePFactory(RCP<const FactoryBase> aggregatesFact, RCP<const FactoryBase> amalgFact, RCP<const FactoryBase> nullspaceFact, RCP<const FactoryBase> AFact, RCP<const FactoryBase> coarseMapFact)
    : aggregatesFact_(aggregatesFact), amalgFact_(amalgFact), nullspaceFact_(nullspaceFact), AFact_(AFact),coarseMapFact_(coarseMapFact),
      QR_(false)//,
      /*domainGidOffset_(0)*/ {

    //stridedBlockId_ = -1; // default: blocked map with constant blocksize "NSDim"
  }
 
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TentativePFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    fineLevel.DeclareInput("A", AFact_.get(), this);
    fineLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
    fineLevel.DeclareInput("Nullspace",  nullspaceFact_.get(), this);
    fineLevel.DeclareInput("UnAmalgamationInfo", amalgFact_.get(), this);
    fineLevel.DeclareInput("CoarseMap", coarseMapFact_.get(), this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TentativeWithQR(bool value) { QR_ = value; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TentativeWithQR() { return QR_; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level & fineLevel, Level & coarseLevel) const {

    FactoryMonitor m(*this, "Tentative prolongator", coarseLevel);

    RCP<Matrix> A = fineLevel.Get< RCP<Matrix> >("A", AFact_.get());

    RCP<Aggregates>  aggregates = fineLevel.Get< RCP<Aggregates> >("Aggregates", aggregatesFact_.get());
    RCP<AmalgamationInfo> amalgInfo = fineLevel.Get< RCP<AmalgamationInfo> >("UnAmalgamationInfo", amalgFact_.get());
    RCP<MultiVector> nullspace  = fineLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());
    RCP<const Map> coarseMap = fineLevel.Get<RCP<const Map> >("CoarseMap",coarseMapFact_.get());

    // Build
    RCP<MultiVector> coarseNullspace; RCP<Matrix> Ptentative; // output of MakeTentative()

    MakeTentative(*A, *aggregates, *amalgInfo, *nullspace, coarseMap, coarseNullspace, Ptentative);

    // Level Set
    coarseLevel.Set("Nullspace", coarseNullspace, this);
    coarseLevel.Set("P", Ptentative, this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MakeTentative(
                     const Matrix& fineA, const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const MultiVector & fineNullspace, RCP<const Map> coarseMap,
                     RCP<MultiVector> & coarseNullspace, RCP<Matrix> & Ptentative) const
  {
    RCP<const Teuchos::Comm<int> > comm = fineA.getRowMap()->getComm();

    // number of aggregates
    GO numAggs = aggregates.GetNumAggregates();

    // Compute array of aggregate sizes (in dofs).
    ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs,0);
    AmalgamationFactory::ComputeUnamalgamatedAggregateSizes(aggregates, amalgInfo, aggSizes);

    // find size of the largest aggregate.
    LO maxAggSize=0;
    for (typename Teuchos::ArrayRCP<LO>::iterator i=aggSizes.begin(); i!=aggSizes.end(); ++i) {
      if (*i > maxAggSize) maxAggSize = *i;
    }

    // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
    // aggToRowMap[i][j] is the jth DOF in aggregate i
    ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
    AmalgamationFactory::UnamalgamateAggregates(aggregates, amalgInfo, aggSizes, aggToRowMap);

    // dimension of fine level nullspace
    const size_t NSDim = fineNullspace.getNumVectors();

    // index base for coarse Dof map (usually 0)
    GO indexBase=fineA.getRowMap()->getIndexBase();

    const RCP<const Map> nonUniqueMap = AmalgamationFactory::ComputeUnamalgamatedImportDofMap(aggregates, amalgInfo, aggSizes);
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

    Ptentative = rcp(new CrsMatrixWrap(rowMapForPtent, NSDim, Xpetra::StaticProfile));

    // prerequisites: rowMapForPtent, NSDim

    // Set up storage for the rows of the local Qs that belong to other processors.
    // FIXME This is inefficient and could be done within the main loop below with std::vector's.
    RCP<const Map> colMap = fineA.getColMap();
    Array<GO> ghostGIDs;
    for (LO j=0; j<numAggs; ++j) {
      for (LO k=0; k<aggSizes[j]; ++k) {
        if (rowMapForPtent->isNodeGlobalElement(aggToRowMap[j][k]) == false) {
          ghostGIDs.push_back(aggToRowMap[j][k]);
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

    QR_Interface<SC,LO> qrWidget(NSDim);

    // reserve enough memory for QR decomposition
    size_t localQRsize = maxAggSize*NSDim;
    if(Teuchos::as<size_t>(maxAggSize) < NSDim) localQRsize = NSDim * NSDim; // make sure, that we always can store the coarse nsp matrix (R in QR decomp.) even if we have only too small aggregates for GEQRF routine.
    ArrayRCP<SC> localQR(localQRsize); // The submatrix of the nullspace to be orthogonalized.

    //Allocate temporary storage for the tentative prolongator.
    GO nFineDofs = nonUniqueMap->getNodeNumElements();
    ArrayRCP<GO> rowPtr(nFineDofs+1);
    for (GO i=0; i<=nFineDofs; ++i)
      rowPtr[i] = i*NSDim;
    ArrayRCP<GO> colPtr(maxAggSize*NSDim,0);
    ArrayRCP<SC> valPtr(maxAggSize*NSDim,0.);

    //*****************************************************************
    //Loop over all aggregates and calculate local QR decompositions.
    //*****************************************************************
    GO qctr=0; //for indexing into Ptent data vectors
    for (LO agg=0; agg<numAggs; ++agg)
    {
      LO myAggSize = aggSizes[agg];
      // For each aggregate, extract the corresponding piece of the nullspace and put it in the flat array,
      // "localQR" (in column major format) for the QR routine.
      for (size_t j=0; j<NSDim; ++j) {
        bool bIsZeroNSColumn = true;
        for (LO k=0; k<myAggSize; ++k) {
          // aggToRowMap[i][k] is the kth DOF in the ith aggregate
          // fineNS[j][n] is the nth entry in the jth NS vector
          try{
            //SC nsVal = fineNS[j][ colMap->getLocalElement(aggToRowMap[agg][k]) ]; // extract information from fine level NS
            SC nsVal = fineNS[j][ nonUniqueMap->getLocalElement(aggToRowMap[agg][k]) ]; // extract information from fine level NS
            localQR[j* myAggSize + k] = nsVal;
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
            std::cout << "aggToRowMap["<<agg<<"][" << k << "] = " << aggToRowMap[agg][k] << std::endl;
            std::cout << "id aggToRowMap[agg][k]=" << aggToRowMap[agg][k] << " is global element in nonUniqueMap = " << nonUniqueMap->isNodeGlobalElement(aggToRowMap[agg][k]) << std::endl;
            std::cout << "colMap local id aggToRowMap[agg][k]=" << nonUniqueMap->getLocalElement(aggToRowMap[agg][k]) << std::endl;
            std::cout << "fineNS...=" << fineNS[j][ nonUniqueMap->getLocalElement(aggToRowMap[agg][k]) ] << std::endl;
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
          std::cout << localQR[ myAggSize*j + i]; std::cout << "\t";
        }
        std::cout << std::endl;
      }
#endif

      Xpetra::global_size_t offset=agg*NSDim;

      if(myAggSize >= Teuchos::as<LocalOrdinal>(NSDim)) {
        // calculate QR decomposition (standard)
         // R is stored in localQR (size: myAggSize x NSDim)
         qrWidget.Compute(myAggSize, localQR);

         // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
         // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
         // This stores the (offset+k)th entry only if it is local according to the coarseMap.
         for (size_t j=0; j<NSDim; ++j) {
           for (size_t k=0; k<=j; ++k) {
             try {
               if (coarseMap->isNodeLocalElement(offset+k)) {
                   coarseNS[j][offset+k] = localQR[ myAggSize*j + k ]; //TODO is offset+k the correct local ID?!
               }
             }
             catch(...) {
               std::cout << "caught error in coarseNS insert, j="<<j<<", offset+k = "<<offset+k<<std::endl;
             }
           }
         }

         // Calculate Q, the tentative prolongator.
         // The Lapack GEQRF call only works for myAggsize >= NSDim
         qrWidget.ExtractQ(myAggSize, localQR);

         // end default case (myAggSize >= NSDim)
      } else {  // sepcial handling for myAggSize < NSDim (i.e. 1pt nodes)
        // construct R by hand, i.e. keep first myAggSize rows untouched
        //GetOStream(Warnings0,0) << "TentativePFactory (WARNING): aggregate with " << myAggSize << " DOFs and nullspace dim " << NSDim << ". special handling of QR decomposition." << std::endl;

        // copy initial localQR values in temporary variable
        ArrayRCP<SC> tmplocalQR(localQRsize);
        for (size_t i=0; i<myAggSize; i++) {
          // loop over cols
          for (size_t j=0; j<NSDim; j++) {
            tmplocalQR[ myAggSize*j + i] = localQR[myAggSize*j + i];
          }
        }

        // copy temporary variables back to correct positions (column size is now NSDim instead of myAggSize
        for (size_t j=0; j<NSDim; j++) {  // loop over cols
          for (size_t i=0; i<myAggSize; i++) {  // loop over rows
            localQR[ NSDim*j + i] = tmplocalQR[myAggSize*j + i];
          }
          // fill NSDim-myAggSize rows with default null space
          for(size_t i=myAggSize; i<NSDim; i++) { // loop over rows
             if(i==j) localQR[NSDim*j+i] = 1.0;
             else localQR[NSDim*j+i] = 0.0;
          }
        }

        // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
        // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
        // This stores the (offset+k)th entry only if it is local according to the coarseMap.

        for (size_t j=0; j<NSDim; ++j) {
          for (size_t k=0; k<=j; ++k) {
            try {
              if (coarseMap->isNodeLocalElement(offset+k)) {
                  coarseNS[j][offset+k] = localQR[ NSDim*j + k ]; // agg has only one node
              }
            }
            catch(...) {
              std::cout << "caught error in coarseNS insert, j="<<j<<", offset+k = "<<offset+k<<std::endl;
            }
          }
        }

#if 0
      std::cout << "R" << std::endl;
      // loop over rows
      for (size_t i=0; i<NSDim; i++) {
        // loop over cols
        for (size_t j=0; j<NSDim; j++) {
          if(j>=i) { std::cout << localQR[ NSDim*j + i]; std::cout << "\t"; }
          else
            std::cout << "00\t";
          //std::cout << localQR[ NSDim*j + i]; std::cout << "\t";
        }
        std::cout << std::endl;
      }
#endif

        // Calculate Q, the tentative prolongator.
        // The Lapack GEQRF call only works for myAggsize >= NSDim
        // special handling for very small aggregates (with myAggsize < NSDim)

        // calculate identity matrix with zero columns for the last NSDim-myAggsize columns.
        for (size_t i=0; i<Teuchos::as<size_t>(myAggSize); i++) {
          // loop over cols
          for (size_t j=0; j<NSDim; j++) {
            if (j==i) localQR[ myAggSize*j + i] = 1.0;
            else localQR[ myAggSize*j + i] = 0.0;
          }
        }
      } // end else (special handling for 1pt aggregates)

#if 0
      std::cout << "Q" << std::endl;
      // loop over rows
      for (size_t i=0; i<myAggSize; i++) {
        // loop over cols
        for (size_t j=0; j<NSDim; j++) {
          std::cout << localQR[ myAggSize*j + i]; std::cout << "\t";
        }
        std::cout << std::endl << std::endl;
      }
#endif

      //Process each row in the local Q factor.  If the row is local to the current processor
      //according to the rowmap, insert it into Ptentative.  Otherwise, save it in ghostQ
      //to be communicated later to the owning processor.
      //FIXME -- what happens if maps are blocked?
      for (GO j=0; j<myAggSize; ++j) {
        //This loop checks whether row associated with current DOF is local, according to rowMapForPtent.
        //If it is, the row is inserted.  If not, the row number, columns, and values are saved in
        //MultiVectors that will be sent to other processors.
        GO globalRow = aggToRowMap[agg][j];
        //TODO is the use of Xpetra::global_size_t below correct?
        if( rowMapForPtent->isNodeGlobalElement(globalRow) == false )
        {
          ghostQrows[qctr] = globalRow;
          for (size_t k=0; k<NSDim; ++k) {
            ghostQcols[k][qctr] = coarseMap->getGlobalElement(agg*NSDim+k);
            ghostQvals[k][qctr] = localQR[k*myAggSize+j];
          }
          ++qctr;
        } else {
          LO nnz=0;
          for (size_t k=0; k<NSDim; ++k) {
            try{
              if (localQR[k*myAggSize+j] != 0.) {
                colPtr[nnz] = coarseMap->getGlobalElement(agg * NSDim + k);
                valPtr[nnz] = localQR[k*myAggSize+j];
                ++nnz;
              }
            }
            catch(...) {
              std::cout << "caught error in colPtr/valPtr insert, current index="<<nnz<<std::endl;
            }
          } //for (size_t k=0; k<NSDim; ++k)

          try{
            Ptentative->insertGlobalValues(globalRow,
                                           colPtr.view(0,nnz),
                                           valPtr.view(0,nnz));
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

    // Import ghost parts of Q factors and insert into Ptentative.
    // First import just the global row numbers.
    RCP<Xpetra::Vector<GO,LO,GO,Node> > targetQrowNums = Xpetra::VectorFactory<GO,LO,GO,Node>::Build(rowMapForPtent);
    targetQrowNums->putScalar(-1);
    targetQrowNums->doImport(*ghostQrowNums,*importer,Xpetra::INSERT);
    ArrayRCP< GO > targetQrows = targetQrowNums->getDataNonConst(0);

    // Now create map based on just the row numbers imported.
    Teuchos::Array<GO> gidsToImport;
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

    valPtr = ArrayRCP<SC>(NSDim,0.);
    colPtr = ArrayRCP<GO>(NSDim,0);
    for (typename Array<GO>::iterator r=gidsToImport.begin(); r!=gidsToImport.end(); ++r) {
      if (targetQvalues->getLocalLength() > 0) {
        for (size_t j=0; j<NSDim; ++j) {
          valPtr[j] = targetQvals[j][reducedMap->getLocalElement(*r)];
          colPtr[j] = targetQcols[j][reducedMap->getLocalElement(*r)];
        }
        Ptentative->insertGlobalValues(*r, colPtr.view(0,NSDim), valPtr.view(0,NSDim));
      } //if (targetQvalues->getLocalLength() > 0)
    }

    Ptentative->fillComplete(coarseMap,fineA.getDomainMap()); //(domain,range) of Ptentative
    
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
