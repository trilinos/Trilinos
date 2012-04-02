#ifndef MUELU_TENTATIVEPFACTORY_DEF_HPP
#define MUELU_TENTATIVEPFACTORY_DEF_HPP

#include <Teuchos_LAPACK.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsOperator.hpp>

#include "MueLu_TentativePFactory_decl.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

#undef OLDCODE

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TentativePFactory(RCP<const FactoryBase> aggregatesFact, RCP<const FactoryBase> nullspaceFact, RCP<const FactoryBase> AFact)
    : aggregatesFact_(aggregatesFact), nullspaceFact_(nullspaceFact), AFact_(AFact),
      QR_(false),
      domainGidOffset_(0) {
  }
 
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TentativePFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {

    fineLevel.DeclareInput("A", AFact_.get(), this);
    fineLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
    fineLevel.DeclareInput("Nullspace",  nullspaceFact_.get(), this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TentativeWithQR(bool value) { QR_ = value; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TentativeWithQR() { return QR_; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const { //TODO
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level & fineLevel, Level & coarseLevel) const {

    // get data from fine level
    RCP<Operator>    A          = fineLevel.Get< RCP<Operator> >("A", AFact_.get());
    RCP<Aggregates>  aggregates = fineLevel.Get< RCP<Aggregates> >("Aggregates", aggregatesFact_.get());
    RCP<MultiVector> nullspace  = fineLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());

    FactoryMonitor m(*this, "Tentative prolongator", coarseLevel);

    // Build
    RCP<MultiVector> coarseNullspace; RCP<Operator> Ptentative; // output of MakeTentative()
    MakeTentative(*A, *aggregates, *nullspace, coarseNullspace, Ptentative);

    // Level Set
    coarseLevel.Set("Nullspace", coarseNullspace, nullspaceFact_.get()); //FIXME !!!!
    //coarseLevel.Set("Nullspace", coarseNullspace, this);
    coarseLevel.Set("P", Ptentative, this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDomainMapOffset(GlobalOrdinal offset) {
    TEUCHOS_TEST_FOR_EXCEPTION(offset < 0, Exceptions::RuntimeError, "MueLu::TentativePFactory::SetDomainMapOffset: domain map offset for coarse gids of tentative prolongator must not be smaller than zero. Error.");
    domainGidOffset_ = offset;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDomainMapOffset() const {
    return domainGidOffset_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MakeTentative(
                     const Operator& fineA, const Aggregates& aggregates, const MultiVector & fineNullspace,
                     RCP<MultiVector> & coarseNullspace, RCP<Operator> & Ptentative) const
  {
    RCP<const Teuchos::Comm<int> > comm = fineA.getRowMap()->getComm();

    // number of aggregates
    GO numAggs = aggregates.GetNumAggregates();

    // dimension of fine level nullspace
    const size_t NSDim = fineNullspace.getNumVectors();

    // number of coarse level dofs (fixed by number of aggregats and nullspace dimension)
    GO nCoarseDofs = numAggs*NSDim;

    // Compute array of aggregate sizes (in dofs).
    ArrayRCP<LO> aggSizes  = aggregates.ComputeAggregateSizes();

    // Calculate total #dofs in local aggregates, find size of the largest aggregate.
    LO maxAggSize=0;
    LO numDofsInLocalAggs=0;
    for (typename Teuchos::ArrayRCP<LO>::iterator i=aggSizes.begin(); i!=aggSizes.end(); ++i) {
      if (*i > maxAggSize) maxAggSize = *i;
      numDofsInLocalAggs += *i;
    }

    // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
    // aggToRowMap[i][j] is the jth DOF in aggregate i
    // TODO: aggToRowMap lives in the column map of A (with overlapping). Note: ComputeAggregateToRowMap
    // returns the local DOFs, that are transformed to global Dofs using the col map later. Wouldn't it be
    // smarter to compute the global dofs in ComputeAggregateToRowMap?
    ArrayRCP< ArrayRCP<LO> > aggToRowMap(numAggs);
    aggregates.ComputeAggregateToRowMap(aggToRowMap);

    // Create the numbering for the new row map for Ptent as follows:
    // convert LIDs in aggToRowmap to GIDs, put them in an ArrayView,
    // and call the Map constructor that takes arbitrary distributions.

    //

    RCP<const Map> colMap = fineA.getColMap();
#ifdef OLDCODE
    // this is not necessary since rowMapForPtent is set to A->RowMap later
    // do we need an overlapping row map for Ptent?? i think communication is done later
    // with the null space and the QR decomposition
    std::vector<GO> globalIdsForPtent;
    for (int i=0; i< aggToRowMap.size(); ++i) {
      for (int k=0; k< aggToRowMap[i].size(); ++k) {
        //just for debugging
        //RCP<const Map> colMap = fineA.getColMap();
        //LO lid = aggToRowMap[i][k];
        //LO mylen = colMap->getNodeNumElements(); //error happens here ... problem with Views?
        //GO gid = colMap->getGlobalElement(lid);
        //globalIdsForPtent.push_back(gid);
        globalIdsForPtent.push_back(colMap()->getGlobalElement(aggToRowMap[i][k]));
      }
    }

    Teuchos::ArrayView<GO> gidsForPtent(&globalIdsForPtent[0],globalIdsForPtent.size());
#endif

    GO indexBase=fineA.getRowMap()->getIndexBase();

#ifdef OLDCODE
    RCP<const Map > rowMapForPtent =
    MapFactory::Build(
                      fineA.getRowMap()->lib(),
                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                      gidsForPtent,
                      indexBase,
                      fineA.getRowMap()->getComm()
                      );
#else
    RCP<const Map > rowMapForPtent = Teuchos::null; // TODO move me to line 230
#endif



    // Allocate workspace for LAPACK QR routines.
    ArrayRCP<SC> localQR(maxAggSize*NSDim); // The submatrix of the nullspace to be orthogonalized.
    LO           workSize = NSDim;          // Length of work. Must be at least dimension of nullspace.
    // LAPACK may calculate better value, returned in work[0].
    ArrayRCP<SC> work(workSize);            // (in/out) work vector
    // FIXME DGEQRF documentation says this should be min(M,N) where B=MxN
    ArrayRCP<SC> tau(NSDim);                // (out) scalar factors of elementary reflectors, input to DORGQR
    LO           info=0;                      // (out) =0: success; =i, i<0: i-th argument has illegal value

    // build coarse level maps (= domain map of transfer operator)
    RCP<const Map > coarseMap = Teuchos::null;
    if (domainGidOffset_ == 0) {
      // default: no offset for domain gids.
      // the gids for the domain Dofs are contiguously numbered starting from zero and equally distributed over all procs
      coarseMap = MapFactory::Build(fineA.getRowMap()->lib(),
                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), nCoarseDofs,
                                                    indexBase, fineA.getRowMap()->getComm()); //JG:Xpetra::global_size_t>?
    } else {

      // TODO: we need something smarter here. It would be nice to have strided maps for the coarse map
      // -> needed for blocked problems

      // the number of domain dof gids starts from domainGidOffset_ > 0
      // first build a tentative coarse Dof map (equally distributed, starting with Gids beginning with 0)
      RCP<const Map > tentativecoarseMap = MapFactory::Build(fineA.getRowMap()->lib(),
                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), nCoarseDofs,
                                                    indexBase, fineA.getRowMap()->getComm()); //JG:Xpetra::global_size_t>?

      // use same distribution as for tentativecoarseMap. shift Gids by adding domain gid offset
      std::vector<GO> globalDomainIdsForPtent(tentativecoarseMap->getNodeNumElements());
      for(GO dolid = 0; dolid < Teuchos::as<GO>(tentativecoarseMap->getNodeNumElements()); dolid++) {
        globalDomainIdsForPtent[dolid] = tentativecoarseMap->getGlobalElement(dolid) + domainGidOffset_;
      }
      Teuchos::ArrayView<GO> domainGidsForPtent(&globalDomainIdsForPtent[0],globalDomainIdsForPtent.size());
      coarseMap = MapFactory::Build(fineA.getRowMap()->lib(),
                                                    tentativecoarseMap->getGlobalNumElements(),//Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                    domainGidsForPtent,
                                                    indexBase,
                                                    fineA.getRowMap()->getComm()); //JG:Xpetra::global_size_t>?
    }

    //Allocate storage for the coarse nullspace.
    coarseNullspace = MultiVectorFactory::Build(coarseMap,NSDim);

    ArrayRCP< ArrayRCP<SC> > coarseNS(NSDim);
    for (size_t i=0; i<NSDim; ++i)
      if (coarseMap->getNodeNumElements() > 0) coarseNS[i] = coarseNullspace->getDataNonConst(i);

    // Builds overlapped nullspace.
    const RCP<const Map> nonUniqueMap = aggregates.GetDofMap(); // fetch overlapping Dofmap from Aggregates structure

    GO nFineDofs = nonUniqueMap->getNodeNumElements();
    const RCP<const Map> uniqueMap    = fineA.getDomainMap(); //FIXME won't work for systems
    RCP<const Import> importer = ImportFactory::Build(uniqueMap, nonUniqueMap);
    RCP<MultiVector> fineNullspaceWithOverlap = MultiVectorFactory::Build(nonUniqueMap,NSDim);
    fineNullspaceWithOverlap->doImport(fineNullspace,*importer,Xpetra::INSERT);

    // Pull out the nullspace vectors so that we can have random access.
    // (Question -- do we have to do this?)
    ArrayRCP< ArrayRCP<const SC> > fineNS(NSDim);
    for (size_t i=0; i<NSDim; ++i)
      fineNS[i] = fineNullspaceWithOverlap->getData(i);
 
    //Allocate temporary storage for the tentative prolongator.
    // TODO Right now, this is stored as a point matrix.
    // TODO we should be able to allocate something other than a CRS matrix
    ArrayRCP<GO> rowPtr(nFineDofs+1);
    for (GO i=0; i<=nFineDofs; ++i)
      rowPtr[i] = i*NSDim;
    ArrayRCP<GO> colPtr(maxAggSize*NSDim,0);
    ArrayRCP<SC> valPtr(maxAggSize*NSDim,0.);

    // Note! DoF's in an aggregate may not be local or consecutive.
    // Use fineA's domain map to FillComplete.

    //This makes the rowmap of Ptent the same as that of fineA.
    //This requires moving some parts of some local Q's to other processors
    //because aggregates can span processors.
    rowMapForPtent = fineA.getRowMap();

    Ptentative = rcp(new CrsOperator(rowMapForPtent, NSDim, Xpetra::StaticProfile));

    // Set up storage for the rows of the local Qs that belong to other processors.
    // FIXME This is inefficient and could be done within the main loop below with std::vector's.
    Array<GO> ghostGIDs;
    for (LO j=0; j<numAggs; ++j)
      {
        for (LO k=0; k<aggSizes[j]; ++k) {
          //TODO is the use of Xpetra::global_size_t below correct?
          if (rowMapForPtent->getGlobalElement(aggToRowMap[j][k]) == (GO) Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid() )
            {
              ghostGIDs.push_back(colMap->getGlobalElement(aggToRowMap[j][k]));
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

    //used in the case of just one nullspace vector
    Teuchos::Array<Magnitude> norms(NSDim);
    fineNullspace.norm2(norms);

    Teuchos::LAPACK<LO,SC> lapack;

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
          //aggToRowMap[i][k] is the kth DOF in the ith aggregate
          //fineNS[j][n] is the nth entry in the jth NS vector
          try{
            SC nsVal = fineNS[j][ aggToRowMap[agg][k] ]; // extract information from fine level NS
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
            std::cout << "fineNS...=" << fineNS[j][ aggToRowMap[agg][k] ] << std::endl;
            std::cerr << "caught an error!" << std::endl;
          }
        } //for (LO k=0 ...
        TEUCHOS_TEST_FOR_EXCEPTION(bIsZeroNSColumn == true, Exceptions::RuntimeError, "MueLu::TentativePFactory::MakeTentative: fine level NS part has a zero column. Error.");
      } //for (LO j=0 ...

      int intFineNSDim = NSDim;

      if (NSDim == 1) {
        //only one nullspace vector, so normalize by hand
        Magnitude dtemp=0;
        //scalar type might be complex, so take absolute value.
        for (LO k=0; k<myAggSize; ++k) {dtemp += std::abs(localQR[k])*std::abs(localQR[k]);}
        dtemp = Teuchos::ScalarTraits<Magnitude>::squareroot(dtemp);
        tau[0] = localQR[0];
        localQR[0] = dtemp;
      } else {
        //Perform the QR.  Upon return, R is stored explicitly, Q is stored implicitly as product
        //of reflection matrices.
        lapack.GEQRF( myAggSize, intFineNSDim, localQR.getRawPtr(), myAggSize,
                      tau.getRawPtr(), work.getRawPtr(), workSize, &info );
      }

      if (info != 0) {
        std::ostringstream tbuf;
        tbuf << info;
        std::string msg = "MakeTentativeWithQR: dgeqrf (LAPACK QR routine) returned error code " + tbuf.str();
        throw(Exceptions::RuntimeError(msg));
      }

      // LAPACK may have determined a better length for the work array.  Returns it in work[0],
      // so we cast to avoid compiler warnings.  Taking a look at the NETLIB reference implementation
      // CGEQRF (complex), work[0] is assigned an integer, so it's safe to take the magnitude.
      // Scalar type might be complex, so take absolute value.
      if ( std::abs(work[0]) > workSize) {
        workSize = (int) std::abs(work[0]);
        work = ArrayRCP<SC>(workSize);
      } else
        workSize = (int) std::abs(work[0]); //TODO huh, think about this -- should it ever shrink?

      // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
      // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
      // This stores the (offset+k)th entry only if it is local according to the coarseMap.
      Xpetra::global_size_t offset=agg*NSDim;
      for (size_t j=0; j<NSDim; ++j) {
        for (size_t k=0; k<=j; ++k) {
          try {
            if (coarseMap->isNodeLocalElement(offset+k))
              coarseNS[j][offset+k] = localQR[ myAggSize*j + k ]; //TODO is offset+k the correct local ID?!
          }
          catch(...) {
            std::cout << "caught error in coarseNS insert, j="<<j<<", offset+k = "<<offset+k<<std::endl;
          }
        }
      }

      // Calculate Q, the tentative prolongator, explicitly.  This requires calling a second LAPACK routine.

      if (NSDim == 1) {
        //again, only one nullspace vector, so calculate Q by hand
        Magnitude dtemp = std::abs(localQR[0]);
        localQR[0] = tau[0];
        dtemp = 1 / dtemp;
        for (LO i=0; i<myAggSize; ++i)
          localQR[i] *= dtemp;
      } else {
        ExtractQ( lapack, myAggSize, intFineNSDim, localQR, tau, work, workSize, info );
      }

      if (info != 0) {
        std::ostringstream tbuf;
        tbuf << info;
        std::string msg = "MakeTentativeWithQR: dorgqr (LAPACK auxiliary QR routine) returned error code "
          + tbuf.str();
        throw(Exceptions::RuntimeError(msg));
      }

      // LAPACK may have determined a better length for the work array.  Returns it in work[0],
      // so we cast to avoid compiler warnings.
      // Scalar type might be complex, so take absolute value.
      if ( std::abs(work[0]) > workSize) {
        workSize = (int) std::abs(work[0]);
        work = ArrayRCP<SC>(workSize);
      } else
        workSize = (int) std::abs(work[0]); //TODO huh, think about this -- should it ever shrink?

      //Process each row in the local Q factor.  If the row is local to the current processor
      //according to the rowmap, insert it into Ptentative.  Otherwise, save it in ghostQ
      //to be communicated later to the owning processor.
      //FIXME -- what happens if maps are blocked?
      for (GO j=0; j<myAggSize; ++j) {
        //This loop checks whether row associated with current DOF is local, according to rowMapForPtent.
        //If it is, the row is inserted.  If not, the row number, columns, and values are saved in
        //MultiVectors that will be sent to other processors.
        LO localRow = aggToRowMap[agg][j];
        //TODO is the use of Xpetra::global_size_t below correct?
        if (rowMapForPtent->getGlobalElement(localRow) == (GO) Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid())
        {
          ghostQrows[qctr] = colMap()->getGlobalElement(localRow);
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
            Ptentative->insertGlobalValues(rowMapForPtent->getGlobalElement(localRow),
                                           colPtr.view(0,nnz),
                                           valPtr.view(0,nnz));
          }
          catch(...) {
            std::cout << "pid " << fineA.getRowMap()->getComm()->getRank()
                      << "caught error during Ptent row insertion, local row "
                      << localRow << ", global row "
                      << rowMapForPtent->getGlobalElement(localRow)
                      << std::endl;
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

  } //MakeTentative()

  //! Non-member templated function to handle extracting Q from QR factorization for different Scalar types.
  template <class Scalar, class LocalOrdinal>
  void ExtractQ(Teuchos::LAPACK<LocalOrdinal,Scalar> &lapack, LocalOrdinal myAggSize,
                int intFineNSDim, ArrayRCP<Scalar> &localQR, ArrayRCP<Scalar> &tau,
                ArrayRCP<Scalar> &work, LocalOrdinal &workSize, LocalOrdinal &info)
  {
    lapack.ORGQR(myAggSize, intFineNSDim, intFineNSDim, localQR.getRawPtr(),
                 myAggSize, tau.getRawPtr(), work.getRawPtr(), workSize, &info );
  }

  //! Non-member specialized function to handle extracting Q from QR factorization for Scalar==complex
  template <class LocalOrdinal>
  void ExtractQ(Teuchos::LAPACK<LocalOrdinal, std::complex<double> > &lapack,
                LocalOrdinal myAggSize, int intFineNSDim, ArrayRCP<std::complex<double> > &localQR,
                ArrayRCP<std::complex<double> > &tau, ArrayRCP<std::complex<double> > &work,
                LocalOrdinal &workSize, LocalOrdinal &info)
  {
    lapack.UNGQR(myAggSize, intFineNSDim, intFineNSDim, localQR.getRawPtr(),
                 myAggSize, tau.getRawPtr(), work.getRawPtr(), workSize, &info );
  }

} //namespace MueLu

//TODO: noQR_

// TODO ReUse: If only P or Nullspace is missing, TentativePFactory can be smart and skip part of the computation.

#define MUELU_TENTATIVEPFACTORY_SHORT
#endif // MUELU_TENTATIVEPFACTORY_DEF_HPP
