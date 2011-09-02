#ifndef MUELU_TENTATIVEPFACTORY_HPP
#define MUELU_TENTATIVEPFACTORY_HPP

#include <iostream>

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_ImportFactory.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp" // TODO: inheritence of TentativePFactory
#include "MueLu_Level.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class TentativePFactory class.
    @brief Factory for building tentative prolongator.

    Factory for creating tentative prolongator.   Nullspace vectors are split across aggregates so that they
    have local support.  The vectors with local support are factored via LAPACK QR.  The Q becomes the
    tentative prolongator, and the R becomes the coarse nullspace. 
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class TentativePFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{
    
    /*! @brief Constructor.
      \param AggregationFact -- (optional) factory that creates aggregates.
    */
    TentativePFactory(RCP<FactoryBase> aggregatesFact = Teuchos::null, RCP<FactoryBase> nullspaceFact = Teuchos::null, RCP<FactoryBase> AFact = Teuchos::null)
      : aggregatesFact_(aggregatesFact), nullspaceFact_(nullspaceFact), AFact_(Teuchos::null),
        QR_(false) { }
    
    //! Destructor.
    virtual ~TentativePFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level & fineLevel, Level & coarseLevel) const {
      fineLevel.Input("A", AFact_());                   // TAW: what's this??
      fineLevel.Input("Aggregates", aggregatesFact_());
      fineLevel.Input("Nullspace",  nullspaceFact_());

      fineLevel.Request("A", AFact_);
      fineLevel.Request("Aggregates", aggregatesFact_);
      fineLevel.Request("Nullspace",  nullspaceFact_);

    }

    //@}

    //! @name Set/Get methods.
    //@{

    void TentativeWithQR(bool value) { QR_ = value; }

    bool TentativeWithQR() { return QR_; }

    //@}

    //! @name Build methods.
    //@{

    bool Build(Level & fineLevel, Level & coarseLevel) const { //TODO
      return BuildP(fineLevel, coarseLevel);
    }

    bool BuildP(Level & fineLevel, Level & coarseLevel) const {
      // Level Get
      RCP<Operator>    A          = fineLevel.Get< RCP<Operator> >("A", AFact_);
      RCP<Aggregates>  aggregates = fineLevel.Get< RCP<Aggregates> >("Aggregates", aggregatesFact_);
      RCP<MultiVector> nullspace  = fineLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_);

      // Build
      std::ostringstream buf; buf << coarseLevel.GetLevelID(); //TODO remove/hide
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TentativePFactory::MakeTentative_"+buf.str()));
      timer->start(true);

      RCP<MultiVector> coarseNullspace; RCP<Operator> Ptentative; // output of MakeTentative()
      MakeTentative(*A, *aggregates, *nullspace, coarseNullspace, Ptentative);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(A->getRowMap()->getComm()));

      // Level Set
      coarseLevel.Set("Nullspace", coarseNullspace, nullspaceFact_);
      coarseLevel.Set("Ptent", Ptentative, this); //TODO: should be P when 'extended needs' implemented

      return true;
    }

    //@}

  private:
    //! @name Static methods.
    //@{
    
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;
    
    /*! @brief Make tentative prolongator with QR.

    - FIXME There is no attempt to detect if the aggregate is too small to support the NS.
    */
     void MakeTentative(const Operator& fineA, const Aggregates& aggregates, const MultiVector & fineNullspace, //-> INPUT
                        RCP<MultiVector> & coarseNullspace, RCP<Operator> & Ptentative) const                   //-> OUTPUT 
    {

      RCP<const Teuchos::Comm<int> > comm = fineA.getRowMap()->getComm();

      GO numAggs = aggregates.GetNumAggregates();

      const size_t NSDim = fineNullspace.getNumVectors();
      GO nCoarseDofs = numAggs*NSDim;
      // Compute array of aggregate sizes.
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
      ArrayRCP< ArrayRCP<LO> > aggToRowMap(numAggs);

      aggregates.ComputeAggregateToRowMap(aggToRowMap);

      // Create the numbering for the new row map for Ptent as follows:
      // convert LIDs in aggToRowmap to GIDs, put them in an ArrayView,
      // and call the Map constructor that takes arbitrary distributions.

      //

      std::vector<GO> globalIdsForPtent;
      for (int i=0; i< aggToRowMap.size(); ++i) {
        for (int k=0; k< aggToRowMap[i].size(); ++k) {
          //just for debugging
          //RCP<const Map> colMap = fineA.getColMap();
          //LO lid = aggToRowMap[i][k];
          //LO mylen = colMap->getNodeNumElements(); //error happens here ... problem with Views?
          //GO gid = colMap->getGlobalElement(lid);
          //globalIdsForPtent.push_back(gid);
          globalIdsForPtent.push_back(fineA.getColMap()->getGlobalElement(aggToRowMap[i][k]));
        }
      }

      Teuchos::ArrayView<GO> gidsForPtent(&globalIdsForPtent[0],globalIdsForPtent.size());

      LO indexBase=fineA.getRowMap()->getIndexBase();

      RCP<const Map > rowMapForPtent =
        MapFactory::Build(
                          fineA.getRowMap()->lib(),
                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                          gidsForPtent,
                          indexBase,
                          fineA.getRowMap()->getComm()
                          );

      // Allocate workspace for LAPACK QR routines.
      ArrayRCP<SC> localQR(maxAggSize*NSDim); // The submatrix of the nullspace to be orthogonalized.
      LO           workSize = NSDim;          // Length of work. Must be at least dimension of nullspace.
                                              // LAPACK may calculate better value, returned in work[0].
      ArrayRCP<SC> work(workSize);            // (in/out) work vector
                                              // FIXME DGEQRF documentation says this should be min(M,N) where B=MxN
      ArrayRCP<SC> tau(NSDim);                // (out) scalar factors of elementary reflectors, input to DORGQR
      LO           info=0;                      // (out) =0: success; =i, i<0: i-th argument has illegal value

      //Allocate storage for the coarse nullspace.

      RCP<const Map > coarseMap = MapFactory::Build(fineA.getRowMap()->lib(),
                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), nCoarseDofs,
                                                    indexBase, fineA.getRowMap()->getComm()); //JG:Xpetra::global_size_t>?

      coarseNullspace = MultiVectorFactory::Build(coarseMap,NSDim);
      ArrayRCP< ArrayRCP<SC> > coarseNS(NSDim);
      for (size_t i=0; i<NSDim; ++i)
        coarseNS[i] = coarseNullspace->getDataNonConst(i);

      // Builds overlapped nullspace.
      const RCP<const Map> nonUniqueMap = aggregates.GetMap();
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

      // Ptentative's row map is such that all DoF's in an aggregate (NO LONGER TRUE)
      // are local and consecutive.  Use fineA's domain map to FillComplete. (STILL TRUE)

      //This makes the rowmap of Ptent the same as that of fineA.
      //This requires moving some parts of some local Q's to other processors
      //because aggregates can span processors.
      rowMapForPtent = fineA.getRowMap();

      Ptentative = rcp(new CrsOperator(rowMapForPtent, NSDim));

      // Set up storage for the rows of the local Qs that belong to other processors.
      // FIXME This is inefficient and could be done within the main loop below with std::vector's.
      RCP<const Map> colMap = fineA.getColMap();
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
      RCP<MultiVector> ghostQcolumns = MultiVectorFactory::Build(ghostQMap,NSDim);
      RCP<MultiVector> ghostQrowNums = MultiVectorFactory::Build(ghostQMap,1);
      ArrayRCP< ArrayRCP<SC> > ghostQvals;
      ArrayRCP< ArrayRCP<SC> > ghostQcols;
      ArrayRCP< SC > ghostQrows;
      if (ghostQvalues->getLocalLength() > 0) {
        ghostQvals.resize(NSDim);
        ghostQcols.resize(NSDim);
        for (size_t i=0; i<NSDim; ++i) {
          ghostQvals[i] = ghostQvalues->getDataNonConst(i);
          ghostQcols[i] = ghostQcolumns->getDataNonConst(i);
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
            for (LO k=0; k<myAggSize; ++k) {
              //aggToRowMap[i][k] is the kth DOF in the ith aggregate
              //fineNS[j][n] is the nth entry in the jth NS vector
              try{
                localQR[j* myAggSize + k] = fineNS[j][ aggToRowMap[agg][k] ];
              }
              catch(...) {
                std::cerr << "caught an error!" << std::endl;
              }
            } //for (LO k=0 ...
          } //for (LO j=0 ...

          int intFineNSDim = NSDim;

          if (NSDim == 1) {
            //only one nullspace vector, so normalize by hand
            Magnitude dtemp=0;
            for (LO k=0; k<myAggSize; ++k) {dtemp += localQR[k]*localQR[k];}
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
          // so we cast to avoid compiler warnings.
          if ( work[0] > workSize) {
            workSize = (int) work[0];
            work = ArrayRCP<SC>(workSize);
          } else
            workSize = (int) work[0]; //TODO huh, think about this -- should it ever shrink?

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
            Magnitude dtemp = localQR[0];
            localQR[0] = tau[0];
            dtemp = 1 / dtemp;
            for (LO i=0; i<myAggSize; ++i)
              localQR[i] *= dtemp;
          } else {
            lapack.ORGQR( myAggSize, intFineNSDim, intFineNSDim, localQR.getRawPtr(),
                          myAggSize, tau.getRawPtr(), work.getRawPtr(), workSize, &info );
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
          if ( work[0] > workSize) {
            workSize = (int) work[0];
            work = ArrayRCP<SC>(workSize);
          } else
            workSize = (int) work[0]; //TODO huh, think about this -- should it ever shrink?

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
                ghostQrows[qctr] = fineA.getColMap()->getGlobalElement(localRow);
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
      RCP<MultiVector> targetQrowNums = MultiVectorFactory::Build(rowMapForPtent,1);
      targetQrowNums->putScalar(-1.0);
      targetQrowNums->doImport(*ghostQrowNums,*importer,Xpetra::INSERT);
      ArrayRCP< SC > targetQrows = targetQrowNums->getDataNonConst(0);

      // Now create map based on just the row numbers imported.
      Teuchos::Array<GO> gidsToImport;
      for (typename ArrayRCP<SC>::iterator r=targetQrows.begin(); r!=targetQrows.end(); ++r) {
        if (*r > -1.0) {
          gidsToImport.push_back((GO)*r);
        }
      }
      RCP<const Map > reducedMap = MapFactory::Build( fineA.getRowMap()->lib(),
                                                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                      gidsToImport, indexBase, fineA.getRowMap()->getComm()    );

      // Import using the row numbers that this processor will receive.
      importer = ImportFactory::Build(ghostQMap, reducedMap);

      RCP<MultiVector> targetQcolumns = MultiVectorFactory::Build(reducedMap,NSDim);
      targetQcolumns->doImport(*ghostQcolumns,*importer,Xpetra::INSERT);
      RCP<MultiVector> targetQvalues = MultiVectorFactory::Build(reducedMap,NSDim);
      targetQvalues->doImport(*ghostQvalues,*importer,Xpetra::INSERT);

      ArrayRCP< ArrayRCP<SC> > targetQvals;
      ArrayRCP<ArrayRCP<SC> > targetQcols;
      if (targetQvalues->getLocalLength() > 0) {
        targetQvals.resize(NSDim);
        targetQcols.resize(NSDim);
        for (size_t i=0; i<NSDim; ++i) {
          targetQvals[i] = targetQvalues->getDataNonConst(i);
          targetQcols[i] = targetQcolumns->getDataNonConst(i);
        }
      }

      for (typename Array<GO>::iterator r=gidsToImport.begin(); r!=gidsToImport.end(); ++r) {
        for (size_t j=0; j<NSDim; ++j) {
          valPtr[j] = targetQvals[j][reducedMap->getLocalElement(*r)];
          colPtr[j] = (GO) targetQcols[j][reducedMap->getLocalElement(*r)];
        }
        Ptentative->insertGlobalValues((GO)*r, colPtr.view(0,NSDim), valPtr.view(0,NSDim));
      }

      Ptentative->fillComplete(coarseMap,fineA.getDomainMap()); //(domain,range) of Ptentative

    } //MakeTentative()

    static void BuildAggregates() {
      throw(Exceptions::NotImplemented("TentativePFactory: BuildAggregates not implemented"));
    }
    static void MakeNoQRTentative() {
      throw(Exceptions::NotImplemented("TentativePFactory: MakeNoQRTentative not implemented"));
    }
    //@}

  private:
    RCP<FactoryBase> aggregatesFact_; //! Factory that creates aggregates
    RCP<FactoryBase> nullspaceFact_;
    RCP<FactoryBase> AFact_;

    bool QR_; //! use QR decomposition for improving nullspace information per default

  }; //class TentativePFactory

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_SHORT

#endif //ifndef MUELU_TENTATIVEPFACTORY_HPP
