#ifndef MUELU_UCAGGREGATIONCOMMHELPER_HPP
#define MUELU_UCAGGREGATIONCOMMHELPER_HPP

#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Aggregates.hpp"

namespace MueLu {

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class UCAggregationCommHelper : public BaseClass {

    typedef double Scalar; // Scalar type only used for weight: always a double.
#include "MueLu_UseShortNames.hpp"

  private:
    RCP<const Import> import_;
    mutable RCP<const Import> winnerImport_; //FIXME get rid of "mutable"
    mutable RCP<Import> pushWinners_; //FIXME get rid of mutable
    RCP<Vector> tempVec_;
    mutable RCP<Vector> perturbWt_;
    mutable RCP<Vector> postComm_;
    mutable RCP<Vector> candidateWinners_;
    mutable ArrayRCP<GO> myWinners_;
    mutable int numMyWinners_;
    mutable RCP<Map> winnerMap_;
    mutable int numCalls_;
    int myPID_;

    //     uniqueMap                A subset of weight.getMap() where each GlobalId 
    //                              has only one unique copy on one processor.
    //                              Normally, weight.getMap() would have both locals
    //                              and ghost elements while uniqueMap would just
    //                              have the locals. It should be possible to
    //                              remove this or make it an optional argument
    //                              and use some existing Epetra/Tpetra capability to 
    //                              make a uniqueMap.
    //
    //     import_                  This corresponds precisely to 
    //                                   Import import_(
    //                                           weight.getMap(), uniqueMap);
    //                              This could also be eliminated and created
    //                              here, but for efficiency user's must pass in.
    //
  public:

    UCAggregationCommHelper(const RCP<const Map> & uniqueMap, const RCP<const Map> & nonUniqueMap)
    {
      import_ = ImportFactory::Build(uniqueMap, nonUniqueMap);
      tempVec_ = VectorFactory::Build(uniqueMap,false); //zeroed out before use
      numMyWinners_ = 0;
      numCalls_ = 0;
      myPID_ = uniqueMap->getComm()->getRank();
    }

    ~UCAggregationCommHelper() { if (!myPID_) std::cout << "ArbitrateAndCommunicate has been called " << numCalls_ << " times" << std::endl;}

    inline void ArbitrateAndCommunicate(Vector &weights, Aggregates &aggregates, const bool perturb) const
    {
      ArbitrateAndCommunicate(weights, *aggregates.GetProcWinner(), &*aggregates.GetVertex2AggId(), perturb);
    }

    //
    // For each GlobalId associated with weight.getMap():
    //
    //      1) find the maximum absolute value of weight[] distributed across all
    //         processors and assign this to all local elements of weight[] (across 
    //         processors) associated with the GlobalId.
    //      2) set procWinner[] to the MyPid() that had the largest element.
    //         procWinner[] is still set if only one processor owns a GlobalId. 
    //
    //         The ONLY CASE when procWinner[i] is NOT set corresponds to when
    //         all local weights associated with a GlobalId are zero. This allows
    //         one to effectively skip the maximum/winner calculation for a subset
    //         of GlobalId's.  This might occur when a processor has already
    //         claimed ownership for a GlobalId and so all local copies have
    //         the same value. We want to skip the maximum calculation with 
    //         tiebreaking to avoid another processor claiming ownership.
    //
    //      3) optionally, set companion[] (across all relevant processors) to the
    //         local companion value associated with the procWinner[] processor.
    //
    //  Input:
    //     weight                   Vector of weights. ASSUMED TO BE nonnegative.
    //
    //     procWinner               Allocated but contents ignored.
    //
    //     companion                Either NULL or allocated but contents ignored.
    //                              If NULL, step 3 above is skipped.
    //

    //     perturb                  Optional arguments that is either true or 
    //                              false (default: true). weight is perturbed
    //                              and the perturbed values are used in step 1)
    //                              above. Returned values reflect the perturbed
    //                              data. This option avoids having lots of
    //                              tiebreaks where the large MyPid() always wins.
    //
    //  Output:
    //     weight                   weight[k] <-- Max(weight[k1],...,weight[kn])
    //                              where weight[kj] live on different processors
    //                              but have the same GlobalId as weight[k] on
    //                              this processor.
    //
    //     procWinner               procWinner[k] <-- MyPid associated with the
    //                              kj yielding the max in 
    //                                    Max(weight[k1],...,weight[kn]) .
    //                              See weight Output comments.
    //                              NOTE: If all input weight[kj]'s are zero,
    //                                    then procWinner[k] is left untouched.
    //
    //     companion                If not null, 
    //                                 companion[k] <-- companion[kj] where
    //                              companion[kj] lives on processor procWinner[k].
    //                              and corresponds to the same GlobalId as k.
    //                              NOTE: If for a particlar GlobalId, no processor
    //                                    has a value of procWinner that matches
    //                                    its MyPid, the corresponding companion
    //                                    is not altered.
    //
    void ArbitrateAndCommunicate(Vector &weight_, LOVector &procWinner_, LOVector *companion, const bool perturb) const
    {
      const RCP<const Map> weightMap = weight_.getMap();
      const size_t nodeNumElements = weightMap->getNodeNumElements();
      const RCP<const Teuchos::Comm<int> > & comm = weightMap->getComm();
      int MyPid = comm->getRank(); // TODO:remove the getMap() step
      ++numCalls_;

      //short-circuit if only one process
      if (comm->getSize() == 1) {
        ArrayRCP<SC> serialWeight = weight_.getDataNonConst(0);
        ArrayRCP<LO> serialProcWinner = procWinner_.getDataNonConst(0);
        for (size_t i=0; i < nodeNumElements; ++i) {
          if (serialWeight[i] > 0) {
            serialWeight[i] = 0;
            serialProcWinner[i] = MyPid;
          }
        }
        //companion doesn't change
        return;
      }

#ifdef COMPARE_IN_OUT_VECTORS
      RCP<Vector> in_weight_ = VectorFactory::Build(weight_.getMap());
      {
        ArrayRCP<SC> in_weight = in_weight_->getDataNonConst(0);
        ArrayRCP<SC> weight = weight_.getDataNonConst(0);
        for (size_t i=0; i < nodeNumElements; ++i) in_weight[i] = weight[i];
      }
      RCP<LOVector> in_procWinner_ = LOVectorFactory::Build(procWinner_.getMap());
      {
        ArrayRCP<LO> in_procWinner = in_procWinner_->getDataNonConst(0);
        ArrayRCP<LO> procWinner = procWinner_.getDataNonConst(0);
        for (size_t i=0; i < nodeNumElements; ++i) in_procWinner[i] = procWinner[i];
      }
      RCP<LOVector> in_companion;
      {
        if (companion != NULL) {
          in_companion = LOVectorFactory::Build(companion->getMap());
          ArrayRCP<LO> in_comp = in_companion->getDataNonConst(0);
          ArrayRCP<LO> comp = companion->getDataNonConst(0);
          for (size_t i=0; i < nodeNumElements; ++i) in_comp[i] = comp[i];
        }
      }
#endif

      if (perturb)
      {
        if (perturbWt_ == Teuchos::null || !perturbWt_->getMap()->isSameAs(*weightMap)) {
          perturbWt_ = VectorFactory::Build(weightMap,false); //no need to zero out because this will be randomized

          // Modify seed of the random algorithm used by perturbWt_->randomize()
          {
            ST::seedrandom(static_cast<unsigned int>(MyPid*2 + (int) (11*ST::random())));
            for (int i = 0; i < 10; ++i) ST::random();
              perturbWt_->setSeed(static_cast<unsigned int>(ST::random()));
          }
          perturbWt_->randomize(); 
          ArrayRCP<SC> lperturbWt = perturbWt_->getDataNonConst(0);
          for (size_t i=0; i < nodeNumElements; ++i)
            lperturbWt[i] = 1e-7*fabs(lperturbWt[i]); //FIXME this won't work for general SC
#ifdef COMPARE_IN_OUT_VECTORS
          ArrayRCP<SC> locperturbWt = perturbWt_->getDataNonConst(0);
          for (size_t i=0; i < nodeNumElements; ++i)
            printf("perturbWt[%d] = %15.10e\n",i,locperturbWt[i]);
#endif
        } //if (perturbWt_ == Teuchos::null || ...

        ArrayRCP<SC> weight = weight_.getDataNonConst(0); // TODO: const?
        ArrayRCP<SC> perturbWt = perturbWt_->getDataNonConst(0);

        // Note: maxValue() not available for Tpetra
        //SC largestGlobalWeight = weight_.maxValue();
        SC largestGlobalWeight = weight_.normInf();
        for (size_t i=0; i < nodeNumElements; ++i) {
          if (weight[i] != 0.) {
            weight[i] += largestGlobalWeight*perturbWt[i];
          }
        }
        //TODO is it necessary to return the *perturbed* weights?
      } //if (perturb)
  
      // Communicate weights and store results in PostComm (which will be copied
      // back into weights later. When multiple processors have different weights
      // for the same GID, we take the largest weight. After this fragment every
      // processor should have the same value for PostComm[] even when multiple
      // copies of the same Gid are involved.

      if (postComm_ == Teuchos::null || !postComm_->getMap()->isSameAs(*weightMap) )
        postComm_ = VectorFactory::Build(weightMap);

      //note: postComm_ is zeroed either in build above, or in loop below upon last touch.

      NonUnique2NonUnique(weight_, *postComm_, Xpetra::ABSMAX);
   
      // Let every processor know who is the procWinner. For nonunique
      // copies of the same Gid, this corresponds to the processor with
      // the highest Wt[]. When several processors have the same positive value
      // for weight[] (which is also the maximum value), the highest proc id
      // is declared the procWinner.
      //
      // Note:This is accomplished by filling a vector with MyPid+1 if weight[k] is
      //      nonzero and PostComm[k]==weight[k]. NonUnique2NonUnique(...,AbsMax)
      //      is invoked to let everyone know the procWinner.
      //      One is then subtracted so that procWinner[i] indicates the 
      //      Pid of the winning processor.
      //      When all weight's for a GID are zero, the associated procWinner's
      //      are left untouched.

      if (candidateWinners_ == Teuchos::null || !candidateWinners_->getMap()->isSameAs(*weightMap) )
        candidateWinners_ = VectorFactory::Build(weightMap,false);
      //note: candidateWinners_ is initialized below

      ArrayRCP<SC> weight = weight_.getDataNonConst(0);

      {
        ArrayRCP<SC> candidateWinners = candidateWinners_->getDataNonConst(0);
        ArrayRCP<SC> postComm = postComm_->getDataNonConst(0);
        for (size_t i=0; i < nodeNumElements; ++i) {
          if (postComm[i] == weight[i]) candidateWinners[i] = (SC) MyPid+1;
          else                          candidateWinners[i] = 0;
          weight[i]=postComm[i]; 
        }
      }
      NonUnique2NonUnique(*candidateWinners_, *postComm_, Xpetra::ABSMAX);

      // Note: 
      //                      associated CandidateWinners[]
      //    weight[i]!=0  ==> on some proc is equal to its ==> postComm[i]!=0
      //                      MyPid+1.
      //          
      int numMyWinners = 0;
      ArrayRCP<LO> procWinner = procWinner_.getDataNonConst(0);
      {
        ArrayRCP<SC> postComm = postComm_->getDataNonConst(0);
        for (size_t i=0; i < nodeNumElements; ++i)  {
          if ( weight[i] != 0.) procWinner[i] = ((int) (postComm[i])) - 1;
          weight[i] = 0.;    //we are done with weight
          postComm[i] = 0.;  //avoids having to initialize postComm_ on next call to ArbitrateAndCommunicate
          if (procWinner[i] == MyPid) ++numMyWinners;
        }
      }

      weight = Teuchos::null; //TODO why do we do this?

      if (companion != NULL) {
        // Now build a new Map, WinnerMap which just consists of procWinners. 
        // This is done by extracting the Gids for Wt, and shoving
        // the subset that correspond to procWinners in MyWinners.
        // WinnerMap is then constructed using MyWinners.
        //
        // In order to avoid regenerating winnerMap_, the following are checked:
        //   1) Do the local number of entries in MyWinners differ?  If so, regenerate/repopulate MyWinners and regenerate winnerMap_.
        //   2) If the local number of entries in MyWinners are the same, do any entries differ?  If so, repopulate MyWinners and
        //      regenerate winnerMap_.
   
        ArrayView<const GO> myGids = weightMap->getNodeElementList(); //== weightMap->MyGlobalElements(myGids);
        bool realloc=false;
        if (numMyWinners != numMyWinners_ || winnerMap_ == Teuchos::null) {
          // The local number of entries in MyWinners_ have changed since the last invocation, so reallocate myWinners_.
          myWinners_ = ArrayRCP<GO>(numMyWinners);
          realloc=true;
          //std::cout << MyPid << ":  numMyWinners has changed : (old) " << numMyWinners_ << ", (new) " << numMyWinners << std::endl;
          numMyWinners_ = numMyWinners;
        }

#ifdef JG_DEBUG
        procWinner = Teuchos::null;
        std::cout << MyPid << ": nodeNumElements=" << nodeNumElements << std::endl;
        std::cout << MyPid << ": procWinner=" << procWinner_ << std::endl;
        procWinner = procWinner_.getDataNonConst(0);
#endif

        if (realloc==true) {
          // The local number of entries in MyWinners have changed since the last invocation, so repopulate MyWinners_.
          numMyWinners = 0;
          for (size_t i = 0; i < nodeNumElements; ++i) {
            if (procWinner[i] == MyPid) {
              myWinners_[numMyWinners++] = myGids[i];
            }
          }
        } else {
          // The local number of entries in MyWinners are the same as the last invocation, but
          // we still must check if any entries differ from the last invocation.
          bool entryMismatch=false;
          numMyWinners = 0;
          for (size_t i = 0; i < nodeNumElements; ++i) {
            if (procWinner[i] == MyPid) {
              if (myWinners_[numMyWinners++] != myGids[i]) {
                entryMismatch=true;
                break;
              }
            }
          }

          if (entryMismatch == true) {
             // Entries differ from last invocation, so repopulate myWinners_.
            realloc=true;
            numMyWinners = 0;
            for (size_t i = 0; i < nodeNumElements; ++i) {
              if (procWinner[i] == MyPid) {
                myWinners_[numMyWinners++] = myGids[i];
              }
            }
          }
        } //if (realloc==true) ... else

        procWinner = Teuchos::null;

#ifdef JG_DEBUG
        std::cout << MyPid << ": numMyWinners=" << numMyWinners << std::endl;
        std::cout << MyPid << ": myWinners_" << myWinners_ << std::endl;
        for(int i=0;i<numMyWinners; i++)
          std::cout << MyPid << ": myWinners_[locId=" << i << "] = " << myWinners_[i] << std::endl;

#endif

#ifdef HAVE_MPI
        //See whether any process has determined that winnerMap_ must be regenerated.
        int irealloc,orealloc;
        if (realloc) irealloc=1;
        else         irealloc=0;
        MPI_Allreduce(&irealloc,&orealloc,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);  //FIXME MPI group should come from comm, may not be WORLD
        if (orealloc == 1) realloc=true;
        else               realloc=false;
#endif

        if (realloc) {
          // Either the number of entries or the value have changed since the last invocation, so reallocation the map.
          const Xpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
          winnerMap_ = MapFactory::Build(weightMap->lib(), GSTI, myWinners_(), 0, weightMap->getComm());
        }
       
        // Pull the Winners out of companion
        //     JustWinners <-- companion[Winners];
   
        RCP<LOVector> justWinners = LOVectorFactory::Build(winnerMap_);

#ifdef JG_DEBUG
        RCP<Teuchos::FancyOStream> out = rcp(new Teuchos::FancyOStream(rcp(&std::cout,false)));
        std::cout << MyPid << ": justWinners(Vector in)=" << *justWinners << std::endl;
        justWinners->describe(*out, Teuchos::VERB_EXTREME);
#endif

        if ( winnerImport_ == Teuchos::null
             || !winnerImport_->getSourceMap()->isSameAs(*weightMap)
             || !winnerImport_->getTargetMap()->isSameAs(*winnerMap_)  )
          winnerImport_ = ImportFactory::Build(weightMap, winnerMap_);
        RCP<const Import> winnerImport = winnerImport_;
        try
          {
            justWinners->doImport(*companion, *winnerImport, Xpetra::INSERT);
          }
        catch(std::exception& e)
          {
            std::cout << MyPid << ": ERR2: An exception occurred." << std::endl;
            throw e;
          }

        // Put the JustWinner values back into companion so that
        // all nonunique copies of the same Gid have the procWinner's
        // version of the companion.
        //#define JG_DEBUG
#ifdef JG_DEBUG
        RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        fos->setOutputToRootOnly(-1);
        if (!weightMap->getComm()->getRank())
          std::cout << "------ winnerMap_ ------" << std::endl;
        winnerMap_->describe(*fos,Teuchos::VERB_EXTREME);
        if (!weightMap->getComm()->getRank())
          std::cout << "------ weightMap ------" << std::endl;
        weightMap->getComm()->barrier();
        weightMap->describe(*fos,Teuchos::VERB_EXTREME);
        //std::cout << *winnerMap_ << std::endl;
        //std::cout << *weightMap << std::endl;
        sleep(5);
        exit(1);
#endif
#ifdef JG_DEBUG
#undef JG_DEBUG
#endif

        if ( pushWinners_ == Teuchos::null
             || !pushWinners_->getSourceMap()->isSameAs(*winnerMap_)
             || !pushWinners_->getTargetMap()->isSameAs(*weightMap)  )
          pushWinners_ = ImportFactory::Build(winnerMap_,weightMap);
        RCP<Import> pushWinners = pushWinners_;
        //RCP<Import> pushWinners = ImportFactory::Build(winnerMap_, weightMap); // VERSION1
        //RCP<Export> pushWinners = ExportFactory::Build(winnerMap_, weightMap); // VERSION4
        try
          {
            companion->doImport(*justWinners, *pushWinners, Xpetra::INSERT);   // VERSION1 Slow
            //companion->doExport(*justWinners, *winnerImport_, Xpetra::INSERT);   // JJH this should work... but exception
//             if (weightMap->lib() == Xpetra::UseEpetra)
//               justWinners->doExport(*companion, *winnerImport, Xpetra::INSERT);  // VERSION2 Tpetra doc is wrong
//             else if (weightMap->lib() == Xpetra::UseTpetra)
//               companion->doExport(*justWinners, *winnerImport, Xpetra::INSERT);     // VERSION3 - TODO: will certainly not work with Epetra? (change Xpetra?)
	    //companion->doExport(*justWinners, *pushWinners, Xpetra::INSERT);     // VERSION4
//             else throw "lib()";
          }
        catch(std::exception& e)
          {
            throw e;
          }
        //#define JG_DEBUG
#ifdef JG_DEBUG
        //            RCP<Teuchos::FancyOStream> out = rcp(new Teuchos::FancyOStream(rcp(&std::cout,false)));
            //->describe(*out, Teuchos::VERB_EXTREME);

            // std::cout << MyPid << ": ERR3: An exception occurred." << std::endl;

            std::cout << MyPid << ": numMyWinners=" << numMyWinners << std::endl;
            
            std::cout << MyPid << ": justWinners(Vector in)=" << std::endl;
            justWinners->describe(*out, Teuchos::VERB_EXTREME);

            std::cout << MyPid << ": companion(Vector out)=" << std::endl;
            companion->describe(*out, Teuchos::VERB_EXTREME);

            // std::cout << MyPid << ": pushWinners(Import(winnerMap_, weight_.Map))=" << *pushWinners << std::endl;
            std::cout << MyPid << ": winnerMap_=" << *winnerMap_ << std::endl;
            std::cout << MyPid << ": weight_.Map=" << *weightMap << std::endl;
#endif
            //  throw e;
            // throw 1;
      }

#ifdef COMPARE_IN_OUT_VECTORS
      if (MyPid == 0) {
        std::cout << "==============" << std::endl;
        std::cout << "call #" << numCalls << " (1-based)" << std::endl;
        std::cout << "==============" << std::endl;
      }
      /*
      bool sameWeight=true;
      bool sameWinner=true;
      bool sameComp=true;
      */
      std::string sameOrDiff;
      {
        ArrayRCP<SC> in_weight = in_weight_->getDataNonConst(0);
        ArrayRCP<SC> weight = weight_.getDataNonConst(0);
        if (MyPid == 0) std::cout << "==============\nweight\n==============\n" << std::endl;
        for (size_t i=0; i < weight_.getLocalLength(); ++i) {
          if (in_weight[i] - weight[i] != 0) sameOrDiff = "  <<<<";
          else                           sameOrDiff = " ";
          std::cout << std::setw(3) << i<<": " << in_weight[i] << "   " << weight[i] << sameOrDiff << in_weight[i] - weight[i] << std::endl;
          /*
          if (in_weight[i] != weight[i]) {
            sameWeight=false;
            std::cout << "\n\nin and out weight DIFFER\n\n" << std::endl;
            std::cout << "i="<<i<<", in=" << in_weight[i] << " , out=" << weight[i] << std::endl;
            break;
          }
          */
        }
      }

      {
        ArrayRCP<LO> in_procWinner = in_procWinner_->getDataNonConst(0);
        ArrayRCP<LO> procWinner = procWinner_.getDataNonConst(0);
        if (MyPid == 0) std::cout << "==============\nprocWinner\n==============\n" << std::endl;
        for (size_t i=0; i < procWinner_.getLocalLength(); ++i) {
          if (in_procWinner[i] != procWinner[i]) sameOrDiff = "  <<<<";
          else                           sameOrDiff = " ";
          std::cout << std::setw(3) << i<<": " << in_procWinner[i] << "   " << procWinner[i] << sameOrDiff << std::endl;
          /*
          if (in_procWinner[i] != procWinner[i]) {
            sameWinner=false;
            std::cout << "\n\nin and out procWinner DIFFER\n\n" << std::endl;
            std::cout << "i="<<i<<", in=" << in_procWinner[i] << ", out=" << procWinner[i] << std::endl;
            break;
          }
          */
        }
      }

      {
        if (companion != NULL) {
          ArrayRCP<LO> in_comp = in_companion->getDataNonConst(0);
          ArrayRCP<LO> comp = companion->getDataNonConst(0);
        if (MyPid == 0) std::cout << "==============\ncompanion\n==============\n" << std::endl;
          for (size_t i=0; i < companion->getLocalLength(); ++i) {
            if (in_comp[i] != comp[i]) sameOrDiff = "  <<<<";
            else                           sameOrDiff = " ";
            std::cout << std::setw(3) << i<<": " << in_comp[i] << "   " << comp[i] << sameOrDiff << std::endl;
            /*
            if (in_comp[i] != comp[i]) {
              sameComp=false;
            std::cout << "\n\nin and out companion DIFFER\n\n" << std::endl;
              std::cout << "i="<<i<<", in=" << in_comp[i] << ", out=" << comp[i] << std::endl;
              break;
            }
            */
          }
        }
      }
#endif
    } //ArbitrateAndCommunicate(Vector&, LOVector &, LOVector *, const bool) const

    // Redistribute data in source to dest where both source and dest might have 
    // multiple copies of the same global id across many processors. The source
    // may not have the same value for all of these multiple copies, but on 
    // termination dest will have a unique value for each global id.  When multiple
    // copies exist in source, 'what' determines how they are combined to make a 
    // unique value in dest (see CombineMode).
    //
    //  Input:
    //     source                   Vector where multiple copies of some GlobalIds
    //                              might exist and might have different values.
    //
    //     dest                     Allocated but contents ignored.
    //
    //     what                     Determines how multiple copies of the same
    //                              GlobalId are combined (see CombineMode).
    //
    //  Output:
    //
    //     dest                     dest has redistributed data from source where
    //                              'what' determines how multiple copies of source
    //                              values associated with the same GlobalId are
    //                              combined into a unique value on all processors.
    //
    void NonUnique2NonUnique(const Vector &source, Vector &dest, const Xpetra::CombineMode what) const
    {
      tempVec_->putScalar(0.);
     
      try
        {
          tempVec_->doExport(source, *import_, what);
          dest.doImport(*tempVec_,   *import_, Xpetra::INSERT);
        }
      catch(std::exception& e)
        {
          int MyPid = tempVec_->getMap()->getComm()->getRank();
          std::cout << MyPid << ": ERR1: An exception occurred." << std::endl;
          throw e;
        }
    }
   
  };


}

#define MUELU_UCAGGREGATIONCOMMHELPER_SHORT
#endif //ifndef MUELU_UCAGGREGATIONCOMMHELPER_HPP

//JG: 
// - procWinner is an array of proc ID -> LocalOrdinal
// - companion == aggregates.GetVertex2AggId() == local aggregate ID -> LocalOrdinal
