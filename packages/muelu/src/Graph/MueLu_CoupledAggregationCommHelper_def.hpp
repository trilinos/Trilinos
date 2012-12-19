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
#ifndef MUELU_COUPLEDAGGREGATIONCOMMHELPER_DEF_HPP
#define MUELU_COUPLEDAGGREGATIONCOMMHELPER_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_ExportFactory.hpp>

#include "MueLu_CoupledAggregationCommHelper_decl.hpp"

#include "MueLu_Utilities.hpp" // maxAll

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoupledAggregationCommHelper<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoupledAggregationCommHelper(const RCP<const Map> & uniqueMap, const RCP<const Map> & nonUniqueMap) {
    import_ = ImportFactory::Build(uniqueMap, nonUniqueMap);
    tempVec_ = VectorFactory::Build(uniqueMap,false); //zeroed out before use
    numMyWinners_ = 0;
    numCalls_ = 0;
    myPID_ = uniqueMap->getComm()->getRank();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoupledAggregationCommHelper<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ArbitrateAndCommunicate(Vector &weight_, LOVector &procWinner_, LOVector *companion, const bool perturb) const {
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
      maxAll(comm,irealloc,orealloc);
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

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoupledAggregationCommHelper<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NonUnique2NonUnique(const Vector &source, Vector &dest, const Xpetra::CombineMode what) const {
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

}

#endif // MUELU_COUPLEDAGGREGATIONCOMMHELPER_DEF_HPP

