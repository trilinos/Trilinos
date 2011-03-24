#ifndef MUELU_AGGALGORITHM2COMM_HPP
#define MUELU_AGGALGORITHM2COMM_HPP

#include <Cthulhu_VectorFactory.hpp>
#include <Cthulhu_ImportFactory.hpp>
#include <Cthulhu_MapFactory.hpp>

class AggAlgorithm2Comm {

private:
  RCP<const Cthulhu::Import<int,int> > import_;

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

  AggAlgorithm2Comm(const RCP<const Cthulhu::Map<int,int> > & uniqueMap, const RCP<const Cthulhu::Map<int,int> > & nonUniqueMap)
  {
    import_ = Cthulhu::ImportFactory<int,int>::Build(uniqueMap, nonUniqueMap);
  }

  ~AggAlgorithm2Comm() {}

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
  int ArbitrateAndCommunicate(Cthulhu::Vector<double> &weight_, Cthulhu::Vector<int> &procWinner_, Cthulhu::Vector<int> *companion, const bool perturb) const
  {
    int MyPid = weight_.getMap()->getComm()->getRank(); // TODO:remove the getMap() step

    if (perturb) {
      RCP<Cthulhu::Vector<double> > perturbWt_ = Cthulhu::VectorFactory<double>::Build(weight_.getMap());

      // Note: maxValue() not available for Tpetra
      //double largestGlobalWeight = weight_.maxValue();
      double largestGlobalWeight = weight_.normInf();

      //TODO
      //    Epetra_Util util;
      //    util.SetSeed( (unsigned int) MyPid*2 + (int) (11*rand()));
      //    for (int i = 0; i < 10; i++) util.SetSeed( (unsigned int) util.RandomInt() );

      // perturbWt.SetSeed( (unsigned int) util.RandomInt() );
      perturbWt_->randomize(); 

      Teuchos::ArrayRCP<double> weight = weight_.getDataNonConst(0);
      Teuchos::ArrayRCP<double> perturbWt = perturbWt_->getDataNonConst(0);

      for (size_t i=0; i < weight_.getMap()->getNodeNumElements(); i++) {
        if (weight[i] == 0.) perturbWt[i] = 0.;
        else perturbWt[i] = weight[i] + 1.0e-7*largestGlobalWeight*fabs(perturbWt[i]);
      }
      for (size_t i=0; i < weight_.getMap()->getNodeNumElements(); i++) weight[i] = perturbWt[i]; 
    }
  
    // Communicate weights and store results in PostComm (which will be copied
    // back into weights later. When multiple processors have different weights
    // for the same GID, we take the largest weight. After this fragment every
    // processor should have the same value for PostComm[] even when multiple
    // copies of the same Gid are involved.

    RCP<Cthulhu::Vector<double> > postComm_ = Cthulhu::VectorFactory<double>::Build(weight_.getMap());

    Teuchos::ArrayRCP<double> postComm = postComm_->getDataNonConst(0);
    postComm_->putScalar(0.0);

    NonUnique2NonUnique(weight_, *postComm_, Cthulhu::ABSMAX);

   
    // Let every processor know who is the procWinner. For nonunique
    // copies of the same Gid, this corresponds to the procesosr with
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

    RCP<Cthulhu::Vector<double> > candidateWinners_ = Cthulhu::VectorFactory<double>::Build(weight_.getMap());

    Teuchos::ArrayRCP<double> candidateWinners = candidateWinners_->getDataNonConst(0);
    candidateWinners_->putScalar(0.0);

    Teuchos::ArrayRCP<double> weight = weight_.getDataNonConst(0);
    for (size_t i=0; i < weight_.getMap()->getNodeNumElements(); i++) {
      if (postComm[i] == weight[i]) candidateWinners[i] = (double) MyPid+1;
    }

    for (size_t i=0; i < weight_.getMap()->getNodeNumElements(); i++) weight[i]=postComm[i]; 

    NonUnique2NonUnique(*candidateWinners_, *postComm_, Cthulhu::ABSMAX);

    // Note: 
    //                      associated CandidateWinners[]
    //    weight[i]!=0  ==> on some proc is equal to its ==> postComm[i]!=0
    //                      MyPid+1.
    //          
    Teuchos::ArrayRCP<int> procWinner = procWinner_.getDataNonConst(0);
    for (size_t i=0; i < weight_.getMap()->getNodeNumElements(); i++)  {
      if ( weight[i] != 0.) procWinner[i] = ((int) (postComm[i])) - 1;
    }

    if (companion != NULL) {
      // Now build a new Map, WinnerMap which just consists of procWinners. 
      // This is done by extracting the Gids for Wt, and shoving
      // the subset that correspond to procWinners in MyWinners.
      // WinnerMap is then constructed using MyWinners.
   
      int numMyWinners = 0;
      for (size_t i = 0; i < weight_.getMap()->getNodeNumElements(); i++) {
        if (procWinner[i] == MyPid) numMyWinners++;
      }
   
      //    int *myGids    = new int[weight_.getMap()->getNodeNumElements()+1];
      //    int *myWinners = new int[numMyWinners+1];
      Teuchos::ArrayView<const int> myGids = weight_.getMap()->getNodeElementList(); //== weight_.getMap()->MyGlobalElements(myGids);
      Teuchos::ArrayRCP<int> myWinners(numMyWinners);

      numMyWinners = 0;
      for (size_t i = 0; i < weight_.getMap()->getNodeNumElements(); i++) {
        if (procWinner[i] == MyPid)
          myWinners[numMyWinners++] = myGids[i];
      }
    
      // Cthulhu::EpetraMap winnerMap(-1, numMyWinners, myWinners, 0, weight_.getMap()->getComm());    
      Cthulhu::global_size_t g = -1; //TODO for Tpetra -1 == ??
      RCP<Cthulhu::Map<int> > winnerMap = Cthulhu::MapFactory<int>::Build(weight_.getMap()->lib(), g, myWinners(), 0, weight_.getMap()->getComm());
       
      // Pull the Winners out of companion
      //     JustWinners <-- companion[Winners];
   
      RCP<Cthulhu::Vector<int> > justWinners = Cthulhu::VectorFactory<int>::Build(winnerMap);
      RCP<const Cthulhu::Import<int,int> > winnerImport = Cthulhu::ImportFactory<int>::Build(winnerMap,weight_.getMap());

      justWinners->doImport(*companion, *winnerImport, Cthulhu::INSERT);
   
      // Put the JustWinner values back into companion so that
      // all nonunique copies of the same Gid have the procWinner's
      // version of the companion.
   
      RCP<Cthulhu::Import<int> > pushWinners = Cthulhu::ImportFactory<int>::Build(weight_.getMap(), winnerMap);
      companion->doImport(*justWinners, *pushWinners, Cthulhu::INSERT);

    }

    return 0; //TODO
  }

  inline int ArbitrateAndCommunicate(Cthulhu::Vector<double> &weights, MueLu::Aggregates<int,int> &aggregates, const bool perturb) const
  {
    return ArbitrateAndCommunicate(weights, *aggregates.GetProcWinner(), &*aggregates.GetVertex2AggId(), perturb);
  }

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
  int NonUnique2NonUnique(const Cthulhu::Vector<double> &source, Cthulhu::Vector<double> &dest, const Cthulhu::CombineMode what) const
  {
    RCP<Cthulhu::Vector<double> > temp = Cthulhu::VectorFactory<double>::Build(import_->getSourceMap());
     
//      cout << *source.getMap() << endl;
//      cout << *import_->getTargetMap() << endl;

    TEST_FOR_EXCEPTION(!source.getMap()->isSameAs(*import_->getTargetMap()), std::runtime_error, "Source Maps don't match.");

//     cout << "EXPORT" << std::endl;

    temp->doExport(source, *import_, what);

//     cout << "IMPORT" << std::endl;

//     std::cout << *dest.getMap() << std::endl;
//     std::cout << *import_->getTargetMap() << std::endl;

    dest.doImport(*temp,   *import_, Cthulhu::INSERT);
     


    return 0;
  }
   
};
#endif //ifndef MUELU_AGGALGORITHM2COMM_HPP
