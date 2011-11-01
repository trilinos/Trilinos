#ifndef MUELU_UCAGGREGATIONCOMMHELPER_DECL_HPP
#define MUELU_UCAGGREGATIONCOMMHELPER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MapFactory.hpp>

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
    ;

    ~UCAggregationCommHelper() ;

    void ArbitrateAndCommunicate(Vector &weights, Aggregates &aggregates, const bool perturb) const
    ;

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
    ; //ArbitrateAndCommunicate(Vector&, LOVector &, LOVector *, const bool) const

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
    ;
   
  };


}

//JG: 
// - procWinner is an array of proc ID -> LocalOrdinal
// - companion == aggregates.GetVertex2AggId() == local aggregate ID -> LocalOrdinal

#define MUELU_UCAGGREGATIONCOMMHELPER_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_UCAGGREGATIONCOMMHELPER_DECL_HPP
