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
#ifndef MUELU_UCAGGREGATIONCOMMHELPER_DECL_HPP
#define MUELU_UCAGGREGATIONCOMMHELPER_DECL_HPP

#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_UCAggregationCommHelper_fwd.hpp"

#include "MueLu_Aggregates.hpp"

namespace MueLu {

  /*!
    @class UCAggregationCommHelper
    @brief Helper class for providing arbitrated communication across processors

    For more details, see the comments for the ArbitrateAndCommunicate methods.
  */

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class UCAggregationCommHelper : public BaseClass {

    typedef double Scalar; // Scalar type only used for weight: always a double.
#undef MUELU_UCAGGREGATIONCOMMHELPER_SHORT
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

    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    UCAggregationCommHelper(const RCP<const Map> & uniqueMap, const RCP<const Map> & nonUniqueMap);

    //! Destructor.
    ~UCAggregationCommHelper() { }

    //@}

    /*!
     @brief This method assigns unknowns to aggregates.

            Tie-breaking is possible is using random weights.

            @param[in] weights vector of weights that help determine ownership.

            @param[in,out] aggregates aggregate data structure

            @param[in] perturb flag indicating whether weights should be randomly perturbed for tie-breaking purposes.
    */
    void ArbitrateAndCommunicate(Vector &weights, Aggregates &aggregates, const bool perturb) const {
      ArbitrateAndCommunicate(weights, *aggregates.GetProcWinner(), &*aggregates.GetVertex2AggId(), perturb);
    }

    /*!
    @brief This class uses a weighted rendezvous algorithm to do a global reduction on a vector that may be based on a non unique map.

    A non-unique map is one that has at least one global ID that occurs on two or more processes.  For each repeated ID \f$i\f$, the
    algorithm finds the maximum value \f$v[i]\f$ in the weight vector \f$v\f$.  This value is communicated to all processors that
    have \f$i\f$ in their local map.  More details are below.


     For each GlobalId \f$K\f$ associated with weight.getMap():

          -# Find the maximum absolute value of \f$weight[K]\f$ across all
             processors and assign this to all local elements of weight[] (across
             processors) that are associated with \f$K\f$.
          -# Set procWinner[] to the MyPid() that had the largest element.
             procWinner[] is still set if only one processor owns a GlobalId.

             The ONLY CASE when procWinner[i] is NOT set corresponds to when
             all local weights associated with a GlobalId are zero. This allows
             one to effectively skip the maximum/winner calculation for a subset
             of GlobalId's.  This might occur when a processor has already
             claimed ownership for a GlobalId and so all local copies have
             the same value. We want to skip the maximum calculation with
             tiebreaking to avoid another processor claiming ownership.

          -# Optionally, set companion[] (across all relevant processors) to the
             local companion value associated with the procWinner[] processor.

         @param weight[in,out]
                                 - On input, vector of NONNEGATIVE weights.
                                 - On output, \f$ \mbox{weight}[k]  \Leftarrow  \max(\mbox{weight}[k_{p1}],\dots,\mbox{weight}[k_{pn}]) \f$
                                  where \f$ \mbox{weight}[k_{pj}] \f$ is processor \f$pj\f$'s value for GID \f$k\f$.

         @param procWinner[in,out]
                                  - On input, allocated but contents ignored.
                                  - On output, \f$\mbox{procWinner}[k] \Leftarrow pj\f$  such that
                                    \f$\mbox{weight}[k_{pj}] = \max(\mbox{weight}[k_{p1}],...,\mbox{weight}[k_{pn}])\f$, where
                                    \f$ \mbox{weight}[k_{pj}] \f$ is processor \f$pj\f$'s value for GID \f$k\f$.
                                  NOTE: If all input \f$\mbox{weight}[k_{pi}]\f$'s are zero, then \f$\mbox{procWinner}[k]\f$ is left untouched.

         @param companion[in,out]
                                  - On input, either NULL or allocated but contents ignored.  If NULL, step 3 above is skipped.
                                  - On output, if not null, \f$\mbox{companion}[k] \Leftarrow \mbox{companion}[k_j]\f$ where
                                  \f$\mbox{companion}[k_j]\f$ lives on processor \f$\mbox{procWinner}[k]\f$.
                                  and corresponds to the same GlobalId as \f$k\f$.
                                  NOTE: If for a particular GlobalId, no processor
                                        has a value of procWinner that matches
                                        its MyPid, the corresponding companion
                                        is not altered.


         @param perturb[in]                  Optional arguments that is either true or
                                             false (default: true). weight is perturbed
                                             and the perturbed values are used in step 1)
                                             above. Returned values reflect the perturbed
                                             data. This option avoids having lots of
                                             tiebreaks where the large MyPid() always wins.

      */
      /*
      Output:
         @param weight            \f$ weight[k]  \Leftarrow  \max(weight[k_1],\dots,weight[k_n]) \f$
                                  where \f$ weight[k_j] \f$ live on different processors
                                  but have the same GlobalId as weight[k] on this processor.

         @param procWinner               procWinner[k] <-- MyPid associated with the
                                  kj yielding the max in
                                        Max(weight[k1],...,weight[kn]) .
                                  See weight Output comments.
                                  NOTE: If all input weight[kj]'s are zero,
                                        then procWinner[k] is left untouched.

         @param companion                If not null,
                                     companion[k] <-- companion[kj] where
                                  companion[kj] lives on processor procWinner[k].
                                  and corresponds to the same GlobalId as k.
                                  NOTE: If for a particlar GlobalId, no processor
                                        has a value of procWinner that matches
                                        its MyPid, the corresponding companion
                                        is not altered.
    */
    void ArbitrateAndCommunicate(Vector &weight, LOVector &procWinner, LOVector *companion, const bool perturb) const; //ArbitrateAndCommunicate(Vector&, LOVector &, LOVector *, const bool) const

    /*!  @brief Redistribute data in source to dest where both source and dest might have multiple copies of the same global id across many processors.

       The source may not have the same value for all of these multiple copies, but on
       termination dest will have a unique value for each global id.  When multiple
       copies exist in source, 'what' determines how they are combined to make a
       unique value in dest (see CombineMode).

        Input:
           @param[in] source        Vector where multiple copies of some GlobalIds
                                    might exist and might have different values.

           @param[in,out] dest      On input, allocated but contents ignored.
                                    On output, contains redistributed data from source where
                                    'what' determines how multiple copies of source
                                    values associated with the same GlobalId are
                                    combined into a unique value on all processors.

           @param[in] what          Determines how multiple copies of the same
                                    GlobalId are combined (see CombineMode).
    */
    void NonUnique2NonUnique(const Vector &source, Vector &dest, const Xpetra::CombineMode what) const;

  };


}

//JG:
// - procWinner is an array of proc ID -> LocalOrdinal
// - companion == aggregates.GetVertex2AggId() == local aggregate ID -> LocalOrdinal

#define MUELU_UCAGGREGATIONCOMMHELPER_SHORT
#endif // MUELU_UCAGGREGATIONCOMMHELPER_DECL_HPP
