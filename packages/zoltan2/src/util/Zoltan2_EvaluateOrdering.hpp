// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_EvaluateOrdering.hpp
 *  \brief Defines the Zoltan2_EvaluateOrdering.hpp class.
 */

#ifndef ZOLTAN2_EVALUATEORDERING_HPP
#define ZOLTAN2_EVALUATEORDERING_HPP

#include <Zoltan2_EvaluateBaseClass.hpp>
#include <Zoltan2_OrderingSolution.hpp>

namespace Zoltan2{

/*! \brief A class that computes and returns quality metrics.
 *  \A base class for the local and global ordering versions.
 */
template <typename Adapter>
class EvaluateOrdering : public EvaluateBaseClass<Adapter> {

private:
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

  // To do - this is only appropriate for the local ordering
  // Need to decide how to organize these classes
  // Do we potentially want local + global in this class
  // Do we want to eliminate the EvaluateLocalOrdering and EvaluateGlobalOrdering
  // derived classes? Or perhaps make them completely independent of each other
  lno_t bandwidth;
  lno_t envelope;
  lno_t separatorSize;

  void sharedConstructor(const Adapter *ia,
                         ParameterList *p,
                         const RCP<const Comm<int> > &problemComm,
                         const LocalOrderingSolution<lno_t> *localSoln,
                         const GlobalOrderingSolution<gno_t> *globalSoln);
public:

  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateOrdering(
    const Adapter *ia,
    ParameterList *p,
    const LocalOrderingSolution<lno_t> *localSoln,
    const GlobalOrderingSolution<gno_t> *globalSoln)
  {
    Teuchos::RCP<const Comm<int> > problemComm = Tpetra::getDefaultComm();
    sharedConstructor(ia, p, problemComm, localSoln, globalSoln);
  }

  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateOrdering(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const LocalOrderingSolution<lno_t> *localSoln,
    const GlobalOrderingSolution<gno_t> *globalSoln)
  {
    sharedConstructor(ia, p, problemComm, localSoln, globalSoln);
  }

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateOrdering(
    const Adapter *ia,
    ParameterList *p,
    MPI_Comm comm,
    const LocalOrderingSolution<lno_t> *localSoln,
    const GlobalOrderingSolution<gno_t> *globalSoln)
  {
    RCP<Teuchos::OpaqueWrapper<MPI_Comm> > wrapper =
        Teuchos::opaqueWrapper(comm);
    RCP<const Comm<int> > problemComm =
        rcp<const Comm<int> >(new Teuchos::MpiComm<int>(wrapper));
    sharedConstructor(ia, p, problemComm, localSoln, globalSoln);
  }
#endif

  lno_t getBandwidth() const { return bandwidth; }
  lno_t getEnvelope() const { return envelope; }
  lno_t getSeparatorSize() const { return separatorSize; }

  /*! \brief Print all metrics of type metricType based on the metric object type
   *  Note that parent class currently suppresses this if the list is empty.
   */
  virtual void printMetrics(std::ostream &os) const {

    // To Do - complete this formatting
    os << "Ordering Metrics" << std::endl;
    os << std::setw(20) << " " << std::setw(11) << "ordered" << std::endl;
    os << std::setw(20) << "envelope" << std::setw(11) << std::setprecision(4)
       << envelope << std::endl;
    os << std::setw(20) << "bandwidth" << std::setw(11) << std::setprecision(4)
       << bandwidth << std::endl;
  }

  void localOrderingMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const Adapter *ia,
    const LocalOrderingSolution<typename Adapter::lno_t> *localSoln)
  {
    env->debug(DETAILED_STATUS, "Entering orderingMetrics"); // begin

    typedef StridedData<lno_t, scalar_t> input_t;

    // get graph
    std::bitset<NUM_MODEL_FLAGS> modelFlags;
    RCP<GraphModel<base_adapter_t> > graph;
    const RCP<const base_adapter_t> bia =
      rcp(dynamic_cast<const base_adapter_t *>(ia), false);
    graph = rcp(new GraphModel<base_adapter_t>(bia,env,comm,modelFlags));
    ArrayView<const gno_t> Ids;
    ArrayView<input_t> vwgts;
    ArrayView<const gno_t> edgeIds;
    ArrayView<const offset_t> offsets;
    ArrayView<input_t> wgts;
    ArrayView<input_t> vtx;
    graph->getEdgeList(edgeIds, offsets, wgts);
    lno_t numVertex = graph->getVertexList(Ids, vwgts);

    lno_t * perm = localSoln->getPermutationView();

    // print as matrix - this was debugging code which can be deleted later
    #define MDM
    #ifdef MDM
    for( int checkRank = 0; checkRank < comm->getSize(); ++checkRank ) {
      comm->barrier();
      if( checkRank == comm->getRank() ) {
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Inspect rank: " << checkRank << std::endl;
        std::cout << std::endl;
        if(numVertex < 30) { // don't spam if it's too many...
          Array<lno_t> oldMatrix(numVertex*numVertex);
          Array<lno_t> newMatrix(numVertex*numVertex);

          // print the solution permutation
          std::cout << std::endl << "perm:  ";
          for(lno_t n = 0; n < numVertex; ++n) {
            std::cout << " " << perm[n] << " ";
          }

          lno_t * iperm = localSoln->getPermutationView(true);
          std::cout << std::endl << "iperm: ";
          for(lno_t n = 0; n < numVertex; ++n) {
            std::cout << " " << iperm[n] << " ";
          }
          std::cout << std::endl;
          // write 1's to old matrix (original form) and new matrix (using solution)
          for (lno_t y = 0; y < numVertex; y++) {
            for (offset_t n = offsets[y]; n < offsets[y+1]; ++n) {
              lno_t x = static_cast<lno_t>(edgeIds[n]); // to resolve
              if (x < numVertex && y < numVertex) { // to develop - for MPI this may not be local
                oldMatrix[x + y*numVertex] = 1;
                newMatrix[perm[x] + perm[y]*numVertex] = 1;
              }
            }
          }

          // print oldMatrix
          std::cout << std::endl << "original graph in matrix form:" << std::endl;
          for(lno_t y = 0; y < numVertex; ++y) {
            for(lno_t x = 0; x < numVertex; ++x) {
              std::cout << " " << oldMatrix[x + y*numVertex];
            }
            std::cout << std::endl;
          }

          // print newMatrix
          std::cout << std::endl << "reordered graph in matrix form:" << std::endl;
          for(lno_t y = 0; y < numVertex; ++y) {
            for(lno_t x = 0; x < numVertex; ++x) {
              std::cout << " " << newMatrix[x + y*numVertex];
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
        }
      }

      comm->barrier();
    }
    #endif // Ends temporary logging which can be deleted later

    // calculate bandwidth and envelope for unsolved and solved case
    lno_t bw_right = 0;
    lno_t bw_left = 0;
    envelope = 0;

    for (lno_t j = 0; j < numVertex; j++) {
      lno_t y = Ids[j];
      for (offset_t n = offsets[j]; n < offsets[j+1]; ++n) {
        lno_t x = static_cast<lno_t>(edgeIds[n]); // to resolve
        if(x < numVertex) {
          lno_t x2 = perm[x];
          lno_t y2 = perm[y];

          // solved bandwidth calculation
          lno_t delta_right = y2 - x2;
          if (delta_right > bw_right) {
            bw_right = delta_right;
          }
          lno_t delta_left = y2 - x2;
          if (delta_left > bw_left) {
            bw_left = delta_left;
          }

          // solved envelope calculation
          if(delta_right > 0) {
            envelope += delta_right;
          }
          if(delta_left > 0) {
            envelope += delta_left;
          }
          envelope += 1; // need to check this - do we count center?
        }
      }
    }

    bandwidth = (bw_left + bw_right + 1);

    // TO DO - No implementation yet for this metric
    separatorSize = 0;

    env->debug(DETAILED_STATUS, "Exiting orderingMetrics"); // end
  }
};

// sharedConstructor
template <typename Adapter>
void EvaluateOrdering<Adapter>::sharedConstructor(
  const Adapter *ia,
  ParameterList *p,
  const RCP<const Comm<int> > &comm,
  const LocalOrderingSolution<lno_t> *localSoln,
  const GlobalOrderingSolution<gno_t> *globalSoln)
{
  RCP<const Comm<int> > problemComm = (comm == Teuchos::null) ?
    Tpetra::getDefaultComm() : comm;

  RCP<Environment> env;
  try{
    env = rcp(new Environment(*p, problemComm));
  }
  Z2_FORWARD_EXCEPTIONS

  env->debug(DETAILED_STATUS, std::string("Entering EvaluateOrdering"));
  env->timerStart(MACRO_TIMERS, "Computing ordering metrics");

  try{
    // May want to move these into the specific derived classes
    // But it depends on whether we eventually may have both types and perhaps
    // want to combine the metrics
    if(localSoln) {
      localOrderingMetrics(env, problemComm, ia, localSoln);
    }

    if(globalSoln) {
      throw std::logic_error("EvaluateOrdering not set up for global ordering.");
    }
  }
  Z2_FORWARD_EXCEPTIONS;
  env->timerStop(MACRO_TIMERS, "Computing ordering metrics");

  env->debug(DETAILED_STATUS, std::string("Exiting EvaluateOrdering"));
}

template <typename Adapter>
class EvaluateLocalOrdering : public EvaluateOrdering<Adapter> {
private:
  typedef typename Adapter::lno_t lno_t;

public:
  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateLocalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const LocalOrderingSolution<lno_t> *localSoln) :
    EvaluateOrdering<Adapter>(ia, p, localSoln, nullptr) {}

  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateLocalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const LocalOrderingSolution<lno_t> *localSoln) :
    EvaluateOrdering<Adapter>(ia, p, problemComm, localSoln, nullptr) {}

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateLocalOrdering(
    const Adapter *ia,
    ParameterList *p,
    MPI_Comm comm,
    const LocalOrderingSolution<lno_t> *localSoln) :
    EvaluateOrdering<Adapter>(ia, p, comm, localSoln, nullptr) {}
#endif
};

template <typename Adapter>
class EvaluateGlobalOrdering : public EvaluateOrdering<Adapter> {
private:
  typedef typename Adapter::gno_t gno_t;

public:
  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateGlobalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const GlobalOrderingSolution<gno_t> *globalSoln) :
    EvaluateOrdering<Adapter>(ia, p, nullptr, globalSoln) {}

  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateGlobalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const GlobalOrderingSolution<gno_t> *globalSoln) :
    EvaluateOrdering<Adapter>(ia, p, problemComm, nullptr, globalSoln) {}

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateGlobalOrdering(
    const Adapter *ia,
    ParameterList *p,
    MPI_Comm comm,
    const GlobalOrderingSolution<gno_t> *globalSoln) :
    EvaluateOrdering<Adapter>(ia, p, comm, nullptr, globalSoln) {}
#endif
};

}   // namespace Zoltan2

#endif
