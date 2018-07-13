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

/*! \file Zoltan2_PartitioningSolution.hpp
    \brief Defines the PartitioningSolution class.
*/

#ifndef _ZOLTAN2_MATRIXPARTITIONINGSOLUTION_HPP_
#define _ZOLTAN2_MATRIXPARTITIONINGSOLUTION_HPP_

// namespace Zoltan2 {
// template <typename Adapter>
// class PartitioningSolution;
// }

// #include <Zoltan2_Environment.hpp>
// #include <Zoltan2_Solution.hpp>
// #include <Zoltan2_GreedyMWM.hpp>
// #include <Zoltan2_Algorithm.hpp>
// #include <Zoltan2_CoordinatePartitioningGraph.hpp>
// #include <cmath>
// #include <algorithm>
// #include <vector>
// #include <limits>


namespace Zoltan2 {

/*! \brief A PartitioningSolution is a solution to a partitioning problem.

    It is initialized by a PartitioningProblem,
    written to by an algorithm, and may be read by the user or by
    a data migration routine in an input adapter.

    \todo Problem computes metrics using the Solution.  Should
  Solution have a pointer to the metrics, since it may persist after
  the Problem is gone?
    \todo save an RCB tree, so it can be used in repartitioning, and
                supplied to the caller.
    \todo doxyfy the comments in this file.
*/

template <typename Adapter>
  class MatrixPartitioningSolution : public Solution
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
#endif

/*! \brief Constructor when part sizes are not supplied.
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param comm the communicator for the problem associated with
 *             this solution
 *    \param algorithm  Algorithm, if any, used to compute the solution.
 *
 *   It is possible that part sizes were supplied on other processes,
 *   so this constructor does do a check to see if part sizes need
 *   to be globally calculated.
 */

  MatrixPartitioningSolution( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const RCP<Algorithm<Adapter> > &algorithm = Teuchos::null);


  ////////////////////////////////////////////////////////////////////
  // Information that the algorithm may wish to query.

  ////////////////////////////////////////////////////////////////////
  // Method used by the algorithm to set results.

  /*! \brief The algorithm uses setIDLists to set the solution.
   *
   *   \param rowIDs  List of row numbers that the nonzeros for this processor
   *      contain.
   *
   *   \param colIDs  List of column numbers that the nonzeros for this processor
   *      contain.
   *
   *   \param domainIDs Domain vector Ids assigned to this process
   *
   *   \param rangeIDs Range vector Ids assigned to this process
   *
   */

  void setIDLists(ArrayRCP<part_t> &rowIDs,ArrayRCP<part_t> &colIDs, 
		  ArrayRCP<part_t> &domainIDs, ArrayRCP<part_t> &rangeIDs);

  ////////////////////////////////////////////////////////////////////

  // TODO: Figure out if I want to do this
  // 
  // /*! \brief Remap a new partition for maximum overlap with an input partition.
  //  *
  //  * Assumptions for this version:
  //  * input part assignment == processor rank for every local object.
  //  * assuming nGlobalParts <= num ranks
  //  * TODO:  Write a version that takes the input part number as input;
  //  *        this change requires input parts in adapters to be provided in
  //  *        the Adapter.
  //  * TODO:  For repartitioning, compare to old remapping results; see Zoltan1.
  //  */

  // void RemapParts();

  // ////////////////////////////////////////////////////////////////////
  // /* Return the weight of objects staying with a given remap.
  //  * If remap is NULL, compute weight of objects staying with given partition
  //  */
  // long measure_stays(part_t *remap, int *idx, part_t *adj, long *wgt,
  //                    part_t nrhs, part_t nlhs)
  // {
  //   long staying = 0;
  //   for (part_t i = 0; i < nrhs; i++) { 
  //     part_t k = (remap ? remap[i] : i);
  //     for (part_t j = idx[k]; j < idx[k+1]; j++) { 
  //       if (i == (adj[j]-nlhs)) {
  //         staying += wgt[j];
  //         break;
  //       }
  //     }
  //   }
  //   return staying;
  // }

  ////////////////////////////////////////////////////////////////////
  // Results that may be queried by the user, by migration methods,
  // or by metric calculation methods.
  // We return raw pointers so users don't have to learn about our
  // pointer wrappers.

  /*! \brief Return the communicator associated with the solution.
   */
  inline const RCP<const Comm<int> > &getCommunicator() const { return comm_;}

  /*! \brief Return the environment associated with the solution.
   */
  inline const RCP<const Environment> &getEnvironment() const { return env_;}


  /*! \brief Returns list of global row IDs of the nonzeros owned by this process
   */
  const gno_t *getRowIdsView() const
  {  
    if (rowIDs_.size() > 0) 
      return rowIDs_.getRawPtr();
    else                   
      return NULL;
  }


  /*! \brief Returns list of global column IDs of the nonzeros owned by this process
   */
  const gno_t *getColIdsView() const
  {  
    if (colIDs_.size() > 0) 
      return colIDs_.getRawPtr();
    else                   
      return NULL;
  }


  // /*! \brief Returns the process list corresponding to the global ID list.
  //     \return The return value is a NULL pointer if part IDs are
  //               synonomous with process IDs.
  //  */
  // const int *getProcListView() const {
  //   if (procs_.size() > 0) return procs_.getRawPtr();
  //   else                   return NULL;
  // }


  // /*! \brief Get the processes containing a part.
  //  *  \param partId a part number from 0 to one less than the global number
  //  *                of parts.
  //  *  \param procMin on return will be set to minimum proc number
  //  *  \param procMax on return will be set to maximum proc number
  //  *
  //  * Normally \c procMin and \c procMax are the same value and a part
  //  * is assigned to one process.  But if there are more processes than
  //  * parts, it's possible that a part will be divided across more than
  //  * one process.
  //  */
  // void getProcsForPart(part_t partId, part_t &procMin, part_t &procMax) const
  // {
  //   env_->localInputAssertion(__FILE__, __LINE__, "invalid part id",
  //     partId >= 0 && partId < nGlobalParts_, BASIC_ASSERTION);

  //   partToProcsMap(partId, procMin, procMax);
  // }

private:


  // void setPartDistribution();

  // void setPartSizes(ArrayView<ArrayRCP<part_t> > reqPartIds,
  //   ArrayView<ArrayRCP<scalar_t> > reqPartSizes);

  // void computePartSizes(int widx, ArrayView<part_t> ids,
  //   ArrayView<scalar_t> sizes);

  // void broadcastPartSizes(int widx);


  RCP<const Environment> env_;             // has application communicator
  const RCP<const Comm<int> > comm_;       // the problem communicator

  //  part_t nGlobalParts_;// target global number of parts

  ////////////////////////////////////////////////////////////////
  // The algorithm sets these values upon completion.


  // There is a decision to made to decide what to store in the solution.
  // 1. One possibility is to store part numbers for each local id.  This
  //    is what was implemented for PartitioninSolution and this was advantageous
  //    for this case since (unlike 2) an extra communication was not needed
  //    after each process assigned its localids a part.  
  // 2. Alternatively, each process can store the set of ids that it owns.  For 1D
  //    this would require an additional communication step to place the ids to the
  //    correct processor.  This would also complicated the case where number of
  //    parts != number of processes.  This is similar to what is needed to build
  //    row or column epetra/tpetra maps.
  //
  // However, things are more complicated for 2D Cartesian partitioning (maybe all
  // of partitioning).  In this case, we need to assign matrix nonzeros to both
  // a process row and a process column.  This means that extra communication will
  // be needed for option #1 in order to determine both of these, which are needed
  // to determine a part assignment for a given nonzero (or alternatively both the
  // process rows and column assignments for a 2D block of matrix rows and columns,
  // although this may complicated part assignment).
  // This eliminates one of the advantages of method 1 vs. method 2.  

  // For method 2 and 2D, each process would store the set of rowIDs and columnIDs
  // for which it owns the nonzeros. Method 2 still
  // has the disadvantage that it becomes more complicated for #parts != #procs. 
  // However, it would have what is needed to build both the row and column ids.
  // For now, we are going with Method #2.

  // Need a plan to handle #parts != #procs

  ArrayRCP<gno_t> rowIDs_;      // Row Ids assigned to this process
  ArrayRCP<gno_t> colIDs_;      // Col Ids assigned to this process

  ArrayRCP<gno_t> domainIDs_;     // Domain vector Ids assigned to this process
  ArrayRCP<gno_t> rangeIDs_;      // Range vector Ids assigned to this process


  bool haveSolution_;


  ////////////////////////////////////////////////////////////////
  // Algorithm used to compute the solution; 
  // Not sure if this is necessary
  const RCP<Algorithm<Adapter> > algorithm_;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
MatrixPartitioningSolution<Adapter>::MatrixPartitioningSolution(const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm, const RCP<Algorithm<Adapter> > &algorithm)
    : env_(env), comm_(comm), rowIDs_(), colIDs_(), haveSolution_(false), 
      algorithm_(algorithm)
{
  env_->memory("After construction of solution");
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void MatrixPartitioningSolution<Adapter>::setIDLists(ArrayRCP<part_t> &rowIDs,ArrayRCP<part_t> &colIDs, 
					             ArrayRCP<part_t> &domainIDs, ArrayRCP<part_t> &rangeIDs)
{
  env_->debug(DETAILED_STATUS, "Entering setParts");

  rowIDs_=rowIDs;
  colIDs_=colIDs;
  domainIDs_=domainIDs;
  rangeIDs_=rangeIDs;

  haveSolution_ = true;

  env_->memory("After Solution has processed algorithm's answer");
  env_->debug(DETAILED_STATUS, "Exiting setParts");
}
////////////////////////////////////////////////////////////////////////////////




}  // namespace Zoltan2

#endif
