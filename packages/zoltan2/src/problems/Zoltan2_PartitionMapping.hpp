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

/*! \file Zoltan2_PartMapping.hpp
    \brief Defines the PartMapping class.
*/

#ifndef _ZOLTAN2_PARTITIONMAPPING_HPP_
#define _ZOLTAN2_PARTITIONMAPPING_HPP_
#include "Zoltan2_Model.hpp"
#include "Zoltan2_PartitioningSolution.hpp"
#include "Teuchos_Comm.hpp"
#include "Zoltan2_Environment.hpp"
#include "Zoltan2_MachineRepresentation.hpp"

namespace Zoltan2 {

/*! \brief PartitionMapping maps a solution or an input distribution to ranks.
*/

template <typename Adapter>
  class PartitionMapping
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::zgid_t zgid_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
#endif

  const Teuchos::Comm<int> *comm;
  const Zoltan2::MachineRepresentation <scalar_t> *machine;
  const Zoltan2::Model<typename Adapter::base_adapter_t> *model;
  const Zoltan2::PartitioningSolution<Adapter> *soln;
  const Environment *env;

/*! \brief Constructor 
 *  Constructor builds the map from parts to ranks.
 *  KDDKDD WILL NEED THE SOLUTION FOR INTELLIGENT MAPPING
 *  KDDKDD BUT MAY WANT TO SET PART SIZES BASED ON CAPABILITY OF A RANK.
 *  KDDKDD SO WHEN SHOULD THE MAP BE CREATED?
 */
  PartitionMapping(
    const Teuchos::Comm<int> *comm_,
    const Zoltan2::MachineRepresentation<scalar_t> *machine_, // If NULL, assume homogeneous
                                                  // Make optional
    const Zoltan2::Model<typename Adapter::base_adapter_t> *model_, // Needed to get information about
                                         // the application data (coords, graph)
    const Zoltan2::PartitioningSolution<Adapter> *soln_, // Needed for mapping a partition
    const Environment *envConst_  // Perhaps envConst should be optional
                                         // so applications can create a mapping
                                         // directly
  ):comm(comm_),
      machine(machine_),
      model(model_),
      soln(soln_),
      env(envConst_)
          {} ;

  PartitionMapping():
          comm(0),
          machine(0),
          model(0),
          soln(0),
          env(0){};

  PartitionMapping(const Environment *envConst_):
          comm(0),
          machine(0),
          model(0),
          soln(0),
          env(envConst_){};

  PartitionMapping(
          const Environment *envConst_,
          const Teuchos::Comm<int> *comm_,
          const MachineRepresentation<scalar_t> *machine_
          ):
          comm(comm_),
          machine(machine_),
          model(0),
          soln(0),
          env(envConst_){};


  virtual ~PartitionMapping(){}

/*! \brief Returns the number of parts to be assigned to this process.
 */
  virtual size_t getLocalNumberOfParts() const = 0;

  /*! \brief Get the parts belonging to a process.
   *  \param procId a process rank
   *  \param numParts on return will be set the number of parts belonging
   *                    to the process.
   *  \param parts on return will be a pointer to the parts assigned to procId
   *
   */
  // TODO:  KDDKDD Decide whether information should be avail for any process
  // TODO:  KDDKDD (requiring more storage or a directory) or only for the 
  // TODO:  KDDKDD local process.
  // TODO:  KDDKDD Could require O(nprocs) storage
  virtual void getPartsForProc(int procId, part_t &numParts, part_t *parts)
    const = 0;

  /*! \brief Get the processes containing a part.
   *  \param partId a part number from 0 to one less than the global number
   *                of parts.
   *  \param numProcs on return will be the number of procs owning partId
   *  \param procs on return will be prointer to the procs owning partId
   */
  // TODO:  KDDKDD Arguments should be count and array, not min and max.
  // TODO:  KDDKDD Could require O(nGlobalParts) storage
  virtual void getProcsForPart(part_t partId, part_t &numProcs, part_t *procs) const = 0;

private:
};

}  // namespace Zoltan2

#endif
