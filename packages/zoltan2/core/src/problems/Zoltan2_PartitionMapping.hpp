// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Zoltan2_Algorithm.hpp"

namespace Zoltan2 {

/*! \brief PartitionMapping maps a solution or an input distribution to ranks.
*/

template <typename Adapter>
  class PartitionMapping :public Algorithm<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::scalar_t pcoord_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
#endif

  const Teuchos::RCP <const Teuchos::Comm<int> >comm;
  const Teuchos::RCP <const Zoltan2::MachineRepresentation <pcoord_t,part_t> > machine;
  const Teuchos::RCP <const Adapter > input_adapter;
  const Teuchos::RCP <const Zoltan2::PartitioningSolution<Adapter> >soln;
  const Teuchos::RCP <const Environment >env;
  const part_t num_parts;
  const part_t *solution_parts;


/*! \brief Constructor 
 *  Constructor builds the map from parts to ranks.
 *  KDDKDD WILL NEED THE SOLUTION FOR INTELLIGENT MAPPING
 *  KDDKDD BUT MAY WANT TO SET PART SIZES BASED ON CAPABILITY OF A RANK.
 *  KDDKDD SO WHEN SHOULD THE MAP BE CREATED?
 */
  PartitionMapping(
      const Teuchos::RCP <const Teuchos::Comm<int> >comm_,
      const Teuchos::RCP <const Zoltan2::MachineRepresentation<pcoord_t,part_t> >machine_, // If NULL, assume homogeneous
                                                    // Make optional
      const Teuchos::RCP <const Adapter> input_adapter_, // Needed to get information about
                                           // the application data (coords, graph)
      const Teuchos::RCP <const Zoltan2::PartitioningSolution<Adapter> >soln_, // Needed for mapping a partition
      const Teuchos::RCP <const Environment > envConst_  // Perhaps envConst should be optional
                                           // so applications can create a mapping
                                           // directly
    ):comm(comm_),
      machine(machine_),
      input_adapter(input_adapter_),
      soln(soln_),
      env(envConst_),num_parts(soln_->getActualGlobalNumberOfParts()),
      solution_parts(soln_->getPartListView())
          {} ;

  PartitionMapping(
      const Teuchos::RCP <const Teuchos::Comm<int> >comm_,
      const Teuchos::RCP <const Zoltan2::MachineRepresentation<pcoord_t,part_t> >machine_, // If NULL, assume homogeneous
                                                    // Make optional
      const Teuchos::RCP <const Adapter> input_adapter_, // Needed to get information about
                                           // the application data (coords, graph)
      const part_t num_parts_,
      const part_t *result_parts,
      const Teuchos::RCP <const Environment > envConst_  // Perhaps envConst should be optional
                                           // so applications can create a mapping
                                           // directly
    ):comm(comm_),
      machine(machine_),
      input_adapter(input_adapter_),
      soln(),
      env(envConst_),num_parts(num_parts_),
      solution_parts(result_parts)
          {} ;

  PartitionMapping(
      const Teuchos::RCP <const Teuchos::Comm<int> >comm_,
      const Teuchos::RCP <const Environment > envConst_  // Perhaps envConst should be optional
                                           // so applications can create a mapping
                                           // directly
    ):comm(comm_),
      machine(),
      input_adapter(),
      soln(),
      env(envConst_),num_parts(0),
      solution_parts(NULL)
          {} ;

  PartitionMapping():
          comm(0),
          machine(0),
          input_adapter(0),
          soln(0),
          env(0),
          solution_parts(NULL){};

  PartitionMapping(const Teuchos::RCP <const Environment >envConst_):
          comm(0),
          machine(0),
          input_adapter(0),
          soln(0),
          env(envConst_),num_parts(0),
          solution_parts(NULL){};

  PartitionMapping(
      const Teuchos::RCP <const Environment > envConst_,
      const Teuchos::RCP <const Teuchos::Comm<int> >comm_,
      const Teuchos::RCP <const MachineRepresentation<pcoord_t,part_t> >machine_
      ):
          comm(comm_),
          machine(machine_),
          input_adapter(0),
          soln(0),
          env(envConst_),num_parts(0),
          solution_parts(NULL){};


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
  virtual void getPartsForProc(int procId, part_t &numParts, part_t *&parts)
    const = 0;

  /*! \brief Get the processes containing a part.
   *  \param partId a part number from 0 to one less than the global number
   *                of parts.
   *  \param numProcs on return will be the number of procs owning partId
   *  \param procs on return will be prointer to the procs owning partId
   */
  // TODO:  KDDKDD Arguments should be count and array, not min and max.
  // TODO:  KDDKDD Could require O(nGlobalParts) storage
  virtual void getProcsForPart(part_t partId, part_t &numProcs, part_t *&procs) const = 0;

private:
};

}  // namespace Zoltan2

#endif
