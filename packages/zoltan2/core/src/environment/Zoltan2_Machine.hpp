// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MACHINE_HPP_
#define _ZOLTAN2_MACHINE_HPP_

#include <Teuchos_Comm.hpp>
#include <Zoltan2_config.h>
namespace Zoltan2{

/*! \brief MachineClass
 *  Base class for representing machine coordinates, networks, etc.
 */
template <typename pcoord_t, typename part_t>
class Machine{

protected:
    int numRanks;
    int myRank;

public:
    /*! \brief Constructor MachineRepresentation Class
     *  \param comm_ Communication object.
     */

    Machine(const Teuchos::Comm<int> &comm) :
      numRanks(comm.getSize()), myRank(comm.getRank())
    { }

    virtual ~Machine(){ }

    /*! \brief indicates whether or not the machine has coordinates
     */
    bool hasMachineCoordinates() const {
      return false;  // Coordinates not available in this machine
    }

    /*! \brief returns the dimension (number of coords per node) in the machine
     */
    int getMachineDim() const {
      return 0;  // Coordinates not available in this machine
    }

    /*! \brief sets the number of unique coordinates in each machine dimension
     *  return true if coordinates are available
     */
    bool getMachineExtent(int *nxyz) const {
      return false;  // Extent not available in this machine
    }

    /*! \brief if the machine has a wrap-around tourus link in each dimension.
     *  return true if the information is available
     */
    bool getMachineExtentWrapArounds(bool *wrap_around) const {
      return false;  // Extent not available in this machine
    }

    /*! \brief getMyCoordinate function
     *  set the machine coordinate xyz of the current process
     *  return true if current process' coordinates are available
     */
    bool getMyMachineCoordinate(pcoord_t *xyz) const { 
      return false;  // Coordinates not available in this machine
    }

    /*! \brief getCoordinate function
     *  set the machine coordinate xyz of any rank process
     *  return true if coordinates are available by rank
     */
    bool getMachineCoordinate(const int rank, pcoord_t *xyz) const {
      return false;  // Coordinates not available by rank
    }

    /*! \brief getCoordinate function
     *  set the machine coordinate xyz of any node by nodename
     *  return true if coordinates are available by nodename
     */
    bool getMachineCoordinate(const char *nodename, pcoord_t *xyz) const {
      return false;  // Coordinates not available by nodename
    }

    /*! \brief getProcDim function
     *  set the coordinates of all ranks 
     *  allCoords[i][j], i=0,...,getMachineDim(), j=0,...,getNumRanks(),
     *  is the i-th dimensional coordinate for rank j.
     *  return true if coordinates are available for all ranks
     */
    bool getAllMachineCoordinatesView(pcoord_t **allCoords) const { 
      return false;  // Coordinates not available in this machine
    }

    /*! \brief getNumRanks function 
     *  return the number of ranks.
     */
    int getNumRanks() const { return numRanks; }

    /*! \brief getHopCount function
     *  set hops between rank1 and rank2
     *  return true if coordinates are available 
     */
    virtual bool getHopCount(int rank1, int rank2, pcoord_t &hops) const {
      return false;
    }

    /*! \brief getNumUniqueGroups function
     *  return the number of unique Dragonfly network groups in provided 
     *  allocation.
     *
     *  Equals the length of group_count member data, if available,
     *  otherwise we consider the whole allocation to be one group.
     */
    virtual part_t getNumUniqueGroups() const {
      return 1;
    }

    /*! \brief getGroupCount function
     *  return the number of ranks in each group (RCA X-dim, e.g. first dim)
     *
     *  Ex, 4 ranks with coord (3, 1, 1) and 8 ranks with coord (5, 2, 4), 
     *  will produce
     *  grp_count = [0, 0, 0, 4, 0, 8, 0, ...] 
     *  which is trimmed and returned as
     *  grp_count = [4, 8]
     *
     *
     *  (Currently only for Zoltan2_MachineDragonflyRCA, and used for 
     *  MultiJagged's first cut in "algorithms/partition/Zoltan2_TaskMapper.hpp"
     *  thru "problems/Zoltan2_MappingProblem.hpp".
     *  return true if group_count is available
     */
    virtual bool getGroupCount(part_t *grp_count) const {
      return false;
    }

    // KDD TODO: Add Graph interface and methods supporting full LDMS interface.

};
}
#endif
