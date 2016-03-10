#ifndef _ZOLTAN2_MACHINEREPRESENTATION_HPP_
#define _ZOLTAN2_MACHINEREPRESENTATION_HPP_

#include <Teuchos_Comm.hpp>
#include <Zoltan2_MachineForTesting.hpp>
#include <Zoltan2_MachineTopoMgr.hpp>
namespace Zoltan2{

/*! \brief MachineRepresentation Class
 *  Base class for representing machine coordinates, networks, etc.
 */
template <typename pcoord_t, typename part_t>
class MachineRepresentation{

public:

#if defined(HAVE_ZOLTAN2_LDMS)
    typedef MachineLDMS<pcoord_t,part_t> machine_t;
#elif defined(HAVE_ZOLTAN2_RCA)
    typedef MachineRCA<pcoord_t,part_t> machine_t;
#elif defined(HAVE_ZOLTAN2_TOPOMANAGER)
    typedef MachineTopoMgr<pcoord_t,part_t> machine_t;
#else
    typedef MachineForTesting<pcoord_t,part_t> machine_t;
#endif

    /*! \brief Constructor MachineRepresentation Class
     *  \param comm_ Communication object.
     */

    MachineRepresentation(const Teuchos::Comm<int> &comm) :
      machine(new machine_t(comm)) { }

    ~MachineRepresentation() {delete machine;}

    // Interface functions follow.
    // They are just wrappers around the specific machine's functions.

    /*! \brief indicates whether or not the machine has coordinates
     */
    inline bool hasMachineCoordinates() const { 
      return machine->hasMachineCoordinates();
    }

    /*! \brief returns the dimension (number of coords per node) in the machine
     */
    inline int getMachineDim() const { return machine->getMachineDim(); }

    /*! \brief sets the number of unique coordinates in each machine dimension
     *  return true if coordinates are available
     */
    inline bool getMachineExtent(int *nxyz) const { 
      return machine->getMachineExtent(nxyz);
    }

    /*! \brief getMyCoordinate function
     *  set the machine coordinate xyz of the current process
     *  return true if current process' coordinates are available
     */
    inline bool getMyMachineCoordinate(pcoord_t *xyz) const { 
      return machine->getMyMachineCoordinate(xyz);
    }

    /*! \brief getCoordinate function
     *  set the machine coordinate xyz of any rank process
     *  return true if coordinates are available by rank
     */
    inline bool getMachineCoordinate(const part_t rank,
                                     pcoord_t *xyz) const {
      return machine->getMachineCoordinate(rank, xyz);
    }

    /*! \brief getCoordinate function
     *  set the machine coordinate xyz of any node by nodename
     *  return true if coordinates are available by nodename
     */
    inline bool getMachineCoordinate(const char *nodename,
                                     pcoord_t *xyz) const {
      return machine->getMachineCoordinate(nodename, xyz);
    }

    /*! \brief getProcDim function
     *  set the coordinates of all ranks 
     *  allCoords[i][j], i=0,...,getMachineDim(), j=0,...,getNumRanks(),
     *  is the i-th dimensional coordinate for rank j.
     *  return true if coordinates are available for all ranks
     */
    inline bool getAllMachineCoordinatesView(pcoord_t **&allCoords) const { 
      return machine->getAllMachineCoordinatesView(allCoords);
    }

    /*! \brief return the number of ranks.
     */
    inline int getNumRanks() const { return machine->getNumRanks(); }

    // KDD TODO: Add Graph interface and methods supporting full LDMS interface.

private:
    machine_t *machine;
};
}
#endif
