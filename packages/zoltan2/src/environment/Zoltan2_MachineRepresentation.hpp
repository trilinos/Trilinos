#ifndef _ZOLTAN2_MACHINEREPRESENTATION_HPP_
#define _ZOLTAN2_MACHINEREPRESENTATION_HPP_

#include <Teuchos_Comm.hpp>
#include <Zoltan2_MachineForTesting.hpp>
#include <Zoltan2_MachineTopoMgr.hpp>
#include <Zoltan2_MachineTopoMgrForTest.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Zoltan2_MachineRCA.hpp>
#include <Zoltan2_MachineRCAForTest.hpp>
#include <Zoltan2_Environment.hpp>

//#define HAVE_ZOLTAN2_BGQTEST
namespace Zoltan2{

/*! \brief MachineRepresentation Class
 *  Base class for representing machine coordinates, networks, etc.
 */
template <typename pcoord_t, typename part_t>
class MachineRepresentation{

public:
    typedef pcoord_t machine_pcoord_t;
    typedef part_t machine_part_t;
#if defined(HAVE_ZOLTAN2_LDMS)
    typedef MachineLDMS<pcoord_t,part_t> machine_t;
#elif defined(HAVE_ZOLTAN2_RCALIB)
    typedef MachineRCA<pcoord_t,part_t> machine_t;
#elif defined(HAVE_ZOLTAN2_TOPOMANAGER)
    typedef MachineTopoMgr<pcoord_t,part_t> machine_t;
#elif defined(HAVE_ZOLTAN2_BGQTEST)
    typedef MachineBGQTest<pcoord_t,part_t> machine_t;
#else
    typedef MachineForTesting<pcoord_t,part_t> machine_t;
    //typedef MachineBGQTest<pcoord_t,part_t> machine_t;
    //typedef MachineRCATest<pcoord_t,part_t> machine_t;
#endif

    /*! \brief Constructor MachineRepresentation Class
     *  \param comm_ Communication object.
     */

    MachineRepresentation(const Teuchos::Comm<int> &comm) :
      machine(new machine_t(comm)) {
   }



    MachineRepresentation(const Teuchos::Comm<int> &comm, const Teuchos::ParameterList &pl) :
          machine(new machine_t(comm, pl)) { }

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

    /*! \brief if the machine has a wrap-around tourus link in each dimension.
     *  return true if the information is available
     */
    bool getMachineExtentWrapArounds(bool *wrap_around) const {
      return machine->getMachineExtentWrapArounds(wrap_around);
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
    inline bool getMachineCoordinate(const int rank,
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


    inline bool getHopCount(int rank1, int rank2, pcoord_t &hops) const {
      return machine->getHopCount(rank1, rank2, hops);
    }

    /*! \brief Set up validators specific to this Problem
    */
    static void getValidParameters(Teuchos::ParameterList & pl)
    {
      //TODO: This should be positive integer validator.
      pl.set("Machine_Optimization_Level", 10,
          "Machine Coordinate Transformation Method",
          Environment::getAnyIntValidator());

      // validator for file does not have to exist
      RCP<Teuchos::FileNameValidator> file_not_required_validator =
        Teuchos::rcp( new Teuchos::FileNameValidator(false) );

      // bool parameter
      pl.set("Input_RCA_Machine_Coords", "", "Input File for input machine coordinates",
          file_not_required_validator);
    }


    // KDD TODO: Add Graph interface and methods supporting full LDMS interface.

private:
    machine_t *machine;
};
}
#endif
