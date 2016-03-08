#ifndef _ZOLTAN2_MACHINE_TOPOMANAGER_HPP_
#define _ZOLTAN2_MACHINE_TOPOMANAGER_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#ifdef HAVE_ZOLTAN2_TOPOMANAGER
#include <TopoManager.h>
#endif

namespace Zoltan2{

/*! \brief A Machine Class for testing only
 *  A more realistic machine should be used for task mapping.
 */

template <typename pcoord_t, typename part_t>
class Zoltan2_MachineTopoMgr : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  Zoltan2_MachineTopoMgr(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
#if defined (CMK_BLUEGENEQ)
    networkDim(6),  tmgr(comm.getSize()),
#elif defined (CMK_BLUEGENEP)
    networkDim(4), tmgr(comm.getSize()),
#else
    networkDim(3),
#endif
    procCoords(NULL)
  {
    //allocate memory for processor coordinates.
    procCoords = new pcoord_t *[networkDim];
    for (int i = 0; i < networkDim; ++i){
      procCoords[i] = new pcoord_t[this->numRanks];
      memset(procCoords[i], 0, sizeof(pcoord_t) * this->numRanks);
    }

    //obtain the coordinate of the processor.
    pcoord_t *xyz = new pcoord_t[networkDim];
    getMyMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; i++)
      procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    //reduceAll the coordinates of each processor.
    gatherMachineCoordinates(comm);
  }

  virtual ~MachineForTesting() {
    for (int i = 0; i < networkDim; i++){
      delete [] procCoords[i];
    }
    delete [] procCoords;
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return networkDim; }

  bool getMachineExtent(part_t *nxyz) const {
#if defined (CMK_BLUEGENEQ)
    nxyz[0] = tmgr.getDimNA();
    nxyz[1] = tmgr.getDimNB();
    nxyz[2] = tmgr.getDimNC();
    nxyz[3] = tmgr.getDimND();
    nxyz[4] = tmgr.getDimNE();
    nxyz[5] = tmgr.getDimNT());
#elif defined (CMK_BLUEGENEP)
    nxyz[0] = tmgr.getDimNX();
    nxyz[1] = tmgr.getDimNY();
    nxyz[2] = tmgr.getDimNZ();
    nxyz[3] = tmgr.getDimNT();
#else
#endif
    return true;
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {

#if defined (CMK_BLUEGENEQ)
    tmgr.rankToCoordinates(i, xyz[0], xyz[1], xyz[2], xyz[3], xyz[4], xyz[5]);
#elif defined (CMK_BLUEGENEP)
    tmgr.rankToCoordinates(i, xyz[0], xyz[1], xyz[2], xyz[3]);
#else
#endif
  }



  bool getMachineCoordinate(const char *nodename, pcoord_t *xyz) {
    return false;  // cannot yet return from nodename
  }

  bool getAllMachineCoordinatesView(pcoord_t **&allCoords) const {
    allCoords = procCoords;
    return true;
  }

private:

  int networkDim;

#ifdef HAVE_ZOLTAN2_TOPOMANAGER
  TopoManager tmgr;
#endif
  pcoord_t **procCoords;   // KDD Maybe should be RCP?


  void gatherMachineCoordinates(const Teuchos::Comm<int> &comm) {
    // reduces and stores all machine coordinates.
    pcoord_t *tmpVect = new pcoord_t [this->numRanks];

    for (int i = 0; i < networkDim; i++) {
      Teuchos::reduceAll<int, pcoord_t>(comm, Teuchos::REDUCE_SUM,
                                        this->numRanks, procCoords[i], tmpVect);
      pcoord_t *tmp = tmpVect;
      tmpVect = procCoords[i];
      procCoords[i] = tmp;
    }
    delete [] tmpVect;
  }
};
}
#endif
