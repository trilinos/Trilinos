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
class MachineTopoMgr : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  MachineTopoMgr(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
#if defined (CMK_BLUEGENEQ)
    networkDim(6),  tmgr(comm.getSize()),
#elif defined (CMK_BLUEGENEP)
    networkDim(4), tmgr(comm.getSize()),
#else
    networkDim(3),
#endif
    procCoords(NULL), machine_extent(NULL)
  {
    machine_extent = new int[networkDim];
    this->getMachineExtent(this->machine_extent);
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

  virtual ~MachineTopoMgr() {
    for (int i = 0; i < networkDim; i++){
      delete [] procCoords[i];
    }
    delete [] procCoords;
    delete [] machine_extent;
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return networkDim; }

  bool getMachineExtent(int *nxyz) const {
#if defined (CMK_BLUEGENEQ)
    nxyz[0] = tmgr.getDimNA();
    nxyz[1] = tmgr.getDimNB();
    nxyz[2] = tmgr.getDimNC();
    nxyz[3] = tmgr.getDimND();
    nxyz[4] = tmgr.getDimNE();
    nxyz[5] = tmgr.getDimNT();
#elif defined (CMK_BLUEGENEP)
    nxyz[0] = tmgr.getDimNX();
    nxyz[1] = tmgr.getDimNY();
    nxyz[2] = tmgr.getDimNZ();
    nxyz[3] = tmgr.getDimNT();
#else
    return false;
#endif
    return true;
  }

  //MD TODO: Not always it has wrap-around links.
  bool getMachineExtentWrapArounds(bool *wrap_around) const {
#if defined (CMK_BLUEGENEQ)
    wrap_around[0] = true;
    wrap_around[1] = true;
    wrap_around[2] = true;
    wrap_around[3] = true;
    wrap_around[4] = true;
    wrap_around[5] = true;
#elif defined (CMK_BLUEGENEP)
    wrap_around[0] = true;
    wrap_around[1] = true;
    wrap_around[2] = true;
    wrap_around[3] = true;
#else
#endif
    return true;
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {
#if defined (CMK_BLUEGENEQ)
    int a,b,c,d,e,t;
    tmgr.rankToCoordinates(this->myRank, a,b,c,d,e,t);
    xyz[0] = a; xyz[1] = b; xyz[2] = c; xyz[3] = d; xyz[4] = e; xyz[5] = t;
    //std::cout << "me:" << this->myRank << " " << a << " " << b << " " << c << " " << d << " " << e << " " << t << std::endl;
#elif defined (CMK_BLUEGENEP)
    int a,b,c,t;
    tmgr.rankToCoordinates(this->myRank, a,b,c,t);
    xyz[0] = a; xyz[1] = b; xyz[2] = c; xyz[3] = t;
#else
    return false;
#endif
   	
    return true;
  }

  bool getMachineExtentWrapArounds(part_t *wrap_around) const {
    wrap_around[0] = true;
    wrap_around[1] = true;
    wrap_around[2] = true;
    wrap_around[3] = true;
    wrap_around[4] = true;
    wrap_around[5] = true;
    return true;
  }
  inline bool getMachineCoordinate(const int rank,
                                   pcoord_t *xyz) const {
    return false;
  }


  bool getMachineCoordinate(const char *nodename, pcoord_t *xyz) {
    return false;  // cannot yet return from nodename
  }

  bool getAllMachineCoordinatesView(pcoord_t **&allCoords) const {
    allCoords = procCoords;
    return true;
  }

  virtual bool getHopCount(int rank1, int rank2, pcoord_t &hops){
    hops = 0;
    for (int i = 0; i < networkDim - 1; ++i){
      pcoord_t distance = procCoords[i][rank1] - procCoords[i][rank2];
      if (distance < 0 ) distance = -distance;
      if (machine_extent[i] - distance < distance) distance = machine_extent[i] - distance;
      hops += distance;
    }
    return true;
  }


private:

  int networkDim;

#ifdef HAVE_ZOLTAN2_TOPOMANAGER
  TopoManager tmgr;
#endif
  pcoord_t **procCoords;   // KDD Maybe should be RCP?

  part_t *machine_extent;

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
