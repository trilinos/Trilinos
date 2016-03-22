#ifndef _ZOLTAN2_MACHINE_TOPOMANAGERTEST_HPP_
#define _ZOLTAN2_MACHINE_TOPOMANAGERTEST_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>


namespace Zoltan2{

/*! \brief A Machine Class for testing only
 *  A more realistic machine should be used for task mapping.
 */

template <typename pcoord_t, typename part_t>
class MachineBGQTest : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  MachineBGQTest(const Teuchos::Comm<int> &comm, const Teuchos::ParameterList &pl ):
      Machine<pcoord_t,part_t>(comm),
      networkDim(6),
      procCoords(NULL){
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

  MachineBGQTest(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
    networkDim(6),
    procCoords(NULL)
  {

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

  virtual ~MachineBGQTest() {
    for (int i = 0; i < networkDim; i++){
      delete [] procCoords[i];
    }
    delete [] procCoords;
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return networkDim; }

  bool getMachineExtentWrapArounds(part_t *wrap_around) const {
    wrap_around[0] = true;
    wrap_around[1] = true;
    wrap_around[2] = true;
    wrap_around[3] = true;
    wrap_around[4] = true;
    wrap_around[5] = true;
    return true;
  }


  bool getMachineExtentWrapArounds(bool *wrap_around) const {
    wrap_around[0] = wrap_around[1] = wrap_around[2] =
        wrap_around[3] = wrap_around[4] = wrap_around[5] = true;
    return true;
  }

  bool getMachineExtent(part_t *nxyz) const {
    const int rank_per_node = 1;
    switch (this->numRanks / rank_per_node){
    case 0:
    case 1:
      nxyz[0] = nxyz[1] = nxyz[2] = nxyz[3] = 1;
      nxyz[4] = 1;
      nxyz[5] = 1;
      break;
    case 2:
      nxyz[0] = nxyz[1] = nxyz[2] = nxyz[3] = 1;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 4:
      nxyz[0] = nxyz[1] = nxyz[2] = 1;
      nxyz[3] = 2;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 8:
      nxyz[0] = nxyz[1] = 1;
      nxyz[2] = 2;
      nxyz[3] = 2;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 16:
      nxyz[0] = 1;
      nxyz[1] = 2;
      nxyz[2] = 2;
      nxyz[3] = 2;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 32:
      nxyz[0] = 1;
      nxyz[1] = 2;
      nxyz[2] = 2;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 64:
      nxyz[0] = 1;
      nxyz[1] = 2;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 128:
      nxyz[0] = 1;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 256:
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 1;
      nxyz[5] = 1;
      break;
    case 512:
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 1024:
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 8;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 2048:
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 16;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    case 4096:
      nxyz[0] = 4;
      nxyz[1] = 8;
      nxyz[2] = 4;
      nxyz[3] = 16;
      nxyz[4] = 2;
      nxyz[5] = 1;
      break;
    default:
      std::cerr << "Too many ranks to test" << std::endl;
      break;
    }

    nxyz[5] = rank_per_node;
    return true;
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {
    int a,b,c,d,e,t;

    int me = this->myRank;
    t = me % machine_extent[5];

    me = me / machine_extent[5];
    e = me % machine_extent[4];

    me = me / machine_extent[4];
    d = me % machine_extent[3];

    me = me / machine_extent[3];
    c = me % machine_extent[2];

    me = me / machine_extent[2];
    b = me % machine_extent[1];

    me = me / machine_extent[1];
    a = me % machine_extent[0];

    xyz[0] = a; xyz[1] = b; xyz[2] = c; xyz[3] = d; xyz[4] = e; xyz[5] = t;

    //std::cout << "me:" << this->myRank << " " << a << " " << b << " " << c << " " << d << " " << e << " " << t << std::endl;
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
  pcoord_t **procCoords;   // KDD Maybe should be RCP?
  part_t machine_extent[6];

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
