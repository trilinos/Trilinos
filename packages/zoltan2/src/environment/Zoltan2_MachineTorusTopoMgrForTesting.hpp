#ifndef _ZOLTAN2_MACHINE_TORUS_TOPOMANAGERTEST_HPP_
#define _ZOLTAN2_MACHINE_TORUS_TOPOMANAGERTEST_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>


namespace Zoltan2{

/*! \brief A TopoManager Machine Class (Torus Networks) for testing only
 *  A more realistic machine should be used for task mapping.
 */

template <typename pcoord_t, typename part_t>
class MachineTorusBGQTest : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  MachineTorusBGQTest(const Teuchos::Comm<int> &comm ):
      Machine<pcoord_t,part_t>(comm),
      networkDim(6),
      procCoords(NULL),machine_extent(NULL),
      delete_transformed_coords(false), 
      transformed_network_dim(0),
      transformed_coordinates (NULL), pl(NULL) 
  {
    transformed_network_dim = networkDim - 1;
    transformed_coordinates = procCoords;

    machine_extent = new int[networkDim];
    this->getMachineExtent(this->machine_extent);
    machine_extent[5] = 1;
    
    // Allocate memory for processor coordinates.
    procCoords = new pcoord_t *[networkDim];
    for (int i = 0; i < networkDim; ++i) {
      procCoords[i] = new pcoord_t[this->numRanks];
      memset(procCoords[i], 0, sizeof(pcoord_t) * this->numRanks);
    }

    // Obtain the coordinate of the processor.
    pcoord_t *xyz = new pcoord_t[networkDim];
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; i++)
      procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    // reduceAll the coordinates of each processor.
    gatherMachineCoordinates(comm);

  }

  MachineTorusBGQTest(const Teuchos::Comm<int> &comm, 
                 const Teuchos::ParameterList &pl_):
    Machine<pcoord_t,part_t>(comm),
    networkDim(6),
    procCoords(NULL),machine_extent(NULL),
    delete_transformed_coords(false), 
    transformed_network_dim(0),
    transformed_coordinates (NULL),
    pl(&pl_)
  {
    transformed_network_dim = networkDim - 1;
    transformed_coordinates = procCoords;
    machine_extent = new int[networkDim];

    this->getMachineExtent(this->machine_extent);
    machine_extent[5] = 1;
    
    // Allocate memory for processor coordinates.
    procCoords = new pcoord_t *[networkDim];
    for (int i = 0; i < networkDim; ++i) {
      procCoords[i] = new pcoord_t[this->numRanks];
      memset(procCoords[i], 0, sizeof(pcoord_t) * this->numRanks);
    }

    // Obtain the coordinate of the processor.
    pcoord_t *xyz = new pcoord_t[networkDim];
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; i++)
      procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    // reduceAll the coordinates of each processor.
    gatherMachineCoordinates(comm);

    const Teuchos::ParameterEntry *pe = 
      this->pl->getEntryPtr("Machine_Optimization_Level");

    if (pe) {
      int optimization_level = 0;

      optimization_level = pe->getValue<int>(&optimization_level);

      if (optimization_level == 0) {
        transformed_network_dim = networkDim - 1;
        transformed_coordinates = procCoords;
      }
      else if (optimization_level >= 1) {
        transformed_network_dim = networkDim - 2;
        transformed_coordinates = procCoords;
      }
    }
  }


  virtual ~MachineTorusBGQTest() {
    for (int i = 0; i < networkDim; i++) {
      delete [] procCoords[i];
    }
    delete [] procCoords;

    delete [] machine_extent;
    if (delete_transformed_coords) {
      for (int i = 0; i < transformed_network_dim; i++) {
        delete [] transformed_coordinates[i];
      }
      delete [] transformed_coordinates;
    }
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return transformed_network_dim; }

  bool getMachineExtentWrapArounds(bool *wrap_around) const {
    int dim = 0;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    return true;
  }


  bool getMachineExtentWrapArounds(part_t *wrap_around) const {
    int dim = 0;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    return true;
  }

  bool getMachineExtent(part_t *extent) const {
    part_t nxyz[6];

    nxyz[0] = 1;
    nxyz[1] = 1;
    nxyz[2] = 1;
    nxyz[3] = 1;
    nxyz[4] = 1;
    nxyz[5] = 1;
    const int rank_per_node = 1;
    const int num_nodes = this->numRanks / rank_per_node;

    if (num_nodes <= 1) {
      nxyz[0] = nxyz[1] = nxyz[2] = nxyz[3] = 1;
      nxyz[4] = 1;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 2) {
      nxyz[0] = nxyz[1] = nxyz[2] = nxyz[3] = 1;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 4) {
      nxyz[0] = nxyz[1] = nxyz[2] = 1;
      nxyz[3] = 2;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 8) {
      nxyz[0] = nxyz[1] = 1;
      nxyz[2] = 2;
      nxyz[3] = 2;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 16) {
      nxyz[0] = 1;
      nxyz[1] = 2;
      nxyz[2] = 2;
      nxyz[3] = 2;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 32) {
      nxyz[0] = 1;
      nxyz[1] = 2;
      nxyz[2] = 2;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 64) {
      nxyz[0] = 1;
      nxyz[1] = 2;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 128) {
      nxyz[0] = 2;
      nxyz[1] = 2;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 256) {
      nxyz[0] = 4;
      nxyz[1] = 2;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 512) {
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 4;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 1024) {
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 8;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 2048) {
      nxyz[0] = 4;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 16;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else if (num_nodes <= 4096) {
      nxyz[0] = 8;
      nxyz[1] = 4;
      nxyz[2] = 4;
      nxyz[3] = 16;
      nxyz[4] = 2;
      nxyz[5] = 1;
    }
    else {
      std::cerr << "Too many ranks to test" << std::endl;
    }
    nxyz[5] = rank_per_node;

    int dim = 0;
    if (dim < transformed_network_dim)
      extent[dim] = nxyz[dim];
    ++dim;
    if (dim < transformed_network_dim)
      extent[dim] = nxyz[dim];
    ++dim;
    if (dim < transformed_network_dim)
      extent[dim] = nxyz[dim];
    ++dim;
    if (dim < transformed_network_dim)
      extent[dim] = nxyz[dim];
    ++dim;
    if (dim < transformed_network_dim)
      extent[dim] = nxyz[dim];
    ++dim;
    if (dim < transformed_network_dim)
      extent[dim] = nxyz[dim];

    return true;
  }


  bool getMyMachineCoordinate(pcoord_t *xyz) {
    for (int i = 0; i < this->transformed_network_dim; ++i) {
      xyz[i] = transformed_coordinates[i][this->myRank];
    }
    return true;
  }

  bool getMyActualMachineCoordinate(pcoord_t *xyz) {
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

//    std::cout << "me:" << this->myRank << " " << a << " " << b 
//      << " " << c << " " << d << " " << e << " " << t << std::endl;
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

  virtual bool getHopCount(int rank1, int rank2, pcoord_t &hops) const override {

    hops = 0;
    for (int i = 0; i < networkDim - 1; ++i) {
      pcoord_t distance = procCoords[i][rank1] - procCoords[i][rank2];
      if (distance < 0 ) 
        distance = -distance;
      if (machine_extent[i] - distance < distance) 
        distance = machine_extent[i] - distance;
      hops += distance;
    }

/*
    if (this->myRank == 0) {
      std::cout << "rank1:" << rank1 << " rank2:" << rank2 
        << " hops:" << hops << std::endl;
    }
*/

    return true;
  }


private:

  int networkDim;
  pcoord_t **procCoords;   // KDD Maybe should be RCP?
  part_t *machine_extent;

  bool delete_transformed_coords;
  int transformed_network_dim;
  pcoord_t **transformed_coordinates;
  const Teuchos::ParameterList *pl;

  void gatherMachineCoordinates(const Teuchos::Comm<int> &comm) {
    // Reduces and stores all machine coordinates.
    pcoord_t *tmpVect = new pcoord_t [this->numRanks];

    for (int i = 0; i < networkDim; i++) {
      Teuchos::reduceAll<int, pcoord_t>(comm, Teuchos::REDUCE_SUM,
                                        this->numRanks, 
                                        procCoords[i], tmpVect);
      pcoord_t *tmp = tmpVect;
      tmpVect = procCoords[i];
      procCoords[i] = tmp;
    }
    delete [] tmpVect;
  }
};
}
#endif
