#ifndef _ZOLTAN2_MACHINE_RCALIB_HPP_
#define _ZOLTAN2_MACHINE_RCALIB_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#ifdef HAVE_ZOLTAN2_RCALIB
extern "C"{
#include <rca_lib.h>
}
#endif


namespace Zoltan2{

/*! \brief A Machine Class for testing only
 *  A more realistic machine should be used for task mapping.
 */

template <typename pcoord_t, typename part_t>
class MachineRCA : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  MachineRCA(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
    networkDim(3),  
    procCoords(NULL), machine_extent(NULL),
    pl(NULL)
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
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; i++)
      procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    //reduceAll the coordinates of each processor.
    gatherMachineCoordinates(comm);

  }

  MachineRCA(const Teuchos::Comm<int> &comm, const Teuchos::ParameterList &pl_ ):
    Machine<pcoord_t,part_t>(comm),
    networkDim(3),
    procCoords(NULL), machine_extent(NULL),
    pl(&pl_)
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
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; i++)
      procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    //reduceAll the coordinates of each processor.
    gatherMachineCoordinates(comm);

    const Teuchos::ParameterEntry *pe = this->pl->getEntryPtr("machine_coord_transformation");
    this->printAllocation();
    if (pe){

      std::string approach;
      approach = pe->getValue<std::string>(&approach);

      if (approach == "EIGNORE"){
        for (int i = 0; i < this->numRanks; ++i){
          this->transform_coords( procCoords[0][i], procCoords[1][i], procCoords[2][i]);
        }
       	pcoord_t mx = this->machine_extent[0];
        pcoord_t my = this->machine_extent[1];
        pcoord_t mz = this->machine_extent[2];

        this->transform_coords(mx, my, mz);
        this->machine_extent[0] = mx;
        this->machine_extent[1] = my;
        this->machine_extent[2] = mz;


        if(this->myRank == 0) std::cout << "Transforming the coordinates" << std::endl;
        this->printAllocation();
      }
    }
  }

  


  virtual ~MachineRCA() {
    for (int i = 0; i < networkDim; i++){
      delete [] procCoords[i];
    }
    delete [] procCoords;
    delete [] machine_extent;
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return this->networkDim;/*transformed_network_dim;*/  }

  bool getMachineExtent(int *nxyz) const {
#if defined (HAVE_ZOLTAN2_RCALIB)
    mesh_coord_t mxyz;
    rca_get_max_dimension(&mxyz);
    int dim = 0;
    nxyz[dim++] = mxyz.mesh_x + 1; //x 
    nxyz[dim++] = mxyz.mesh_y + 1; //y
    nxyz[dim++] = mxyz.mesh_z + 1; //z
    return true;
#else
    return false;
#endif
  }

  void printAllocation(){
    if(this->myRank == 0){
      for (int i = 0; i < this->numRanks; ++i){ 
        std::cout << "Rank:" << i << " " << procCoords[0][i] << " " << procCoords[1][i] << " " << procCoords[2][i] << std::endl;
      } 
      std::cout << "Machine Extent:" << " " << this->machine_extent[0] << " " << this->machine_extent[1] << " " << this->machine_extent[2] << std::endl;
    }
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {
    for (int i = 0; i < this->networkDim; ++i){
      xyz[i] = procCoords[i][this->myRank];
    }
    return true;
  }

  bool getMyActualMachineCoordinate(pcoord_t *xyz) {
#if defined (HAVE_ZOLTAN2_RCALIB)
    rs_node_t nodeInfo;  /* Cray node info for node running this function */
    rca_get_nodeid(&nodeInfo);
    int NIDs = (int)nodeInfo.rs_node_s._node_id;  /* its node ID */

    mesh_coord_t node_coord;
    int returnval = rca_get_meshcoord((uint16_t)NIDs, &node_coord);
    if (returnval == -1){
      return false;
    }
    xyz[0] = node_coord.mesh_x;
    xyz[1] = node_coord.mesh_y;
    xyz[2] = node_coord.mesh_z;
    return true;
#else
    return false;
#endif
  }

  inline bool getMachineCoordinate(const int rank,
                                   pcoord_t *xyz) const {
    for (int i = 0; i < this->networkDim; ++i){
      xyz[i] = procCoords[i][rank];
    }
    return true;
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
    for (int i = 0; i < networkDim; ++i){
      pcoord_t distance = procCoords[i][rank1] - procCoords[i][rank2];
      if (distance < 0 ) distance = -distance;
      if (machine_extent[i] - distance < distance) distance = machine_extent[i] - distance;
      hops += distance;
    }
    return true;
  }


private:

  int networkDim;

  pcoord_t **procCoords; //, **transformed_proc_coords;   // KDD Maybe should be RCP?
  part_t *machine_extent; //, **transformed_machine_extent;
  const Teuchos::ParameterList *pl;
  //bool delete_tranformed_coords;

/*
  bool delete_transformed_coords;
  int transformed_network_dim;
  pcoord_t **transformed_coordinates;
*/
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

  //on titan X links are all uniform. Their bandwidth is 75 GB/s. 
  //distance is related to 1/bandwidth, so we choose to scale it with 8. (75x8 = 600)
  pcoord_t _get_transformed_x(pcoord_t x){
    return x * 8;
  }
  
  //on titan Y links are 2 kinds. The ones between [0-1] [2-3] [4-5] [6-7] and so on, are mezzanine links.
  //their bandwidth is 75 GB/s. They are scaled with 8.
  //The distances betwen [1-2] [3-4] [5-6] and so on are Y cables, with bandwidth 37.5 GB/s. They are scaled with 16.
  pcoord_t _get_transformed_y(pcoord_t y){
    size_t yordinal = y;
    return yordinal/2 * (16 + 8) + (yordinal %2) * 8;
  }
  
  //on titan X links are also 2 kinds. The ones between [0-7] [8-15] [16-23]  and so on, are backplane links.
  //their bandwidth is 120 GB/s. They are scaled with 5. (120 * 5 = 600)
  //The distances betwen [7-8] [15-16] [23-0] and so on are Y cables, with bandwidth 75 GB/s. They are scaled with 8.
  pcoord_t _get_transformed_z(pcoord_t z){
    size_t zordinal = z;
    //return zordinal  /8 * (8)   + (zordinal/8 * 7 + zordinal % 8) * 5;
    return zordinal * 5 + (zordinal / 8) * 3;
  }
  
  void transform_coords( pcoord_t &x,  pcoord_t &y,  pcoord_t &z){
    x = this->_get_transformed_x(x);
    y = this->_get_transformed_y(y);
    z = this->_get_transformed_z(z);
  }
 
};
}
#endif
