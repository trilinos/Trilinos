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

/*! \brief A Dragonfly (e.g. Cori & Trinity) Machine Class for 
 *  testing only. A more realistic machine should be used for 
 *  task mapping.
 *
 *  Nodes in Cori are divided into groups(RCA x_dim) of 96 nodes
 *  and all groups are connected with an all-to-all connection.
 *  Within a group, nodes are arranged into 6 rows (RCA y_dim) 
 *  and 16 columns (RCA z_dim).
 *  All nodes within a row are connected with an all-to-all. Same 
 *  for columns. Therefore: 
 *  (3, 2, 1) -> (3, 4, 11) will take 2 hops.
 *  (5, 1, 1) -> (5, 1, 5)  will take 1 hop.
 *
 *  We represent this "nearness" by transformation using a high 
 *  dimensional coord-system of size (1 + N_y + N_z). The first 
 *  element represents the group, the next N_y elements represent
 *  the row, the next N_z elements represent the columns.
 *  
 *  (9, 2, 1) in the high-dim space:
 *  __ __
 *  | 9 |______ Group, 1 element
 *  | 0 |
 *  | 0 |
 *  | 1 |
 *  | 0 |
 *  | 0 |
 *  | 0 |______ Row, N_y elements
 *  | 0 |
 *  | 1 |
 *  | 0 |
 *  | 0 |
 *  |...|
 *  | 0 |
 *  -- --______ Col, N_z elements 
 *
 *
 *  To assist with MultiJagged coord partitioning we stretch the 
 *  dimensions. If RCA coords are (3, 2, 14), we first transform 
 *  the X, Y, Z by 
 *  
 *  X_new = 2 * X * N_Y * N_Z;
 *  Y_new = Y + X * N_Y * N_Z;
 *  Z_new = Z + X * N_Y * N_Z;
 *
 *  Then transformed coords are (576, 290, 302) and in high-dim
 *  space:
 *
 *  (3,2,10) -> (576, 290, 296) -> 
 *
 *  (576,| 288, 288, 289, 288, 288, 288,| 288, ..., 288, 289, 288)
 *
 *  Now Coordinates are distance sqrt(2) apart if 1 hop, and
 *  distance 2 apart if 2 hops.  
 *
 *  NOTE: Does not account for dragonfly's dynamic routing 
 */

template <typename pcoord_t, typename part_t>
class MachineDragonfly : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: Dragonfly (e.g. Cori & Trinity) network 
   *  machine description;
   *
   *  Does not do coord transformation.
   *  
   *  \param comm Communication object.
   */

  MachineDragonfly(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
    transformed_networkDim(3), 
    actual_networkDim(3),
    transformed_procCoords(NULL), 
    actual_procCoords(NULL),
    transformed_machine_extent(NULL),
    actual_machine_extent(NULL),
    is_transformed(false), 
    pl(NULL)
  {
    actual_machine_extent = new int[actual_networkDim];
    this->getActualMachineExtent(this->actual_machine_extent);
    
    // transformed dims = 1 + N_y + N_z
    transformed_networkDim = 1 + actual_machine_extent[1] + 
      actual_machine_extent[2];
    actual_machine_extent = new int[transformed_networkDim];

    // allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];
    
    for (int i = 0; i < networkDim; ++i){
      actual_procCoords[i] = new pcoord_t[this->numRanks];
      memset(actual_procCoords[i], 0, 
          sizeof(pcoord_t) * this->numRanks);
    }

    pcoord_t *xyz = new pcoord_t[transform]
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; ++i)
      actual_procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    // reduceAll the coordinates of each processor.
    gatherActualMachineCoordinates(this->actual_procCoords,
        this->actual_networkDim, comm);    

/*    // 3v 
    machine_extent = new int[networkDim];
    this->getActualMachineExtent(this->machine_extent);
    // actual is transformed or not?
    actual_machine_extent = machine_extent;

    //allocate memory for processor coordinates.
    actual_procCoords = procCoords = new pcoord_t *[networkDim];
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
*/
  }

  // No necessary wrap arounds for dragonfly networks. Groups
  // have wrap around, but group all-to-all connection makes unneccessary.
  virtual bool getMachineExtentWrapArounds(bool *wrap_around) const {
    return false;
  }

  MachineRCA(const Teuchos::Comm<int> &comm, 
             const Teuchos::ParameterList &pl_ ):
    Machine<pcoord_t,part_t>(comm),
    transformed_networkDim(3), 
    actual_networkDim(3),
    transformed_procCoords(NULL), 
    actual_procCoords(NULL),
    transformed_machine_extent(NULL),
    actual_machine_extent(NULL),
    is_transformed(false), 
    pl(&pl_)
  {

    actual_machine_extent = new int[actual_networkDim];
    this->getActualMachineExtent(this->actual_machine_extent);
    
    // transformed dims = 1 + N_y + N_z
    transformed_networkDim = 1 + actual_machine_extent[1] + 
      actual_machine_extent[2];
    actual_machine_extent = new int[transformed_networkDim];

    // allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];
    
    for (int i = 0; i < networkDim; ++i){
      transformed_procCoords[i] = new pcoord_t[this->numRanks];
      memset(transformed_procCoords[i], 0, 
          sizeof(pcoord_t) * this->numRanks);
    }

    pcoord_t *xyz = new pcoord_t[transform]
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < networkDim; ++i)
      actual_procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

   // Transform with mach opt level
    if (pe2){
      int optimization_level;
      optimization_level = pe2->getValue<int>(&optimization_level);

      if (optimization_level == 1){
        is_transformed = true;
        // +1 = z coordinate
        this->networkDim = machine_extent[1] * machine_extent[2] + 1;  



        pcoord_t *xyz = new pcoord_t[transform]
        getMyActualMachineCoordinate(xyz);
        for (int i = 0; i < networkDim; ++i)
          actual_procCoords[i][this->myRank] = xyz[i];
        delete [] xyz;

        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i){
          procCoords[i] = new pcoord_t[this->numRanks] ;//this->proc_coords[permutation[i]];
        }
        for (int i = 0; i < this->numRanks; ++i){
          procCoords[0][i] = this->actual_procCoords[0][i] * 8;
          int yordinal = this->actual_procCoords[1][i];
          procCoords[1][i] = yordinal/2 * (16 + 8) + (yordinal %2) * 8;
          int zordinal = this->actual_procCoords[2][i];
          procCoords[2][i] = zordinal * 5 + (zordinal / 8) * 3;
        }
        int mx = this->machine_extent[0];
        int my = this->machine_extent[1];
        int mz = this->machine_extent[2];


        this->machine_extent = new int[networkDim];
        this->machine_extent[0] = mx * 8;
        this->machine_extent[1] = my/2 * (16 + 8) + (my %2) * 8;
        this->machine_extent[2] = mz * 5 + (mz / 8) * 3;
        if(this->myRank == 0) std::cout << "Transforming the coordinates" << std::endl;
        //this->printAllocation();
      }

      // reduceAll the transformed coordinates of each processor.
      gatherMachineCoordinates(this->transformed_procCoords, 
          this->transformed_networkDim, comm);    
    }
    
/*    
    actual_machine_extent = machine_extent = new int[networkDim];
    this->getActualMachineExtent(this->machine_extent);
    actual_machine_extent = machine_extent;

    //allocate memory for processor coordinates.
    actual_procCoords = procCoords = new pcoord_t *[networkDim];
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

    const Teuchos::ParameterEntry *pe2 = this->pl->getEntryPtr("Machine_Optimization_Level");
    //this->printAllocation();
    
    // Transform with mach opt level
    if (pe2){
      int optimization_level;
      optimization_level = pe2->getValue<int>(&optimization_level);

      if (optimization_level == 1){
        is_transformed = true;
        // +1 = z coordinate
        this->networkDim = machine_extent[1] * machine_extent[2] + 1;  
        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i){
          procCoords[i] = new pcoord_t[this->numRanks] ;//this->proc_coords[permutation[i]];
        }
        for (int i = 0; i < this->numRanks; ++i){
          procCoords[0][i] = this->actual_procCoords[0][i] * 8;
          int yordinal = this->actual_procCoords[1][i];
          procCoords[1][i] = yordinal/2 * (16 + 8) + (yordinal %2) * 8;
          int zordinal = this->actual_procCoords[2][i];
          procCoords[2][i] = zordinal * 5 + (zordinal / 8) * 3;
        }
        int mx = this->machine_extent[0];
        int my = this->machine_extent[1];
        int mz = this->machine_extent[2];


        this->machine_extent = new int[networkDim];
        this->machine_extent[0] = mx * 8;
        this->machine_extent[1] = my/2 * (16 + 8) + (my %2) * 8;
        this->machine_extent[2] = mz * 5 + (mz / 8) * 3;
        if(this->myRank == 0) std::cout << "Transforming the coordinates" << std::endl;
        //this->printAllocation();
      }
      else if(optimization_level >= 3){
        is_transformed = true;
        this->networkDim = 6;
        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i){
          procCoords[i] = new pcoord_t[this->numRanks] ;//this->proc_coords[permutation[i]];
        }

        //this->machine_extent[0] = this->actual_machine_extent
        this->machine_extent = new int[networkDim];

        this->machine_extent[0] = ceil (int (this->actual_machine_extent[0]) / 2.0) * 64 ;
        this->machine_extent[3] = 2 * 8 ;
        this->machine_extent[1] = ceil(int (this->actual_machine_extent[1])  / 2.0) * 8 * 2400;
        this->machine_extent[4] = 2 * 8;
        this->machine_extent[2] = ceil((int (this->actual_machine_extent[2])) / 8.0) * 160;
        this->machine_extent[5] = 8 * 5;

        for (int k = 0; k < this->numRanks ; k++){
          //This part is for titan.
          //But it holds for other 3D torus machines such as Bluewaters.

          //Bandwitdh along
          // X = 75
          // Y = 37.5 or 75 --- everyother has 37.5 --- Y[0-1] =75 but Y[1-2]=37.5
          // Z = 75 or 120 ---- Y[0-1-2-3-4-5-6-7] = 120, Y[7-8] = 75

          //Along X we make groups of 2. Then scale the distance with 64.
          //First dimension is represents x/2
          procCoords[0][k] = (int (this->actual_procCoords[0][k]) / 2) * 64;
          //Then the 3rd dimension is x%2. distance is scaled with 8, reversely proportional with bw=75
          procCoords[3][k] = (int (this->actual_procCoords[0][k]) % 2) * 8 ;

          //Along Y. Every other one has the slowest link. So we want distances between Y/2 huge.
          //We scale Y/2 with 2400 so that we make sure that it is the first one we divie.
          procCoords[1][k] = (int (this->actual_procCoords[1][k])  / 2) * 8 * 2400;
          //The other one is scaled with 8 as in X.
          procCoords[4][k] = (int (this->actual_procCoords[1][k])  % 2) * 8;

          //We make groups of 8 along Z. Then distances between these groups are scaled with 160.
          //So that it is more than 2x distance than the distance with X grouping.
          //That is we scale the groups of Zs with 160. Groups of X with 64.
          //Zs has 8 processors connecting them, while X has only one. We want to divide along
          //Z twice before dividing along X.
          procCoords[2][k] = ((int (this->actual_procCoords[2][k])) / 8) * 160;
          //In the second group everything is scaled with 5, as bw=120
          procCoords[5][k] = ((int (this->actual_procCoords[2][k])) % 8) * 5;
        }
      }
      else if(optimization_level == 2){
        //This is as above case. but we make groups of 3 along X instead.
        is_transformed = true;
        this->networkDim = 6;
        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i){
          procCoords[i] = new pcoord_t[this->numRanks] ;//this->proc_coords[permutation[i]];
        }

        //this->machine_extent[0] = this->actual_machine_extent
        this->machine_extent = new int[networkDim];

        this->machine_extent[0] = ceil(int (this->actual_machine_extent[0]) / 3.0) * 128 ;
        this->machine_extent[3] = 3 * 8 ;
        this->machine_extent[1] = ceil(int (this->actual_machine_extent[1])  / 2.0) * 8 * 2400;
        this->machine_extent[4] = 2 * 8;
        this->machine_extent[2] = ceil((int (this->actual_machine_extent[2])) / 8.0) * 160;
        this->machine_extent[5] = 8 * 5;


        for (int k = 0; k < this->numRanks ; k++){
          //This part is for titan.
          //But it holds for other 3D torus machines such as Bluewaters.

          //Bandwitdh along
          // X = 75
          // Y = 37.5 or 75 --- everyother has 37.5 --- Y[0-1] =75 but Y[1-2]=37.5
          // Z = 75 or 120 ---- Y[0-1-2-3-4-5-6-7] = 120, Y[7-8] = 75

          //In this case we make groups of 3. along X.
          procCoords[0][k] = (int (this->actual_procCoords[0][k]) / 3) * 128;
          //Then the 3rd dimension is x%2. distance is scaled with 8, reversely proportional with bw=75
          procCoords[3][k] = (int (this->actual_procCoords[0][k]) % 3) * 8 ;

          //Along Y. Every other one has the slowest link. So we want distances between Y/2 huge.
          //We scale Y/2 with 2400 so that we make sure that it is the first one we divie.
          procCoords[1][k] = (int (this->actual_procCoords[1][k])  / 2) * 8 * 2400;
          //The other one is scaled with 8 as in X.
          procCoords[4][k] = (int (this->actual_procCoords[1][k])  % 2) * 8;


          procCoords[2][k] = ((int (this->actual_procCoords[2][k])) / 8) * 160;
          //In the second group everything is scaled with 5, as bw=120
          procCoords[5][k] = ((int (this->actual_procCoords[2][k])) % 8) * 5;
        }
      }
    }
*/

  }

  virtual ~MachineRCA() {
    if (is_transformed){
      is_transformed = false;
      for (int i = 0; i < actual_networkDim; ++i){
        delete [] actual_procCoords[i];
      }
      delete [] actual_procCoords;
      delete [] actual_machine_extent;
    }
    for (int i = 0; i < networkDim; ++i){
      delete [] procCoords[i];
    }
    delete [] procCoords;
    delete [] machine_extent;
  }

  bool hasMachineCoordinates() const { return true; }
  int getTransformedMachineDim() const { return this->transformed_networkDim; }
  int getActualMachineDim() const { return this->actual_networkDim;  }

  bool getTransformedMachineExtent(int *nxyz) const {
    if (is_transformed) {
      for (int dim = 0; dim < transformed_networkDim; ++dim)
        nxyz[dim] = this->transformed_machine_extent[dim]

      return true;
    }
    else 
      return false;
  }

  // RCA Machine Extents
  bool getActualMachineExtent(int *nxyz) const {
#if defined (HAVE_ZOLTAN2_RCALIB)
    mesh_coord_t mxyz;
    rca_get_max_dimension(&mxyz);

    int dim = 0;
    nxyz[dim++] = mxyz.mesh_x + 1; // X - group [0, ~100]
    nxyz[dim++] = mxyz.mesh_y + 1; // Y - row within group [0, 5]
    nxyz[dim++] = mxyz.mesh_z + 1; // Z - col within row [0, 15]
    return true;
#else
    return false;
#endif
  }

  void printAllocation(){
    if(this->myRank == 0){
      for (int i = 0; i < this->numRanks; ++i){ 
        std::cout << "Rank:" << i 
          << " " << procCoords[0][i] 
          << " " << procCoords[1][i] 
          << " " << procCoords[2][i] << std::endl;
      } 
      std::cout << "Machine Extent:" 
        << " " << this->machine_extent[0] 
        << " " << this->machine_extent[1] 
        << " " << this->machine_extent[2] << std::endl;
    }
  }

  bool getMyTransformedMachineCoordinate(pcoord_t *xyz) {
    for (int i = 0; i < this->networkDim; ++i){
      xyz[i] = procCoords[i][this->myRank];
    }
    return true;
  }

  bool getMyActualMachineCoordinate(pcoord_t *xyz) {
#if defined (HAVE_ZOLTAN2_RCALIB)
    // Cray node info for current node
    rs_node_t nodeInfo; 
    rca_get_nodeid(&nodeInfo);

    // Current node ID
    int NIDs = (int)nodeInfo.rs_node_s._node_id;

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

  inline bool getTransformedMachineCoordinate(const int rank,
                                        pcoord_t *xyz) const{
    for (int i = 0; i < this->transformed_networkDim; ++i){
      xyz[i] = procCoords[i][rank];
    }
    return true;
  }

  inline bool getActualMachineCoordinate(const int rank,
                                   pcoord_t *xyz) const {
    for (int i = 0; i < this->actual_networkDim; ++i){
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

  int transformed_networkDim;
  int actual_networkDim;

  pcoord_t **transformed_procCoords;
  pcoord_t **actual_procCoords;

  part_t *transformed_machine_extent;
  part_t *actual_machine_extent;
  bool is_transformed;


  const Teuchos::ParameterList *pl;
  //bool delete_tranformed_coords;

/*
  bool delete_transformed_coords;
  int transformed_network_dim;
  pcoord_t **transformed_coordinates;
*/
  void gatherMachineCoordinates(pcoord_t *coords, int netDim, 
      const Teuchos::Comm<int> &comm) {
    // reduces and stores all machine coordinates.
    pcoord_t *tmpVect = new pcoord_t [this->numRanks];

    for (int i = 0; i < netDim; ++i) {
      Teuchos::reduceAll<int, pcoord_t>(comm, Teuchos::REDUCE_SUM,
                                        this->numRanks, coords[i], tmpVect);
      pcoord_t *tmp = tmpVect;
      tmpVect = coords[i];
      coords[i] = tmp;
    }
    delete [] tmpVect;
  }

};
}
#endif
