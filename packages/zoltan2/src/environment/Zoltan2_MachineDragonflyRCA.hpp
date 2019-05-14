#ifndef _ZOLTAN2_MACHINE_DRAGONFLY_RCALIB_HPP_
#define _ZOLTAN2_MACHINE_DRAGONFLY_RCALIB_HPP_

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
 *  task mapping
 *
 *  Requires RCA library to run and $ZOLTAN2_MACHINE_DRAGONFLY = TRUE. 
 */

/*! \brief A Dragonfly (e.g. Cori & Trinity) RCA Machine Class
 *
 *  NOTE: REQUIRES RCA library to run
 *
 *  Nodes in Cori are divided into groups(RCA x_dim) of 384 nodes
 *  and all groups are connected with an all-to-all connection.
 *  Within a group, clusters of 4 nodes are arranged into 6 
 *  rows (RCA y_dim) and 16 columns (RCA z_dim).
 *  All nodes within a row are connected with an all-to-all. Same 
 *  for columns. Therefore: 
 *  (3, 2, 1) -> (3, 4, 11) will take 2 hops.
 *  (5, 1, 1) -> (5, 1, 5)  will take 1 hop.
 *
 *  We represent the "nearness" by transformation using a high 
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
 *  the X by 
 *
 *  X_new = 2 * X * N_Y * N_Z;
 *
 *  Then transformed coords are (576, 2, 14) and in high-dim
 *  space:
 *
 *  (3, 2, 14) -> (576, 2, 14) -> 
 *
 *  (576,| 0, 0, 1, 0, 0, 0,| 0, ..., 0, 1, 0)
 *  
 *  Now Coordinates are distance sqrt(2) apart if 1 hop, and
 *  distance 2 apart if 2 hops.  
 *
 *  NOTE: Does not account for dragonfly's dynamic routing 
 */
template <typename pcoord_t, typename part_t>
class MachineDragonflyRCA : public Machine <pcoord_t, part_t> {

public:

  /*! \brief Constructor: Dragonfly (e.g. Cori & Trinity) network 
   *  machine description;
   *
   *  Does not do coord transformation.
   *  
   *  \param comm Communication object.
   */
  MachineDragonflyRCA(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
    transformed_networkDim(3), 
    actual_networkDim(3),
    transformed_procCoords(NULL), 
    actual_procCoords(NULL),
    transformed_machine_extent(NULL),
    actual_machine_extent(NULL),
    is_transformed(false), 
    pl(NULL) {

    actual_machine_extent = new int[actual_networkDim];
    this->getActualMachineExtent(this->actual_machine_extent);
   
    // Number of parts in each Group (i.e. RCA's X coord == Grp g)
    group_count = new int[actual_machine_extent[0]];

    // Transformed dims = 1 + N_y + N_z
    transformed_networkDim = 1 + actual_machine_extent[1] + 
      actual_machine_extent[2];
    actual_machine_extent = new int[transformed_networkDim];

    // Allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];
        
    for (int i = 0; i < actual_networkDim; ++i) {
      actual_procCoords[i] = new pcoord_t[this->numRanks];
      memset(actual_procCoords[i], 0, 
          sizeof(pcoord_t) * this->numRanks);
    }

    pcoord_t *xyz = new pcoord_t[transformed_networkDim];
    getMyActualMachineCoordinate(xyz);
    for (int i = 0; i < actual_networkDim; ++i)
      actual_procCoords[i][this->myRank] = xyz[i];
    delete [] xyz;

    // reduceAll the coordinates of each processor.
    gatherMachineCoordinates(this->actual_procCoords,
        this->actual_networkDim, comm); 
  }

  // No necessary wrap arounds for dragonfly networks. Groups
  // have wrap around, but group all-to-all connection makes unneccessary.
  virtual bool getMachineExtentWrapArounds(bool *wrap_around) const {
    return false;
  }


 /*! \brief Constructor: Dragonfly (e.g. Cori & Trinity) network 
   *  machine description;
   *
   *  Does coord transformation if parameter list has a "Machine 
   *  Optimization Level" parameter set. 
   *  
   *  \param comm Communication object.
   *  \param pl   Parameter List
   */

  MachineDragonflyRCA(const Teuchos::Comm<int> &comm, 
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
     
    // Number of parts in each Group (i.e. RCA's X coord == Grp g)
    group_count = new int[actual_machine_extent[0]];
    
    // Allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];
    
    pcoord_t *xyz = new pcoord_t[actual_networkDim];
    getMyActualMachineCoordinate(xyz);

    const Teuchos::ParameterEntry *pe2 = 
      this->pl->getEntryPtr("Machine_Optimization_Level");

    // Transform with mach opt level
    if (pe2) {
      int optimization_level;
      optimization_level = pe2->getValue<int>(&optimization_level);

      if (optimization_level > 0) {
        is_transformed = true;

        if (this->myRank == 0) 
          std::cout << "Transforming the coordinates" << std::endl;
       
        // Transformed dims = 1 + N_y + N_z
        transformed_networkDim = 1 + actual_machine_extent[1] + 
          actual_machine_extent[2];
        transformed_machine_extent = new int[transformed_networkDim];

        transformed_procCoords = new pcoord_t *[transformed_networkDim];

        // Allocate memory for transformed coordinates
        for (int i = 0; i < transformed_networkDim; ++i) {
          transformed_procCoords[i] = new pcoord_t[this->numRanks];
          memset(transformed_procCoords[i], 0,
                 sizeof(pcoord_t) * this->numRanks);
        }
        
        // Calculate transformed coordinates and machine extents
        int nx = this->actual_machine_extent[0];
        int ny = this->actual_machine_extent[1];
        int nz = this->actual_machine_extent[2];
        
        transformed_procCoords[0][this->myRank] = 2 * xyz[0] * ny * nz;

        for (int i = 1; i < 1 + ny; ++i) {
          // Shift y-coord given a group, xyz[0];
          transformed_procCoords[i][this->myRank] = 0;
          // Increment in the dim where y-coord present  
          if (xyz[1] == i - 1)
            transformed_procCoords[i][this->myRank] = 2;
        }
        for (int i = 1 + ny; i < transformed_networkDim; ++i) {
          // Shift z-coord given a group, xyz[0];
          transformed_procCoords[i][this->myRank] = 0;
          // Increment in the dim where z-coord present
          if (xyz[2] == i - (1 + ny))
            transformed_procCoords[i][this->myRank] = 1;
        }

        this->transformed_machine_extent = new int[transformed_networkDim];
        
        // Max shifted high dim coordinate system
        this->transformed_machine_extent[0] = 2 * (nx - 1) * ny * nz;
        for (int i = 1; i < 1 + ny; ++i) {
          this->transformed_machine_extent[i] = 2;
        }
        for (int i = 1 + ny; i < transformed_networkDim; ++i) {
          this->transformed_machine_extent[i] = 1;
        }

        // reduceAll the transformed coordinates of each processor.
        gatherMachineCoordinates(this->transformed_procCoords, 
          this->transformed_networkDim, comm);    

        this->printAllocation();
      }
    }
    // If no coordinate transformation, gather actual coords
    if (!is_transformed) {

      for (int i = 0; i < actual_networkDim; ++i) {
        actual_procCoords[i] = new pcoord_t[this->numRanks];
        memset(actual_procCoords[i], 0, 
               sizeof(pcoord_t) * this->numRanks);
      }

      for (int i = 0; i < actual_networkDim; ++i)
        actual_procCoords[i][this->myRank] = xyz[i];

      // reduceAll the actual coordinates of each processor
      gatherMachineCoordinates(this->actual_procCoords,
          this->actual_networkDim, comm);

      this->printAllocation();
    }
    delete [] xyz;
  }

  virtual ~MachineDragonflyRCA() {
    if (is_transformed) {
      is_transformed = false;
      for (int i = 0; i < transformed_networkDim; ++i) {
        delete [] transformed_procCoords[i];
      }
      delete [] transformed_machine_extent; 
    }
    else {
      for (int i = 0; i < actual_networkDim; ++i) {
        delete [] actual_procCoords[i];
      }
    }
    delete [] actual_procCoords;
    delete [] actual_machine_extent;
    delete [] transformed_procCoords; 
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const {
    if (is_transformed) 
      return this->transformed_networkDim;
    else
      return this->actual_networkDim;
  }

  bool getTransformedMachineExtent(int *nxyz) const {
    if (is_transformed) {
      for (int dim = 0; dim < transformed_networkDim; ++dim)
        nxyz[dim] = this->transformed_machine_extent[dim];

      return true;
    }
    else 
      return false;
  }

  // Return RCA machine extents
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

  bool getMachineExtent(int *nxyz) const {
    if (is_transformed) 
      this->getTransformedMachineExtent(nxyz);
    else
      this->getActualMachineExtent(nxyz);
    
    return true;
  }

  int getGroupCount(int *grp_count) const override {
    
    grp_count = group_count;

    return actual_machine_extent[0];
  }

  void printAllocation() {
    if (this->myRank == 0) {
      // Print transformed coordinates and extents
      if (is_transformed) {
        for (int i = 0; i < this->numRanks; ++i) { 
          std::cout << "Rank:" << i;
            for (int j = 0; j < this->transformed_networkDim; ++j) {
              std::cout << " " << transformed_procCoords[j][i]; 
            }
            std::cout << std::endl;  
        } 

        std::cout << std::endl << "Transformed Machine Extent: ";
        for (int i = 0; i < this->transformed_networkDim; ++i) {
          std::cout << " " << this->transformed_machine_extent[i];
        }
        std::cout << std::endl;
      }
      // Print actual coordinates and extents
      else {
        for (int i = 0; i < this->numRanks; ++i) { 
          std::cout << "Rank:" << i;
            for (int j = 0; j < this->actual_networkDim; ++j) {
              std::cout << " " << actual_procCoords[j][i]; 
            }
            std::cout << std::endl;  
        } 

        std::cout << std::endl << "Actual Machine Extent: ";
        for (int i = 0; i < this->actual_networkDim; ++i) {
          std::cout << " " << this->actual_machine_extent[i];
        }
        std::cout << std::endl;
      }
    }
  }

  bool getMyTransformedMachineCoordinate(pcoord_t *xyz) {
    if (is_transformed) {
      for (int i = 0; i < this->transformed_networkDim; ++i) {
        xyz[i] = transformed_procCoords[i][this->myRank];
      }

      return true;
    }
    else
      return false;
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
    if (returnval == -1) {
      return false;
    }
    xyz[0] = node_coord.mesh_x;
    xyz[1] = node_coord.mesh_y;
    xyz[2] = node_coord.mesh_z;

    group_count[(int)xyz[0]]++;

    return true;
#else
    return false;
#endif
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {
    if (is_transformed) 
      this->getMyTransformedMachineCoordinate(xyz);
    else
      this->getMyActualMachineCoordinate(xyz);
  
    return true;
  }

  inline bool getMachineCoordinate(const int rank,
                                   pcoord_t *xyz) const {
    if (is_transformed) {
      for (int i = 0; i < this->transformed_networkDim; ++i) {
        xyz[i] = transformed_procCoords[i][rank];
      }
    }
    else {
      for (int i = 0; i < this->actual_networkDim; ++i) {
        xyz[i] = actual_procCoords[i][rank];
      }
    }

    return true;   
  }

  bool getMachineCoordinate(const char *nodename, pcoord_t *xyz) {
    return false;  // cannot yet return from nodename
  }
 
  bool getAllMachineCoordinatesView(pcoord_t **&allCoords) const {
    if (is_transformed) {
      allCoords = transformed_procCoords; 
    }
    else {
      allCoords = actual_procCoords;
    }

    return true;
  }

  // Return (approx) hop count from rank1 to rank2. Does not account for 
  // dynamic routing.
  bool getHopCount(int rank1, int rank2, pcoord_t &hops) override {
    
    hops = 0;

    if (is_transformed) {     
      // Case: ranks in different groups
      // Does not account for location of group to group connection. 
      // (Most group to group messages will take 5 hops)
      if (transformed_procCoords[0][rank1] != 
          transformed_procCoords[0][rank2]) 
      {
        hops = 5;
        return true;
      }

      // Case: ranks in same group
      // For each 2 differences in transformed_coordinates then
      // 1 hop
      for (int i = 1; i < transformed_networkDim; ++i) {
        if (transformed_procCoords[i][rank1] != 
            transformed_procCoords[i][rank2]) 
          ++hops;
      }
      hops /= 2;
    }
    else {

      // Case: ranks in different groups
      // Does not account for location of group to group connection. 
      // (Nearly all group to group messages will take 5 hops)
      if (actual_procCoords[0][rank1] != 
          actual_procCoords[0][rank2]) 
      {
        hops = 5;
        return true;
      }
      
      // Case: ranks in same group
      // For each difference in actual_coordinates then
      // 1 hop
      for (int i = 1; i < actual_networkDim; ++i) {
        if (actual_procCoords[i][rank1] != 
            actual_procCoords[i][rank2]) 
          ++hops;
      } 
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
  part_t *group_count;
  bool is_transformed;

  const Teuchos::ParameterList *pl;
//  bool delete_tranformed_coords;

/*
  bool delete_transformed_coords;
  int transformed_network_dim;
  pcoord_t **transformed_coordinates;
*/

  void gatherMachineCoordinates(pcoord_t **&coords, int netDim, 
      const Teuchos::Comm<int> &comm) {
    // Reduces and stores all machine coordinates.
    pcoord_t *tmpVect = new pcoord_t [this->numRanks];

    for (int i = 0; i < netDim; ++i) {
      Teuchos::reduceAll<int, pcoord_t>(comm, Teuchos::REDUCE_SUM,
                                        this->numRanks, 
                                        coords[i], tmpVect);
      pcoord_t *tmp = tmpVect;
      tmpVect = coords[i];
      coords[i] = tmp;
    }
    delete [] tmpVect;
  }

};

} // namespace Zoltan2

#endif
