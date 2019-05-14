#ifndef _ZOLTAN2_MACHINE_DRAGONFLY_RCALIBTEST_HPP_
#define _ZOLTAN2_MACHINE_DRAGONFLY_RCALIBTEST_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#include <cstdlib>    /* srand, rand */
#include <fstream>
#include <string>

namespace Zoltan2{

/*! \brief A Dragonfly (e.g. Cori & Trinity) Machine Class for 
 *  testing only. A more realistic machine should be used for 
 *  task mapping.
 *
 *  Does NOT require RCA library to run. 
 */

template <typename pcoord_t, typename part_t>
class MachineDragonflyRCAForTesting : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: Dragonfly (e.g. Cori & Trinity) RCA network 
   *  machine description;
   *
   *  Does not do coord transformation.
   *  
   *  \param comm Communication object.
   */

  MachineDragonflyRCAForTesting(const Teuchos::Comm<int> &comm):
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

  MachineDragonflyRCAForTesting(const Teuchos::Comm<int> &comm, 
             const Teuchos::ParameterList &pl_ ):
    Machine<pcoord_t,part_t>(comm),
    transformed_networkDim(3), 
    actual_networkDim(3),
    transformed_procCoords(NULL), 
    actual_procCoords(NULL),
    transformed_machine_extent(NULL),
    actual_machine_extent(NULL),
    is_transformed(false), 
    pl(&pl_) {

    actual_machine_extent = new int[actual_networkDim];
    this->getActualMachineExtent(this->actual_machine_extent);
    
    // Number of parts in each Group (i.e. RCA's X coord == Grp g)
    group_count = new int[actual_machine_extent[0]];

    memset(group_count, 0, sizeof(int) * actual_machine_extent[0]);
    
    std::cout << "\nMAKE DRAGONFLYRCAFORTESTING\n";

    // Allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];
    
    pcoord_t *xyz = new pcoord_t[actual_networkDim];
    getMyActualMachineCoordinate(xyz);

    for (int i = 0; i < actual_machine_extent[0]; ++i)
      std::cout << "\nRank: " << this->myRank << " 1group_count[" << i << "]: " << group_count[i] << "\n";

//    Teuchos::reduce<int, int>(group_count, group_count, 
//                              this->numRanks, 
//                              Teuchos::REDUCE_SUM,
//                              0, comm);

    std::cout << "\nRank: " << this->myRank << " Comm size: " << comm.getSize() << "\n";


    int *tmp_vec = new int[actual_machine_extent[0]];

    memset(tmp_vec, 0, sizeof(int) * actual_machine_extent[0]);

//    comm.barrier();

    Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_SUM,
                                        actual_machine_extent[0], 
                                        group_count, tmp_vec);

//    comm.barrier();

    group_count = tmp_vec;

    for (int i = 0; i < actual_machine_extent[0]; ++i)
      std::cout << "\nRank: " << this->myRank << " 2group_count[" << i << "]: " << group_count[i] << "\n";

    for (int i = 0; i < actual_machine_extent[0]; ++i)
      std::cout << "\nRank: " << this->myRank << " 2tmp_vec[" << i << "]: " << tmp_vec[i] << "\n";
    
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

        const Teuchos::ParameterEntry *pe_x =
          this->pl->getEntryPtr("Machine_X_Stretch"); 
        const Teuchos::ParameterEntry *pe_y =
          this->pl->getEntryPtr("Machine_Y_Stretch"); 
        const Teuchos::ParameterEntry *pe_z =
          this->pl->getEntryPtr("Machine_Z_Stretch"); 

        // Default X,Y,Z stretches
        int x_stretch = 3;
        int y_stretch = 2;
        int z_stretch = 1;

        if (pe_x)  
          x_stretch = pe_x->getValue<int>(&x_stretch);
        if (pe_y)  
          y_stretch = pe_y->getValue<int>(&y_stretch);
        if (pe_x)  
          z_stretch = pe_z->getValue<int>(&z_stretch);
         
        transformed_procCoords[0][this->myRank] = x_stretch * xyz[0] * ny * nz;
        
        for (int i = 1; i < 1 + ny; ++i) {
          // Shift y-coord given a group, xyz[0];
          transformed_procCoords[i][this->myRank] = 0;
          // Increment in the dim where y-coord present  
          if (xyz[1] == i - 1) {
            transformed_procCoords[i][this->myRank] = y_stretch; 
          }
        }
        for (int i = 1 + ny; i < transformed_networkDim; ++i) {
          // Shift z-coord given a group, xyz[0];
          transformed_procCoords[i][this->myRank] = 0;
          // Increment in the dim where z-coord present
          if (xyz[2] == i - (1 + ny))
            transformed_procCoords[i][this->myRank] = z_stretch;
        }

        this->transformed_machine_extent = new int[transformed_networkDim];
        
        // Max shifted high dim coordinate system
        this->transformed_machine_extent[0] = x_stretch * (nx - 1) * ny * nz;
        for (int i = 1; i < 1 + ny; ++i) {
          this->transformed_machine_extent[i] = y_stretch; 
        }
        for (int i = 1 + ny; i < transformed_networkDim; ++i) {
          this->transformed_machine_extent[i] = z_stretch;
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

  virtual ~MachineDragonflyRCAForTesting() {
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
/*
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
*/
    
    nxyz[0] = 11; // X - group
    nxyz[1] = 6;  // Y - row within group
    nxyz[2] = 16; // Z - col within group
  
    // Needed for test/unit_test/Machine.cpp PASS
//    nxyz[0] = 4;
//    nxyz[1] = 8;
//    nxyz[2] = 12;

    return true; 
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
          std::cout << "Rank:" << i << "  ";
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
/*
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
    return true;
#else
    return false;
#endif
*/
    srand(this->myRank);

    int x = rand() % 11;
    int y = rand() % 6;
    int z = rand() % 16;

    xyz[0] = x; 
    xyz[1] = y; 
    xyz[2] = z; 

      // Needed for test/unit_test/Machine.cpp PASS
//    xyz[0] = this->myRank;
//    xyz[1] = this->numRanks;
//    xyz[2] = this->numRanks + 1;

    std::cout << "\nGROUP: " << xyz[0] << "\n"; 

    group_count[x]++;
    
    return true; 
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
  virtual bool getHopCount(int rank1, int rank2, pcoord_t &hops) {
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
  int *group_count;
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
