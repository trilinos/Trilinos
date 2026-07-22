// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MACHINE_DRAGONFLY_RCALIBTEST_HPP_
#define _ZOLTAN2_MACHINE_DRAGONFLY_RCALIBTEST_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#include <cstdlib>    /* srand, rand */
#include <fstream>
#include <string>

namespace Zoltan2{

/*! \brief A Dragonfly (e.g. Cori, Trinity, Theta) Machine Class for 
 *  testing only. A more realistic machine should be used for 
 *  task mapping (i.e. see Zoltan2_MachineDragonflyRCA.hpp).
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
    num_unique_groups(0),
    group_count(NULL),
    is_transformed(false),
    pl(NULL) {
   
    actual_machine_extent = new int[actual_networkDim];
    this->getActualMachineExtent(this->actual_machine_extent);
    
    // Number of ranks in each Dragonfly network group 
    // (i.e. RCA's X coord == Grp g)
    group_count = new part_t[actual_machine_extent[0]];

    memset(group_count, 0, sizeof(part_t) * actual_machine_extent[0]);
    
    // Transformed dims = 1 + N_y + N_z
    transformed_networkDim = 1 + actual_machine_extent[1] + 
      actual_machine_extent[2];
    transformed_machine_extent = new int[transformed_networkDim];

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

    // Gather number of ranks in each Dragonfly network group 
    // from across all ranks
    part_t *tmp_vec = new part_t[actual_machine_extent[0]];
    memset(tmp_vec, 0, sizeof(part_t) * actual_machine_extent[0]);

    Teuchos::reduceAll<int, part_t>(comm, Teuchos::REDUCE_SUM,
                                    actual_machine_extent[0], 
                                    group_count, 
                                    tmp_vec);

    // remove zero entries from reduced array
    num_unique_groups = 0;

    for (int i = 0; i < actual_machine_extent[0]; ++i) {
      if (tmp_vec[i] > 0) {
        ++num_unique_groups;
      }
    }

    // Reset group_count array to new size
    delete[] group_count;
    group_count = new part_t[num_unique_groups];
   
    int pos = 0;
    for (int i = 0; i < actual_machine_extent[0]; ++i) {
      if (tmp_vec[i] > 0) {
        group_count[pos] = tmp_vec[i];
        ++pos;
      }
    } 

    delete[] tmp_vec;

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
   *  Optimization Level > 0" parameter set. 
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
    num_unique_groups(0),
    group_count(NULL), 
    is_transformed(false),
    pl(&pl_) { 

    actual_machine_extent = new int[actual_networkDim];
    this->getActualMachineExtent(this->actual_machine_extent);
    
    // Number of parts in each Dragonfly network group 
    // (i.e. RCA's X coord == Grp g)
    group_count = new part_t[actual_machine_extent[0]];

    memset(group_count, 0, sizeof(part_t) * actual_machine_extent[0]);
    
    // Allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];
    
    pcoord_t *xyz = new pcoord_t[actual_networkDim];
    getMyActualMachineCoordinate(xyz);

    // Gather number of ranks in each Dragonfly network group 
    // from across all ranks
    part_t *tmp_vec = new part_t[actual_machine_extent[0]];
    memset(tmp_vec, 0, sizeof(part_t) * actual_machine_extent[0]);

    Teuchos::reduceAll<int, part_t>(comm, Teuchos::REDUCE_SUM,
                                    actual_machine_extent[0], 
                                    group_count, 
                                    tmp_vec);

    // remove zero entries from reduced vector
    num_unique_groups = 0;

    for (int i = 0; i < actual_machine_extent[0]; ++i) {
      if (tmp_vec[i] > 0) {
        ++num_unique_groups;
      }
    }

    // Reset group_count array to new size (# of nonzero groups)
    delete[] group_count;
    group_count = new part_t[num_unique_groups];
   
    int pos = 0;
    for (int i = 0; i < actual_machine_extent[0]; ++i)
    {
      if (tmp_vec[i] > 0) {
        group_count[pos] = tmp_vec[i];
        ++pos;
      }
    } 
    delete[] tmp_vec;
     
    const Teuchos::ParameterEntry *pe2 = 
      this->pl->getEntryPtr("Machine_Optimization_Level");

    // Transform with mach opt level
    if (pe2) {
      int optimization_level;
      optimization_level = pe2->getValue<int>(&optimization_level);

      if (optimization_level > 0) {
        is_transformed = true;
       
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
        
        // Transform X coords 
        transformed_procCoords[0][this->myRank] = 
          x_stretch * xyz[0] * ny * nz;
        
        // Transform Y coords
        for (int i = 1; i < 1 + ny; ++i) {
          // Shift y-coord given a group, xyz[0];
          transformed_procCoords[i][this->myRank] = 0;
          // Increment in the dim where y-coord present  
          if (xyz[1] == i - 1) {
            transformed_procCoords[i][this->myRank] = y_stretch; 
          }
        }
        // Transform Z coords
        for (int i = 1 + ny; i < transformed_networkDim; ++i) {
          // Shift z-coord given a group, xyz[0];
          transformed_procCoords[i][this->myRank] = 0;
          // Increment in the dim where z-coord present
          if (xyz[2] == i - (1 + ny))
            transformed_procCoords[i][this->myRank] = z_stretch;
        }

        this->transformed_machine_extent = new int[transformed_networkDim];
        
        // Maximum extents in shifted high dim coordinate system
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

  // Destructor
  virtual ~MachineDragonflyRCAForTesting() {
    if (is_transformed) {
      is_transformed = false;
      if (this->numRanks > 1) {
        for (int i = 0; i < transformed_networkDim; ++i) {
          delete [] transformed_procCoords[i];
        }
      } 
      delete [] transformed_machine_extent;
    }
    else {
      if (this->numRanks > 1) {
        for (int i = 0; i < actual_networkDim; ++i) {
          delete [] actual_procCoords[i];
        }
      } 
    }
 
    delete [] actual_procCoords;
    delete [] transformed_procCoords;  
   
    delete [] actual_machine_extent;
    delete [] group_count;
  }

  bool hasMachineCoordinates() const { return true; }

  // Return dimensions of coords, transformed or actual
  int getMachineDim() const {
    if (is_transformed) 
      return this->transformed_networkDim;
    else
      return this->actual_networkDim;
  }

  // Return the transformed maximum machine extents
  bool getTransformedMachineExtent(int *nxyz) const {
    if (is_transformed) {
      for (int dim = 0; dim < transformed_networkDim; ++dim)
        nxyz[dim] = this->transformed_machine_extent[dim];

      return true;
    }
    else 
      return false;
  }

  // Return the fake "RCA" machine extents for testing
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

  // Return machine extents, transformed or actual
  bool getMachineExtent(int *nxyz) const {
    if (is_transformed) 
      this->getTransformedMachineExtent(nxyz);
    else
      this->getActualMachineExtent(nxyz);
    
    return true;
  }

  // Return number of groups (RCA X-dim) with allocated nodes
  part_t getNumUniqueGroups() const override{
    return this->num_unique_groups;
  }

  // Return number of ranks in each group (RCA X-dim) in an allocation
  bool getGroupCount(part_t *grp_count) const override {
        
    if (group_count != NULL) {
      for (int i = 0; i < num_unique_groups; ++i) { 
        grp_count[i] = this->group_count[i];
      }

      return true;
    }
    else 
      return false;
  }

  // Print allocation coords and extents on rank 0, transformed or actual
  void printAllocation() {
    if (this->myRank >= 0) {
      // Print transformed coordinates and extents
      if (is_transformed) {
        for (int i = 0; i < this->numRanks; ++i) { 
          std::cout << "Rank:" << i << "  ";
            for (int j = 0; j < this->transformed_networkDim; ++j) {
              std::cout << " " << this->transformed_procCoords[j][i]; 
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

  // Return transformed coord for this rank
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

  // Return the fake "RCA" coord for this rank for testing
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

    group_count[x]++;
    
    return true; 
  }

  // Return machine coordinate for this rank, transformed or actual
  bool getMyMachineCoordinate(pcoord_t *xyz) {
    if (is_transformed) 
      this->getMyTransformedMachineCoordinate(xyz);
    else
      this->getMyActualMachineCoordinate(xyz);
  
    return true;
  }

  // Return machine coord of given rank, transformed or actual
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

 // Return view of all machine coords, transformed or actual 
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
  // Dragonfly's dynamic routing.
  virtual bool getHopCount(int rank1, int rank2, pcoord_t &hops) const override {
    hops = 0;
      
    if (rank1 == rank2)
      return true;
    if (rank1 >= this->numRanks || rank2 >= this->numRanks) {
      std::cerr << "Rank outside bounds for the machine ranks";
      exit(1);
    }

    if (this->is_transformed) {    
      // Case: ranks in different groups (i.e. different RCA x-coords)
      // Does not account for location of group to group connection. 
      // (Most group to group messages will take 5 hops)
      if (this->transformed_procCoords[0][rank1] != 
          this->transformed_procCoords[0][rank2]) 
      {
        hops = 5;
    
        return true;
      }

      // Case: ranks in same group
      // For each 2 differences in transformed_coordinates then
      // 1 hop
      for (int i = 1; i < this->transformed_networkDim; ++i) {
        if (this->transformed_procCoords[i][rank1] != 
            this->transformed_procCoords[i][rank2]) 
          ++hops;
      }
      hops /= 2;    
    }
    else {
      // Case: ranks in different groups
      // Does not account for location of group to group connection. 
      // (Nearly all group to group messages will take 5 hops)
      if (this->actual_procCoords[0][rank1] != 
          this->actual_procCoords[0][rank2]) 
      {
        hops = 5;
        return true;
      }
      
      // Case: ranks in same group
      // For each difference in actual_coordinates then
      // 1 hop
      for (int i = 1; i < this->actual_networkDim; ++i) {
        if (this->actual_procCoords[i][rank1] != 
            this->actual_procCoords[i][rank2]) 
          ++hops;
      }
    }

    return true;
  }

private:

  // # of dimensions in the stored coordinates, transformed or actual
  int transformed_networkDim;
  int actual_networkDim;

  // Machine Coordinates
  pcoord_t **transformed_procCoords;
  pcoord_t **actual_procCoords;

  // Maximum extents for each dimension, transformed or actual
  part_t *transformed_machine_extent;
  part_t *actual_machine_extent;

  // Number of groups (RCA X-dim) with nonzero nodes allocated
  part_t num_unique_groups; 
  // Distribution of nodes in each group (zero node groups have been trimmed)
  part_t *group_count;

  // Are out coordinates transformed?
  bool is_transformed;

  const Teuchos::ParameterList *pl;


  // reduceAll the machine coordinates
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
