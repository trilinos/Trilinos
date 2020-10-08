#ifndef _ZOLTAN2_MACHINE_FATTREE_LIBTEST_HPP_
#define _ZOLTAN2_MACHINE_FATTREE_LIBTEST_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#include <fstream>
#include <string>
#include <unistd.h>

namespace Zoltan2{

/*! \brief A FatTree (e.g. Astra, Summit, & Sierra) Machine Class for
 *  production.
 */

template <typename pcoord_t, typename part_t>
class MachineFatTree : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: FatTree (e.g. Astra, Summit, & Sierra) network
   *  machine description;
   *
   *  Does not do coord transformation.
   *
   *  \param comm Communication object.
   */

  MachineFatTree(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
    transformed_networkDim(3),
    actual_networkDim(3),
    nborhoods_per_row(2),
    transformed_procCoords(NULL),
    actual_procCoords(NULL),
    transformed_machine_extent(3),
    actual_machine_extent(3),
    num_unique_groups(1),
    group_count(),
    num_unique_subgroups(),
    subgroup_counts(),
    is_transformed(false),
    pl(NULL) {


    if (this->myRank == 0)
        std::cout << "\nCreating FatTree Machine!. Entered no PL FatTree\n" << std::endl;

//    actual_machine_extent = new int[actual_networkDim];
    getActualMachineExtent(actual_machine_extent);

//    transformed_machine_extent = new int[transformed_networkDim];

    transformed_machine_extent =
        { nborhoods_per_row * actual_machine_extent[0],
          actual_machine_extent[1] / nborhoods_per_row,
          actual_machine_extent[2]};



    // Number of parts in each FatTree network neighborhood
    // a neighborhood is defined as links connected by the L1 switch
    // (below director level switches)
    // (i.e. Row * nborhoods_per_row +
    // ( Col / (num_cols / nborhoods_per_row))  == Grp g)


    group_count.resize(actual_machine_extent[0]);
    subgroup_counts.resize(actual_machine_extent[0]);

    for (size_t i = 0; i < subgroup_counts.size(); ++i) {
      subgroup_counts[i].resize(actual_machine_extent[1]);
    }

    // Allocate memory for processor coords, Use arrays for speed
    actual_procCoords = new pcoord_t *[actual_networkDim];
    transformed_procCoords = new pcoord_t *[transformed_networkDim];

    for (int i = 0; i < actual_networkDim; ++i) {
      actual_procCoords[i] = new pcoord_t[this->numRanks];
      memset(actual_procCoords[i], 0,
          sizeof(pcoord_t) * this->numRanks);
    }

//    pcoord_t *xyz = new pcoord_t[actual_networkDim];

    std::vector<pcoord_t> xyz(actual_networkDim);

    if (this->myRank == 0)
      std::cout << "About to get xyz" << std::endl;

    getMyActualMachineCoordinate(xyz);

    if (this->myRank == 0)
      std::cout << "Got xyz" << std::endl;

    for (size_t i = 0; i < xyz.size(); ++i)
      actual_procCoords[i][this->myRank] = xyz[i];
//    delete [] xyz;


    group_count[int(xyz[0])] = 1;
    subgroup_counts[int(xyz[0])][int(xyz[1])] = 1;

    calc_groups(comm);

    // reduceAll the coordinates of each processor.
    gatherMachineCoordinates(actual_procCoords,
                             actual_networkDim, comm);

    printAllocation();
  }

 /*! \brief Constructor: FatTree (e.g. Cori & Trinity) network
   *  machine description;
   *
   *  Does coord transformation if parameter list has a "Machine
   *  Optimization Level > 0" parameter set.
   *
   *  \param comm Communication object.
   *  \param pl   Parameter List
   */
  MachineFatTree(const Teuchos::Comm<int> &comm,
             const Teuchos::ParameterList &pl_ ):
    Machine<pcoord_t,part_t>(comm),
    transformed_networkDim(3),
    actual_networkDim(3),
    nborhoods_per_row(2),
    transformed_procCoords(NULL),
    actual_procCoords(NULL),
    transformed_machine_extent(3),
    actual_machine_extent(3),
    num_unique_groups(1),
    group_count(),
    num_unique_subgroups(),
    subgroup_counts(),
    is_transformed(false),
    pl(&pl_) {

    if (this->myRank == 0)
      std::cout << "\nCreating FatTree Machine!. Entered  PL FatTree for transformation\n" << std::endl;

//    actual_machine_extent = new int[actual_networkDim];
    getActualMachineExtent(actual_machine_extent);

    transformed_machine_extent =
      {nborhoods_per_row * actual_machine_extent[0],
       actual_machine_extent[1] / nborhoods_per_row,
       actual_machine_extent[2]};

    // Allocate memory for processor coords
    actual_procCoords = new pcoord_t *[actual_networkDim];

    // Get my real network coordinate
    std::vector<pcoord_t> xyz(actual_networkDim);

    if (this->myRank == 0)
      std::cout << "About to get xyz" << std::endl;

    getMyActualMachineCoordinate(xyz);

    if (this->myRank == 0)
      std::cout << "Got xyz" << std::endl;


    if(xyz[0] > actual_machine_extent[0] ||
       xyz[1] > actual_machine_extent[1] ||
       xyz[2] > actual_machine_extent[2]) {

      std::cout << "Rank: " << this->myRank
        << " Coord: (" << xyz[0]
        << ", " << xyz[1] << ", " << xyz[2]
        << ") is beyond extents ("
        << actual_machine_extent[0] << ", "
        << actual_machine_extent[0] << ", "
        << actual_machine_extent[0] << ")! Exiting. "
        << std::endl;

        exit(0);
    }

    const Teuchos::ParameterEntry *pe2 =
      this->pl->getEntryPtr("Machine_Optimization_Level");

    if (pe2 || 1) {

      int optimization_level = pe2->getValue<int>(&optimization_level);

      if (optimization_level > 0 || 1) {
        is_transformed = true;

        if (this->myRank == 0)
          std::cout << "\nOptimizing!" << std::endl;

        transformed_procCoords = new pcoord_t *[transformed_networkDim];

        // Allocate memory for transformed coordinates
        for (int i = 0; i < transformed_networkDim; ++i) {
          transformed_procCoords[i] = new pcoord_t[this->numRanks];
          memset(transformed_procCoords[i], 0,
                 sizeof(pcoord_t) * this->numRanks);
        }
        // Number of parts in each FatTree network neighborhood
        // a neighborhood is defined as links connected by the L1 switch
        // (below director level switches)
        // (i.e. Row_idx * nborhoods_per_row +
        // ( Col_idx / (num_cols / nborhoods_per_row))  == Grp g)
        group_count.resize(transformed_machine_extent[0]);
        subgroup_counts.resize(transformed_machine_extent[0]);

        for (int i = 0; i < transformed_machine_extent[0]; ++i) {
          subgroup_counts[i].resize(transformed_machine_extent[1]);
        }


        int transformed_group =
            nborhoods_per_row * xyz[0] +
            int(xyz[1] / transformed_machine_extent[1]);
        int transformed_subgroup =
            int(xyz[1]) % transformed_machine_extent[1];
        int transformed_node = xyz[2];

        // A1-A18 is one neighborhood, A19-A36 is another neighborhood
        transformed_procCoords[0][this->myRank] = transformed_group;
        // Rack within neighborhood
        transformed_procCoords[1][this->myRank] = transformed_subgroup;
        // Node within rack
        transformed_procCoords[2][this->myRank] = transformed_node;

        // Group Count calculations

        group_count[transformed_group] = 1;
        subgroup_counts[transformed_group][transformed_subgroup] = 1;

        calc_groups(comm);

        // reduceAll the transformed coordinates of each processor.
        gatherMachineCoordinates(transformed_procCoords,
                                 transformed_networkDim, comm);

        printAllocation();
      }
    }


    // If no coordinate transformation, gather actual coords
    if (!is_transformed) {

      if (this->myRank == 0)
        std::cout << "\nNot Transforming" << std::endl;

      for (int i = 0; i < actual_networkDim; ++i) {
        actual_procCoords[i] = new pcoord_t[this->numRanks];
        memset(actual_procCoords[i], 0,
               sizeof(pcoord_t) * this->numRanks);
      }

      for (int i = 0; i < actual_networkDim; ++i)
        actual_procCoords[i][this->myRank] = xyz[i];

      // Number of parts in each FatTree network neighborhood
      // a neighborhood is defined as links connected by the L1 switch
      // (below director level switches)
      // (i.e. Row * nborhoods_per_row +
      // ( Col / (num_cols / nborhoods_per_row))  == Grp g)
      group_count.resize(actual_machine_extent[0]);
      subgroup_counts.resize(actual_machine_extent[0]);

      for (int i = 0; i < actual_machine_extent[0]; ++i) {
        subgroup_counts[i].resize(actual_machine_extent[1]);
      }

      group_count[int(xyz[0])] = 1;
      subgroup_counts[int(xyz[0])][int(xyz[1])] = 1;

      calc_groups(comm);

      // reduceAll the actual coordinates of each processor
      gatherMachineCoordinates(actual_procCoords,
                               actual_networkDim, comm);

      printAllocation();
    }

  }

  // Destructor
  virtual ~MachineFatTree() {

    if (is_transformed) {
      is_transformed = false;
      if (this->numRanks > 1) {
        for (int i = 0; i < transformed_networkDim; ++i) {
          delete [] transformed_procCoords[i];
        }
      }
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
  }

  bool hasMachineCoordinates() const { return true; }

  // Return dimensions of coords, transformed or actual
  int getMachineDim() const {
    if (is_transformed)
      return transformed_networkDim;
    else
      return actual_networkDim;
  }

  // Return the transformed maximum machine extents
  bool getTransformedMachineExtent(std::vector<int> &nxyz) const {
    if (is_transformed) {

      nxyz.assign(transformed_machine_extent.begin(),
                  transformed_machine_extent.end());

      return true;
    }
    else
      return false;
  }

  // Return the fake "RCA" machine extents for testing
  bool getActualMachineExtent(std::vector<int> &nxyz) const {

   // Actual FatTree Coordinate Extents
    nxyz[0] = 8;  // X - Row of rack on Summit floor, [A-H]
    nxyz[1] = 36; // Y - Col of rack on Summit floor, [01-36]
    nxyz[2] = 18; // Z - Node within rack,            [01-18]

    // Transformed FatTree Coordinate Extents
//    nxyz[0] = 16; // X - L1 switch neighborhood (not director level switch)
//    nxyz[1] = 18; // Y - Rack within L1 switch neighborhood
//    nxyz[2] = 18; // Z - Node within rack

    // Needed for test/unit_test/Machine.cpp PASS
//    nxyz[0] = 4;
//    nxyz[1] = 8;
//    nxyz[2] = 12;

    return true;
  }

  // Return machine extents, transformed or actual
  bool getMachineExtent(int *nxyz) const {

    std::vector<int> extent;

    if (is_transformed) {
      extent.resize(transformed_networkDim);
      getTransformedMachineExtent(extent);
    }
    else {
      extent.resize(actual_networkDim);
      getActualMachineExtent(extent);
    }

    std::copy(extent.begin(), extent.end(), nxyz);

    return true;
  }

  // Return number of groups or neighborhoods with allocated nodes
  part_t getNumUniqueGroups() const override {
    return num_unique_groups;
  }

  // Return number of ranks in each group (RCA X-dim) in an allocation
//  bool getGroupCount(part_t *grp_count) const override {
  bool getGroupCount2(std::vector<part_t> &grp_count) const override {

    if (group_count.size() > 0) {
//      std::copy(group_count.begin(), group_count.end(), grp_count);

      grp_count = group_count;

      return true;
    }
    else
      return false;
  }

  // Return the number of unique subgroups within each group
  // that contain allocated nodes.
  // Ex. num_unique_groups = 4, and group_count = [16, 8, 8, 4]
  // and let subgroup_counts =
  //      [[4, 0, 0, 0, ..., 8, 4],    // (3 unique subgroups)
  //       [2, 2, 0, 0, ..., 2, 2],    // (4 unique subgroups)
  //       [0, 0, 2, 2, ..., 2, 2],    // (4 unique subgroups)
  //       [2, 0, 0, 0, ..., 2, 0]]    // (2 unique subgroups)
  //
  // then num_unique subgroups = [3, 4, 4, 2]
//  bool getNumUniqueSubgroups(part_t *num_unique_subgrps) const override {
  bool getNumUniqueSubgroups(std::vector<part_t> &num_unique_subgrps) const override {

    if (num_unique_subgroups.size() > 0) {

//      std::copy(num_unique_subgroups.begin(), num_unique_subgroups.end(), num_unique_subgrps);
      num_unique_subgrps = num_unique_subgroups;

      return true;
    }
    else
      return false;
  }

  // Return the allocated node counts subgroup counts within each group.
  //
  // Ex. num_unique_groups = 4, and group_count = [16, 8, 8, 4]
  // and let the original subgroup_counts =
  //      [[4, 0, 0, 0, ..., 8, 4],    // (3 unique subgroups)
  //       [2, 2, 0, 0, ..., 2, 2],    // (4 unique subgroups)
  //       [0, 0, 2, 2, ..., 2, 2],    // (4 unique subgroups)
  //       [2, 0, 0, 0, ..., 2, 0]]    // (2 unique subgroups)
  //
  // then num_unique subgroups = [3, 4, 4, 2]
  // now the returned subgroup_counts =
  //    [[4, 8, 4, ...],    // (3 unique subgroups)
  //     [2, 2, 2, 2, ...],    // (4 unique subgroups)
  //     [2, 2, 2, 2, ...],    // (4 unique subgroups)
  //     [2, 2, ...]]    // (2 unique subgroups)
  //
  // then num_unique subgroups = [3, 4, 4, 2]
//  bool getSubgroupCounts(part_t **subgrp_counts) const override {
  bool getSubgroupCounts(std::vector<std::vector<part_t>> &subgrp_counts) const override {

    if (subgroup_counts.size() > 0 && subgroup_counts[0].size() > 0) {
//      for (int i = 0; i < num_unique_groups; ++i) {

//        std::copy(subgroup_counts[i].begin(),
//                  subgroup_counts[i].end(),
//                  subgrp_counts[i]);

//      }

      subgrp_counts = subgroup_counts;

      return true;
    }
    else
      return false;
  }



  // Print allocation coords and extents on rank 0, transformed or actual
  void printAllocation() {
    if (this->myRank == 0) {

      std::cout << "\nPrinting Allocation:\n" << std::endl;

      // Print transformed coordinates and extents
      if (is_transformed) {
        std::cout << "Transformed Machine Coordinates:\n";
        for (int i = 0; i < this->numRanks; ++i) {
          std::cout << "Rank: " << i << "  ";
            for (int j = 0; j < transformed_networkDim; ++j) {
              std::cout << " " << transformed_procCoords[j][i];
            }
            std::cout << std::endl;
        }

        std::cout << std::endl << "Transformed Machine Extent: ";
        for (int i = 0; i < transformed_networkDim; ++i) {
          std::cout << " " << transformed_machine_extent[i];
        }
        std::cout << std::endl << std::endl;
      }
      // Print actual coordinates and extents
      else {
        std::cout << "Actual Machine Coordinates:\n";
        for (int i = 0; i < this->numRanks; ++i) {
          std::cout << "Rank: " << i;
            for (int j = 0; j < actual_networkDim; ++j) {
              std::cout << " " << actual_procCoords[j][i];
            }
            std::cout << std::endl;
        }

        std::cout << std::endl << "Actual Machine Extent: ";
        for (int i = 0; i < actual_networkDim; ++i) {
          std::cout << " " << actual_machine_extent[i];
        }
        std::cout << std::endl << std::endl;
      }
    }
  }

  // Return transformed coord for this rank
  bool getMyTransformedMachineCoordinate(std::vector<pcoord_t> &xyz) {

    if (is_transformed) {
      for (int i = 0; i < transformed_networkDim; ++i) {
        xyz[i] = transformed_procCoords[i][this->myRank];
      }

      return true;
    }
    else
      return false;
  }

  // Return the fake "RCA" coord for this rank for testing
  bool getMyActualMachineCoordinate(std::vector<pcoord_t> &xyz) {

    char hostname[7];
    int rc = gethostname(hostname, sizeof(hostname));

//    if (rc != 0) {
//      std::cout << "\nrc: " << rc << " Error reading hostname: " << hostname << ". Done!" << std::endl;
//      exit(1);
//    }

    convertHostnameToCoordinate(hostname, xyz);

    return true;
  }

  // Return machine coordinate for this rank, transformed or actual
  bool getMyMachineCoordinate(pcoord_t *xyz) {

    std::vector<pcoord_t> coord;

    if (is_transformed) {
      coord.resize(transformed_networkDim);
      this->getMyTransformedMachineCoordinate(coord);
    }
    else {
      coord.resize(actual_networkDim);
      this->getMyActualMachineCoordinate(coord);
    }

    std::copy(coord.begin(), coord.end(), xyz);

    return true;
  }

  // Return machine coord of given rank, transformed or actual
  inline bool getMachineCoordinate(const int rank,
                                   pcoord_t *xyz) const {
    if (is_transformed) {
      for (int i = 0; i < transformed_networkDim; ++i) {
        xyz[i] = transformed_procCoords[i][rank];
      }
    }
    else {
      for (int i = 0; i < actual_networkDim; ++i) {
        xyz[i] = actual_procCoords[i][rank];
      }
    }

    return true;
  }

  bool getMachineCoordinate(const char *nodename, pcoord_t *xyz) {
    return false; // cannot yet return from nodename
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

  // No necessary wrap arounds for dragonfly networks. Groups
  // have wrap around, but group all-to-all connection makes unneccessary.
  virtual bool getMachineExtentWrapArounds(bool *wrap_around) const {
    return false;
  }

  // FatTree machine requires two levels of nonuniform partitioning
  // 1.) Group level    (nbohood switches)
  // 2.) Subgroup level (racks with switch nborhood).
  virtual int getNumNonuniformLevels() const override {
    return 2;
  }

  // Return (approx) hop count from rank1 to rank2. Does not account for
  // FatTree's dynamic routing.
  bool getHopCount(int rank1, int rank2, pcoord_t &hops) const {

    if (rank1 == rank2)
      return true;
    if (rank1 >= this->numRanks || rank2 >= this->numRanks) {
      std::cerr << "Rank " << rank1 << " " << rank2 << " are outside bounds for the machine ranks: " << this->numRanks << std::endl;
      exit(1);
    }

    if (this->is_transformed) {
      // Case: ranks in different groups (i.e. different RCA x-coords)
      // Does not account for location of group to group connection.
      // (Most group to group messages will take 5 hops)
      if (transformed_procCoords[0][rank1] !=
          transformed_procCoords[0][rank2]) {
        hops = 8;
      }
      else if (transformed_procCoords[1][rank1] !=
               transformed_procCoords[1][rank2]) {
        hops = 4;
      }
      else
        hops = 2;
    }
    else {
      // Case: ranks in different groups
      // Does not account for location of group to group connection.
      // (Nearly all group to group messages will take 5 hops)
      if (actual_procCoords[0][rank1] !=
          actual_procCoords[0][rank2]) {
        hops = 8;
      }
      // Same row but different nborhoods
      else if ((actual_procCoords[0][rank1] ==
                actual_procCoords[0][rank2]) &&
                actual_procCoords[1][rank1] /
                  (actual_machine_extent[1] / nborhoods_per_row) !=
                actual_procCoords[1][rank2] /
                  (actual_machine_extent[1] / nborhoods_per_row)) {
        hops = 8;
      }
      // Same nborhood but different rack
      else if (actual_procCoords[1][rank1] !=
               actual_procCoords[1][rank2]) {
        hops = 4;
      }
      else
        hops = 2;
    }

    return true;
  }

private:

  // # of dimensions in the stored coordinates, transformed or actual
  int transformed_networkDim;
  int actual_networkDim;

  // # of neighborhoods per row of racks
  int nborhoods_per_row;

  // Machine Coordinates
  pcoord_t **transformed_procCoords;
  pcoord_t **actual_procCoords;

  // Maximum extents for each dimension, transformed or actual
//  part_t *transformed_machine_extent;
//  part_t *actual_machine_extent;

  std::vector<int> transformed_machine_extent;
  std::vector<int> actual_machine_extent;

  // Number of groups (FatTree neighborhoods) with nonzero
  // nodes allocated
  part_t num_unique_groups;
  // Distribution of nodes in each group (zero node groups
  // have been trimmed)
//  part_t *group_count;

  std::vector<part_t> group_count;

  std::vector<part_t> num_unique_subgroups;
  std::vector<std::vector<part_t>> subgroup_counts;

  // Number of subgroup
//  part_t *num_unique_subgroups;
//  part_t **subgroup_counts;

  // Are out coordinates transformed?
  bool is_transformed;

  const Teuchos::ParameterList *pl;


  // reduceAll the machine coordinates
  void gatherMachineCoordinates(pcoord_t** coords, int netDim,
      const Teuchos::Comm<int> &comm) {
    // Reduces and stores all machine coordinates.
    pcoord_t *tmpVect = new pcoord_t [this->numRanks];

    for (int i = 0; i < netDim; ++i) {
      Teuchos::reduceAll<int, pcoord_t>(comm, Teuchos::REDUCE_SUM,
                                        this->numRanks,
                                        coords[i], tmpVect);

      // Memory Leak on tmpVect's initial memory location
      pcoord_t *tmp = tmpVect;
      tmpVect = coords[i];
      coords[i] = tmp;
    }

    delete [] tmpVect;
  }

  // Calculate group and subgroup member values.
  // Requires that group_count and subgroup_counts have been incremented for the current rank
  void calc_groups(const Teuchos::Comm<int> &comm) {

    std::vector<int> extents;

    if (is_transformed)
      extents = transformed_machine_extent;
    else
      extents = actual_machine_extent;

    // Gather number of ranks in each FatTree network group
    // from across all ranks
    std::vector<part_t> tmp_vec(group_count.size());

    Teuchos::reduceAll<int, part_t>(comm, Teuchos::REDUCE_SUM,
                                    extents[0],
                                    group_count.data(),
                                    tmp_vec.data());

    // Remote zeros from tmp group vector
    std::vector<int>::iterator nonzeros_iter =
        std::remove_if(tmp_vec.begin(),
                       tmp_vec.end(),
                       [](part_t x) { return x == 0; });

    tmp_vec.resize( nonzeros_iter -  tmp_vec.begin() );

    // remove zero entries from reduced vector
    num_unique_groups = tmp_vec.size();

    group_count = tmp_vec;

    // Gather number of ranks in each FatTree network subgroup
    // from across all ranks
    std::vector<part_t> tmp_subgrp(extents[1]);

    int row_idx = 0;

    // Trim subgroup_counts 2d array
    for (size_t i = 0; i < subgroup_counts.size(); ++i) {

      // Find global subgroup counts for the current group
      Teuchos::reduceAll<int, part_t>(comm, Teuchos::REDUCE_SUM,
                                      extents[1],
                                      subgroup_counts[i].data(),
                                      tmp_subgrp.data());

      int num_nonzeros_subgrp =
          std::count_if(tmp_subgrp.begin(),
                        tmp_subgrp.end(),
                        [](part_t x){ return x > 0; });

      // Shift nonzeros to the right and trim empty rows
      if (num_nonzeros_subgrp > 0) {
        num_unique_subgroups.push_back(num_nonzeros_subgrp);

        size_t pos_j = 0;
        for (size_t j = 0; j < tmp_subgrp.size(); ++j) {
          if (tmp_subgrp[j] != 0) {

            tmp_subgrp[pos_j] = tmp_subgrp[j];

            if (pos_j != j)
              tmp_subgrp[j] = 0;
            pos_j++;

          }
        }

        subgroup_counts[row_idx] = tmp_subgrp;

        row_idx++;
      }
    }

    subgroup_counts.resize(num_unique_groups);

/*
    if (this->myRank == 0) {
      std::cout << "\nTransformed: " << is_transformed << std::endl;
      std::cout << "Num_Uniques_Groups: " << num_unique_groups << std::endl;

      std::cout << "GroupCount: ";
      for (int i = 0; i < num_unique_groups; ++i) {
        std::cout << " " << group_count[i];
      }
      std::cout << std::endl;

      std::cout << "Num_Unique_Subgroups: ";
      for (size_t i = 0; i < num_unique_subgroups.size(); ++i) {
        std::cout << " " << num_unique_subgroups[i];
      }
      std::cout << std::endl;

      std::cout << "\nSubgroup_Counts: " << std::endl;
      for (int i = 0; i < num_unique_groups; ++i) {
        std::cout << "\nGroup " << i << ":  ";
        for (int j = 0; j < extents[1]; ++j) {
          std::cout << " " << subgroup_counts[i][j];
        }
      }
      std::cout << std::endl;
    }
*/
  }

  // Convert hostname to coordinate
  bool convertHostnameToCoordinate(const char *nodename, std::vector<pcoord_t> &xyz) {

    // [A-H]
    int x = nodename[0] - 'a';

    // [0-35]
    int y10 = nodename[1] - '0';
    int y1 = nodename[2] - '0';
    int y = y10 * 10 + y1 - 1;

    // [0-17]
    int z10 = nodename[4] - '0';
    int z1 = nodename[5] - '0';
    int z = z10 * 10 + z1 - 1;


    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return true;
  }

};

} // namespace Zoltan2

#endif
