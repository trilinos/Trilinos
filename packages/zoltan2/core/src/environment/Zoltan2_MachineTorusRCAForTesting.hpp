// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MACHINE_TORUS_RCALIBTEST_HPP_
#define _ZOLTAN2_MACHINE_TORUS_RCALIBTEST_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#include <cstdlib>     /* srand, rand */
#include <fstream>
#include <string>

namespace Zoltan2{

/*! \brief An RCA Machine Class (Torus Networks) for testing only
 *  A more realistic machine should be used for task mapping.
 */

template <typename pcoord_t, typename part_t>
class MachineTorusRCAForTesting : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  MachineTorusRCAForTesting(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
    networkDim(3), actual_networkDim(3),
    procCoords(NULL), actual_procCoords(NULL),
    machine_extent(NULL),actual_machine_extent(NULL),
    is_transformed(false), pl(NULL)
  {
    actual_machine_extent = machine_extent = new int[networkDim];
    this->getRealMachineExtent(this->machine_extent);
    actual_machine_extent = machine_extent;

    // Allocate memory for processor coordinates.
    actual_procCoords = procCoords = new pcoord_t *[networkDim];
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

  virtual bool getMachineExtentWrapArounds(bool *wrap_around) const {
    int dim = 0;
    int transformed_network_dim = networkDim;

    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    return true;
  }

  MachineTorusRCAForTesting(const Teuchos::Comm<int> &comm, 
                       const Teuchos::ParameterList &pl_):
    Machine<pcoord_t,part_t>(comm),
    networkDim(3), actual_networkDim(3),
    procCoords(NULL), actual_procCoords(NULL),
    machine_extent(NULL),actual_machine_extent(NULL),
    is_transformed(false), pl(&pl_)
  {

    actual_machine_extent = machine_extent = new int[networkDim];
    this->getRealMachineExtent(this->machine_extent);
    actual_machine_extent = machine_extent;

    // Allocate memory for processor coordinates.
    actual_procCoords = procCoords = new pcoord_t *[networkDim];


    const Teuchos::ParameterEntry *pe1 = 
      this->pl->getEntryPtr("Input_RCA_Machine_Coords");
    if (pe1) {
      std::string input_coord_file;
      input_coord_file = pe1->getValue<std::string>(&input_coord_file);
      if (input_coord_file != "") {

        if (this->myRank == 0) {
          std::vector < std::vector <pcoord_t> > proc_coords(networkDim);
          std::fstream machine_coord_file(input_coord_file.c_str());

          part_t i = 0;
          pcoord_t a,b, c;
          machine_coord_file >> a >> b >> c;
          while(!machine_coord_file.eof()) {
            proc_coords[0].push_back(a);
            proc_coords[1].push_back(b);
            proc_coords[2].push_back(c);
            ++i;
            machine_coord_file >> a >> b >> c;
          }

          machine_coord_file.close();
          std::cout << "Rewriting numprocs from:" 
            << this->numRanks << " to:" << i << std::endl;
          this->numRanks = i;

          for(int ii = 0; ii < networkDim; ++ii) {
            procCoords[ii] = new pcoord_t[this->numRanks];
            for (int j = 0; j < this->numRanks; ++j) {
              procCoords[ii][j] = proc_coords[ii][j];
            }
          }
        }
        comm.broadcast(0, sizeof(int), (char *) &(this->numRanks));

        if (this->myRank != 0) {
          for (int i = 0; i < networkDim; ++i) {
            procCoords[i] = new pcoord_t[this->numRanks];
            memset(procCoords[i], 0, sizeof(pcoord_t) * this->numRanks);
          }
        }
      }
    }
    else {
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
    }

    // reduceAll the coordinates of each processor.
    gatherMachineCoordinates(comm);

    const Teuchos::ParameterEntry *pe2 = 
      this->pl->getEntryPtr("Machine_Optimization_Level");
//    this->printAllocation();
    if (pe2) {
      int optimization_level;
      optimization_level = pe2->getValue<int>(&optimization_level);

      if (optimization_level == 1) {
        is_transformed = true;
        this->networkDim = 3;
        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i) {
          procCoords[i] = new pcoord_t[this->numRanks] ;
//          this->proc_coords[permutation[i]];
        }
        for (int i = 0; i < this->numRanks; ++i) {
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
        if(this->myRank == 0) 
          std::cout << "Transforming the coordinates" << std::endl;
//        this->printAllocation();
      }
      else if(optimization_level >= 3) {
        is_transformed = true;
        this->networkDim = 6;
        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i) {
          procCoords[i] = new pcoord_t[this->numRanks] ;
//        this->proc_coords[permutation[i]];
        }

//        this->machine_extent[0] = this->actual_machine_extent
        this->machine_extent = new int[networkDim];

        this->machine_extent[0] = 
          ceil (int (this->actual_machine_extent[0]) / 2.0) * 64 ;
        this->machine_extent[3] = 2 * 8 ;
        this->machine_extent[1] = 
          ceil(int (this->actual_machine_extent[1])  / 2.0) * 8 * 2400;
        this->machine_extent[4] = 2 * 8;
        this->machine_extent[2] = 
          ceil((int (this->actual_machine_extent[2])) / 8.0) * 160;
        this->machine_extent[5] = 8 * 5;

        for (int k = 0; k < this->numRanks ; k++) {
          // This part is for titan.
          // But it holds for other 3D torus machines such as Bluewaters.

          // Bandwitdh along
          // X = 75
          // Y = 37.5 or 75 --- everyother has 37.5 
          // --- Y[0-1] =75 but Y[1-2]=37.5
          // Z = 75 or 120 ---- Y[0-1-2-3-4-5-6-7] = 120, Y[7-8] = 75

          // Along X we make groups of 2. Then scale the distance with 64.
          // First dimension is represents x/2
          procCoords[0][k] = (int (this->actual_procCoords[0][k]) / 2) * 64;
          // Then the 3rd dimension is x%2. distance is scaled with 8, 
          // reversely proportional with bw=75
          procCoords[3][k] = (int (this->actual_procCoords[0][k]) % 2) * 8 ;

          // Along Y. Every other one has the slowest link. So we want 
          // distances between Y/2 huge.
          // We scale Y/2 with 2400 so that we make sure that it is the 
          // first one we divie.
          procCoords[1][k] = 
            (int (this->actual_procCoords[1][k])  / 2) * 8 * 2400;
          // The other one is scaled with 8 as in X.
          procCoords[4][k] = (int (this->actual_procCoords[1][k])  % 2) * 8;

          // We make groups of 8 along Z. Then distances between these 
          // groups are scaled with 160.
          // So that it is more than 2x distance than the distance with X 
          // grouping.
          // That is we scale the groups of Zs with 160. Groups of X with 64.
          // Zs has 8 processors connecting them, while X has only one. We 
          // want to divide along Z twice before dividing along X.
          procCoords[2][k] = 
            ((int (this->actual_procCoords[2][k])) / 8) * 160;
          // In the second group everything is scaled with 5, as bw=120
          procCoords[5][k] = ((int (this->actual_procCoords[2][k])) % 8) * 5;
        }
      }
      else if(optimization_level == 2) {
        //  This is as above case. but we make groups of 3 along X instead.
        is_transformed = true;
        this->networkDim = 6;
        procCoords = new pcoord_t * [networkDim];
        for(int i = 0; i < networkDim; ++i) {
          procCoords[i] = new pcoord_t[this->numRanks] ;
//          this->proc_coords[permutation[i]];
        }

//        this->machine_extent[0] = this->actual_machine_extent
        this->machine_extent = new int[networkDim];

        this->machine_extent[0] = 
          ceil(int (this->actual_machine_extent[0]) / 3.0) * 128 ;
        this->machine_extent[3] = 3 * 8 ;
        this->machine_extent[1] = 
          ceil(int (this->actual_machine_extent[1])  / 2.0) * 8 * 2400;
        this->machine_extent[4] = 2 * 8;
        this->machine_extent[2] = 
          ceil((int (this->actual_machine_extent[2])) / 8.0) * 160;
        this->machine_extent[5] = 8 * 5;


        for (int k = 0; k < this->numRanks ; k++) {
          // This part is for titan.
          // But it holds for other 3D torus machines such as Bluewaters.

          // Bandwitdh along
          // X = 75
          // Y = 37.5 or 75 --- everyother has 37.5 
          // --- Y[0-1] =75 but Y[1-2]=37.5
          // Z = 75 or 120 ---- Y[0-1-2-3-4-5-6-7] = 120, Y[7-8] = 75

          // In this case we make groups of 3. along X.
          procCoords[0][k] = (int (this->actual_procCoords[0][k]) / 3) * 128;
          // Then the 3rd dimension is x%2. distance is scaled with 8, 
          // reversely proportional with bw=75
          procCoords[3][k] = (int (this->actual_procCoords[0][k]) % 3) * 8 ;

          // Along Y. Every other one has the slowest link. So we want 
          // distances between Y/2 huge.
          // We scale Y/2 with 2400 so that we make sure that it is the 
          // first one we divie.
          procCoords[1][k] = 
            (int (this->actual_procCoords[1][k])  / 2) * 8 * 2400;
          // The other one is scaled with 8 as in X.
          procCoords[4][k] = (int (this->actual_procCoords[1][k])  % 2) * 8;


          procCoords[2][k] = 
            ((int (this->actual_procCoords[2][k])) / 8) * 160;
          // In the second group everything is scaled with 5, as bw=120
          procCoords[5][k] = ((int (this->actual_procCoords[2][k])) % 8) * 5;
        }
      }
    }
  }

  virtual ~MachineTorusRCAForTesting() {
    if (is_transformed) {
      is_transformed = false;
      for (int i = 0; i < actual_networkDim; i++) {
        delete [] actual_procCoords[i];
      }
      delete [] actual_procCoords;
      delete [] actual_machine_extent;
    }
    for (int i = 0; i < networkDim; i++) {
      delete [] procCoords[i];
    }
    delete [] procCoords;
    delete [] machine_extent;
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return this->networkDim; }
  int getRealMachineDim() const { return this->actual_networkDim; }

  bool getMachineExtent(int *nxyz) const {
    if (is_transformed) {
      return false;
    }
    else {
      int dim = 0;
      nxyz[dim++] = this->machine_extent[0]; //x
      nxyz[dim++] = this->machine_extent[1]; //y
      nxyz[dim++] = this->machine_extent[2]; //z
      return true;
    }
  }

  bool getRealMachineExtent(int *nxyz) const {
    int dim = 0;
    nxyz[dim++] = 25; //x
    nxyz[dim++] = 16; //y
    nxyz[dim++] = 24; //z
    return true;
  }


  void printAllocation() {
    if(this->myRank == 0) {
      for (int i = 0; i < this->numRanks; ++i) {
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

  bool getMyMachineCoordinate(pcoord_t *xyz) {
    for (int i = 0; i < this->networkDim; ++i) {
      xyz[i] = procCoords[i][this->myRank];
    }
    return true;
  }

  bool getMyActualMachineCoordinate(pcoord_t *xyz) {
    xyz[0] = rand() % 25;
    xyz[1] = rand() % 16;
    xyz[2] = rand() % 24;
    return true;
  }

  inline bool getMachineCoordinate(const int rank,
                                   pcoord_t *xyz) const {
    for (int i = 0; i < this->networkDim; ++i) {
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

  virtual bool getHopCount(int rank1, int rank2, pcoord_t &hops) const override {
    hops = 0;
    for (int i = 0; i < networkDim; ++i) {
      pcoord_t distance = procCoords[i][rank1] - procCoords[i][rank2];
      if (distance < 0) 
        distance = -distance;
      if (machine_extent[i] - distance < distance) 
        distance = machine_extent[i] - distance;
      hops += distance;
    }
    return true;
  }


private:

  int networkDim;
  int actual_networkDim;

  pcoord_t **procCoords;
  pcoord_t **actual_procCoords;

  part_t *machine_extent;
  part_t *actual_machine_extent;
  bool is_transformed;


  const Teuchos::ParameterList *pl;

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
                                        this->numRanks, 
                                        procCoords[i], tmpVect);
      pcoord_t *tmp = tmpVect;
      tmpVect = procCoords[i];
      procCoords[i] = tmp;
    }
    delete [] tmpVect;
  }

};

} // namespace Zoltan2
#endif
