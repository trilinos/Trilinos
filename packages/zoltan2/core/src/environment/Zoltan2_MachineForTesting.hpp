// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MACHINEFORTESTING_HPP_
#define _ZOLTAN2_MACHINEFORTESTING_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

namespace Zoltan2{

/*! \brief A Machine Class for testing only
 *  A more realistic machine should be used for task mapping.
 */

template <typename pcoord_t, typename part_t>
class MachineForTesting : public Machine<pcoord_t, part_t> {

public:
  /*! \brief Constructor: A default machine description used only for testing; 
   * it does not contain actual machine information.
   *  \param comm Communication object.
   */

  MachineForTesting(const Teuchos::Comm<int> &comm): 
    Machine<pcoord_t,part_t>(comm),
    networkDim(3),
    procCoords(NULL)
  {
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

  MachineForTesting(const Teuchos::Comm<int> &comm, const Teuchos::ParameterList &pl):
    Machine<pcoord_t,part_t>(comm),
    networkDim(3),
    procCoords(NULL)
  {
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

  virtual ~MachineForTesting() {
    for (int i = 0; i < networkDim; i++){
      delete [] procCoords[i];
    }
    delete [] procCoords;
  }

  bool hasMachineCoordinates() const { return true; }

  int getMachineDim() const { return networkDim; }

  bool getMachineExtent(int *nxyz) const { 
    // Ficticious machine extent
    nxyz[0] = this->numRanks;
    nxyz[1] = 2*this->numRanks;
    nxyz[2] = 3*this->numRanks;
    return true; 
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {
    return getMachineCoordinate(this->myRank, xyz);
  }

  bool getMachineCoordinate(const int rank, pcoord_t *xyz) {
    // Ficticious machine coordinates
    // part_t slice = part_t(pow(double(this->numRanks), double(1.0/networkDim))
    //                           + 0.5);
    // part_t m = rank;
    // for (int i = 0; i < networkDim; ++i){
    //   xyz[i] = m / part_t(pow(slice, double(networkDim - i - 1)));
    //   m = m % part_t(pow(double(slice), double(networkDim - i - 1)));
    // }

    xyz[0] = rank;
    xyz[1] = this->numRanks;
    xyz[2] = this->numRanks+1;
    return true;
  }

  bool getMachineCoordinate(const char *nodename, pcoord_t *xyz) {
    return false;  // cannot yet return from nodename
  }

  bool getAllMachineCoordinatesView(pcoord_t **&allCoords) const {
    allCoords = procCoords;
    return true;
  }

private:

  int networkDim;

  pcoord_t **procCoords;   // KDD Maybe should be RCP?

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
