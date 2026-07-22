// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MACHINE_TORUS_TOPOMANAGER_HPP_
#define _ZOLTAN2_MACHINE_TORUS_TOPOMANAGER_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Machine.hpp>

#ifdef HAVE_ZOLTAN2_TOPOMANAGER
#include <TopoManager.h>
#endif

namespace Zoltan2{

/*! \brief A TopoManager Machine Class on Torus Networks
 */

template <typename pcoord_t, typename part_t>
class MachineTorusTopoMgr : public Machine <pcoord_t, part_t> {

public:
  /*! \brief Constructor: A BlueGeneQ network machine description;
   *  \param comm Communication object.
   */

  MachineTorusTopoMgr(const Teuchos::Comm<int> &comm):
    Machine<pcoord_t,part_t>(comm),
#if defined (CMK_BLUEGENEQ)
    networkDim(6),  tmgr(comm.getSize()),
#elif defined (CMK_BLUEGENEP)
    networkDim(4), tmgr(comm.getSize()),
#else
    networkDim(3),
#endif
    procCoords(NULL), machine_extent(NULL),
    delete_transformed_coords(false), 
    transformed_network_dim(0),
    transformed_coordinates (NULL), pl(NULL)
  {
    transformed_network_dim = networkDim - 1;
    transformed_coordinates = procCoords;
    machine_extent = new int[networkDim];
    this->getMachineExtent(this->machine_extent);
    
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

  MachineTorusTopoMgr(const Teuchos::Comm<int> &comm, 
                 const Teuchos::ParameterList &pl_ ):
    Machine<pcoord_t,part_t>(comm),
#if defined (CMK_BLUEGENEQ)
    networkDim(6),  tmgr(comm.getSize()),
#elif defined (CMK_BLUEGENEP)
    networkDim(4), tmgr(comm.getSize()),
#else
    networkDim(3),
#endif
    procCoords(NULL), machine_extent(NULL),
    delete_transformed_coords(false), 
    transformed_network_dim(0),
    transformed_coordinates (NULL),
    pl(&pl_)
  {
    transformed_network_dim = networkDim - 1;
    transformed_coordinates = procCoords;
    machine_extent = new int[networkDim];
    this->getMachineExtent(this->machine_extent);
    
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

  virtual ~MachineTorusTopoMgr() {
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

  int getRealMachineDim() const { return networkDim; }

  bool getMachineExtent(int *nxyz) const {
#if defined (CMK_BLUEGENEQ)
    int dim = 0;
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNA();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNB();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNC();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimND();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNE();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNT();
    return true;
#elif defined (CMK_BLUEGENEP)
    int dim = 0;
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNX();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNY();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNZ();
    if (dim < transformed_network_dim)
      nxyz[dim++] = tmgr.getDimNT();
    return true;
#else
    return false;
#endif
  }

  // MD TODO: Not always it has wrap-around links.
  bool getMachineExtentWrapArounds(bool *wrap_around) const {
#if defined (CMK_BLUEGENEQ)
    // Leave it as this for now, figure out if there is a way to 
    // determine tourus from topomanager.
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
#elif defined (CMK_BLUEGENEP)
    int dim = 0;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
    if (dim < transformed_network_dim)
      wrap_around[dim++] = true;
#else
#endif
    return true;
  }

  bool getMyMachineCoordinate(pcoord_t *xyz) {
    for (int i = 0; i < this->transformed_network_dim; ++i) {
      xyz[i] = transformed_coordinates[i][this->myRank];
    }
    return true;
  }

  bool getMyActualMachineCoordinate(pcoord_t *xyz) {
#if defined (CMK_BLUEGENEQ)
    int a,b,c,d,e,t;
    tmgr.rankToCoordinates(this->myRank, a,b,c,d,e,t);
    xyz[0] = a; xyz[1] = b; xyz[2] = c; xyz[3] = d; xyz[4] = e; xyz[5] = t;
    //std::cout << "me:" << this->myRank 
    //  << " " << a << " " << b << " " << c << " " << d 
    //  << " " << e << " " << t << std::endl;
    return true;
#elif defined (CMK_BLUEGENEP)
    int a,b,c,t;
    tmgr.rankToCoordinates(this->myRank, a,b,c,t);
    xyz[0] = a; xyz[1] = b; xyz[2] = c; xyz[3] = t;
    return true;
#else
    return false;
#endif
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
    return true;
  }


private:

  int networkDim;

#ifdef HAVE_ZOLTAN2_TOPOMANAGER
  TopoManager tmgr;
#endif
  pcoord_t **procCoords;   // KDD Maybe should be RCP?
  part_t *machine_extent;
  const Teuchos::ParameterList *pl;


  bool delete_transformed_coords;
  int transformed_network_dim;
  pcoord_t **transformed_coordinates;

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
}
#endif
