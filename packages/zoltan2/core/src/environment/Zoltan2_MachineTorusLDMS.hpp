// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MACHINEDEFAULT_HPP_
#define _ZOLTAN2_MACHINEDEFAULT_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Zoltan2{

/*! \brief A Default MachineRepresentation Class
 *
 *  Work In Progress, use another Machine Type.
 */

template <typename nNo_t, typename nCoord_t>
class DefaultMachine : public MachineRepresentation<nNo_t, nCoord_t> {

private:
  int networkDim;
  int numProcs;
  int myRank;

  nCoord_t **procCoords;   // KDD Maybe should be RCP?

public:
  /*! \brief Constructor MachineRepresentation Class
   *  \param comm_ Communication object.
   */

  MachineRepresentation(const Comm<int> &comm):
    networkDim(0), 
    numProcs(comm.getSize()), 
    myRank(comm.getRank()),
    procCoords(NULL)
  {
    // Will need this constructor to be specific to RAAMP (MD).
    // Will need a default constructor using, e.g., GeometricGenerator
    // or nothing at all, for when RAAMP is not available as TPL.
    //
    // (AG) In addition, need to be able to run without special
    // privileges in system (e.g., on hopper).
    // Notes:  For now, all cores connected to same NIC will get the
    // same coordinates; later, we could add extra coordinate dimensions
    // to represent nodes or dies (using hwloc info through RAAMP
    // data object).

    // (MD) will modify mapping test to use machine representation
    // #ifdef HAVE_ZOLTAN2_OVIS

    // Call initializer for RAAMP data object (AG)

    // get network dimension.
    // TODO change.
    // Call RAAMP Data Object to get the network dimension (AG)
    networkDim = 3;

    //allocate memory for processor coordinates.
    procCoords = new nCoord_t *[networkDim];
    for (int i = 0; i < networkDim; ++i) {
      procCoords[i] = new nCoord_t [numProcs];
      memset (procCoords[i], 0, sizeof(nCoord_t) * numProcs);
    }
    // Obtain the coordinate of the processor.
    this->getMyCoordinate(/*nCoord_t &xyz[networkDim]*/);
    // Copy xyz into appropriate spot in procCoords. (MD)  
    // KDD I agree with this

    // reduceAll the coordinates of each processor.
    this->gatherMachineCoordinates();
  }


  /*! \brief Constructor MachineRepresentation Class
   *  \param comm_ RCP Communication object.
   */
  MachineRepresentation(const RCP<Comm<int> > &comm_):
    networkDim(0), 
    numProcs(comm_->getSize()), 
    procCoords(0), 
    comm(comm_) 
  {
    // Will need this constructor to be specific to RAAMP (MD).
    // Will need a default constructor using, e.g., GeometricGenerator
    // or nothing at all, for when RAAMP is not available as TPL.
    //
    // (AG) In addition, need to be able to run without special
    // privileges in system (e.g., on hopper).  
    // Notes:  For now, all cores connected to same NIC will get the
    // same coordinates; later, we could add extra coordinate dimensions
    // to represent nodes or dies (using hwloc info through RAAMP
    // data object).

    // (MD) will modify mapping test to use machine representation
    // #ifdef HAVE_ZOLTAN2_OVIS

    // Call initializer for RAAMP data object (AG)

    // get network dimension.
    // TODO change.
    // Call RAAMP Data Object to get the network dimension (AG)
    networkDim = 3;

    // Allocate memory for processor coordinates.
    procCoords = new nCoord_t *[networkDim];
    for (int i = 0; i < networkDim; ++i) {
      procCoords[i] = new nCoord_t [numProcs];
      memset (procCoords[i], 0, sizeof(nCoord_t) * numProcs);
    }
    // Obtain the coordinate of the processor.
    this->getMyCoordinate(/*nCoord_t &xyz[networkDim]*/);
    // Copy xyz into appropriate spot in procCoords. (MD)  // KDD I Agree.

    // reduceAll the coordinates of each processor.
    this->gatherMachineCoordinates();
  }


  /*! \brief getMyCoordinate function
   *  stores the coordinate of the current processor in procCoords[*][rank]
   */
  void getMyCoordinate(/* nCoord_t &xyz[networkDim]*/) {  
    // KDD Enable the argument rather
    // KDD than writing into array here

    // Call RAAMP system to get coordinates and store in xyz (MD)
    // What is the RAAMP call?  (AG)
    // AG will return a view (pointer) to RAAMP's data.
    // We will copy it into xyz.

//KDD #if defined(HAVE_ZOLTAN2_LDMS)
//KDD #elif defined(HAVE_ZOLTAN2_TOPOMGR)
//KDD #elif defined(HAVE_ZOLTAN2_RCA)
//KDD #else

    // The code below may be good for the default constructor, perhaps,
    // but it should copy the data into xyz instead of the procCoords.
    int myRank = comm->getRank();

    int slice = int (pow( double(numProcs), double(1.0 / networkDim)) + 0.5 );

    int m = myRank;
    for (int i = 0; i < networkDim; ++i) {
      procCoords[i][myRank] = m / int(pow(slice, double(networkDim - i - 1)));
      m = m % int(pow(double(slice), double(networkDim - i - 1)));
    }
//KDD #endif
  }

  // KDD Need to return coordinate of any rank?
  // void getCoordinate(partId_t rank, nCoord_t &xyz[networkDim]) { }

  /*! \brief gatherMachineCoordinates function
   *  reduces and stores all machine coordinates.
   */
  void gatherMachineCoordinates() {  // KDD Should be private
    nCoord_t *tmpVect = new nCoord_t [numProcs];

    for (int i = 0; i < networkDim; ++i) {
      reduceAll<int, nCoord_t>(
          *comm,
          Teuchos::REDUCE_SUM,
          numProcs,
          procCoords[i],
          tmpVect);
      nCoord_t *tmp = tmpVect;
      tmpVect = procCoords[i];
      procCoords[i] = tmp;
    }
    delete [] tmpVect;
  }

  /*! \brief destructor of the class
   * free memory in procCoords.
   */
  virtual ~MachineRepresentation() {
    for (int i = 0; i < networkDim; ++i) {
      delete [] procCoords[i];
    }
    delete [] procCoords;
    // Free/release THE RAAMP Data Object.
    // Deinitialize/finalize/whatever (AG)
  }

  /*! \brief getProcDim function
   * returns the dimension of the physical processor layout.
   */
  int getProcDim() const{  // KDD Maybe getNetworkDim or getProcCoordDim
    return networkDim;
  }

  /*! \brief getProcDim function
   * returns the coordinates of processors in two dimensional array.
   */
  nCoord_t** getProcCoords() const{  
    // KDD Make clear that returning a View; maybe return ArrayView
    return procCoords;
  }

  /*! \brief getNumProcs function
   * returns the number of processors.
   */
  int getNumProcs() const{
    return numProcs;
  }

  // KDD TODO:  Need more for full LDMS interface.

};
}
#endif
