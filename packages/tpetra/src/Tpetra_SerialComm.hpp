// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_SERIALCOMM_HPP
#define TPETRA_SERIALCOMM_HPP

#include "Tpetra_Comm.hpp"
#include "Tpetra_Object.hpp"

namespace Tpetra {

//! Tpetra::SerialComm:  The Tpetra Serial Communication Class.
/*! The SerialComm class is an implementation of Tpetra::Comm, providing the general
  information and services needed for other Tpetra classes to run on a serial computer.
*/

template<typename PacketType, typename OrdinalType>
class SerialComm : public Object, public virtual Comm<PacketType, OrdinalType> {
public:
  //@{ \name Constructor/Destructor Methods
  //! Constructor
  /*! Builds an instance of a serial communicator.  Even
    if the application is running in parallel via MPI, this communicator
    will execute in serial.  The access functions return the number of
    memory images to be 1 and the image ID to be 0.
  */
  SerialComm();
  
  //! Copy constructor
  /*! Makes an exact copy of an existing SerialComm instance.
  */
  SerialComm(SerialComm<PacketType, OrdinalType> const& comm);
  
  //! Destructor.
  /*! Completely deletes a SerialComm object.  
    \warning Note:  All objects that depend on a SerialComm instance 
    should be destroyed prior to calling this function.
  */
  ~SerialComm();
  //@}

	//@{ \name Image Info Methods
	//! getMyImageID - In serial mode, always returns 0.
	int getMyImageID() const {return(0);};
	//! getNumImages - In serial mode, always returns 1.
	int getNumImages() const {return(1);};
	//@}

  //@{ \name Barrier Methods
  //! Barrier function. 
  /*! A no-op for a serial communicator.
   */
  void barrier() const {};
  //@}

  //@{ \name Broadcast Methods
  //! SerialComm Broadcast function.
  /*! A no-op for a serial communicator.
    \param myVals InOut
           On entry, the root image contains the list of values.  On exit,
	   all images will have the same list of values.  Note that values must be
	   allocated on all images before the broadcast.
    \param count In
           On entry, contains the length of the list of myVals.
    \param root In
           On entry, contains the imageID from which all images will receive a copy of myVals.
  */
  void broadcast(PacketType* myVals, OrdinalType const count, int const root) const {};
  //@}

  //@{ \name Gather Methods
  //! SerialComm All Gather function.
  /*! A copy for a serial communicator.
    \param myVals In
           On entry, contains the list of values, to be sent to all images.
    \param allVals Out
           On exit, contains the list of values from all images. Must by of size numImages*count.
    \param count In
           On entry, contains the length of the list of myVals.
  */
  void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const;
  //@}

  //@{ \name Sum Methods
  //! SerialComm Global Sum function.
  /*! A copy for a serial communicator.
    \param partialSums In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all images.
    \param globalSums Out
           On exit, contains the list of values summed across all images.
    \param count In
           On entry, contains the length of the list of values.
  */
  void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const;
  //@}
	
  //@{ \name Max/Min Methods
  //! SerialComm Global Max function.
  /*! A copy for a serial communicator.
    \param partialMaxs In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all images.
    \param globalMaxs Out
           On exit, contains the list of values summed across all images.
    \param count In
           On entry, contains the length of the list of values.
  */
  void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const;
  //! SerialComm Global Min function.
  /*! A copy for a serial communicator.
    \param partialMins In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all images.
    \param globalMins Out
           On exit, contains the list of values summed across all images.
    \param count In
           On entry, contains the length of the list of values.
  */
  void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const;
  //@}

  //@{ \name Parallel Prefix Methods
  //! SerialComm Scan Sum function.
  /*! A copy for a serial communicator.
    \param myVals In
           On entry, contains the list of values to be summed across all images.
    \param scanSums Out
           On exit, contains the list of values summed across images 0 through i.
    \param count In
           On entry, contains the length of the list of values.
  */
  void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const;
  //@}

	//@{ \name I/O Methods
	//! Print methods
	void print(ostream& os) const {os << "::Memory Image " << getMyImageID() << " of " << getNumImages() << " total images" << endl;};
	void printInfo(ostream& os) const {print(os);};
	//@}

}; // class SerialComm


// begin Tpetra_SerialComm.cpp
//=============================================================================

  // default constructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm() 
    : Object("Tpetra::Comm[Serial]") {}
  
  // copy constructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm(SerialComm<PacketType, OrdinalType> const& comm) 
    : Object(comm.label()) {}
  
  // destructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::~SerialComm() {}
  
  // gather all
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      allVals[i] = myVals[i];
  }
  
  // sum
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      globalSums[i] = partialSums[i];
  }
  
  // min/max
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      globalMaxs[i] = partialMaxs[i];
  }
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      globalMins[i] = partialMins[i];
  }
  
  // scansum
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      scanSums[i] = myVals[i];
  }

	// print

//=============================================================================
// end Tpetra_SerialComm.cpp

} // namespace Tpetra

#endif // TPETRA_SERIALCOMM_HPP
