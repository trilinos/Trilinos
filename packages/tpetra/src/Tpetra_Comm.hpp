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

#ifndef TPETRA_COMM_HPP
#define TPETRA_COMM_HPP

#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {

//! Tpetra::Comm:  The Tpetra Communication Abstract Base Class.
/*! The Tpetra Comm class is an interface that encapsulates the general
  information and services needed for other Tpetra classes to run on a parallel computer.
  
  Comm currently has one default implementation, via SerialComm, for serial execution.
  (A second default implementation is planned. MpiComm will be for MPI
  distributed memory execution.)  It is meant to insulate the user from
  the specifics of communication that are not required for normal
  manipulation of linear algebra objects.  Most Comm interfaces are similar to MPI
  interfaces, except that the type of data is not required as an argument since C++ can bind
  to the appropriate interface based on argument typing.
*/

template<typename PacketType, typename OrdinalType>
class Comm {
  public:

  //@{ \name Constructor/Destructor Methods
  //! Destructor
  virtual ~Comm() {};
  //@}
  
  //@{ \name Barrier Methods
  //! Barrier. 
  /*! Each image must stop until all images have reached the barrier.
   */
  virtual void barrier() const = 0;
  //@}

  //@{ \name Broadcast Methods
  //! Broadcast
  /*!Take list of input values from the root image and sends to all other images.
    \param myVals InOut
           On entry, the root image contains the list of values.  On exit,
	   all images will have the same list of values.  Note that values must be
	   allocated on all images before the broadcast.
    \param count In
           On entry, contains the length of myVals.
    \param root In
           On entry, contains the imageID from which all images will receive a copy of myVals.
  */
  virtual void broadcast(PacketType* myVals, OrdinalType const count, int const root) const = 0;
  //@}

  //@{ \name Gather Methods
  //! Gather All function.
  /*! Take list of input values from all images in the communicator and creates an ordered contiguous list of
    those values on each image.
    \param myVals In
           On entry, contains the list of values to be sent to all images.
    \param allVals Out
           On exit, contains the list of values from all images. Must be of size numImages*count.
    \param count In
           On entry, contains the length of myVals.
  */
  virtual void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const = 0;
  //@}

  //@{ \name Sum Methods
  //! Global Sum function.
  /*!Take list of input values from all images in the communicator, computes the sum and returns the
    sum to all images.
    \param partialSums In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all images.
    \param globalSums Out
           On exit, contains the list of values summed across all images.
    \param count In
           On entry, contains the length of partialSums.
  */
  virtual void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const = 0;
  //@}
	
  //@{ \name Max/Min Methods
  //! Global Max function.
  /*! Take list of input values from all images in the communicator, computes the max and returns the
    max to all images.
    \param partialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
	   using these Partial Maxs, the max across all images will be computed.
    \param globalMaxs Out
           On exit, contains the list of maxs computed across all images.
    \param count In
           On entry, contains the length of partialMaxs.
  */
  virtual void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const = 0;
  //! Global Min function.
  /*! Take list of input values from all images in the communicator, computes the min and returns the
    min to all images.
    \param partialMins In
           On entry, contains the list of values, usually partial mins computed locally;
	   using these Partial Mins, the min across all images will be computed.
    \param globalMins Out
           On exit, contains the list of mins computed across all images.
    \param count In
           On entry, contains the length of partialMins.
  */
  virtual void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const = 0;
  //@}

  //@{ \name Parallel Prefix Methods
  //! Scan Sum function.
  /*! Take list of input values from all images in the communicator, computes the scan sum and returns it 
    to all images such that image i contains the sum of values from image 0 up to and including
   image i.
    \param myVals In
           On entry, contains the list of values to be summed across all images.
    \param scanSums Out
           On exit, contains the list of values summed across images 0 through i.
    \param count In
           On entry, contains the length of myVals.
  */
  virtual void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const = 0;
  //@}

	//@{ \name I/O Methods
	//! printInfo
	virtual void printInfo(ostream& os) const = 0;
	//@}
	
}; // class Comm

} // namespace Tpetra

#endif // TPETRA_COMM_HPP
