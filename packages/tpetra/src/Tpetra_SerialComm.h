/*Paul
04-Aug-2002 Status: Templated for class T. All Epetra methods except Distributor. Fully documented. Switched to images.
03-Sept-2002 Took out Directory and ImageID methods. Templated for PacketType, OrdinalType.
*/

#ifndef _TPETRA_SERIALCOMM_H_
#define _TPETRA_SERIALCOMM_H_

#include "Tpetra_Comm.h"
#include "Tpetra_Object.h"

namespace Tpetra {

//! Tpetra::SerialComm:  The Tpetra Serial Communication Class.
/*! The SerialComm class is an implementation of Tpetra::Comm, providing the general
  information and services needed for other Tpetra classes to run on a serial computer.
*/

template<class PacketType, class OrdinalType>
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
  SerialComm(const SerialComm<PacketType, OrdinalType>& Comm);
  
  //! Destructor.
  /*! Completely deletes a SerialComm object.  
    \warning Note:  All objects that depend on a SerialComm instance 
    should be destroyed prior to calling this function.
  */
  ~SerialComm();
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
  void broadcast(PacketType* myVals, OrdinalType count, int root) const {};
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
  void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType count) const;
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
  void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType count) const;
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
  void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType count) const;
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
  void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType count) const;
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
  void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType count) const;
  //@}

}; // class SerialComm

} // namespace Tpetra

#include "Tpetra_SerialComm.cpp"
#endif // _TPETRA_SERIALCOMM_H_
