/*Paul
04-Aug-2002 Status: Templated for class T. All Epetra methods except Distributor. Fully documented. Switched to images.
03-Sept-2002 Took out Directory and ImageID methods. Templated for PacketType, OrdinalType.
12-Oct-2002 Added some consts (still some left).
12-Nov-2002 Changed remaining template<class...> to template<typename...>
19-Nov-2002 myImageID and numImages moved back from Platform, print method updated
*/

#ifndef _TPETRA_SERIALCOMM_HPP_
#define _TPETRA_SERIALCOMM_HPP_

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
  SerialComm(const SerialComm<PacketType, OrdinalType>& Comm);
  
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
  void broadcast(PacketType* myVals, const OrdinalType count, const int root) const {};
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
  void gatherAll(PacketType* myVals, PacketType* allVals, const OrdinalType count) const;
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
  void sumAll(PacketType* partialSums, PacketType* globalSums, const OrdinalType count) const;
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
  void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, const OrdinalType count) const;
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
  void minAll(PacketType* partialMins, PacketType* globalMins, const OrdinalType count) const;
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
  void scanSum(PacketType* myVals, PacketType* scanSums, const OrdinalType count) const;
  //@}

	//@{ \name I/O Methods
	//! Print methods
	void print(ostream& os) const {os << "::Memory Image " << getMyImageID() << " of " << getNumImages() << " total images" << endl;};
	void printInfo(ostream& os) const {print(os);};
	//@}

}; // class SerialComm

} // namespace Tpetra

#include "Tpetra_SerialComm.cpp"
#endif // _TPETRA_SERIALCOMM_HPP_
