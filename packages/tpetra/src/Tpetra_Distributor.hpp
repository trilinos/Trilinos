/*Paul
18-Nov-2002 Copied from Epetra_Distributor.h
19-Nov-2002 Templated for <PT,OT>, modified slightly.
23-Nov-2002 do methods named Do temporarily
25-Nov-2002 do methods fixed
*/

#ifndef _TPETRA_DISTRIBUTOR_HPP_
#define _TPETRA_DISTRIBUTOR_HPP_

#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {

//! Tpetra::Distributor:  The Tpetra Gather/Scatter Setup Abstract Base Class.
/*! The Distributor class is an interface that encapsulates the general
  information and services needed for other Tpetra classes to perform gather/scatter
  operations on a parallel computer.
  A Distributor object is actually produced by calling a method in the Tpetra::Platform class.
  
  Distributor currently has one default implementations, Tpetra::SerialDistributor, 
	for serial execution, although support for MPI distributed memory execution is planned.
  It is meant to insulate the user from the specifics of communication 
	that are not required for normal manipulation of linear algebra objects.
*/

template<typename PacketType, typename OrdinalType>
class Distributor {
	
  public:
  //@{ \name Constructor and Destructor
  //! clone constructor.
  virtual Distributor<PacketType, OrdinalType>* clone() = 0;
  //! Destructor.
  virtual ~Distributor() {};
  //@}

  //@{ \name Gather/Scatter Constructors
  //! Create Distributor object using list of Image IDs to which we export
  /*! Take a list of Image IDs and construct a plan for efficiently scattering to these images.
      Return the number of IDs being sent to me.
    \param numExportIDs In
           Number of IDs that need to be sent from this image.
    \param exportImageIDs In
           List of images that will get the exported IDs.
    \param deterministic In
           If set to true, communication will be deterministic (repeatable) from call to call.
    \param numRemoteIDs Out
           Number of IDs this image will be receiving.
  */
  virtual void createFromSends(const OrdinalType& numExportIDs, const OrdinalType* exportImageIDs,
															 const bool& deterministic, OrdinalType& numRemoteIDs ) = 0;

  //! Create Distributor object using list of Remote global IDs and corresponding ImageIDs
  /*! Take a list of global IDs and construct a plan for efficiently scattering to these images.
      Return the number and list of IDs being sent by me.
    \param numRemoteIDs In
           Number of IDs this image will be receiving.
    \param remoteGIDs In
           List of IDs that this image wants.
    \param remoteImageIDs In
           List of images that will send the remote IDs.
    \param deterministic In
           If set to true, communication will be deterministic (repeatable) from call to call.
    \param numExportIDs Out
           Number of IDs that need to be sent from this image.
    \param exportImageIDs Out
           List of images that will get the exported IDs.
  */
  virtual void createFromRecvs(const OrdinalType& numRemoteIDs, const OrdinalType* remoteGIDs, 
															 const OrdinalType* remoteImageIDs, const bool& deterministic, 
															 OrdinalType& numExportIDs, OrdinalType*& exportGIDs, 
															 OrdinalType*& exportImageIDs) = 0;
  //@}

  //@{ \name Execute Gather/Scatter Operations (Constant size objects)

  //! Execute plan on buffer of export objects in a single step
  virtual void doPostsAndWaits(PacketType* export_objs, const OrdinalType& obj_size, PacketType* import_objs) = 0;

  //! Execute reverse of plan on buffer of export objects in a single step
  virtual void doReversePostsAndWaits(PacketType* export_objs, const OrdinalType& obj_size, PacketType* import_objs) = 0;

  //! Post buffer of export objects (can do other local work before executing Waits)
  virtual void doPosts(PacketType* export_objs, const OrdinalType& obj_size, PacketType* import_objs) = 0;
  //! Wait on a set of posts
  virtual void doWaits(PacketType* export_objs, const OrdinalType& obj_size, PacketType* import_objs) = 0;

  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  virtual void doReversePosts(PacketType* export_objs, const OrdinalType& obj_size, PacketType* import_objs) = 0;

  //! Wait on a reverse set of posts
  virtual void doReverseWaits(PacketType* export_objs, const OrdinalType& obj_size, PacketType* import_objs) = 0;
  //@}

  //@{ \name Execute Gather/Scatter Operations (Non-constant size objects)

  //! Execute plan on buffer of export objects in a single step (object size may vary)
  virtual void doPostsAndWaits(PacketType* export_objs, const OrdinalType*& obj_size, PacketType* import_objs) = 0;

  //! Execute reverse of plan on buffer of export objects in a single step (object size may vary)
  virtual void doReversePostsAndWaits(PacketType* export_objs, const OrdinalType*& obj_size, PacketType* import_objs) = 0;

  //! Post buffer of export objects (can do other local work before executing Waits)
  virtual void doPosts(PacketType* export_objs, const OrdinalType*& obj_size, PacketType* import_objs) = 0;
  //! Wait on a set of posts
  virtual void doWaits(PacketType* export_objs, const OrdinalType*& obj_size, PacketType* import_objs) = 0;

  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  virtual void doReversePosts(PacketType* export_objs, const OrdinalType*& obj_size, PacketType* import_objs) = 0;

  //! Wait on a reverse set of posts
  virtual void doReverseWaits(PacketType* export_objs, const OrdinalType*& obj_size, PacketType* import_objs) = 0;
  //@}

  //@{ \name I/O Methods
	//! printInfo
  virtual void printInfo(ostream& os) const = 0;
  //@}

}; // class Distributor

} // namespace Tpetra

#endif // _TPETRA_DISTRIBUTOR_HPP_
