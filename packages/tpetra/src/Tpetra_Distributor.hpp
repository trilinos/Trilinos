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

#ifndef TPETRA_DISTRIBUTOR_HPP
#define TPETRA_DISTRIBUTOR_HPP

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
  
  template<typename OrdinalType>
  class Distributor {
    
  public:
    //@{ \name Constructor and Destructor

    //! clone constructor.
    virtual Distributor<OrdinalType>* clone() = 0;

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
    virtual void createFromSends(OrdinalType const& numExportIDs, 
                                 std::vector<OrdinalType> const& exportImageIDs,
                                 bool const& deterministic,
                                 OrdinalType& numRemoteIDs) = 0;
    
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
    virtual void createFromRecvs(OrdinalType const& numRemoteIDs, 
                                 std::vector<OrdinalType> const& remoteGIDs, 
                                 std::vector<OrdinalType> const& remoteImageIDs, 
                                 bool const& deterministic, 
                                 OrdinalType& numExportIDs, 
                                 std::vector<OrdinalType>& exportGIDs, 
                                 std::vector<OrdinalType>& exportImageIDs) = 0;

    //@}
    
    //@{ \name Execute Gather/Scatter Operations (Constant size objects)
    
    //! Execute plan on buffer of export objects in a single step
    virtual void doPostsAndWaits(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) = 0;
    
    //! Post buffer of export objects (can do other local work before executing Waits)
    virtual void doPosts(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) = 0;
    
    //! Wait on a set of posts
    virtual void doWaits() = 0;
    
    //! Execute reverse of plan on buffer of export objects in a single step
    virtual void doReversePostsAndWaits(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) = 0;
    
    //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
    virtual void doReversePosts(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) = 0;
    
    //! Wait on a reverse set of posts
    virtual void doReverseWaits() = 0;

    //@}
    
    //@{ \name I/O Methods

    //! printInfo
    virtual void printInfo(ostream& os) const = 0;

    //@}
    
  }; // class Distributor
  
} // namespace Tpetra

#endif // TPETRA_DISTRIBUTOR_HPP
