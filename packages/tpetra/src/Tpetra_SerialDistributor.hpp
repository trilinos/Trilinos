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

#ifndef TPETRA_SERIALDISTRIBUTOR_HPP
#define TPETRA_SERIALDISTRIBUTOR_HPP

#include "Tpetra_Object.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

//! Tpetra::SerialDistributor:  The Tpetra Serial implementation of the Tpetra::Distributor Gather/Scatter Setup Class.
/*! The SerialDistributor class is an Serial implement of Tpetra::Distributor that is essentially a trivial class
    since a serial machine is a trivial parallel machine.
		An SerialDistributor object is actually produced by calling a method in the Tpetra::SerialPlatform class.

		Most SerialDistributor methods throw an error of -1, since they should never be called.
*/

	template<typename OrdinalType>
	class SerialDistributor : public Object, public virtual Distributor<OrdinalType> {
  public:

    //@{ \name Constructor/Destructor
    
    //! Default Constructor.
    SerialDistributor() : Object("Tpetra::SerialDistributor") {};
    
    //! Copy Constructor
    SerialDistributor(SerialDistributor<OrdinalType> const& plan) : Object(plan.label()) {};
    
    //! Clone method
    Distributor<OrdinalType>* clone() {
      SerialDistributor<OrdinalType>* distributor = new SerialDistributor<OrdinalType>(*this); 
      return(distributor); 
    };
    
    //! Destructor.
    virtual ~SerialDistributor() {};

    //@}
    
    //@{ \name Gather/Scatter Constructors

    //! Create Distributor object using list of Image IDs to send to
    void createFromSends(OrdinalType const& numExportIDs, 
                         std::vector<OrdinalType> const& exportImageIDs,
                         bool const& deterministic,
                         OrdinalType& numRemoteIDs) 
		{
      throw reportError("This method should never be called.", -1);
    };

    //! Create Distributor object using list of Image IDs to receive from
    void createFromRecvs(OrdinalType const& numRemoteIDs, 
                         std::vector<OrdinalType> const& remoteGIDs, 
                         std::vector<OrdinalType> const& remoteImageIDs, 
                         bool const& deterministic, 
                         OrdinalType& numExportIDs, 
                         std::vector<OrdinalType>& exportGIDs, 
                         std::vector<OrdinalType>& exportImageIDs)
		{
      throw reportError("This method should never be called.", -1);
    };

    //@}
    
    //@{ \name Execute Gather/Scatter Operations (Constant size objects)

    //! doPostsAndWaits
    void doPostsAndWaits(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) {
      throw reportError("This method should never be called.", -1);
    };
    
    //! doPosts
    void doPosts(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) {
      throw reportError("This method should never be called.", -1);
    };

    //! doWaits
    void doWaits() {
      throw reportError("This method should never be called.", -1);
    };

    //! doReversePostsAndWaits
    void doReversePostsAndWaits(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) {
      throw reportError("This method should never be called.", -1);
    };
    
    //! doReversePosts
    void doReversePosts(char* export_objs, OrdinalType const& obj_size, OrdinalType& len_import_objs, char*& import_objs) {
      throw reportError("This method should never be called.", -1);
    };

    //! doReverseWaits
    void doReverseWaits() {
      throw reportError("This method should never be called.", -1);
    };

    //@}
    
    //@{ \name I/O Methods

    //! print method inherited from Object
    void print(ostream& os) const {};

    //! printInfo method inherited from Distributor
    void printInfo(ostream& os) const {os << *this;};

    //@}
    
  }; // class SerialDistributor
  
} // namespace Tpetra

#endif // TPETRA_SERIALDISTRIBUTOR_HPP
