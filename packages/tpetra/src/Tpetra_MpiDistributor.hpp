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

#ifndef TPETRA_MPIDISTRIBUTOR_HPP
#define TPETRA_MPIDISTRIBUTOR_HPP

#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

  // forward declaration of MpiData, needed to prevent circular inclusions
  // actual #include statement at the end of this file
  class MpiData;
  
  //! Tpetra::MpiDistributor:  The Tpetra MPI implementation of the Tpetra::Distributor Gather/Scatter Setup Class.
  /*! The MpiDistributor class is an MPI implement of Tpetra::Distributor.
		An MpiDistributor object is actually produced by calling a method in the Tpetra::MpiPlatform class.
  */
  
	template<typename PacketType, typename OrdinalType>
	class MpiDistributor : public Object, public virtual Distributor<PacketType, OrdinalType> {
  public:
    
    //@{ \name Constructor/Destructor
    
    //! default constuctor
    /*! This is used when the user wants to create an MpiDistributor instance directly.
      \param In
      Comm - MPI_Comm communicator we will use.
    */
		MpiDistributor(MPI_Comm Comm) 
      : Object("Tpetra::MpiDistributor") 
      , MpiData_()
    {
      MpiData_ = Teuchos::rcp(new MpiData(Comm));
    };
    
    //! platform constructor
    /*! This is used by MpiPlatform to create an MpiDistributor instance. It should
      not be called directly by the user.
      \param In
      MpiData - MpiData inner data class passed in by MpiPlatform.
    */
    MpiDistributor(Teuchos::RefCountPtr<MpiData> const& data)
      : Object("Tpetra::MpiDistributor")
      , MpiData_(data)
    {};
    
    //! Copy Constructor
    MpiDistributor(MpiDistributor<PacketType, OrdinalType> const& distributor) 
      : Object(distributor.label())
      , MpiData_(distributor.MpiData_)
    {};
    
    //! Clone method
    Distributor<PacketType, OrdinalType>* clone() {
      MpiDistributor<PacketType, OrdinalType>* distributor = new MpiDistributor<PacketType, OrdinalType>(*this); 
      return(distributor); 
    };
    
    //! Destructor.
    ~MpiDistributor() {};
    //@}
    

    //@{ \name Gather/Scatter Constructors

    //! Create Distributor object using list of Image IDs to send to
    void createFromSends(OrdinalType const& numExportIDs, 
                         OrdinalType const* exportImageIDs,
                         bool const& deterministic, 
                         OrdinalType& numRemoteIDs ) 
		{
      throw reportError("This method is not implemented yet.", -1);
    };

    //! Create Distributor object using list of Image IDs to receive from
    void createFromRecvs(OrdinalType const& numRemoteIDs, 
                         OrdinalType const* remoteGIDs, 
                         OrdinalType const* remoteImageIDs, 
                         bool const& deterministic, 
                         OrdinalType& numExportIDs, 
                         OrdinalType*& exportGIDs, 
                         OrdinalType*& exportImageIDs)
		{
      throw reportError("This method is not implemented yet.", -1);
    };

    //@}
    

    //@{ \name Constant Size

    //! do
    void doPostsAndWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //! doReverse
    void doReversePostsAndWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };
    
    //! doPosts
    void doPosts(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //! doWaits
    void doWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };
    
    //! doReversePosts
    void doReversePosts(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //! doReverseWaits
    void doReverseWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //@}
    

    //@{ \name Non-Constant Size
 
   //! do
    void doPostsAndWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //! doReverse
    void doReversePostsAndWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };
    
    //! doPosts
    void doPosts(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //! doWaits
    void doWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };
    
    //! doReversePosts
    void doReversePosts(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //! doReverseWaits
    void doReverseWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs) {
      throw reportError("This method is not implemented yet.", -1);
    };

    //@}
    

    //@{ \name I/O Methods

    //! print method inherited from Object
    void print(ostream& os) const {};

    //! printInfo method inherited from Distributor
    void printInfo(ostream& os) const {os << *this;};

    //@}
    
  private:
    
    Teuchos::RefCountPtr<MpiData> MpiData_;
    
    // convenience functions for returning inner data class, both const and nonconst versions.
    MpiData& data() {return(*MpiData_);}
    MpiData const& data() const {return(*MpiData_);}


    
  }; // class MpiDistributor
  
} // namespace Tpetra

#include "Tpetra_MpiData.hpp"

#endif // TPETRA_MPIDISTRIBUTOR_HPP
