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

#ifndef TPETRA_DIRECTORY_HPP
#define TPETRA_DIRECTORY_HPP

namespace Tpetra {

  //! Tpetra::Directory - Directory abstract base class
  
  /*! A Directory object must be created by a call to the Platform createDirectory method. 
      The Directory is needed to allow referencing of non-local elements.
  */
  
  template<typename OrdinalType>
  class Directory {
  public:
    
    //@{ \name Constructors/Destructor.
    //! Directory destructor.
    virtual ~Directory() {};
    //@}
    
    //@{ \name Query methods.
    
    //! getDirectoryEntries : Returns image info for non-local ElementSpace entries
    /*! Given a list of Global Entry IDs, this function returns the list of
        IDs of the owning memory image that correspond to the list of entries.
      \param In
             globalEntries - List of Global IDs being passed in.
      \param Out
             imageIDs - On return contains list of Image IDs owning the Global IDs in question.
    */
    virtual void getDirectoryEntries(std::vector<OrdinalType> const globalEntries, std::vector<OrdinalType> images) const = 0;
    
    //! getDirectoryEntries : Returns image and local id info for non-local ElementSpace entries
    /*! Given a list of Global Entry IDs, this function returns the list of
        image IDs and local IDs on the owning memory image that correspond
        to the list of entries.  If LocalEntries is 0, then local IDs are 
        not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
        element sizes for the requested global entries.
      \param In
             globalEntries - List of Global IDs being passed in.
      \param Out
             imageIDs - On return contains list of Image IDs owning the Global IDs in question.
      \param Out
             localEntries - On return contains the local ID of the global on the owning image. 
    */
    virtual void getDirectoryEntries(std::vector<OrdinalType> const globalEntries, std::vector<OrdinalType> images, std::vector<OrdinalType> localEntries) const = 0;

    //@}
    
  }; // class Directory
  
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_HPP
