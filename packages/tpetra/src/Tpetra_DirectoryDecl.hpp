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

#ifndef TPETRA_DIRECTORY_DECL_HPP
#define TPETRA_DIRECTORY_DECL_HPP

#include <Teuchos_RCP.hpp>
#include "Tpetra_MapDecl.hpp"
#include <Teuchos_Object.hpp>

namespace Tpetra {

  //! Tpetra::Directory
  
  /*! For Map objects, a Directory object must be created to allow referencing
      of non-local elements. Tpetra::Directory produces and contains a uniform linear
      Map and a list of imageIDs allowing non-local elements to be accessed
      by dereferencing throught the Directory.
      
      This class currently has one constructor, taking an Map object.
  */
  
  template<typename OrdinalType>
  class Directory : public Teuchos::Object {
  public:
    
    //@{ \name Constructors/Destructor.
    
    //! constructor
    Directory(const Map<OrdinalType> & map);
    
 private:
    //! copy constructor
    Directory(const Directory<OrdinalType> & Directory);
 public:
    //! destructor.
    ~Directory();
    
    //@}
    
    //@{ \name Query methods.
    
    //! getDirectoryEntries : Returns image info for non-local Map entries
    /*! Given a list of Global Entry IDs, this function returns the list of
      IDs of the owning memory image that correspond to the list of entries.
      \param In
      globalEntries - List of Global IDs being passed in.
      \param Out
      images - On return contains list of Image IDs owning the Global IDs in question.
    */
    void getDirectoryEntries(const std::vector<OrdinalType> & globalEntries, 
                                   std::vector<OrdinalType>& images) const;
    
    //! getDirectoryEntries : Returns image and local id info for non-local Map entries
    /*! Given a list of Global Entry IDs, this function returns the list of
      image IDs and local IDs on the owning memory image that correspond
      to the list of entries.  If LocalEntries is 0, then local IDs are 
      not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
      element sizes for the requested global entries.
      \param In
      globalEntries - List of Global IDs being passed in.
      \param Out
      images - On return contains list of Image IDs owning the Global IDs in question.
      \param Out
      localEntries - On return contains the local ID of the global on the owning image. 
    */
    void getDirectoryEntries(const std::vector<OrdinalType> & globalEntries, 
                                   std::vector<OrdinalType>& images, 
                                   std::vector<OrdinalType>& localEntries) const;
    //@}
    
  private:
    const Map<OrdinalType> map_;
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm_;
    std::vector<OrdinalType> allMinGIDs_;
    std::vector<OrdinalType> imageIDs_;
    std::vector<OrdinalType> LIDs_;
    Teuchos::RCP< Map<OrdinalType> > directoryMap_;
    
    //! Assignment operator (declared but not defined, do not use)
    Directory<OrdinalType>& operator = (const Directory<OrdinalType> & Source);

    // common code for both versions of getDirectoryEntries
    void getEntries(const std::vector<OrdinalType> & globalEntries, 
                          std::vector<OrdinalType>& images, 
                          std::vector<OrdinalType>& localEntries, 
                          bool computeLIDs) const;
    
    // directory setup for non-contiguous ES
    void generateDirectory();
    
  }; // class Directory
  
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_DECL_HPP

