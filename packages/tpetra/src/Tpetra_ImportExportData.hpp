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

#ifndef TPETRA_IMPORTEXPORTDATA_HPP
#define TPETRA_IMPORTEXPORTDATA_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Object.hpp>

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Import,Export, needed to prevent circular inclusions
  template<typename Ordinal> class Import;
  template<typename Ordinal> class Export;
#endif

  template<typename Ordinal>
  class ImportExportData : public Teuchos::Object {
    friend class Import<Ordinal>;
    friend class Export<Ordinal>;
  public:
    ImportExportData(const Map<Ordinal> & source, const Map<Ordinal> & target);
    ~ImportExportData();

  protected:
    // OT vectors
    Teuchos::Array<Ordinal> permuteToLIDs_;
    Teuchos::Array<Ordinal> permuteFromLIDs_;
    Teuchos::Array<Ordinal> remoteLIDs_;
    Teuchos::Array<Ordinal> exportGIDs_;
    // These are ArrayRCP because in the construction of an Import object, they are allocated and returned by a call to 
    Teuchos::ArrayRCP<Ordinal> exportLIDs_;
    Teuchos::ArrayRCP<Ordinal> exportImageIDs_;

    // OTs
    Teuchos_Ordinal numSameIDs_;

    // Maps
    const Map<Ordinal> source_;
    const Map<Ordinal> target_;

    // Platform, Comm, Distributor, etc.
    Teuchos::RCP<const Platform<Ordinal> > platform_;
    Distributor<Ordinal> distributor_;

  private:
    //! Copy constructor (declared but not defined, do not use)
    ImportExportData(ImportExportData<Ordinal> const& rhs);
    //! Assignment operator (declared but not defined, do not use)
    ImportExportData<Ordinal>& operator = (ImportExportData<Ordinal> const& rhs);
  }; // class ImportExportData


  template <typename Ordinal>
  ImportExportData<Ordinal>::ImportExportData(const Map<Ordinal> & source, const Map<Ordinal> & target)
  : Teuchos::Object("Tpetra::ImportExportData")
  , numSameIDs_(0)
  , source_(source)
  , target_(target)
  , platform_(source.getPlatform()->clone())
  , distributor_(platform_->createComm())
  {}

  template <typename Ordinal>
  ImportExportData<Ordinal>::~ImportExportData() 
  {}

} // namespace Tpetra

#endif // TPETRA_IMPORTEXPORTDATA_HPP
