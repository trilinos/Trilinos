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

#ifndef TPETRA_IMPORTDATA_HPP
#define TPETRA_IMPORTDATA_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Object.hpp>
#include "Tpetra_Import.hpp"

namespace Tpetra {

  template<typename Ordinal>
  class ImportData : public Teuchos::Object {
    friend class Import<Ordinal>;
  public:
    // default constructor
    ImportData(const Map<Ordinal> & source, const Map<Ordinal> & target);

    // destructor. no heap-data, so no need to override
    ~ImportData();

  protected:
    // OT vectors
    std::vector<Ordinal> permuteToLIDs_;
    std::vector<Ordinal> permuteFromLIDs_;
    std::vector<Ordinal> remoteLIDs_;
    std::vector<Ordinal> remoteGIDs_;
    std::vector<Ordinal> exportLIDs_;
    std::vector<Ordinal> exportImageIDs_;

    // OTs
    Ordinal numSameIDs_;
    Ordinal numPermuteIDs_;
    Ordinal numRemoteIDs_;
    Ordinal numExportIDs_;

    // Maps
    const Map<Ordinal> source_;
    const Map<Ordinal> target_;

    // Platform, Comm, Distributor, etc.
    Teuchos::RCP<const Platform<Ordinal> > platform_;
    Distributor<Ordinal> distributor_;

  private:
    //! Copy constructor (declared but not defined, do not use)
    ImportData(ImportData<Ordinal> const& rhs);
    //! Assignment operator (declared but not defined, do not use)
    ImportData<Ordinal>& operator = (ImportData<Ordinal> const& rhs);

  }; // class ImportData

  template <typename Ordinal>
  ImportData<Ordinal>::ImportData(const Map<Ordinal> & source, const Map<Ordinal> & target)
  : Teuchos::Object("Tpetra::ImportData")
  , numSameIDs_(Teuchos::OrdinalTraits<Ordinal>::zero())
  , numPermuteIDs_(Teuchos::OrdinalTraits<Ordinal>::zero())
  , numRemoteIDs_(Teuchos::OrdinalTraits<Ordinal>::zero())
  , numExportIDs_(Teuchos::OrdinalTraits<Ordinal>::zero())
  , source_(source)
  , target_(target)
  , platform_(source.getPlatform()->clone())
  , distributor_(platform_->createComm())
  {}

  template <typename Ordinal>
  ImportData<Ordinal>::~ImportData() 
  {}

} // namespace Tpetra

#endif // TPETRA_IMPORTDATA_HPP
