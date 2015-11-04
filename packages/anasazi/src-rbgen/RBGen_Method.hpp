// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef RBGEN_METHOD_HPP
#define RBGEN_METHOD_HPP

#include "RBGen_FileIOHandler.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

/*! \file RBGen_Method.h
    \brief Provides abstract base class for reduced basis methods.
*/

namespace RBGen {

  //! Abstract base class for reduced basis methods.
  template<class DataSetType, class OperatorType>
  class Method {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    Method() {};


    //! Destructor.
    virtual ~Method() {};
    //@}

    //! @name Computation Methods
    //@{

    //! Compute a basis for the provided snapshots.
    virtual void computeBasis() = 0;

    //! Append new snapshots to the set, and update the basis.
    virtual void updateBasis( const Teuchos::RCP< DataSetType >& update_ss ) = 0;

    //@}

    //! @name Get Methods
    //@{

    //! Get the basis computed by the reduced basis method.
    virtual Teuchos::RCP<const DataSetType> getBasis() const = 0;

    //! Returns the computational time taken to compute the reduced basis.
    virtual double getCompTime() const = 0;

    //@}

    //! @name Set Methods
    //@{

    //! Initialize the method with the given parameter list and snapshot set.
    virtual void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                             const Teuchos::RCP< const DataSetType >& ss,
			     const Teuchos::RCP< RBGen::FileIOHandler< OperatorType > >& fileio = Teuchos::null ) = 0;
    
    //! Reset the snapshot set used to compute the reduced basis.
    virtual void Reset( const Teuchos::RCP< DataSetType >& new_ss ) = 0;

    //@}

    //! @name Status Methods
    //@{

    virtual bool isInitialized() = 0;

    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_METHOD_HPP
