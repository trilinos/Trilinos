// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
