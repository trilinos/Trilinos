// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef MATRIX_MARKET_FILE_IO_HANDLER_H
#define MATRIX_MARKET_FILE_IO_HANDLER_H


#include "RBGen_FileIOHandler.hpp"

// Forward declaration of Epetra_MultiVector class
class Epetra_MultiVector;

namespace RBGen {
  
  //! FileIOHandler for reading an Epetra_MultiVector from a Matrix Market file.
  class MatrixMarketFileIOHandler : public virtual FileIOHandler< Epetra_MultiVector > {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    MatrixMarketFileIOHandler();

    //! Destructor.
    virtual ~MatrixMarketFileIOHandler() {};

    //@}

    //! @name Initialization/Reset Methods
    //@{

    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params );

    void Reset() { num_nodes = 0; isInit = false; };

    //@}

    //! @name File Reading Methods
    //@{

    //! Method for reading multiple files and putting them into an Epetra_MultiVector.
    Teuchos::RCP<Epetra_MultiVector> Read( const std::vector<std::string>& filenames );

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one Epetra_MultiVector into a file.
    void Write( const Teuchos::RCP<const Epetra_MultiVector>& MV, const std::string& filename );

    //@}
    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    bool isInitialized() const { return isInit; };

    //@}

  private:
    // Number of nodes. 
    int num_nodes;

    // Whether or not we know the file format.
    bool isInit;

    // File input / output paths
    std::string in_path, out_path;

    // ParameterList that this file handler was initialized with.
    Teuchos::RCP< Teuchos::ParameterList > params_;

  };
  
} // namespace RBGen

#endif // MATRIX_MARKET_FILE_IO_HANDLER_H
