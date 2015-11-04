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


#ifndef RBGEN_BASIC_POD_H
#define RBGEN_BASIC_POD_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

namespace RBGen {
  
  //! Class for producing a basis using LAPACK
  class LapackPOD : public virtual Method<Epetra_MultiVector,Epetra_Operator>, public virtual PODMethod<double> {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    LapackPOD();

    //! Destructor.
    virtual ~LapackPOD() {};
    //@}

    //! @name Computation Methods
    //@{

    void computeBasis();
   
    void updateBasis( const Teuchos::RCP< Epetra_MultiVector >& update_ss ) {};
 
    //@}

    //! @name Get Methods
    //@{
    
    Teuchos::RCP<const Epetra_MultiVector> getBasis() const { return basis_; }

    std::vector<double> getSingularValues() const { return sv_; }

    double getCompTime() const { return comp_time_; }
    //@}
    
    //! @name Set Methods
    //@{
    
        //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                     const Teuchos::RCP< const Epetra_MultiVector >& ss,
		     const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio = Teuchos::null );

    void Reset( const Teuchos::RCP<Epetra_MultiVector>& new_ss ) { ss_ = new_ss;  }

    //@}

    //! @name Status Methods
    //@{

    bool isInitialized() { return isInitialized_; }

    //@}

  private:

    // Is this object initialized.
    bool isInitialized_;

    // Size of the basis that this method will compute.
    int basis_size_;

    // Computational time (using wall clock).
    double comp_time_;

    // Pointers to the snapshots and reduced basis.
    Teuchos::RCP<const Epetra_MultiVector> ss_;
    Teuchos::RCP<Epetra_MultiVector> basis_;

    // Vector holding singular values.
    std::vector<double> sv_;
  };
  
} // end of RBGen namespace

#endif // RBGEN_BASIC_POD_H
