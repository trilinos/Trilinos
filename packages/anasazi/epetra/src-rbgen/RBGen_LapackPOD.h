// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
