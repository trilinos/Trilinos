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

#ifndef RBGEN_ISVDMULTISDB_H
#define RBGEN_ISVDMULTISDB_H

#include "RBGen_IncSVDPOD.h"

//
// computeBasis()
//       |
//   makePass()     ___
//   while () {    /    expand()
//      incStep()  ---- SVD()
//   }             \___ shrink()
//
// makePass(), expand() and shrink() are pure virtual
//
// makePass() is implemented in a base class that decides 
// if a method is Multipass, and if so, in what manner
// makePass() puts the data for the next pass into the proper columns
// of U_, then calls incstep (telling it how many new columns)
//
// incstep calls expand to construct expanded U_ and V_, and B_
// it computes the SVD of B_
// it passes this data to shrink(), which shrink according to the 
// base class.
//
// expand(),shrink() are implemented in a base class that 
// decides the representation: UDV, QRW, QBW
// 
//                                   IncSVDPOD
//                                       |
//      -------------------------------------------------------------------
//      |       |        |           |            |            |          |
//  ISVDUDV  ISVDQRW  ISVDQBW   ISVDMultiCD ISVDMultiSDB ISVDMultiSDB ISVDSingle
//      |       |        |           |            |            |          |
//      ------------------           --------------------------------------
//              \                                       /
//               \                                     /
//                \                                   /
//                 \                                 /
//                  \                               /
//                   \---  Concrete Base Class ----/
//
// Then a concrete base class (one of 3x4==12 varieties) is formed simply through
// inheritence. This is the Template Pattern type of Design Pattern.
//

namespace RBGen {

  //! Class for producing a basis using the Incremental SVD in a single pass
  class ISVDMultiSDB : public virtual IncSVDPOD {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    ISVDMultiSDB();

    //! Destructor.
    virtual ~ISVDMultiSDB() {};
    //@}

    //! @name Computation Methods
    //@{

    //! Update the current basis by appending new snapshots.
    void updateBasis( const Teuchos::RCP< Epetra_MultiVector >& update_ss );

    //@}

    //! @name Set Methods
    //@{

    //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                     const Teuchos::RCP< const Epetra_MultiVector >& init,
                     const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio = Teuchos::null );

    //@}

  protected:

    // private member for performing inc steps
    virtual void makePass();

    // will need workspace for 
    // - W12 = [W1 W2], 2*maxBasisSize  (local)
    // - A W, 2*maxBasisSize  (dist)
    // - W12*V, maxBasisSize (local)
    Teuchos::RCP<Epetra_MultiVector> workAW_, workW_, workWV_;
  };

} // end of RBGen namespace

#endif // RBGEN_ISVDMULTISDB_H
