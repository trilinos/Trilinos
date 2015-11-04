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

#ifndef RBGEN_ISVD_MULTISDAUDV_H
#define RBGEN_ISVD_MULTISDAUDV_H

#include "RBGen_ISVDMultiSDA.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  //! IncSVD method implementing UDV storage with multiple steepest descent (variant A) passes.
  class ISVD_MultiSDAUDV : public virtual ISVDUDV, public virtual ISVDMultiSDA {
    public:
      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_MultiSDAUDV();

      //! Destructor.
      virtual ~ISVD_MultiSDAUDV() {};
      //@}

      //! @name Set Methods
      //@{

      //! Initialize the method with the given parameter list and snapshot set.
      void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
          const Teuchos::RCP< const Epetra_MultiVector >& init,
          const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio = Teuchos::null );

      //@}
  };

} // end of RBGen namespace

#endif // RBGEN_ISVD_MULTISDAUDV_H
