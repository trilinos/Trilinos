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

#include "RBGen_EpetraMVMethodFactory.h"
#include "RBGen_LapackPOD.h"
#include "Teuchos_Assert.hpp"

#include "RBGen_AnasaziPOD.h"
#include "RBGen_ISVD_SingleUDV.h"
#include "RBGen_ISVD_MultiCDUDV.h"
#include "RBGen_ISVD_MultiSDAUDV.h"
#include "RBGen_ISVD_MultiSDBUDV.h"
#include "RBGen_StSVD_RTR.h"

namespace RBGen {
  
  Teuchos::RCP<Method<Epetra_MultiVector,Epetra_Operator> > EpetraMVMethodFactory::create( const Teuchos::ParameterList& params )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist( "Reduced Basis Method" ), std::invalid_argument, "Reduced Basis Method sublist does not exist!");

    // Get the "Reduced Basis Method" sublist.
    const Teuchos::ParameterList& rbmethod_params = params.sublist( "Reduced Basis Method" );

    // Get the file format type
    std::string method = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(rbmethod_params),
                                                             "Method" );

    Teuchos::RCP< Method<Epetra_MultiVector,Epetra_Operator> > RBMethod;

    // POD computed using exact SVD through LAPACK
    if ( method == "Lapack POD" ) {
      RBMethod = Teuchos::rcp( new LapackPOD() );
    } 
    // Inexact POD computed using inexact SVD through Anasazi
    // IncSVDPOD uses Anasazi utility classes, while AnasaziPOD uses Anasazi for the solution
    else if ( method == "IncSVD POD" ) {
      std::string incsvdmethod = rbmethod_params.get<std::string>("IncSVD Method");
      if ( incsvdmethod == "Single/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_SingleUDV() );
      }
      else if ( incsvdmethod == "MultiCD/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_MultiCDUDV() );
      }
      else if ( incsvdmethod == "MultiSDA/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_MultiSDAUDV() );
      }
      else if ( incsvdmethod == "MultiSDB/UDV" ) {
        RBMethod = Teuchos::rcp( new ISVD_MultiSDBUDV() );
      }
    } 
    else if ( method == "StSVD/RTR") {
      RBMethod = Teuchos::rcp( new StSVDRTR() );
    }
    else if ( method == "Anasazi POD" ) {
      RBMethod = Teuchos::rcp( new AnasaziPOD() );
    } else 
    {
      std::string err_str = "Reduced basis method, 'Method = " + method + "', is not recognized!";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, err_str);
    }
    //
    // Return the method created
    //
    return RBMethod;
  }  

} // end of RBGen namespace

