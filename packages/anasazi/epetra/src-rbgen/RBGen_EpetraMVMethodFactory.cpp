// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

