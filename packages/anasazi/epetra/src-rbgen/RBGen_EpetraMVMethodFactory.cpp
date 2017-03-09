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

