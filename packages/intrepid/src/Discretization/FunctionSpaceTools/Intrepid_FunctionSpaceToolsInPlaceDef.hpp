// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_FunctionSpaceToolsInPlaceDef.hpp
    \brief  Definition file for the Intrepid::FunctionSpaceToolsInPlace class.
    \author Created by R. Kirby
*/


namespace Intrepid {

  template<class Scalar, class ArrayType>
  void FunctionSpaceToolsInPlace::HGRADtransformVALUE(ArrayType & inOutVals )
  {
    return;
  }						      

  template<class Scalar, class ArrayType>
  void FunctionSpaceToolsInPlace::HGRADtransformVALUEDual(ArrayType & inOutVals )
  {
    return;
  }						      

  template<class Scalar, class ArrayType, class ArrayTypeJac>
  void FunctionSpaceToolsInPlace::HGRADtransformGRAD(ArrayType   & inOutVals,
						     const ArrayTypeJac & jacobianInverse,
						     const char           transpose) 
  {
    FieldContainer<Scalar> tmp(inOutVals.dimension(3));

    // test for transpose direction, one loop nest for each direction
    if (transpose == 'T')
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int f=0;f<inOutVals.dimension(1);f++)
	      {
		for (int p=0;p<inOutVals.dimension(2);p++)
		  {
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			tmp(d1) = 0.0;
			for (int d2=0;d2<inOutVals.dimension(3);d2++)
			  {
			    tmp(d1) += jacobianInverse(c,p,d2,d1) * inOutVals(c,f,p,d2);
			  }
		      }
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			inOutVals(c,f,p,d1) = tmp(d1);
		      }
		  }
	      }
	  }
      }
    else if (transpose == 'N')
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int f=0;f<inOutVals.dimension(1);f++)
	      {
		for (int p=0;p<inOutVals.dimension(2);p++)
		  {
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			tmp(d1) = 0.0;
			for (int d2=0;d2<inOutVals.dimension(3);d2++)
			  {
			    tmp(d1) += jacobianInverse(c,p,d1,d2) * inOutVals(c,f,p,d2);
			  }
		      }
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			inOutVals(c,f,p,d1) = tmp(d1);
		      }
		  }
	      }
	  }
      }
    else
      {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				    "Intrepid:FunctionSpaceToolsInPlace::HGRADtransformGRAD::Unknown transpose type" );
      }
  }

  template<class Scalar, class ArrayType, class ArrayTypeJac>
  void FunctionSpaceToolsInPlace::HGRADtransformGRADDual(ArrayType   & inOutVals,
							 const ArrayTypeJac & jacobianInverse,
							 const char           transpose) 
  {
    char t_loc;
    if (transpose == 'T') 
      {
	t_loc = 'N';
      }
    else if (transpose == 'N')
      {
	t_loc = 'T';
      }
    else
      {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				    "Intrepid:FunctionSpaceToolsInPlace::HGRADtransformGRADDual::Unknown transpose type" );
      }
    FunctionSpaceToolsInPlace::HGRADtransformGRAD<Scalar,ArrayType,ArrayTypeJac>( inOutVals,
										  jacobianInverse,
										  t_loc); 
  }


  template<class Scalar, class ArrayType, class ArrayTypeJac>
  void FunctionSpaceToolsInPlace::HCURLtransformVALUE(ArrayType           & inOutVals,
						      const ArrayTypeJac  & jacobianInverse,
						      const char            transpose) 
  {
    FunctionSpaceToolsInPlace::HGRADtransformGRAD<Scalar,ArrayType,ArrayTypeJac>( inOutVals , jacobianInverse, transpose );
  }

  template<class Scalar, class ArrayType, class ArrayTypeJac>
  void FunctionSpaceToolsInPlace::HCURLtransformVALUEDual(ArrayType           & inOutVals,
							  const ArrayTypeJac  & jacobianInverse,
							  const char            transpose) 
  {
    char t_loc;
    if (transpose == 'T') 
      {
	t_loc = 'N';
      }
    else if (transpose == 'N')
      {
	t_loc = 'T';
      }
    else
      {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				    "Intrepid:FunctionSpaceToolsInPlace::HCURLtransformVALUEDual::Unknown transpose type" );
      }
    FunctionSpaceToolsInPlace::HCURLtransformVALUEDual<Scalar,ArrayType,ArrayTypeJac>(inOutVals,
										     jacobianInverse,
										     t_loc); 

  }

  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HCURLtransformCURL(ArrayType  & inOutVals,
						     const ArrayTypeJac  & jacobian,
						     const ArrayTypeDet  & jacobianDet,
						     const char            transpose) 
  {
    FieldContainer<Scalar> tmp(inOutVals.dimension(3));

    // test for transpose direction, one loop nest for each direction
    if (transpose == 'T')
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int f=0;f<inOutVals.dimension(1);f++)
	      {
		for (int p=0;p<inOutVals.dimension(2);p++)
		  {
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			tmp(d1) = 0.0;
			for (int d2=0;d2<inOutVals.dimension(3);d2++)
			  {
			    tmp(d1) += jacobian(c,p,d2,d1) * inOutVals(c,f,p,d2);
			  }
		      }
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			inOutVals(c,f,p,d1) = tmp(d1) / jacobianDet(c,p);
		      }
		  }
	      }
	  }
      }
    else if (transpose == 'N')
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int f=0;f<inOutVals.dimension(1);f++)
	      {
		for (int p=0;p<inOutVals.dimension(2);p++)
		  {
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			tmp(d1) = 0.0;
			for (int d2=0;d2<inOutVals.dimension(3);d2++)
			  {
			    tmp(d1) += jacobian(c,p,d1,d2) * inOutVals(c,f,p,d2);
			  }
		      }
		    for (int d1=0;d1<inOutVals.dimension(3);d1++)
		      {
			inOutVals(c,f,p,d1) = tmp(d1) / jacobianDet(c,p);
		      }
		  }
	      }
	  }
      }
    else
      {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				    "Intrepid:FunctionSpaceToolsInPlace::HCURLtransformCURL::Unknown transpose type" );
      }

  }

  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HCURLtransformCURLDual(ArrayType  & inOutVals,
							 const ArrayTypeJac  & jacobian,
							 const ArrayTypeDet  & jacobianDet,
							 const char            transpose) 
  {
    char t_loc;
    if (transpose == 'T') 
      {
	t_loc = 'N';
      }
    else if (transpose == 'N')
      {
	t_loc = 'T';
      }
    else
      {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				    "Intrepid:FunctionSpaceToolsInPlace::HCURLtransformCURLDual::Unknown transpose type" );
      }
    FunctionSpaceToolsInPlace::HCURLtransformCURLDual<Scalar,ArrayType,ArrayTypeJac>(inOutVals,
										    jacobian,
										    jacobianDet,
										    t_loc); 
  }

  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HDIVtransformVALUE(ArrayType           & inOutVals,
						     const ArrayTypeJac  & jacobian,
						     const ArrayTypeDet  & jacobianDet,
						     const char            transpose) 
  {
    FunctionSpaceToolsInPlace::HCURLtransformCURL<Scalar,ArrayType,ArrayTypeJac,ArrayTypeDet>( inOutVals,jacobian,jacobianDet,transpose);
  }

  template<class Scalar, class ArrayType, class ArrayTypeJac, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HDIVtransformVALUEDual(ArrayType           & inOutVals,
							 const ArrayTypeJac  & jacobian,
							 const ArrayTypeDet  & jacobianDet,
							 const char            transpose) 
  {
    char t_loc;
    if (transpose == 'T')
      {
	t_loc = 'N';
      }
    else if (transpose == 'N')
      {
	t_loc = 'T';
      }
    else
      {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				    "FunctionSpaceToolsInPlace::HDIVtransformVALUEDual: invalid transpose character");
      }
    FunctionSpaceToolsInPlace::HDIVtransformVALUE<Scalar,ArrayType,ArrayTypeJac,ArrayTypeDet>( inOutVals,
											       jacobian,
											       jacobianDet ,
											       t_loc );
  }

  template<class Scalar, class ArrayType, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HDIVtransformDIV(ArrayType           & inOutVals,
						   const ArrayTypeDet  & jacobianDet)
  {
    for (int c=0;c<inOutVals.dimension(0);c++)
      {
	for (int f=0;f<inOutVals.dimension(1);f++)
	  {
	    for (int p=0;p<inOutVals.dimension(2);p++)
	      {
		inOutVals(c,f,p) /= jacobianDet(c,p);
	      }
	  }
      }
  }

  template<class Scalar, class ArrayType, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HDIVtransformDIVDual(ArrayType           & inOutVals,
						       const ArrayTypeDet  & jacobianDet)
  {
    FunctionSpaceToolsInPlace::HDIVtransformDIV<Scalar,ArrayType,ArrayTypeDet>( inOutVals ,
										jacobianDet );
  }

  template<class Scalar, class ArrayType, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HVOLtransformVALUE(ArrayType           & inOutVals,
						     const ArrayTypeDet  & jacobianDet)
  {
    FunctionSpaceToolsInPlace::HDIVtransformDIV<Scalar,ArrayType,ArrayTypeDet>( inOutVals,jacobianDet);
  }

  template<class Scalar, class ArrayType, class ArrayTypeDet>
  void FunctionSpaceToolsInPlace::HVOLtransformVALUEDual(ArrayType           & inOutVals,
							 const ArrayTypeDet  & jacobianDet)
  {
    FunctionSpaceToolsInPlace::HVOLtransformVALUEDual( inOutVals ,jacobianDet );
  }



  template<class Scalar, class ArrayType, class ArrayTypeMeasure>
  void FunctionSpaceToolsInPlace::multiplyMeasure(ArrayType       & inOutVals,
						  const ArrayTypeMeasure   & inMeasure)
  {
    if (inOutVals.rank() == 2)  // inOutVals is (C,P)
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int p=0;p<inOutVals.dimension(0);p++)
	      {
		inOutVals(c,p) *= inMeasure(c,p);
	      }
	  }
      }
    else if (inOutVals.rank() == 3) // inOutVals is (C,F,P)
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int f=0;f<inOutVals.dimension(1);f++)
	      {
		for (int p=0;p<inOutVals.dimension(2);p++)
		  {
		    inOutVals(c,f,p) *= inMeasure(c,p);
		  }
	      }
	  }
      }
    else if (inOutVals.rank() == 4) // inOutVals is (C,F,P)
      {
	for (int c=0;c<inOutVals.dimension(0);c++)
	  {
	    for (int f=0;f<inOutVals.dimension(1);f++)
	      {
		for (int p=0;p<inOutVals.dimension(2);p++)
		  {
		    for (int d=0;d<inOutVals.dimension(3);d++)
		      {
			inOutVals(c,f,p,d) *= inMeasure(c,p);
		      }
		  }
	      }
	  }
      }
  } 

} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

