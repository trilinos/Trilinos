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

/** \file   Intrepid_RealSpaceToolsDef.hpp
    \brief  Definition file for utility classes providing basic linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


namespace Intrepid {



template<class Scalar>
void RealSpaceTools<Scalar>::absval(Scalar* absArray, const Scalar* inArray, const int size) {
  for (size_t i=0; i<size; i++) {
    absArray[i] = std::abs(inArray[i]);
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::absval(Scalar* inoutAbsArray, const int size) {
  for (size_t i=0; i<size; i++) {
    inoutAbsArray[i] = std::abs(inoutAbsArray[i]);
  }
}



template<class Scalar>
template<class ArrayAbs, class ArrayIn>
void RealSpaceTools<Scalar>::absval(ArrayAbs & absArray, const ArrayIn & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( ( getrank(inArray) != getrank(absArray) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::absval): Array arguments must have identical ranks!");
    for (size_t i=0; i<getrank(inArray); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inArray.dimension(i)) != static_cast<size_t>(absArray.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::absval): Dimensions of array arguments do not agree!");
    }
#endif 
  
   ArrayWrapper<Scalar,ArrayAbs, Rank<ArrayAbs >::value, false>absArrayWrap(absArray);
   ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inArrayWrap(inArray);
 
   int inArrayRank=getrank(inArray);


   if(inArrayRank==5){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++)
          for (size_t m=0; m<static_cast<size_t>(static_cast<size_t>(inArray.dimension(4))); m++){
         absArrayWrap(i,j,k,l,m) = std::abs(inArrayWrap(i,j,k,l,m));
          }
	}else if(inArrayRank==4){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++){
            absArrayWrap(i,j,k,l) = std::abs(inArrayWrap(i,j,k,l));
          }
	}else if(inArrayRank==3){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++){
         absArrayWrap(i,j,k) = std::abs(inArrayWrap(i,j,k));
          }
	}else if(inArrayRank==2){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++){
         absArrayWrap(i,j) = std::abs(inArrayWrap(i,j));
          }
	}else if(inArrayRank==1){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++){
        absArrayWrap(i) = std::abs(inArrayWrap(i));

          }
	}  
}



template<class Scalar>
template<class ArrayInOut>
void RealSpaceTools<Scalar>::absval(ArrayInOut & inoutAbsArray) {
  for (size_t i=0; i<(size_t)inoutAbsArray.size(); i++) {
    inoutAbsArray[i] = std::abs(inoutAbsArray[i]);
  }
}



template<class Scalar>
Scalar RealSpaceTools<Scalar>::vectorNorm(const Scalar* inVec, const size_t dim, const ENorm normType) {
  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(size_t i = 0; i < dim; i++){
        temp += inVec[i]*inVec[i];
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(inVec[0]);
      for(size_t i = 1; i < dim; i++){
        Scalar absData = std::abs(inVec[i]);
        if (temp < absData) temp = absData;
      }
      break;
    case NORM_ONE:
      for(size_t i = 0; i < dim; i++){
        temp += std::abs(inVec[i]);
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}

template<class Scalar>
template<class ArrayIn>
Scalar RealSpaceTools<Scalar>::vectorNorm(const ArrayIn & inVec, const ENorm normType) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( !(getrank(inVec) >= 1 && getrank(inVec) <= 5) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::vectorNorm): Vector argument must have rank 1!");
#endif
  ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inVecWrap(inVec);
  int inVecRank=getrank(inVecWrap);
  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:{
   if(inVecRank==5){ 
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inVec.dimension(3)); l++)
          for (size_t m=0; m<static_cast<size_t>(inVec.dimension(4)); m++)
      temp += inVecWrap(i,j,k,l,m)*inVecWrap(i,j,k,l,m);
     }else if(inVecRank==4){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inVec.dimension(3)); l++)
      temp += inVecWrap(i,j,k,l)*inVecWrap(i,j,k,l); 	 
	 }else if(inVecRank==3){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
      temp += inVecWrap(i,j,k)*inVecWrap(i,j,k); 	 
	 }else if(inVecRank==2){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      temp += inVecWrap(i,j)*inVecWrap(i,j); 	 
	 }else if(inVecRank==1){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
      temp += inVecWrap(i)*inVecWrap(i); 	 
	 }         
      temp = std::sqrt(temp);
  }
      break;
    case NORM_INF:{

     if(inVecRank==5){
   temp = std::abs(inVecWrap(0,0,0,0,0));
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inVec.dimension(3)); l++)
          for (size_t m=1; m<static_cast<size_t>(inVec.dimension(4)); m++){
         Scalar absData = std::abs(inVecWrap(i,j,k,l,m));
         if (temp < absData) temp = absData;
	    }
	}else if(inVecRank==4){
 temp = std::abs(inVecWrap(0,0,0,0));		
  for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
        for (size_t l=1; l<static_cast<size_t>(inVec.dimension(3)); l++){
         Scalar absData = std::abs(inVecWrap(i,j,k,l));
         if (temp < absData) temp = absData;
	    }	
	}else if(inVecRank==3){
  temp = std::abs(inVecWrap(0,0,0));		
  for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++){
         Scalar absData = std::abs(inVecWrap(i,j,k));
         if (temp < absData) temp = absData;
	    }	
	}else if(inVecRank==2){
  temp = std::abs(inVecWrap(0,0));		
  for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++){
         Scalar absData = std::abs(inVecWrap(i,j));
         if (temp < absData) temp = absData;
	    }	
	}else if(inVecRank==1){
  temp = std::abs(inVecWrap(0));		
  for (size_t i=1; i<static_cast<size_t>(inVec.dimension(0)); i++){
         Scalar absData = std::abs(inVecWrap(i));
         if (temp < absData) temp = absData;
	    }	
	}
}  
      break;
    case NORM_ONE:{
        if(inVecRank==5){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inVec.dimension(3)); l++)
          for (size_t m=0; m<static_cast<size_t>(inVec.dimension(4)); m++){
          temp += std::abs(inVecWrap(i,j,k,l,m));
          }
	}else if(inVecRank==4){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inVec.dimension(3)); l++){
          temp += std::abs(inVecWrap(i,j,k,l));
          }
	}else if(inVecRank==3){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inVec.dimension(2)); k++){
          temp += std::abs(inVecWrap(i,j,k));
          }
	}else if(inVecRank==2){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inVec.dimension(1)); j++){
          temp += std::abs(inVecWrap(i,j));
          }
	}else if(inVecRank==1){
   for (size_t i=0; i<static_cast<size_t>(inVec.dimension(0)); i++){
          temp += std::abs(inVecWrap(i));
          }
	}	
}
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}
/*
template<class Scalar>
template<class ArrayIn>
Scalar RealSpaceTools<Scalar>::vectorNorm(const ArrayIn & inVec, const ENorm normType) {

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( ( inVec.rank() != 1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Vector argument must have rank 1!");
#endif

  int size = inVec.size();

  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(size_t i = 0; i < size; i++){
        temp += inVec[i]*inVec[i];
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(inVec[0]);
      for(size_t i = 1; i < size; i++){
        Scalar absData = std::abs(inVec[i]);
        if (temp < absData) temp = absData;
      }
      break;
    case NORM_ONE:
      for(size_t i = 0; i < size; i++){
        temp += std::abs(inVec[i]);
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}*/
template<class Scalar>
template<class ArrayNorm, class ArrayIn>
void RealSpaceTools<Scalar>::vectorNorm(ArrayNorm & normArray, const ArrayIn & inVecs, const ENorm normType) {

  ArrayWrapper<Scalar,ArrayNorm, Rank<ArrayNorm >::value, false>normArrayWrap(normArray);
  ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inVecsWrap(inVecs);

  size_t arrayRank = getrank(inVecs);
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( arrayRank != getrank(normArray)+1 ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::vectorNorm): Ranks of norm and vector array arguments are incompatible!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::vectorNorm): Rank of vector array must be 2 or 3!");
  for (size_t i=0; i<arrayRank-1; i++) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( inVecs.dimension(i) != normArray.dimension(i) ),
				std::invalid_argument,
				">>> ERROR (RealSpaceTools::vectorNorm): Dimensions of norm and vector arguments do not agree!");
  }
#endif

  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = static_cast<size_t>(inVecs.dimension(arrayRank-1)); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 3:
      dim_i0 = static_cast<size_t>(inVecs.dimension(0));
      dim_i1 = static_cast<size_t>(inVecs.dimension(1));
        switch(normType) {
    case NORM_TWO: {
      for (size_t i0=0; i0<dim_i0; i0++) {
        for (size_t i1=0; i1<dim_i1; i1++) {
          Scalar temp = (Scalar)0;
          for(size_t i = 0; i < dim; i++){
            temp += inVecsWrap(i0,i1,i)*inVecsWrap(i0,i1,i);
          }
          normArrayWrap(i0,i1) = std::sqrt(temp);
        }
      }
      break;
    } // case NORM_TWO

    case NORM_INF: {
      for (size_t i0=0; i0<dim_i0; i0++) {
        for (size_t i1=0; i1<dim_i1; i1++) {
          Scalar temp = (Scalar)0;
          temp = std::abs(inVecsWrap(i0,i1,0));
          for(size_t i = 1; i < dim; i++){
            Scalar absData = std::abs(inVecsWrap(i0,i1,i));
            if (temp < absData) temp = absData;
          }
          normArrayWrap(i0,i1) = temp;
        }
      }
      break;
    } // case NORM_INF

    case NORM_ONE: {
      for (size_t i0=0; i0<dim_i0; i0++) {
        for (size_t i1=0; i1<dim_i1; i1++) {
          Scalar temp = (Scalar)0;
          for(size_t i = 0; i < dim; i++){
            temp += std::abs(inVecsWrap(i0,i1,i));
          }
          normArrayWrap(i0,i1) = temp;
        }
      }
      break;
    } // case NORM_ONE

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
      
      
      
      break;
    case 2:
      dim_i1 = static_cast<size_t>(inVecs.dimension(0));
        switch(normType) {
    case NORM_TWO: {

        for (size_t i1=0; i1<dim_i1; i1++) {
          Scalar temp = (Scalar)0;
          for(size_t i = 0; i < dim; i++){
            temp += inVecsWrap(i1,i)*inVecsWrap(i1,i);
          }
          normArrayWrap(i1) = std::sqrt(temp);
        }
      
      break;
    } // case NORM_TWO

    case NORM_INF: {
        for (size_t i1=0; i1<dim_i1; i1++) {
          Scalar temp = (Scalar)0;
          temp = std::abs(inVecsWrap(i1,0));
          for(size_t i = 1; i < dim; i++){
            Scalar absData = std::abs(inVecsWrap(i1,i));
            if (temp < absData) temp = absData;
          }
          normArrayWrap(i1) = temp;
        }
      break;
    } // case NORM_INF

    case NORM_ONE: {
        for (size_t i1=0; i1<dim_i1; i1++) {
          Scalar temp = (Scalar)0;
          for(size_t i = 0; i < dim; i++){
            temp += std::abs(inVecsWrap(i1,i));
          }
          normArrayWrap(i1) = temp;
        }
      break;
    } // case NORM_ONE

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
      
      
      
      break;
  }


}
/*

template<class Scalar>
template<class ArrayNorm, class ArrayIn>
void RealSpaceTools<Scalar>::vectorNorm(ArrayNorm & normArray, const ArrayIn & inVecs, const ENorm normType) {

  int arrayRank = inVecs.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( ( arrayRank != normArray.rank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Ranks of norm and vector array arguments are incompatible!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Rank of vector array must be 2 or 3!");
    for (int i=0; i<arrayRank-1; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( inVecs.dimension(i) != normArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Dimensions of norm and vector arguments do not agree!");
    }
#endif

  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = inVecs.dimension(arrayRank-1); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 3:
      dim_i0 = inVecs.dimension(0);
      dim_i1 = inVecs.dimension(1);
      break;
    case 2:
      dim_i1 = inVecs.dimension(0);
      break;
  }

  switch(normType) {
    case NORM_TWO: {
      int offset_i0, offset, normOffset;
      for (size_t i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (size_t i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          for(size_t i = 0; i < dim; i++){
            temp += inVecs[offset+i]*inVecs[offset+i];
          }
          normArray[normOffset] = std::sqrt(temp);
        }
      }
      break;
    } // case NORM_TWO

    case NORM_INF: {
      int offset_i0, offset, normOffset;
      for (size_t i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (size_t i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          temp = std::abs(inVecs[offset]);
          for(size_t i = 1; i < dim; i++){
            Scalar absData = std::abs(inVecs[offset+i]);
            if (temp < absData) temp = absData;
          }
          normArray[normOffset] = temp;
        }
      }
      break;
    } // case NORM_INF

    case NORM_ONE: {
      int offset_i0, offset, normOffset;
      for (size_t i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (size_t i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          for(size_t i = 0; i < dim; i++){
            temp += std::abs(inVecs[offset+i]);
          }
          normArray[normOffset] = temp;
        }
      }
      break;
    } // case NORM_ONE

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
}


*/

template<class Scalar>
void RealSpaceTools<Scalar>::transpose(Scalar* transposeMat, const Scalar* inMat, const size_t dim) {
  for(size_t i=0; i < dim; i++){
    transposeMat[i*dim+i]=inMat[i*dim+i];    // Set diagonal elements
    for(size_t j=i+1; j < dim; j++){
      transposeMat[i*dim+j]=inMat[j*dim+i];  // Set off-diagonal elements
      transposeMat[j*dim+i]=inMat[i*dim+j];
    }
  }
}



template<class Scalar>
template<class ArrayTranspose, class ArrayIn>
void RealSpaceTools<Scalar>::transpose(ArrayTranspose & transposeMats, const ArrayIn & inMats) {
  size_t arrayRank = getrank(inMats);
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( arrayRank != getrank(transposeMats) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::transpose): Matrix array arguments do not have identical ranks!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 4) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::transpose): Rank of matrix array must be 2, 3, or 4!");
  for (size_t i=0; i<arrayRank; i++) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(i)) != static_cast<size_t>(transposeMats.dimension(i)) ),
				std::invalid_argument,
				">>> ERROR (RealSpaceTools::transpose): Dimensions of matrix arguments do not agree!");
  }
  TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(arrayRank-2)) != static_cast<size_t>(inMats.dimension(arrayRank-1)) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::transpose): Matrices are not square!");
#endif
  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = static_cast<size_t>(inMats.dimension(arrayRank-2)); // spatial dimension


		ArrayWrapper<Scalar,ArrayTranspose,Rank<ArrayTranspose>::value,false>transposeArr(transposeMats);
		ArrayWrapper<Scalar,ArrayIn,Rank<ArrayIn>::value,true>inputArr(inMats);
  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = static_cast<size_t>(inMats.dimension(0));
      dim_i1 = static_cast<size_t>(inMats.dimension(1));
      
       for (size_t i0=0; i0<dim_i0; i0++) {
    for (size_t i1=0; i1<dim_i1; i1++) {
      for(size_t i=0; i < dim; i++){
		transposeArr(i0,i1,i,i)=inputArr(i0,i1,i,i);         
        for(size_t j=i+1; j < dim; j++){
		  transposeArr(i0,i1,i,j)=inputArr(i0,i1,j,i);	
		  transposeArr(i0,i1,j,i)=inputArr(i0,i1,i,j);	  
        }
      }

    } // i1
  } // i0
      break;
    case 3:
      dim_i1 = static_cast<size_t>(inMats.dimension(0));
    for (size_t i1=0; i1<dim_i1; i1++) {
      for(size_t i=0; i < dim; i++){
		transposeArr(i1,i,i)=inputArr(i1,i,i);         
        for(size_t j=i+1; j < dim; j++){
		  transposeArr(i1,i,j)=inputArr(i1,j,i);	
		  transposeArr(i1,j,i)=inputArr(i1,i,j);
        }
    } // i1
 } 
      break;
  }



}

template<class Scalar>
void RealSpaceTools<Scalar>::inverse(Scalar* inverseMat, const Scalar* inMat, const size_t dim) {

  switch(dim) {
    case 3: {
      size_t i, j, rowID = 0, colID = 0;
      int rowperm[3]={0,1,2};
      int colperm[3]={0,1,2}; // Complete pivoting
      Scalar emax(0);

      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( inMat[i*3+j] ) >  emax){
            rowID = i;  colID = j; emax = std::abs( inMat[i*3+j] );
          }
        }
      }
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
      if( rowID ){
        rowperm[0] = rowID;
        rowperm[rowID] = 0;
      }
      if( colID ){
        colperm[0] = colID;
        colperm[colID] = 0;
      }
      Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          B[i][j] = inMat[rowperm[i]*3+colperm[j]];
        }
      }
      B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
      for(i=0; i < 2; ++i){
        for(j=0; j < 2; ++j){
          S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
        }
      }
      Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif

      Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
      Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

      for(j=0; j<2;j++)
        Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
      for(i=0; i<2;i++)
        Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

      Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
      Bi[1][1] =  Si[0][0];
      Bi[1][2] =  Si[0][1];
      Bi[2][1] =  Si[1][0];
      Bi[2][2] =  Si[1][1];
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          inverseMat[i*3+j] = Bi[colperm[i]][rowperm[j]]; // set inverse
        }
      }
      break;
    } // case 3

    case 2: {

      Scalar determinant    = inMat[0]*inMat[3]-inMat[1]*inMat[2];;
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( (inMat[0]==(Scalar)0) && (inMat[1]==(Scalar)0) &&
                            (inMat[2]==(Scalar)0) && (inMat[3]==(Scalar)0) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
      TEUCHOS_TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
      inverseMat[0] =   inMat[3] / determinant;
      inverseMat[1] = - inMat[1] / determinant;
      //
      inverseMat[2] = - inMat[2] / determinant;
      inverseMat[3] =   inMat[0] / determinant;
      break;
    } // case 2

    case 1: {
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( inMat[0] == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
      inverseMat[0] = (Scalar)1 / inMat[0];
      break;
    } // case 1

  } // switch (dim)
}

template<class Scalar>
template<class ArrayInverse, class ArrayIn>
void RealSpaceTools<Scalar>::inverse(ArrayInverse & inverseMats, const ArrayIn & inMats) {
	
 ArrayWrapper<Scalar,ArrayInverse, Rank<ArrayInverse >::value, false>inverseMatsWrap(inverseMats);
 ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inMatsWrap(inMats);

  size_t arrayRank = getrank(inMats);

#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( arrayRank != getrank(inverseMats) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::inverse): Matrix array arguments do not have identical ranks!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 4) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::inverse): Rank of matrix array must be 2, 3, or 4!");
  for (size_t i=0; i<arrayRank; i++) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(i)) != static_cast<size_t>(inverseMats.dimension(i)) ),
				std::invalid_argument,
				">>> ERROR (RealSpaceTools::inverse): Dimensions of matrix arguments do not agree!");
  }
  TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(arrayRank-2)) != static_cast<size_t>(inMats.dimension(arrayRank-1)) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::inverse): Matrices are not square!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (static_cast<size_t>(inMats.dimension(arrayRank-2)) < 1) || (static_cast<size_t>(inMats.dimension(arrayRank-2)) > 3) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::inverse): Spatial dimension must be 1, 2, or 3!");
#endif

  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = static_cast<size_t>(inMats.dimension(arrayRank-2)); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = static_cast<size_t>(inMats.dimension(0));
      dim_i1 = static_cast<size_t>(inMats.dimension(1));
       switch(dim) {
    case 3: {
     

      for (size_t i0=0; i0<dim_i0; i0++) {
  
        for (size_t i1=0; i1<dim_i1; i1++) {
         

          size_t i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i0,i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i0,i1,i,j) );
              }
            }
          }
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
				      std::invalid_argument,
				      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif
          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMatsWrap(i0,i1,rowperm[i],colperm[j]);
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
				      std::invalid_argument,
				      ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif

          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMatsWrap(i0,i1,i,j) = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
        } // for i1
      } // for i0
      break;
    } // case 3

    case 2: {

      for (size_t i0=0; i0<dim_i0; i0++) {

        for (size_t i1=0; i1<dim_i1; i1++) {
 

          Scalar determinant    = inMatsWrap(i0,i1,0,0)*inMatsWrap(i0,i1,1,1)-inMatsWrap(i0,i1,0,1)*inMatsWrap(i0,i1,1,0);
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( (inMatsWrap(i0,i1,0,0)==(Scalar)0)   && (inMatsWrap(i0,i1,0,1)==(Scalar)0) &&
					(inMatsWrap(i0,i1,1,0)==(Scalar)0) && (inMatsWrap(i0,i1,1,1)==(Scalar)0) ),
				      std::invalid_argument,
				      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
          TEUCHOS_TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
				      std::invalid_argument,
				      ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif
          inverseMatsWrap(i0,i1,0,0)   = inMatsWrap(i0,i1,1,1) / determinant;
          inverseMatsWrap(i0,i1,0,1) = - inMatsWrap(i0,i1,0,1) / determinant;
          //
          inverseMatsWrap(i0,i1,1,0) = - inMatsWrap(i0,i1,1,0) / determinant;
          inverseMatsWrap(i0,i1,1,1) =   inMatsWrap(i0,i1,0,0) / determinant;
        } // for i1
      } // for i0
      break;
    } // case 2

    case 1: {
       for (size_t i0=0; i0<dim_i0; i0++) {
        for (size_t i1=0; i1<dim_i1; i1++) {
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( inMatsWrap(i0,i1,0,0) == (Scalar)0 ),
				      std::invalid_argument,
				      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif

          inverseMatsWrap(i0,i1,0,0) = (Scalar)1 / inMatsWrap(i0,i1,0,0);
        } // for i1
      } // for i2
          
      
      break;
    } // case 1

  } // switch (dim)	
      break;
    case 3:
      dim_i1 = static_cast<size_t>(inMats.dimension(0));
       switch(dim) {
    case 3: {

        for (size_t i1=0; i1<dim_i1; i1++) {

          size_t i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i1,i,j) );
              }
            }
          }
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif

          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMatsWrap(i1,rowperm[i],colperm[j]);
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
				      std::invalid_argument,
				      ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif


          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMatsWrap(i1,i,j) = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
        } // for i1
      
      break;
    } // case 3

    case 2: {

        for (size_t i1=0; i1<dim_i1; i1++) {
         

          Scalar determinant    = inMatsWrap(i1,0,0)*inMatsWrap(i1,1,1)-inMatsWrap(i1,0,1)*inMatsWrap(i1,1,0);
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( (inMatsWrap(i1,0,0)==(Scalar)0)   && (inMatsWrap(i1,0,1)==(Scalar)0) &&
                                        (inMatsWrap(i1,1,0)==(Scalar)0) && (inMatsWrap(i1,1,1)==(Scalar)0) ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
          TEUCHOS_TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif

          inverseMatsWrap(i1,0,0)   = inMatsWrap(i1,1,1) / determinant;
          inverseMatsWrap(i1,0,1) = - inMatsWrap(i1,0,1) / determinant;
          //
          inverseMatsWrap(i1,1,0) = - inMatsWrap(i1,1,0) / determinant;
          inverseMatsWrap(i1,1,1) =   inMatsWrap(i1,0,0) / determinant;
        } // for i1
 
      break;
    } // case 2

    case 1: {
   
      for (size_t i1=0; i1<dim_i1; i1++) {
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
	TEUCHOS_TEST_FOR_EXCEPTION( ( inMatsWrap(i1,0,0) == (Scalar)0 ),
				    std::invalid_argument,
				    ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif


        inverseMatsWrap(i1,0,0) = (Scalar)1 / inMatsWrap(i1,0,0); 
        } 
    
          
      
       break;
      } // case 1

    } // switch (dim)	
      break;
    case 2:
     switch(dim) {
    case 3: {
          size_t i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i,j) );
              }
            }
          }
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif

          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMatsWrap(rowperm[i],colperm[j]);
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif

          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMatsWrap(i,j) = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
   
      break;
    } // case 3

    case 2: {
          Scalar determinant    = inMatsWrap(0,0)*inMatsWrap(1,1)-inMatsWrap(0,1)*inMatsWrap(1,0);
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEUCHOS_TEST_FOR_EXCEPTION( ( (inMatsWrap(0,0)==(Scalar)0)   && (inMatsWrap(0,1)==(Scalar)0) &&
                                        (inMatsWrap(1,0)==(Scalar)0) && (inMatsWrap(1,1)==(Scalar)0) ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
          TEUCHOS_TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
                                      std::invalid_argument,
                                      ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif
          inverseMatsWrap(0,0)   = inMatsWrap(1,1) / determinant;
          inverseMatsWrap(0,1) = - inMatsWrap(0,1) / determinant;
          //
          inverseMatsWrap(1,0) = - inMatsWrap(1,0) / determinant;
          inverseMatsWrap(1,1) =   inMatsWrap(0,0) / determinant;
      
      break;
    } // case 2

    case 1: {
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION( ( inMatsWrap(0,0) == (Scalar)0 ),
				  std::invalid_argument,
				  ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif
           inverseMatsWrap(0,0) = (Scalar)1 / inMatsWrap(0,0);       
      break;
    } // case 1

  }
    break;  
  }

  
}
	



template<class Scalar>
Scalar RealSpaceTools<Scalar>::det(const Scalar* inMat, const size_t dim) {
  Scalar determinant(0);

  switch (dim) {
    case 3: {
      int i,j,rowID = 0;
      int colID = 0;
      int rowperm[3]={0,1,2};
      int colperm[3]={0,1,2}; // Complete pivoting
      Scalar emax(0);

      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( inMat[i*dim+j] ) >  emax){
            rowID = i;  colID = j; emax = std::abs( inMat[i*dim+j] );
          }
        }
      }
      if( emax > 0 ){
        if( rowID ){
          rowperm[0] = rowID;
          rowperm[rowID] = 0;
        }
        if( colID ){
          colperm[0] = colID;
          colperm[colID] = 0;
        }
        Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
        for(i=0; i < 3; ++i){
          for(j=0; j < 3; ++j){
            B[i][j] = inMat[rowperm[i]*dim+colperm[j]];
          }
        }
        B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
        for(i=0; i < 2; ++i){
          for(j=0; j < 2; ++j){
            S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
          }
        }
        determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
        if( rowID ) determinant = -determinant;
        if( colID ) determinant = -determinant;
      }
      break;
    } // case 3

    case 2:
      determinant = inMat[0]*inMat[3]-
                    inMat[1]*inMat[2];
      break;

    case 1:
      determinant = inMat[0];
      break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  } // switch (dim)

  return determinant;
}



template<class Scalar>
template<class ArrayIn>
Scalar RealSpaceTools<Scalar>::det(const ArrayIn & inMat) {

#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inMat) != 2),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Rank of matrix argument must be 2!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMat.dimension(0)) != static_cast<size_t>(inMat.dimension(1)) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Matrix is not square!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (static_cast<size_t>(inMat.dimension(0)) < 1) || (static_cast<size_t>(inMat.dimension(0)) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif
    ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inMatWrap(inMat);

  size_t dim = static_cast<size_t>(inMat.dimension(0));
  Scalar determinant(0);

  switch (dim) {
    case 3: {
      int i,j,rowID = 0;
      int colID = 0;
      int rowperm[3]={0,1,2};
      int colperm[3]={0,1,2}; // Complete pivoting
      Scalar emax(0);

      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( inMatWrap(i,j) ) >  emax){
            rowID = i;  colID = j; emax = std::abs( inMatWrap(i,j) );
          }
        }
      }
      if( emax > 0 ){
        if( rowID ){
          rowperm[0] = rowID;
          rowperm[rowID] = 0;
        }
        if( colID ){
          colperm[0] = colID;
          colperm[colID] = 0;
        }
        Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
        for(i=0; i < 3; ++i){
          for(j=0; j < 3; ++j){
            B[i][j] = inMatWrap(rowperm[i],colperm[j]);
          }
        }
        B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
        for(i=0; i < 2; ++i){
          for(j=0; j < 2; ++j){
            S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
          }
        }
        determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
        if( rowID ) determinant = -determinant;
        if( colID ) determinant = -determinant;
      }
      break;
    } // case 3

    case 2:
      determinant = inMatWrap(0,0)*inMatWrap(1,1)-
                    inMatWrap(0,1)*inMatWrap(1,0);
      break;

    case 1:
      determinant = inMatWrap(0,0);
      break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  } // switch (dim)

  return determinant;
}

template<class Scalar>
template<class ArrayDet, class ArrayIn>
void RealSpaceTools<Scalar>::det(ArrayDet & detArray, const ArrayIn & inMats) {
    ArrayWrapper<Scalar,ArrayDet, Rank<ArrayDet >::value, false>detArrayWrap(detArray);
    ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inMatsWrap(inMats);

  size_t matArrayRank = getrank(inMats);
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( matArrayRank != getrank(detArray)+2 ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::det): Determinant and matrix array arguments do not have compatible ranks!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (matArrayRank < 3) || (matArrayRank > 4) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::det): Rank of matrix array must be 3 or 4!");
  for (size_t i=0; i<matArrayRank-2; i++) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(i)) != static_cast<size_t>(detArray.dimension(i)) ),
				std::invalid_argument,
				">>> ERROR (RealSpaceTools::det): Dimensions of determinant and matrix array arguments do not agree!");
  }
  TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(matArrayRank-2)) != static_cast<size_t>(inMats.dimension(matArrayRank-1)) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::det): Matrices are not square!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (static_cast<size_t>(inMats.dimension(matArrayRank-2)) < 1) || (static_cast<size_t>(inMats.dimension(matArrayRank-2)) > 3) ),
			      std::invalid_argument,
			      ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif

  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = inMats.dimension(matArrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(matArrayRank) {
    case 4:
      dim_i0 = static_cast<size_t>(inMats.dimension(0));
      dim_i1 = static_cast<size_t>(inMats.dimension(1));
        switch(dim) {
    case 3: {
   

      for (size_t i0=0; i0<dim_i0; i0++) {
       
        for (size_t i1=0; i1<dim_i1; i1++) {
      

          int i,j,rowID = 0;
          int colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0), determinant(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i0,i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i0,i1,i,j) );
              }
            }
          }
          if( emax > 0 ){
            if( rowID ){
              rowperm[0] = rowID;
              rowperm[rowID] = 0;
            }
            if( colID ){
              colperm[0] = colID;
              colperm[colID] = 0;
            }
            Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
            for(i=0; i < 3; ++i){
              for(j=0; j < 3; ++j){
                B[i][j] = inMatsWrap(i0,i1,rowperm[i],colperm[j]);
              }
            }
            B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
            for(i=0; i < 2; ++i){
              for(j=0; j < 2; ++j){
                S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
              }
            }
            determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
            if( rowID ) determinant = -determinant;
            if( colID ) determinant = -determinant;
          }
          detArrayWrap(i0,i1)= determinant;
        } // for i1
      } // for i0
      break;
    } // case 3

    case 2: {   

      for (size_t i0=0; i0<dim_i0; i0++) {
      
        for (size_t i1=0; i1<dim_i1; i1++) {

          detArrayWrap(i0,i1) = inMatsWrap(i0,i1,0,0)*inMatsWrap(i0,i1,1,1)-inMatsWrap(i0,i1,0,1)*inMatsWrap(i0,i1,1,0);
        } // for i1
      } // for i0
      break;
    } // case 2

    case 1: {

      for (size_t i0=0; i0<dim_i0; i0++) {
      
        for (size_t i1=0; i1<dim_i1; i1++) {
          detArrayWrap(i0,i1) = inMatsWrap(i0,i1,0,0);
        } // for i1
      } // for i2
      break;
    } // case 1

  } // switch (dim)
      break;
    case 3:
      dim_i1 = static_cast<size_t>(inMats.dimension(0));
        switch(dim) {
    case 3: {   
      for (size_t i1=0; i1<dim_i1; i1++) {

          int i,j,rowID = 0;
          int colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0), determinant(0);
	

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i1,i,j) );
              }
            }
          }
          if( emax > 0 ){
            if( rowID ){
              rowperm[0] = rowID;
              rowperm[rowID] = 0;
            }
            if( colID ){
              colperm[0] = colID;
              colperm[colID] = 0;
            }
            Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
            for(i=0; i < 3; ++i){
              for(j=0; j < 3; ++j){
                B[i][j] = inMatsWrap(i1,rowperm[i],colperm[j]);
              }
            }
            B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
            for(i=0; i < 2; ++i){
              for(j=0; j < 2; ++j){
                S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
              }
            }
            determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
            if( rowID ) determinant = -determinant;
            if( colID ) determinant = -determinant;
          }
          detArrayWrap(i1) = determinant;
	  }       
      break;
    } // case 3

    case 2: {
      for (size_t i1=0; i1<dim_i1; i1++) {
      detArrayWrap(i1) = inMatsWrap(i1,0,0)*inMatsWrap(i1,1,1)-inMatsWrap(i1,0,1)*inMatsWrap(i1,1,0);
      }
      break;
    } // case 2

    case 1: {
      for (size_t i1=0; i1<dim_i1; i1++) {
	detArrayWrap(i1) = inMatsWrap(i1,0,0);
      }
      break;
    } // case 1
}
break;

  } // switch (dim)
		
		}
		


template<class Scalar>
void RealSpaceTools<Scalar>::add(Scalar* sumArray, const Scalar* inArray1, const Scalar* inArray2, const int size) {
  for (size_t i=0; i<size; i++) {
    sumArray[i] = inArray1[i] + inArray2[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::add(Scalar* inoutSumArray, const Scalar* inArray, const int size) {
  for (size_t i=0; i<size; i++) {
    inoutSumArray[i] += inArray[i];
  }
}



template<class Scalar>
template<class ArraySum, class ArrayIn1, class ArrayIn2>
void RealSpaceTools<Scalar>::add(ArraySum & sumArray, const ArrayIn1 & inArray1, const ArrayIn2 & inArray2) {
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inArray1) != getrank(inArray2)) || (getrank(inArray1) != getrank(sumArray)) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (size_t i=0; i<getrank(inArray1); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( (static_cast<size_t>(inArray1.dimension(i)) != static_cast<size_t>(inArray2.dimension(i))) || (static_cast<size_t>(inArray1.dimension(i)) != static_cast<size_t>(sumArray.dimension(i))) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif


  
  
   ArrayWrapper<Scalar,ArraySum, Rank<ArraySum >::value, false>sumArrayWrap(sumArray);
   ArrayWrapper<Scalar,ArrayIn1, Rank<ArrayIn1 >::value, true>inArray1Wrap(inArray1);
   ArrayWrapper<Scalar,ArrayIn2, Rank<ArrayIn2 >::value, true>inArray2Wrap(inArray2);
   int inArrayRank=getrank(inArray1);


   if(inArrayRank==5){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray1.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inArray1.dimension(3)); l++)
          for (size_t m=0; m<static_cast<size_t>(inArray1.dimension(4)); m++){
    sumArrayWrap(i,j,k,l,m) = inArray1Wrap(i,j,k,l,m)+inArray2Wrap(i,j,k,l,m);
          }
	}else if(inArrayRank==4){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray1.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inArray1.dimension(3)); l++){
           sumArrayWrap(i,j,k,l) = inArray1Wrap(i,j,k,l)+inArray2Wrap(i,j,k,l);
          }
	}else if(inArrayRank==3){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray1.dimension(2)); k++){
        sumArrayWrap(i,j,k) = inArray1Wrap(i,j,k)+inArray2Wrap(i,j,k);
          }
	}else if(inArrayRank==2){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++){
         sumArrayWrap(i,j) = inArray1Wrap(i,j)+inArray2Wrap(i,j);
          }
	}else if(inArrayRank==1){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++){
       sumArrayWrap(i) = inArray1Wrap(i)+inArray2Wrap(i);

          }
	}
}



template<class Scalar>
template<class ArraySum, class ArrayIn>
void RealSpaceTools<Scalar>::add(ArraySum & inoutSumArray, const ArrayIn & inArray) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( getrank(inArray) != getrank(inoutSumArray) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (size_t i=0; i<getrank(inArray); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inArray.dimension(i)) != static_cast<size_t>(inoutSumArray.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif
   
  	 ArrayWrapper<Scalar,ArraySum, Rank<ArraySum >::value, false>inoutSumArrayWrap(inoutSumArray);
	 ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inArrayWrap(inArray);
   int inArrayRank=getrank(inArray);


   if(inArrayRank==5){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++)
          for (size_t m=0; m<static_cast<size_t>(static_cast<size_t>(inArray.dimension(4))); m++){
    inoutSumArrayWrap(i,j,k,l,m) += inArrayWrap(i,j,k,l,m);
          }
	}else if(inArrayRank==4){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++){
          inoutSumArrayWrap(i,j,k,l) += inArrayWrap(i,j,k,l);
          }
	}else if(inArrayRank==3){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray.dimension(2)); k++){
          inoutSumArrayWrap(i,j,k) += inArrayWrap(i,j,k);
          }
	}else if(inArrayRank==2){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++){
          inoutSumArrayWrap(i,j) += inArrayWrap(i,j);
          }
	}else if(inArrayRank==1){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++){
          inoutSumArrayWrap(i) += inArrayWrap(i);

          }
	}
}



template<class Scalar>
void RealSpaceTools<Scalar>::subtract(Scalar* diffArray, const Scalar* inArray1, const Scalar* inArray2, const int size) {
  for (size_t i=0; i<size; i++) {
    diffArray[i] = inArray1[i] - inArray2[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::subtract(Scalar* inoutDiffArray, const Scalar* inArray, const int size) {
  for (int i=0; i<size; i++) {
    inoutDiffArray[i] -= inArray[i];
  }
}



template<class Scalar>
template<class ArrayDiff, class ArrayIn1, class ArrayIn2>
void RealSpaceTools<Scalar>::subtract(ArrayDiff & diffArray, const ArrayIn1 & inArray1, const ArrayIn2 & inArray2) {
	 ArrayWrapper<Scalar,ArrayDiff, Rank<ArrayDiff >::value, false>diffArrayWrap(diffArray);
	 ArrayWrapper<Scalar,ArrayIn1, Rank<ArrayIn1 >::value, true>inArray1Wrap(inArray1);
	 ArrayWrapper<Scalar,ArrayIn2, Rank<ArrayIn2 >::value, true>inArray2Wrap(inArray2);	 
	 size_t inArray1Rank=getrank(inArray1);
#ifdef HAVE_INTREPID_DEBUG
	 TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inArray1) != getrank(inArray2)) || (getrank(inArray1) != getrank(diffArray)) ),
				     std::invalid_argument,
				     ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
	 for (size_t i=0; i<getrank(inArray1); i++) {
	   TEUCHOS_TEST_FOR_EXCEPTION( ( (static_cast<size_t>(inArray1.dimension(i)) != static_cast<size_t>(inArray2.dimension(i))) || (static_cast<size_t>(inArray1.dimension(i)) != static_cast<size_t>(diffArray.dimension(i))) ),
				       std::invalid_argument,
				       ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
	 }
#endif
       if(inArray1Rank==5){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray1.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inArray1.dimension(3)); l++)
          for (size_t m=0; m<static_cast<size_t>(inArray1.dimension(4)); m++){
    diffArrayWrap(i,j,k,l,m) = inArray1Wrap(i,j,k,l,m)-inArray2Wrap(i,j,k,l,m);
          }
	}else if(inArray1Rank==4){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray1.dimension(2)); k++)
        for (size_t l=0; l<static_cast<size_t>(inArray1.dimension(3)); l++){
    diffArrayWrap(i,j,k,l) = inArray1Wrap(i,j,k,l)-inArray2Wrap(i,j,k,l);
          }
	}else if(inArray1Rank==3){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++)
      for (size_t k=0; k<static_cast<size_t>(inArray1.dimension(2)); k++){
    diffArrayWrap(i,j,k) = inArray1Wrap(i,j,k)-inArray2Wrap(i,j,k);
          }
	}else if(inArray1Rank==2){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++)
    for (size_t j=0; j<static_cast<size_t>(inArray1.dimension(1)); j++){
    diffArrayWrap(i,j) = inArray1Wrap(i,j)-inArray2Wrap(i,j);
          }
	}else if(inArray1Rank==1){
   for (size_t i=0; i<static_cast<size_t>(inArray1.dimension(0)); i++){
    diffArrayWrap(i) = inArray1Wrap(i)-inArray2Wrap(i);

          }
	}

}


template<class Scalar>
template<class ArrayDiff, class ArrayIn>
void RealSpaceTools<Scalar>::subtract(ArrayDiff & inoutDiffArray, const ArrayIn & inArray) {
	 ArrayWrapper<Scalar,ArrayDiff, Rank<ArrayDiff >::value, false>inoutDiffArrayWrap(inoutDiffArray);
	 ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inArrayWrap(inArray);
   int inArrayRank=getrank(inArray);
#ifdef HAVE_INTREPID_DEBUG
   TEUCHOS_TEST_FOR_EXCEPTION( ( getrank(inArray) != getrank(inoutDiffArray) ),
			       std::invalid_argument,
			       ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
   for (size_t i=0; i<getrank(inArray); i++) {
     TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inArray.dimension(i)) != static_cast<size_t>(inoutDiffArray.dimension(i)) ),
				 std::invalid_argument,
				 ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
   }
#endif

   if(inArrayRank==5){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++)
          for (size_t m=0; m<static_cast<size_t>(static_cast<size_t>(inArray.dimension(4))); m++){
    inoutDiffArrayWrap(i,j,k,l,m) -= inArrayWrap(i,j,k,l,m);
          }
	}else if(inArrayRank==4){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++){
          inoutDiffArrayWrap(i,j,k,l) -= inArrayWrap(i,j,k,l);
          }
	}else if(inArrayRank==3){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++){
          inoutDiffArrayWrap(i,j,k) -= inArrayWrap(i,j,k);
          }
	}else if(inArrayRank==2){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++){
          inoutDiffArrayWrap(i,j) -= inArrayWrap(i,j);
          }
	}else if(inArrayRank==1){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++){
          inoutDiffArrayWrap(i) -= inArrayWrap(i);

          }
	}
  
}


template<class Scalar>
void RealSpaceTools<Scalar>::scale(Scalar* scaledArray, const Scalar* inArray, const int size, const Scalar scalar) {
  for (size_t i=0; i<size; i++) {
    scaledArray[i] = scalar*inArray[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::scale(Scalar* inoutScaledArray, const int size, const Scalar scalar) {
  for (size_t i=0; i<size; i++) {
    inoutScaledArray[i] *= scalar;
  }
}



template<class Scalar>
template<class ArrayScaled, class ArrayIn>
void RealSpaceTools<Scalar>::scale(ArrayScaled & scaledArray, const ArrayIn & inArray, const Scalar scalar) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( getrank(inArray) != getrank(scaledArray) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::scale): Array arguments must have identical ranks!");
    for (size_t i=0; i<getrank(inArray); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inArray.dimension(i)) != static_cast<size_t>(scaledArray.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::scale): Dimensions of array arguments do not agree!");
    }
#endif


  int inArrayRank=getrank(inArray);
   ArrayWrapper<Scalar,ArrayScaled, Rank<ArrayScaled >::value, false>scaledArrayWrap(scaledArray);
   ArrayWrapper<Scalar,ArrayIn, Rank<ArrayIn >::value, true>inArrayWrap(inArray);
         if(inArrayRank==5){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++)
          for (size_t m=0; m<static_cast<size_t>(static_cast<size_t>(inArray.dimension(4))); m++){
    scaledArrayWrap(i,j,k,l,m) = scalar*inArrayWrap(i,j,k,l,m);
          }
	}else if(inArrayRank==4){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++)
        for (size_t l=0; l<static_cast<size_t>(static_cast<size_t>(inArray.dimension(3))); l++){
    scaledArrayWrap(i,j,k,l) = scalar*inArrayWrap(i,j,k,l);
          }
	}else if(inArrayRank==3){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++)
      for (size_t k=0; k<static_cast<size_t>(static_cast<size_t>(inArray.dimension(2))); k++){
    scaledArrayWrap(i,j,k) = scalar*inArrayWrap(i,j,k);
          }
	}else if(inArrayRank==2){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++)
    for (size_t j=0; j<static_cast<size_t>(static_cast<size_t>(inArray.dimension(1))); j++){
    scaledArrayWrap(i,j) = scalar*inArrayWrap(i,j);
          }
	}else if(inArrayRank==1){
   for (size_t i=0; i<static_cast<size_t>(static_cast<size_t>(inArray.dimension(0))); i++){
     scaledArrayWrap(i) = scalar*inArrayWrap(i);

          }
	}
}



template<class Scalar>
template<class ArrayScaled>
void RealSpaceTools<Scalar>::scale(ArrayScaled & inoutScaledArray, const Scalar scalar) {
  // Intrepid::FieldContainer has size type int
  const int theSize = (int) inoutScaledArray.size();
  for (int i=0; i<theSize; i++) {
    inoutScaledArray[i] *= scalar;
  }
}




template<class Scalar>
Scalar RealSpaceTools<Scalar>::dot(const Scalar* inArray1, const Scalar* inArray2, const int size) {
  Scalar dot(0);
  for (int i=0; i<size; i++) {
    dot += inArray1[i]*inArray2[i];
  }
  return dot;  
}



template<class Scalar>
template<class ArrayVec1, class ArrayVec2>
Scalar RealSpaceTools<Scalar>::dot(const ArrayVec1 & inVec1, const ArrayVec2 & inVec2) {
#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inVec1) != 1) || (getrank(inVec2) != 1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inVec1.dimension(0)) != static_cast<size_t>(inVec2.dimension(0)) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
#endif
   ArrayWrapper<Scalar,ArrayVec1, Rank<ArrayVec1 >::value, true>inVec1Wrap(inVec1);
   ArrayWrapper<Scalar,ArrayVec2, Rank<ArrayVec2 >::value, true>inVec2Wrap(inVec2);
  Scalar dot(0);
  for (size_t i=0; i<static_cast<size_t>(inVec1.dimension(0)); i++) {
    dot += inVec1Wrap(i)*inVec2Wrap(i);
  }
  return dot;  

}



template<class Scalar>
template<class ArrayDot, class ArrayVec1, class ArrayVec2>
void RealSpaceTools<Scalar>::dot(ArrayDot & dotArray, const ArrayVec1 & inVecs1, const ArrayVec2 & inVecs2) {

  size_t arrayRank = getrank(inVecs1);

#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( arrayRank != getrank(dotArray)+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Ranks of norm and vector array arguments are incompatible!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( arrayRank != getrank(inVecs2) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Ranks of input vector arguments must be identical!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Rank of input vector arguments must be 2 or 3!");
    for (size_t i=0; i<arrayRank; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inVecs1.dimension(i)) != static_cast<size_t>(inVecs2.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::dot): Dimensions of input vector arguments do not agree!");
    }
    for (size_t i=0; i<arrayRank-1; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inVecs1.dimension(i)) != static_cast<size_t>(dotArray.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::dot): Dimensions of dot-product and vector arrays do not agree!");
    }
#endif
   ArrayWrapper<Scalar,ArrayDot, Rank<ArrayDot >::value, false>dotArrayWrap(dotArray);
   ArrayWrapper<Scalar,ArrayVec1, Rank<ArrayVec1 >::value, true>inVecs1Wrap(inVecs1);
   ArrayWrapper<Scalar,ArrayVec2, Rank<ArrayVec2 >::value, true>inVecs2Wrap(inVecs2);
  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = static_cast<size_t>(inVecs1.dimension(arrayRank-1)); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 3:
      dim_i0 = static_cast<size_t>(inVecs1.dimension(0));
      dim_i1 = static_cast<size_t>(inVecs1.dimension(1));
   for (size_t i0=0; i0<dim_i0; i0++) {
    for (size_t i1=0; i1<dim_i1; i1++) {
      Scalar dot(0);
      for (size_t i=0; i<dim; i++) {
        dot += inVecs1Wrap(i0,i1,i)*inVecs2Wrap(i0,i1,i);
      }
      dotArrayWrap(i0,i1) = dot;
    }
  }
      break;
    case 2:
      dim_i1 = static_cast<size_t>(inVecs1.dimension(0));
     for (size_t i1=0; i1<dim_i1; i1++) {
      Scalar dot(0);
      for (size_t i=0; i<dim; i++) {
        dot += inVecs1Wrap(i1,i)*inVecs2Wrap(i1,i);
      }
      dotArrayWrap(i1) = dot;
    }
      break;
    case 1:
    Scalar dot(0);
     for (size_t i=0; i<dim; i++) {
        dot += inVecs1Wrap(i)*inVecs2Wrap(i);
      }
      dotArrayWrap(0) = dot;
    
    break;      
  }  
  
}



template<class Scalar>
void RealSpaceTools<Scalar>::matvec(Scalar* matVec, const Scalar* inMat, const Scalar* inVec, const size_t dim) {
  for (size_t i=0; i<dim; i++) {
    Scalar sumdot(0);
    for (size_t j=0; j<dim; j++) {
      sumdot += inMat[i*dim+j]*inVec[j];
    }
    matVec[i] = sumdot; 
  }
}



template<class Scalar>
template<class ArrayMatVec, class ArrayMat, class ArrayVec>
void RealSpaceTools<Scalar>::matvec(ArrayMatVec & matVecs, const ArrayMat & inMats, const ArrayVec & inVecs) {
  size_t matArrayRank = getrank(inMats);

#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( matArrayRank != getrank(inVecs)+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Vector and matrix array arguments do not have compatible ranks!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (matArrayRank < 3) || (matArrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Rank of matrix array must be 3 or 4!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( getrank(matVecs) != getrank(inVecs) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Vector arrays must be have the same rank!");
    for (size_t i=0; i<matArrayRank-1; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(i)) != static_cast<size_t>(inVecs.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector and matrix array arguments do not agree!");
    }
    for (size_t i=0; i<getrank(inVecs); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(matVecs.dimension(i)) != static_cast<size_t>(inVecs.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector array arguments do not agree!");
    }
    TEUCHOS_TEST_FOR_EXCEPTION( ( static_cast<size_t>(inMats.dimension(matArrayRank-2)) != static_cast<size_t>(inMats.dimension(matArrayRank-1)) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Matrices are not square!");
#endif
    ArrayWrapper<Scalar,ArrayMatVec, Rank<ArrayMatVec >::value, false>matVecsWrap(matVecs);
    ArrayWrapper<Scalar,ArrayMat, Rank<ArrayMat >::value, true>inMatsWrap(inMats);
    ArrayWrapper<Scalar,ArrayVec, Rank<ArrayVec >::value, true>inVecsWrap(inVecs);
  size_t dim_i0 = 1; // first  index dimension (e.g. cell index)
  size_t dim_i1 = 1; // second index dimension (e.g. point index)
  size_t dim    = static_cast<size_t>(inMats.dimension(matArrayRank-2)); // spatial dimension

  // determine i0 and i1 dimensions
  switch(matArrayRank) {
    case 4:
      dim_i0 = static_cast<size_t>(inMats.dimension(0));
      dim_i1 = static_cast<size_t>(inMats.dimension(1));
   for (size_t i0=0; i0<dim_i0; i0++) {
    for (size_t i1=0; i1<dim_i1; i1++) {
      for (size_t i=0; i<dim; i++) {
        Scalar sumdot(0);
        for (size_t j=0; j<dim; j++) {
          sumdot += inMatsWrap(i0,i1,i,j)*inVecsWrap(i0,i1,j);
        }
        matVecsWrap(i0,i1,i) = sumdot;
      }
    }
  }
      break;
    case 3:
      dim_i1 = static_cast<size_t>(inMats.dimension(0));
  
    for (size_t i1=0; i1<dim_i1; i1++) {
      for (size_t i=0; i<dim; i++) {
        Scalar sumdot(0);
        for (size_t j=0; j<dim; j++) {
          sumdot += inMatsWrap(i1,i,j)*inVecsWrap(i1,j);
        }
        matVecsWrap(i1,i) = sumdot;
      }
    }
        
      break;
          
  }
  
}
template<class Scalar>
template<class ArrayVecProd, class ArrayIn1, class ArrayIn2>
void RealSpaceTools<Scalar>::vecprod(ArrayVecProd & vecProd, const ArrayIn1 & inLeft, const ArrayIn2 & inRight) {
    ArrayWrapper<Scalar,ArrayVecProd, Rank<ArrayVecProd >::value, false>vecProdWrap(vecProd);
    ArrayWrapper<Scalar,ArrayIn1, Rank<ArrayIn1 >::value, true>inLeftWrap(inLeft);
    ArrayWrapper<Scalar,ArrayIn2, Rank<ArrayIn2 >::value, true>inRightWrap(inRight);
#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) all array arguments are required to have matching dimensions and rank: (D), (I0,D) or (I0,I1,D)
     *      (2) spatial dimension should be 2 or 3
     */
    std::string errmsg = ">>> ERROR (RealSpaceTools::vecprod):";
  
    // (1) check rank range on inLeft and then compare the other arrays with inLeft
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inLeft,  1,3), std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inLeft, inRight), std::invalid_argument, errmsg);    
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inLeft, vecProd), std::invalid_argument, errmsg);   
  
    // (2) spatial dimension ordinal = array rank - 1. Suffices to check only one array because we just
    //     checked whether or not the arrays have matching dimensions. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inLeft, getrank(inLeft) - 1,  2,3), std::invalid_argument, errmsg);
  
#endif 

    int spaceDim = static_cast<size_t>(inLeft.dimension(getrank(inLeft) - 1));

    switch(getrank(inLeft) ){
    
    case 1:
      {        
        vecProdWrap(0) = inLeftWrap(1)*inRightWrap(2) - inLeftWrap(2)*inRightWrap(1);
        vecProdWrap(1) = inLeftWrap(2)*inRightWrap(0) - inLeftWrap(0)*inRightWrap(2);              
        vecProdWrap(2) = inLeftWrap(0)*inRightWrap(1) - inLeftWrap(1)*inRightWrap(0);    
      }
      break;
      
    case 2:
      {
        size_t dim0 = static_cast<size_t>(inLeft.dimension(0));
        if(spaceDim == 3) {
          for(size_t i0 = 0; i0 < dim0; i0++){
            vecProdWrap(i0, 0) = inLeftWrap(i0, 1)*inRightWrap(i0, 2) - inLeftWrap(i0, 2)*inRightWrap(i0, 1);
            vecProdWrap(i0, 1) = inLeftWrap(i0, 2)*inRightWrap(i0, 0) - inLeftWrap(i0, 0)*inRightWrap(i0, 2);              
            vecProdWrap(i0, 2) = inLeftWrap(i0, 0)*inRightWrap(i0, 1) - inLeftWrap(i0, 1)*inRightWrap(i0, 0);
          }// i0
        } //spaceDim == 3
        else if(spaceDim == 2){
          for(size_t i0 = 0; i0 < dim0; i0++){
            // vecprod is scalar - do we still want result to be (i0,i1,D)?
            vecProdWrap(i0, 0) = inLeftWrap(i0, 0)*inRightWrap(i0, 1) - inLeftWrap(i0, 1)*inRightWrap(i0, 0);
          }// i0
        }// spaceDim == 2
      }// case 2
      break;
      
    case 3:
      {
        size_t dim0 = static_cast<size_t>(inLeft.dimension(0));
        size_t dim1 = static_cast<size_t>(inLeft.dimension(1));
        if(spaceDim == 3) {
          for(size_t i0 = 0; i0 < dim0; i0++){
            for(size_t i1 = 0; i1 < dim1; i1++){
              vecProdWrap(i0, i1, 0) = inLeftWrap(i0, i1, 1)*inRightWrap(i0, i1, 2) - inLeftWrap(i0, i1, 2)*inRightWrap(i0, i1, 1);
              vecProdWrap(i0, i1, 1) = inLeftWrap(i0, i1, 2)*inRightWrap(i0, i1, 0) - inLeftWrap(i0, i1, 0)*inRightWrap(i0, i1, 2);              
              vecProdWrap(i0, i1, 2) = inLeftWrap(i0, i1, 0)*inRightWrap(i0, i1, 1) - inLeftWrap(i0, i1, 1)*inRightWrap(i0, i1, 0);
            }// i1
          }// i0
        } //spaceDim == 3
        else if(spaceDim == 2){
          for(size_t i0 = 0; i0 < dim0; i0++){
            for(size_t i1 = 0; i1 < dim1; i1++){
              // vecprod is scalar - do we still want result to be (i0,i1,D)?
              vecProdWrap(i0, i1, 0) = inLeftWrap(i0, i1, 0)*inRightWrap(i0, i1, 1) - inLeftWrap(i0, i1, 1)*inRightWrap(i0, i1, 0);
            }// i1
          }// i0
        }// spaceDim == 2
      } // case 3
      break;
      
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
                         ">>> ERROR (RealSpaceTools::vecprod): rank-1,2,3 arrays required");      
  }
  
}


} // namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

