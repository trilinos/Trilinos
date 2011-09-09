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
  for (int i=0; i<size; i++) {
    absArray[i] = std::abs(inArray[i]);
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::absval(Scalar* inoutAbsArray, const int size) {
  for (int i=0; i<size; i++) {
    inoutAbsArray[i] = std::abs(inoutAbsArray[i]);
  }
}



template<class Scalar>
template<class ArrayAbs, class ArrayIn>
void RealSpaceTools<Scalar>::absval(ArrayAbs & absArray, const ArrayIn & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.rank() != absArray.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::absval): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.rank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.dimension(i) != absArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::absval): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.size(); i++) {
    absArray[i] = std::abs(inArray[i]);
  }
}



template<class Scalar>
template<class ArrayInOut>
void RealSpaceTools<Scalar>::absval(ArrayInOut & inoutAbsArray) {
  for (int i=0; i<inoutAbsArray.size(); i++) {
    inoutAbsArray[i] = std::abs(inoutAbsArray[i]);
  }
}



template<class Scalar>
Scalar RealSpaceTools<Scalar>::vectorNorm(const Scalar* inVec, const int dim, const ENorm normType) {
  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(int i = 0; i < dim; i++){
        temp += inVec[i]*inVec[i];
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(inVec[0]);
      for(int i = 1; i < dim; i++){
        Scalar absData = std::abs(inVec[i]);
        if (temp < absData) temp = absData;
      }
      break;
    case NORM_ONE:
      for(int i = 0; i < dim; i++){
        temp += std::abs(inVec[i]);
      }
      break;
    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}



template<class Scalar>
template<class ArrayIn>
Scalar RealSpaceTools<Scalar>::vectorNorm(const ArrayIn & inVec, const ENorm normType) {

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inVec.rank() != 1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Vector argument must have rank 1!");
#endif

  int size = inVec.size();

  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(int i = 0; i < size; i++){
        temp += inVec[i]*inVec[i];
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(inVec[0]);
      for(int i = 1; i < size; i++){
        Scalar absData = std::abs(inVec[i]);
        if (temp < absData) temp = absData;
      }
      break;
    case NORM_ONE:
      for(int i = 0; i < size; i++){
        temp += std::abs(inVec[i]);
      }
      break;
    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}



template<class Scalar>
template<class ArrayNorm, class ArrayIn>
void RealSpaceTools<Scalar>::vectorNorm(ArrayNorm & normArray, const ArrayIn & inVecs, const ENorm normType) {

  int arrayRank = inVecs.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != normArray.rank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Ranks of norm and vector array arguments are incompatible!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Rank of vector array must be 2 or 3!");
    for (int i=0; i<arrayRank-1; i++) {
      TEST_FOR_EXCEPTION( ( inVecs.dimension(i) != normArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Dimensions of norm and vector arguments do not agree!");
    }
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inVecs.dimension(arrayRank-1); // spatial dimension

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
      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          for(int i = 0; i < dim; i++){
            temp += inVecs[offset+i]*inVecs[offset+i];
          }
          normArray[normOffset] = std::sqrt(temp);
        }
      }
      break;
    } // case NORM_TWO

    case NORM_INF: {
      int offset_i0, offset, normOffset;
      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          temp = std::abs(inVecs[offset]);
          for(int i = 1; i < dim; i++){
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
      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          for(int i = 0; i < dim; i++){
            temp += std::abs(inVecs[offset+i]);
          }
          normArray[normOffset] = temp;
        }
      }
      break;
    } // case NORM_ONE

    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::transpose(Scalar* transposeMat, const Scalar* inMat, const int dim) {
  for(int i=0; i < dim; i++){
    transposeMat[i*dim+i]=inMat[i*dim+i];    // Set diagonal elements
    for(int j=i+1; j < dim; j++){
      transposeMat[i*dim+j]=inMat[j*dim+i];  // Set off-diagonal elements
      transposeMat[j*dim+i]=inMat[i*dim+j];
    }
  }
}



template<class Scalar>
template<class ArrayTranspose, class ArrayIn>
void RealSpaceTools<Scalar>::transpose(ArrayTranspose & transposeMats, const ArrayIn & inMats) {
  int arrayRank = inMats.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != transposeMats.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::transpose): Matrix array arguments do not have identical ranks!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::transpose): Rank of matrix array must be 2, 3, or 4!");
    for (int i=0; i<arrayRank; i++) {
      TEST_FOR_EXCEPTION( ( inMats.dimension(i) != transposeMats.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::transpose): Dimensions of matrix arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.dimension(arrayRank-2) != inMats.dimension(arrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::transpose): Matrices are not square!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.dimension(arrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = inMats.dimension(0);
      dim_i1 = inMats.dimension(1);
      break;
    case 3:
      dim_i1 = inMats.dimension(0);
      break;
  }

  int offset_i0, offset;

  for (int i0=0; i0<dim_i0; i0++) {
    offset_i0 = i0*dim_i1;
    for (int i1=0; i1<dim_i1; i1++) {
      offset  = offset_i0 + i1;
      offset *= (dim*dim);

      for(int i=0; i < dim; i++){
        transposeMats[offset+i*dim+i]=inMats[offset+i*dim+i];    // Set diagonal elements
        for(int j=i+1; j < dim; j++){
          transposeMats[offset+i*dim+j]=inMats[offset+j*dim+i];  // Set off-diagonal elements
          transposeMats[offset+j*dim+i]=inMats[offset+i*dim+j];
        }
      }

    } // i1
  } // i0

}



template<class Scalar>
void RealSpaceTools<Scalar>::inverse(Scalar* inverseMat, const Scalar* inMat, const int dim) {

  switch(dim) {
    case 3: {
      int i, j, rowID = 0, colID = 0;
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
      TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
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
      TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
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
      TEST_FOR_EXCEPTION( ( (inMat[0]==(Scalar)0) && (inMat[1]==(Scalar)0) &&
                            (inMat[2]==(Scalar)0) && (inMat[3]==(Scalar)0) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
      TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
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
      TEST_FOR_EXCEPTION( ( inMat[0] == (Scalar)0 ),
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

  int arrayRank = inMats.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != inverseMats.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Matrix array arguments do not have identical ranks!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Rank of matrix array must be 2, 3, or 4!");
    for (int i=0; i<arrayRank; i++) {
      TEST_FOR_EXCEPTION( ( inMats.dimension(i) != inverseMats.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::inverse): Dimensions of matrix arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.dimension(arrayRank-2) != inMats.dimension(arrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Matrices are not square!");
    TEST_FOR_EXCEPTION( ( (inMats.dimension(arrayRank-2) < 1) || (inMats.dimension(arrayRank-2) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Spatial dimension must be 1, 2, or 3!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.dimension(arrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = inMats.dimension(0);
      dim_i1 = inMats.dimension(1);
      break;
    case 3:
      dim_i1 = inMats.dimension(0);
      break;
  }

  switch(dim) {
    case 3: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;
          offset *= 9;

          int i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMats[offset+i*3+j] ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMats[offset+i*3+j] );
              }
            }
          }
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
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
              B[i][j] = inMats[offset+rowperm[i]*3+colperm[j]];
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
          TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
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
              inverseMats[offset+i*3+j] = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
        } // for i1
      } // for i0
      break;
    } // case 3

    case 2: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;;
          offset *= 4;

          Scalar determinant    = inMats[offset]*inMats[offset+3]-inMats[offset+1]*inMats[offset+2];;
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEST_FOR_EXCEPTION( ( (inMats[offset]==(Scalar)0)   && (inMats[offset+1]==(Scalar)0) &&
                                (inMats[offset+2]==(Scalar)0) && (inMats[offset+3]==(Scalar)0) ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
          TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif
          inverseMats[offset]   = inMats[offset+3] / determinant;
          inverseMats[offset+1] = - inMats[offset+1] / determinant;
          //
          inverseMats[offset+2] = - inMats[offset+2] / determinant;
          inverseMats[offset+3] =   inMats[offset] / determinant;
        } // for i1
      } // for i0
      break;
    } // case 2

    case 1: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;;
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
          TEST_FOR_EXCEPTION( ( inMats[offset] == (Scalar)0 ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
#endif
          inverseMats[offset] = (Scalar)1 / inMats[offset];
        } // for i1
      } // for i2
      break;
    } // case 1

  } // switch (dim)
}



template<class Scalar>
Scalar RealSpaceTools<Scalar>::det(const Scalar* inMat, const int dim) {
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
      TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  } // switch (dim)

  return determinant;
}



template<class Scalar>
template<class ArrayIn>
Scalar RealSpaceTools<Scalar>::det(const ArrayIn & inMat) {

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( (inMat.rank() != 2),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Rank of matrix argument must be 2!");
    TEST_FOR_EXCEPTION( ( inMat.dimension(0) != inMat.dimension(1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Matrix is not square!");
    TEST_FOR_EXCEPTION( ( (inMat.dimension(0) < 1) || (inMat.dimension(0) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif

  int dim = inMat.dimension(0);
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
      TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  } // switch (dim)

  return determinant;
}




template<class Scalar>
template<class ArrayDet, class ArrayIn>
void RealSpaceTools<Scalar>::det(ArrayDet & detArray, const ArrayIn & inMats) {

  int matArrayRank = inMats.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( matArrayRank != detArray.rank()+2 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Determinant and matrix array arguments do not have compatible ranks!");
    TEST_FOR_EXCEPTION( ( (matArrayRank < 3) || (matArrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Rank of matrix array must be 3 or 4!");
    for (int i=0; i<matArrayRank-2; i++) {
      TEST_FOR_EXCEPTION( ( inMats.dimension(i) != detArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::det): Dimensions of determinant and matrix array arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.dimension(matArrayRank-2) != inMats.dimension(matArrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Matrices are not square!");
    TEST_FOR_EXCEPTION( ( (inMats.dimension(matArrayRank-2) < 1) || (inMats.dimension(matArrayRank-2) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.dimension(matArrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(matArrayRank) {
    case 4:
      dim_i0 = inMats.dimension(0);
      dim_i1 = inMats.dimension(1);
      break;
    case 3:
      dim_i1 = inMats.dimension(0);
      break;
  }

  switch(dim) {
    case 3: {
      int offset_i0, offset, detOffset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset     = offset_i0 + i1;
          detOffset  = offset;
          offset    *= 9;

          int i,j,rowID = 0;
          int colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0), determinant(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMats[offset+i*3+j] ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMats[offset+i*3+j] );
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
                B[i][j] = inMats[offset+rowperm[i]*3+colperm[j]];
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
          detArray[detOffset] = determinant;
        } // for i1
      } // for i0
      break;
    } // case 3

    case 2: {
      int offset_i0, offset, detOffset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset     = offset_i0 + i1;
          detOffset  = offset;
          offset    *= 4;

          detArray[detOffset] = inMats[offset]*inMats[offset+3]-inMats[offset+1]*inMats[offset+2];;
        } // for i1
      } // for i0
      break;
    } // case 2

    case 1: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;;
          detArray[offset] = inMats[offset];
        } // for i1
      } // for i2
      break;
    } // case 1

  } // switch (dim)
}



template<class Scalar>
void RealSpaceTools<Scalar>::add(Scalar* sumArray, const Scalar* inArray1, const Scalar* inArray2, const int size) {
  for (int i=0; i<size; i++) {
    sumArray[i] = inArray1[i] + inArray2[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::add(Scalar* inoutSumArray, const Scalar* inArray, const int size) {
  for (int i=0; i<size; i++) {
    inoutSumArray[i] += inArray[i];
  }
}



template<class Scalar>
template<class ArraySum, class ArrayIn1, class ArrayIn2>
void RealSpaceTools<Scalar>::add(ArraySum & sumArray, const ArrayIn1 & inArray1, const ArrayIn2 & inArray2) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( (inArray1.rank() != inArray2.rank()) || (inArray1.rank() != sumArray.rank()) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (int i=0; i<inArray1.rank(); i++) {
      TEST_FOR_EXCEPTION( ( (inArray1.dimension(i) != inArray2.dimension(i)) || (inArray1.dimension(i) != sumArray.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray1.size(); i++) {
    sumArray[i] = inArray1[i] + inArray2[i];
  }
}



template<class Scalar>
template<class ArraySum, class ArrayIn>
void RealSpaceTools<Scalar>::add(ArraySum & inoutSumArray, const ArrayIn & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.rank() != inoutSumArray.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.rank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.dimension(i) != inoutSumArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.size(); i++) {
    inoutSumArray[i] += inArray[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::subtract(Scalar* diffArray, const Scalar* inArray1, const Scalar* inArray2, const int size) {
  for (int i=0; i<size; i++) {
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
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( (inArray1.rank() != inArray2.rank()) || (inArray1.rank() != diffArray.rank()) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
    for (int i=0; i<inArray1.rank(); i++) {
      TEST_FOR_EXCEPTION( ( (inArray1.dimension(i) != inArray2.dimension(i)) || (inArray1.dimension(i) != diffArray.dimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray1.size(); i++) {
    diffArray[i] = inArray1[i] - inArray2[i];
  }
}



template<class Scalar>
template<class ArrayDiff, class ArrayIn>
void RealSpaceTools<Scalar>::subtract(ArrayDiff & inoutDiffArray, const ArrayIn & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.rank() != inoutDiffArray.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.rank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.dimension(i) != inoutDiffArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.size(); i++) {
    inoutDiffArray[i] -= inArray[i];
  }
}




template<class Scalar>
void RealSpaceTools<Scalar>::scale(Scalar* scaledArray, const Scalar* inArray, const int size, const Scalar scalar) {
  for (int i=0; i<size; i++) {
    scaledArray[i] = scalar*inArray[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::scale(Scalar* inoutScaledArray, const int size, const Scalar scalar) {
  for (int i=0; i<size; i++) {
    inoutScaledArray[i] *= scalar;
  }
}



template<class Scalar>
template<class ArrayScaled, class ArrayIn>
void RealSpaceTools<Scalar>::scale(ArrayScaled & scaledArray, const ArrayIn & inArray, const Scalar scalar) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.rank() != scaledArray.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::scale): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.rank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.dimension(i) != scaledArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::scale): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.size(); i++) {
    scaledArray[i] = scalar*inArray[i];
  }
}



template<class Scalar>
template<class ArrayScaled>
void RealSpaceTools<Scalar>::scale(ArrayScaled & inoutScaledArray, const Scalar scalar) {
  for (int i=0; i<inoutScaledArray.size(); i++) {
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
    TEST_FOR_EXCEPTION( ( (inVec1.rank() != 1) || (inVec2.rank() != 1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
    TEST_FOR_EXCEPTION( ( inVec1.dimension(0) != inVec2.dimension(0) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
#endif

  Scalar dot(0);
  for (int i=0; i<inVec1.size(); i++) {
    dot += inVec1[i]*inVec2[i];
  }
  return dot;  

}



template<class Scalar>
template<class ArrayDot, class ArrayVec1, class ArrayVec2>
void RealSpaceTools<Scalar>::dot(ArrayDot & dotArray, const ArrayVec1 & inVecs1, const ArrayVec2 & inVecs2) {

  int arrayRank = inVecs1.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != dotArray.rank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Ranks of norm and vector array arguments are incompatible!");
    TEST_FOR_EXCEPTION( ( arrayRank != inVecs2.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Ranks of input vector arguments must be identical!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Rank of input vector arguments must be 2 or 3!");
    for (int i=0; i<arrayRank; i++) {
      TEST_FOR_EXCEPTION( ( inVecs1.dimension(i) != inVecs2.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::dot): Dimensions of input vector arguments do not agree!");
    }
    for (int i=0; i<arrayRank-1; i++) {
      TEST_FOR_EXCEPTION( ( inVecs1.dimension(i) != dotArray.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::dot): Dimensions of dot-product and vector arrays do not agree!");
    }
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inVecs1.dimension(arrayRank-1); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 3:
      dim_i0 = inVecs1.dimension(0);
      dim_i1 = inVecs1.dimension(1);
      break;
    case 2:
      dim_i1 = inVecs1.dimension(0);
      break;
  }

  int offset_i0, offset, dotOffset;
  for (int i0=0; i0<dim_i0; i0++) {
    offset_i0 = i0*dim_i1;
    for (int i1=0; i1<dim_i1; i1++) {
      offset      = offset_i0 + i1;
      dotOffset   = offset;
      offset     *= dim;
      Scalar dot(0);
      for (int i=0; i<dim; i++) {
        dot += inVecs1[offset+i]*inVecs2[offset+i];
      }
      dotArray[dotOffset] = dot;
    }
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::matvec(Scalar* matVec, const Scalar* inMat, const Scalar* inVec, const int dim) {
  for (int i=0; i<dim; i++) {
    Scalar sumdot(0);
    for (int j=0; j<dim; j++) {
      sumdot += inMat[i*dim+j]*inVec[j];
    }
    matVec[i] = sumdot; 
  }
}



template<class Scalar>
template<class ArrayMatVec, class ArrayMat, class ArrayVec>
void RealSpaceTools<Scalar>::matvec(ArrayMatVec & matVecs, const ArrayMat & inMats, const ArrayVec & inVecs) {
  int matArrayRank = inMats.rank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( matArrayRank != inVecs.rank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Vector and matrix array arguments do not have compatible ranks!");
    TEST_FOR_EXCEPTION( ( (matArrayRank < 3) || (matArrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Rank of matrix array must be 3 or 4!");
    TEST_FOR_EXCEPTION( ( matVecs.rank() != inVecs.rank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Vector arrays must be have the same rank!");
    for (int i=0; i<matArrayRank-1; i++) {
      TEST_FOR_EXCEPTION( ( inMats.dimension(i) != inVecs.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector and matrix array arguments do not agree!");
    }
    for (int i=0; i<inVecs.rank(); i++) {
      TEST_FOR_EXCEPTION( ( matVecs.dimension(i) != inVecs.dimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector array arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.dimension(matArrayRank-2) != inMats.dimension(matArrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Matrices are not square!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.dimension(matArrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(matArrayRank) {
    case 4:
      dim_i0 = inMats.dimension(0);
      dim_i1 = inMats.dimension(1);
      break;
    case 3:
      dim_i1 = inMats.dimension(0);
      break;
  }

  int offset_i0, offset, vecOffset;

  for (int i0=0; i0<dim_i0; i0++) {
    offset_i0 = i0*dim_i1;
    for (int i1=0; i1<dim_i1; i1++) {
      offset     = offset_i0 + i1;
      vecOffset  = offset*dim;
      offset     = vecOffset*dim;

      for (int i=0; i<dim; i++) {
        Scalar sumdot(0);
        for (int j=0; j<dim; j++) {
          sumdot += inMats[offset+i*dim+j]*inVecs[vecOffset+j];
        }
        matVecs[vecOffset+i] = sumdot;
      }
    }
  }
}


template<class Scalar>
template<class ArrayVecProd, class ArrayIn1, class ArrayIn2>
void RealSpaceTools<Scalar>::vecprod(ArrayVecProd & vecProd, const ArrayIn1 & inLeft, const ArrayIn2 & inRight) {
  
#ifdef HAVE_INTREPID_DEBUG
  /*
   *   Check array rank and spatial dimension range (if applicable)
   *      (1) all array arguments are required to have matching dimensions and rank: (D), (I0,D) or (I0,I1,D)
   *      (2) spatial dimension should be 2 or 3
   */
  std::string errmsg = ">>> ERROR (RealSpaceTools::vecprod):";
  
  // (1) check rank range on inLeft and then compare the other arrays with inLeft
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inLeft,  1,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inLeft, inRight), std::invalid_argument, errmsg);    
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inLeft, vecProd), std::invalid_argument, errmsg);   
  
  // (2) spatial dimension ordinal = array rank - 1. Suffices to check only one array because we just
  //     checked whether or not the arrays have matching dimensions. 
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inLeft, inLeft.rank() - 1,  2,3), std::invalid_argument, errmsg);
  
#endif

 int spaceDim = inLeft.dimension(inLeft.rank() - 1);

  switch(inLeft.rank() ){
    
    case 1:
      {        
        vecProd(0) = inLeft(1)*inRight(2) - inLeft(2)*inRight(1);
        vecProd(1) = inLeft(2)*inRight(0) - inLeft(0)*inRight(2);              
        vecProd(2) = inLeft(0)*inRight(1) - inLeft(1)*inRight(0);    
      }
      break;
      
    case 2:
      {
        int dim0 = inLeft.dimension(0);
        if(spaceDim == 3) {
          for(int i0 = 0; i0 < dim0; i0++){
            vecProd(i0, 0) = inLeft(i0, 1)*inRight(i0, 2) - inLeft(i0, 2)*inRight(i0, 1);
            vecProd(i0, 1) = inLeft(i0, 2)*inRight(i0, 0) - inLeft(i0, 0)*inRight(i0, 2);              
            vecProd(i0, 2) = inLeft(i0, 0)*inRight(i0, 1) - inLeft(i0, 1)*inRight(i0, 0);
          }// i0
        } //spaceDim == 3
        else if(spaceDim == 2){
          for(int i0 = 0; i0 < dim0; i0++){
            // vecprod is scalar - do we still want result to be (i0,i1,D)?
            vecProd(i0, 0) = inLeft(i0, 0)*inRight(i0, 1) - inLeft(i0, 1)*inRight(i0, 0);
          }// i0
        }// spaceDim == 2
      }// case 2
      break;
      
    case 3:
      {
        int dim0 = inLeft.dimension(0);
        int dim1 = inLeft.dimension(1);
        if(spaceDim == 3) {
          for(int i0 = 0; i0 < dim0; i0++){
            for(int i1 = 0; i1 < dim1; i1++){
              vecProd(i0, i1, 0) = inLeft(i0, i1, 1)*inRight(i0, i1, 2) - inLeft(i0, i1, 2)*inRight(i0, i1, 1);
              vecProd(i0, i1, 1) = inLeft(i0, i1, 2)*inRight(i0, i1, 0) - inLeft(i0, i1, 0)*inRight(i0, i1, 2);              
              vecProd(i0, i1, 2) = inLeft(i0, i1, 0)*inRight(i0, i1, 1) - inLeft(i0, i1, 1)*inRight(i0, i1, 0);
            }// i1
          }// i0
        } //spaceDim == 3
        else if(spaceDim == 2){
          for(int i0 = 0; i0 < dim0; i0++){
            for(int i1 = 0; i1 < dim1; i1++){
              // vecprod is scalar - do we still want result to be (i0,i1,D)?
              vecProd(i0, i1, 0) = inLeft(i0, i1, 0)*inRight(i0, i1, 1) - inLeft(i0, i1, 1)*inRight(i0, i1, 0);
            }// i1
          }// i0
        }// spaceDim == 2
      } // case 3
      break;
      
    default:
      TEST_FOR_EXCEPTION(true, std::invalid_argument, 
                         ">>> ERROR (RealSpaceTools::vecprod): rank-1,2,3 arrays required");      
  }
  
}

} // namespace Intrepid
