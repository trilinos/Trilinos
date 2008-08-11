// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_RealSpaceDef.hpp
    \brief  Definition file for utility classes providing basic linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


namespace Intrepid {

//===========================================================================//
//                                                                           //
//              Member function definitions of the class Point.              //
//                                                                           //
//===========================================================================//

template<class Scalar>
Point<Scalar>::Point(const Point<Scalar>& right) {
  data_.resize(right.getDim());
  data_ = right.data_;
  frameKind_ = right.frameKind_;
}


  
template<class Scalar>
Point<Scalar>::Point(const int dim,
                     const EFrame frameKind): frameKind_(frameKind) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (dim < 1) || (dim > INTREPID_MAX_DIMENSION) ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid dim argument.");
#endif
  data_.assign(dim, (Scalar)0);
}



template<class Scalar>
Point<Scalar>::Point(const Scalar x,
                     const EFrame frameKind) : frameKind_(frameKind) {
  data_.resize(1);
  data_[0] = x;
}


  
template<class Scalar>
Point<Scalar>::Point(const Scalar x, 
                     const Scalar y, 
                     const EFrame frameKind) : frameKind_(frameKind)  {
  data_.resize(2);
  data_[0] = x;
  data_[1] = y;
}


  
template<class Scalar>
Point<Scalar>::Point(const int x, 
                     const int y, 
                     const EFrame frameKind) : frameKind_(frameKind)  {
  data_.resize(2);
  data_[0] = (Scalar)x;
  data_[1] = (Scalar)y;
}


  
template<class Scalar>
Point<Scalar>::Point(const Scalar x, 
                     const Scalar y, 
                     const Scalar z, 
                     const EFrame frameKind) : frameKind_(frameKind)  {
  data_.resize(3);
  data_[0] = x;
  data_[1] = y;
  data_[2] = z;
}


  
template<class Scalar>
Point<Scalar>::Point(const Scalar* dataPtr, 
                     const int dim, 
                     const EFrame frameKind) : frameKind_(frameKind)  {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (dim < 1) || (dim > INTREPID_MAX_DIMENSION) ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid dim argument.");
#endif
  data_.assign(dataPtr, dataPtr + dim);
}



template<class Scalar>
Point<Scalar>::Point(const Teuchos::Array<Scalar>& dataArray,
                     const EFrame frameKind) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (dataArray.size() < 1) || (dataArray.size() > INTREPID_MAX_DIMENSION) ),
                      std::invalid_argument,
                      ">>> ERROR (Point): The data array does not specify a point with valid dimension.");
#endif
  data_ = dataArray;
}

  
template<class Scalar>
inline int Point<Scalar>::getDim() const {
  return data_.size();
}


  
template<class Scalar>
inline EFrame Point<Scalar>::getFrameKind() const {
  return frameKind_;
}


  
template<class Scalar>
inline std::string Point<Scalar>::getFrameName() const {
  return EFrameToString(frameKind_);
}



template<class Scalar>
inline const Teuchos::Array<Scalar> & Point<Scalar>::getCoordinates() const {
  return data_;
}


  
template<class Scalar>
inline void Point<Scalar>::setCoordinates(const Scalar* dataPtr, 
                                          const int targetDim) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this->getDim() != targetDim ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid targetDim argument.");
#endif
  for(int dim = 0; dim < targetDim; dim++){
    data_[dim] = dataPtr[dim];
  }
}


  
template<class Scalar>
inline void Point<Scalar>::setFrameKind(const EFrame newFrameKind) {
  frameKind_ = newFrameKind;
}

  
  
template<class Scalar>
Scalar Point<Scalar>::distance(const Point& endPoint) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this->getDim() != endPoint.getDim() ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid endPoint dimension.");
#endif
  Scalar temp = (Scalar)0.0;
  for(unsigned int i = 0; i < data_.size(); i++) 
  {
    Scalar diff = data_[i] - endPoint.data_[i];
    temp += diff*diff; 
  }
  return (Scalar)std::sqrt(temp);
}


  
template<class Scalar>
Scalar Point<Scalar>::norm(ENorm normType) const{
  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(int i = 0; i < this -> getDim(); i++){
        temp += data_[i]*data_[i]; 
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(data_[0]);
      for(int i = 1; i < this -> getDim(); i++){
        Scalar absData = std::abs(data_[i]);
        if (temp < absData) temp = absData; 
      }
      break;
    case NORM_ONE:
      for(int i = 0; i < this ->getDim(); i++){
        temp += std::abs(data_[i]); 
      }
      break;
    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (Point): Invalid argument normType.");
  }
  return temp;
}


  
template<class Scalar>
const Scalar & Point<Scalar>::operator [] (const int coordinateId) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (coordinateId < 0) || (coordinateId >= this->getDim()) ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid coordinateId argument.");
#endif
  return data_[coordinateId];
}


  
template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator = (const Point<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid right-hand side to '='. Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( this->getDim() != right.getDim() ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid dimension of right-hand side argument to '='.");
#endif
  data_ = right.data_;
  frameKind_ = right.frameKind_;      // changes point type of left!
  return *this;
}



template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator += (const Point<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid right-hand side to '+='. Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( this->getDim() != right.getDim() ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid dimension of right-hand side argument to '+='.");
#endif
  switch (getDim()) {      // does not change point type of left!
    case 3:
      data_[0] += right.data_[0];
      data_[1] += right.data_[1];
      data_[2] += right.data_[2];
      break;
    case 2:
      data_[0] += right.data_[0];
      data_[1] += right.data_[1];
      break;
    case 1:
      data_[0] += right.data_[0];
      break;
    default:
      TEST_FOR_EXCEPTION( ( (getDim() != 1) && (getDim() != 2) && (getDim() != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Point): Invalid point dimension.");
  }
  return *this;
}



template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator ^= (const Point<Scalar>& right)
{                                       // Does not change frameKind_!
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid right-hand side to '^='. Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( (this->getDim() != 3) || (right.getDim() != 3) ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid dimension of operands. Cross product defined only in 3D.");
#endif
  Teuchos::Array<Scalar> tmp(3);
  tmp[0] = data_[1]*right.data_[2]-data_[2]*right.data_[1];
  tmp[1] = data_[2]*right.data_[0]-data_[0]*right.data_[2];
  tmp[2] = data_[0]*right.data_[1]-data_[1]*right.data_[0];
  data_  = tmp;
  return *this;
}



template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator -= (const Point<Scalar>& right)
{                                         // Does not change frameKind_!
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( this == &right ),
                        std::invalid_argument,
                        ">>> ERROR (Point): Invalid right-hand side to '-='. Self-assignment prohibited.");
    TEST_FOR_EXCEPTION( ( this->getDim() != right.getDim() ),
                        std::invalid_argument,
                        ">>> ERROR (Point): Invalid dimension of right-hand side argument to '-='.");
#endif
  switch (getDim()) {
    case 3:
      data_[0] -= right.data_[0];
      data_[1] -= right.data_[1];
      data_[2] -= right.data_[2];
      break;
    case 2:
      data_[0] -= right.data_[0];
      data_[1] -= right.data_[1];
      break;
    case 1:
      data_[0] -= right.data_[0];
      break;
    default:
      TEST_FOR_EXCEPTION( ( (getDim() != 1) && (getDim() != 2) && (getDim() != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Point): Invalid point dimension.");
  }
  return *this;
}



template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator *= (const Scalar left)
{                                         // Does not change frameKind_!
  switch (getDim()) {
    case 3:
      data_[0] *= left;
      data_[1] *= left;
      data_[2] *= left;
      break;
    case 2:
      data_[0] *= left;
      data_[1] *= left;
      break;
    case 1:
      data_[0] *= left;
      break;
    default:
      TEST_FOR_EXCEPTION( ( (getDim() != 1) && (getDim() != 2) && (getDim() != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Point): Invalid point dimension.");
  }
  return *this;
}

//===========================================================================//
//                                                                           //
//           END of member definitions; START friends and related            //
//                                                                           //
//===========================================================================//

template<class Scalar>
const Point<Scalar> operator ^ (const Point<Scalar>& left, const Point<Scalar>& right) {
  Point<Scalar> result(left);
  result ^= right;
  return result;
}



template<class Scalar>
const Point<Scalar> operator + (const Point<Scalar>& left, const Point<Scalar>& right) {
  Point<Scalar> result(left);
  result += right;
  return result;
}



template<class Scalar>
const Point<Scalar> operator - (const Point<Scalar>& left, const Point<Scalar>& right) {
  Point<Scalar> result(left);
  result -= right;
  return result;
}



template<class Scalar>
const Point<Scalar> operator * (const Scalar left, const Point<Scalar>& right) {
  Point<Scalar> result(right);
  result *= left;
  return result;
}



template<class Scalar>
const Scalar operator * (const Point<Scalar>& left, const Point<Scalar>& right) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( left.getDim() != right.getDim() ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Nonmatching dimensions. Points must be of the same size.");
#endif
  Scalar result = (Scalar)0;
  switch (left.getDim()) {
    case 3:
      result = left.getCoordinates()[0]*right.getCoordinates()[0] +
      left.getCoordinates()[1]*right.getCoordinates()[1] +
      left.getCoordinates()[2]*right.getCoordinates()[2];
      break;
    case 2:
      result = left.getCoordinates()[0]*right.getCoordinates()[0] +
      left.getCoordinates()[1]*right.getCoordinates()[1];
      break;
    case 1:
      result = left.getCoordinates()[0]*right.getCoordinates()[0];
      break;
    default:
      TEST_FOR_EXCEPTION( ( (left.getDim() != 1) && (left.getDim() != 2) && (left.getDim() != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Point): Invalid point dimension.");
  }
  return result;
}



template<class Scalar>
std::ostream& operator << (std::ostream& os, const Point<Scalar>& point) {
  // Save the format state of the original ostream os.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os);  

  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.setf(std::ios_base::right);
  int myprec = os.precision();
  
  int dim = point.getDim();
  os << "  " <<  dim << "D "<< point.getFrameName() <<" Point (";
  for (int j=0; j < dim; j++) {
    os << std::setw(myprec+8) << point.getCoordinates()[j];
  }
  os << ")";

  // reset format state of os
  os.copyfmt(oldFormatState);

  return os;
}

// End member, friend, and related function definitions of class Point.



//===========================================================================//
//                                                                           //
//              Member function definitions of the class Matrix.             //
//                                                                           //
//===========================================================================//

template<class Scalar>
Matrix<Scalar>::Matrix(const Matrix<Scalar>& right){
  dim_ = right.dim_;
  elements_.assign(right.elements_.begin(), right.elements_.end());
}



template<class Scalar>
Matrix<Scalar>::Matrix(const int dim){
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (dim < 1) || (dim > INTREPID_MAX_DIMENSION) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Invalid dim argument.");
#endif
  dim_ = dim;
  elements_.assign(dim*dim,Scalar(0));
}



template<class Scalar>
Matrix<Scalar>::Matrix(const Scalar* dataPtr, const int dim){
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (dim < 1) || (dim > INTREPID_MAX_DIMENSION) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Invalid dim argument.");
#endif
  dim_ = dim;
  elements_.resize(dim_*dim_);
  for(int i=0; i < dim_; ++i){
    for(int j=0; j < dim_; ++j){
      elements_[i*dim_+j]=*(dataPtr+i*dim_+j);
    }
  }
}



template<class Scalar>
void Matrix<Scalar>::setElements(const Scalar* dataPtr, const int dim){
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( dim_ != dim ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Invalid dim argument.");
#endif
  for(int i=0; i < dim_; ++i){
    for(int j=0; j < dim_; ++j){
      elements_[i*dim_+j]=*(dataPtr+i*dim_+j);
    }
  }
}



template<class Scalar>
void Matrix<Scalar>::resize(const int dim){
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (dim < 1) || (dim > INTREPID_MAX_DIMENSION) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Invalid dim argument.");
#endif
  dim_ = dim;
  elements_.assign(dim_*dim_,Scalar(0));
}



template<class Scalar>
int Matrix<Scalar>::getDim() const {
  return dim_;
}



template<class Scalar>
Scalar Matrix<Scalar>::getElement(int rowID, int colID) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (rowID < 0) || (rowID >= dim_) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): row ID out of range.");
  TEST_FOR_EXCEPTION( ( (colID < 0) || (colID >= dim_) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): column ID out of range.");
#endif
  return elements_[rowID*dim_+colID];
}



template<class Scalar>
Point<Scalar> Matrix<Scalar>::getColumn(int colID) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (colID < 0) || (colID >= dim_) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): column ID out of range.");
#endif
  Scalar elVec[INTREPID_MAX_DIMENSION];
  for (int i=0; i<dim_; i++)
    elVec[i] = elements_[i*dim_+colID];
  return Point<Scalar>(elVec, dim_);
}



template<class Scalar>
Point<Scalar> Matrix<Scalar>::getRow(int rowID) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (rowID < 0) || (rowID >= dim_) ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): row ID out of range.");
#endif
  Scalar elVec[INTREPID_MAX_DIMENSION];
  for (int i=0; i<dim_; i++)
    elVec[i] = elements_[rowID*dim_+i];
  return Point<Scalar>(elVec, dim_);
}



template<class Scalar>
Matrix<Scalar> Matrix<Scalar>::getTranspose() const {
  Matrix<Scalar> transpose(dim_);
  for(int i=0; i < dim_; i++){
    transpose.elements_[i*dim_+i]=elements_[i*dim_+i];    // Set diagonal elements
    for(int j=i+1; j < dim_; j++){
      transpose.elements_[i*dim_+j]=elements_[j*dim_+i];  // Set off-diagonal elements
      transpose.elements_[j*dim_+i]=elements_[i*dim_+j];
    }
  }
  return transpose;
}



template<class Scalar>
Matrix<Scalar> Matrix<Scalar>::getInverse() const {
  int i,j,rowID = 0, colID = 0;
  int rowperm[3]={0,1,2};
  int colperm[3]={0,1,2}; // Complete pivoting
  Scalar emax(0), determinant(0);
  
  if (dim_==3) {
    for(i=0; i < 3; ++i){
      for(j=0; j < 3; ++j){
        if( std::abs( elements_[i*dim_+j] ) >  emax){
          rowID = i;  colID = j; emax = std::abs( elements_[i*dim_+j] );
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
    Scalar Inverse[3][3], B[3][3], S[2][2], Bi[3][3]; // B=rowperm elements_ colperm, S=Schur complement(Boo)
    for(i=0; i < 3; ++i){
      for(j=0; j < 3; ++j){
        B[i][j] = elements_[rowperm[i]*dim_+colperm[j]];
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
          Inverse[i][j] = Bi[colperm[i]][rowperm[j]]; // return inverse in place
        }
      }
      return Matrix<Scalar>(&Inverse[0][0],dim_);
  }
  else if (dim_==2) {
    Scalar Inverse[2][2];
    determinant= det();
    Inverse[0][0] =   elements_[3]/determinant;
    Inverse[0][1] = - elements_[1]/determinant;
    //
    Inverse[1][0] = - elements_[2]/determinant;
    Inverse[1][1] =   elements_[0]/determinant;
    return Matrix<Scalar>(&Inverse[0][0],dim_);
  }
  else { // dim==1
    Scalar Inverse[1][1];
    determinant= det();
    Inverse[0][0] = (Scalar)1 / elements_[0];
    return Matrix<Scalar>(&Inverse[0][0],dim_);
  }
}

template<class Scalar>
void Matrix<Scalar>::transpose() {
  Scalar temp(0);
  for(int i=0; i < dim_; ++i){
    for(int j=i+1; j < dim_; ++j){
      temp=elements_[i*dim_+j];
      elements_[i*dim_+j]=elements_[j*dim_+i];            // Leave diagonal elements alone!
      elements_[j*dim_+i]=temp;
    }
  }
}



template<class Scalar>
Scalar Matrix<Scalar>::det() const {
  int i,j,rowID = 0;
  int colID = 0;
  int rowperm[3]={0,1,2};
  int colperm[3]={0,1,2}; // Complete pivoting
  Scalar emax(0), determinant(0);

  switch (dim_) {
    case 3:
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( elements_[i*dim_+j] ) >  emax){
            rowID = i;  colID = j; emax = std::abs( elements_[i*dim_+j] );
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
        Scalar B[3][3], S[2][2]; // B=rowperm elements_ colperm, S=Schur complement(Boo)
        for(i=0; i < 3; ++i){
          for(j=0; j < 3; ++j){
            B[i][j] = elements_[rowperm[i]*dim_+colperm[j]];
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
    case 2:
      determinant = elements_[0]*elements_[3]-
                    elements_[1]*elements_[2];
      break;
    case 1:
      determinant = elements_[0];
      break;
    default:
      TEST_FOR_EXCEPTION( ( (dim_ != 1) && (dim_ != 2) && (dim_ != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  }
  return determinant;
}



template<class Scalar>
Scalar Matrix<Scalar>::norm(ENorm normType) const {
  Scalar result(0);
  Scalar temp(0);

  switch(normType) {
    case NORM_ONE:  // std::max column sum of the absolute values
      // Set result equal to 1-norm of the first column
      result = this -> getColumn(0).norm(NORM_ONE);
      
      // If dim_ > 1 compare 1-norm of first column with 1-norm of second column
      if(dim_ > 1){
        temp = getColumn(1).norm(NORM_ONE);
        result = (temp > result) ? temp : result;
      }
      
      // Now result holds the larger of the 1-norms of columns 1 and 2. If dim_=3 compare
      // this number with the 1-norm of the 3rd column:
      if(dim_ == 3) {
        temp = getColumn(2).norm(NORM_ONE);
        result = (temp > result) ? temp : result; 
      }
      break;
    case NORM_INF:  // std::max row sum of absolute values: apply above algorithm to rows
      // Set result equal to 1-norm of the first row
      result = this -> getRow(0).norm(NORM_ONE);
      
      // If dim_ > 1 compare 1-norm of first row with 1-norm of second row
      if(dim_ > 1){
        temp = getRow(1).norm(NORM_ONE);
        result = (temp > result) ? temp : result;
      }
      
      // Now result holds the larger of the 1-norms of rows 1 and 2. If dim_=3 compare
      // this number with the 1-norm of the 3rd row:
      if(dim_ == 3) {
        temp = getRow(2).norm(NORM_ONE);
        result = (temp > result) ? temp : result; 
      }
      break;
    case NORM_FRO:  // square root of sum of all element squares 
      for(int i = 0; i < dim_; i++){
        for(int j = 0; j < dim_; j++){
          Scalar tmp = elements_[i*dim_+j];
          result += tmp*tmp;
        }
      }
      result = std::sqrt(result);
      break;
    case NORM_TWO:
    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_ONE) && (normType != NORM_INF) && (normType != NORM_FRO) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Matrix norm not implemented.");
  }
  return result;
}



template<class Scalar>
void Matrix<Scalar>::invert() {
  Matrix<Scalar> tempMatrix = this -> getInverse();
  *this = tempMatrix;
}


template<class Scalar>
void Matrix<Scalar>::scalarMultiply(const Scalar& scalar) {
  for(int i = 0; i < dim_; i++) {
    for(int j = 0; j < dim_; j++){
      elements_[i*dim_+j] *= scalar; 
    }
  }
}


template<class Scalar>
void Matrix<Scalar>::multiplyLeft(Scalar*        matVec,
                                  const Scalar*  vec) const {
  for (int i = 0; i < dim_; i++) {
    matVec[i] = 0.0;
    for (int j = 0; j < dim_; j++) {
      matVec[i] += elements_[i*dim_+j]*vec[j];
    }
  }
}


template<class Scalar>
void Matrix<Scalar>::rowMultiply(Scalar &                       rowVec,
                                 const int                      rowId,
                                 const FieldContainer<Scalar>&  vec,
                                 const int                      vecBegin) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (rowId < 0) || (dim_ < rowId) ), std::invalid_argument,
                      ">>> ERROR (Matrix): Invalid row number.");
  // vecIndex must be at least 0 and such that there are at least dim more values left in the container
  TEST_FOR_EXCEPTION( ( (vecBegin < 0) || ( (vec.getSize() - dim_) < vecBegin) ), std::invalid_argument,
                      ">>> ERROR (Matrix): Invalid index for the first vector element.");
#endif
  rowVec = 0;
  for (int j = 0; j < dim_; j++) {
    rowVec += elements_[rowId*dim_ + j]*vec[vecBegin + j];
  }
}



template<class Scalar>
void Matrix<Scalar>::multiplyLeft(Teuchos::Array<Scalar> &        matVec,
                                  const Teuchos::Array<Scalar> &  vec) const {
  int dim = vec.size();
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( dim_ != dim ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Matrix and vector dimensions do not match.");
#endif
  for (int i = 0; i < dim_; i++) {
    matVec[i] = 0.0;
    for (int j = 0; j < dim_; j++) {
      matVec[i] += elements_[i*dim_+j]*vec[j];
    }
  }
}



template<class Scalar> 
const Scalar& Matrix<Scalar>::operator () (const int rowId, const int colId) const{
  return elements_[rowId*dim_+colId]; 
}



template<class Scalar>
inline Matrix<Scalar> & Matrix<Scalar>::operator = (const Matrix<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( dim_ != right.dim_ ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Dimensions do not match in assignment statement.");
#endif
  elements_ = right.elements_;
  return *this;
}



template<class Scalar>
inline Matrix<Scalar> & Matrix<Scalar>::operator += (const Matrix<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( dim_ != right.dim_ ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Dimensions do not match in in-place addition.");
#endif
  switch (dim_) {      
    case 3:
      elements_[0] += right.elements_[0];
      elements_[1] += right.elements_[1];
      elements_[2] += right.elements_[2];
      //
      elements_[3] += right.elements_[3];
      elements_[4] += right.elements_[4];
      elements_[5] += right.elements_[5];
      //
      elements_[6] += right.elements_[6];
      elements_[7] += right.elements_[7];
      elements_[8] += right.elements_[8];
      break;
    case 2:
      elements_[0] += right.elements_[0];
      elements_[1] += right.elements_[1];
      //
      elements_[2] += right.elements_[2];
      elements_[3] += right.elements_[3];
      break;
    case 1:
      elements_[0] += right.elements_[0];
      break;
    default:
      TEST_FOR_EXCEPTION( ( (dim_ != 1) && (dim_ != 2) && (dim_ != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  }

  return *this;
}



template<class Scalar>
inline Matrix<Scalar> & Matrix<Scalar>::operator -= (const Matrix<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( dim_ != right.dim_ ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Dimensions do not match in in-place subtraction.");
#endif
  switch (dim_) {      
    case 3:
      elements_[0] -= right.elements_[0];
      elements_[1] -= right.elements_[1];
      elements_[2] -= right.elements_[2];
      //
      elements_[3] -= right.elements_[3];
      elements_[4] -= right.elements_[4];
      elements_[5] -= right.elements_[5];
      //
      elements_[6] -= right.elements_[6];
      elements_[7] -= right.elements_[7];
      elements_[8] -= right.elements_[8];
      break;
    case 2:
      elements_[0] -= right.elements_[0];
      elements_[1] -= right.elements_[1];
      //
      elements_[2] -= right.elements_[2];
      elements_[3] -= right.elements_[3];
      break;
    case 1:
      elements_[0] -= right.elements_[0];
      break;
    default:
      TEST_FOR_EXCEPTION( ( (dim_ != 1) && (dim_ != 2) && (dim_ != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  }

  return *this;
}



template<class Scalar>
inline Matrix<Scalar> & Matrix<Scalar>::operator *= (const Scalar left)
{
  switch (dim_) {      
    case 3:
      elements_[0] *= left;
      elements_[1] *= left;
      elements_[2] *= left;
      //
      elements_[3] *= left;
      elements_[4] *= left;
      elements_[5] *= left;
      //
      elements_[6] *= left;
      elements_[7] *= left;
      elements_[8] *= left;
      break;
    case 2:
      elements_[0] *= left;
      elements_[1] *= left;
      //
      elements_[2] *= left;
      elements_[3] *= left;
      break;
    case 1:
      elements_[0] *= left;
      break;
    default:
      TEST_FOR_EXCEPTION( ( (dim_ != 1) && (dim_ != 2) && (dim_ != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  }

  return *this;
}

//===========================================================================//
//                                                                           //
//  END of member definitions for class Matrix; START friends and related    //
//                                                                           //
//===========================================================================//

template<class Scalar>
const Matrix<Scalar> operator + (const Matrix<Scalar>& left, const Matrix<Scalar>& right){
  Matrix<Scalar> result(left);
  result += right;
  return result;
}



template<class Scalar>
const Matrix<Scalar> operator - (const Matrix<Scalar>& left, const Matrix<Scalar>& right){
  Matrix<Scalar> result(left);
  result -= right;
  return result;
}



template<class Scalar>
const Matrix<Scalar> operator * (const Scalar left, const Matrix<Scalar>& right) {
  Matrix<Scalar> result(right);
  result *= left;
  return result;
}



template<class Scalar>
const Point<Scalar> operator * (const Matrix<Scalar>& mat, const Point<Scalar>& vec) {
  Scalar Product[3];
  int dim = vec.getDim();
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( mat.getDim() != dim ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Matrix and vector dimensions do not match.");
#endif
  for (int i = 0; i < dim; ++i) {
    Product[i] = 0.0;
    for (int j = 0; j < dim; ++j) {
      Product[i] += mat.getElement(i,j)*vec.getCoordinates()[j];
    }
  }
  return Point<Scalar>(Product, dim);
}



template<class Scalar>
const Point<Scalar> operator * ( const Point<Scalar>& vec, const Matrix<Scalar>& mat) {
  Scalar Product[3];
  int dim = vec.getDim();
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( mat.getDim() != dim ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Matrix and vector dimensions do not match.");
#endif
  for (int i = 0; i < dim; ++i)
  {
    Product[i] = 0.0;
    for (int j = 0; j < dim; ++j)
      Product[i] += mat.getElement(j,i)*vec.getCoordinates()[j];
  }
  return Point<Scalar>(Product, dim);
}



template<class Scalar>
const Matrix<Scalar> operator * (const Matrix<Scalar>& lmat, const Matrix<Scalar>& rmat) {
  int dim = lmat.getDim();
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( rmat.getDim() != dim ),
                      std::invalid_argument,
                      ">>> ERROR (Matrix): Matrix dimensions do not match.");
#endif
  Scalar Product[9];
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j){
      Product[i*dim+j] = (Scalar)0;
      for (int k = 0; k < dim; ++k){
        Product[i*dim+j] += lmat.getElement(i,k)*rmat.getElement(k,j);
      }
    }
  }
  return Matrix<Scalar>(Product, dim);
}



template<class Scalar>
std::ostream& operator << (std::ostream& os, const Matrix<Scalar>& matrix) {
  // Save the format state of the original ostream os.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os);

  short dim = matrix.getDim();

  int myprec = os.precision();
  int outWidth = myprec+10;

  os << " " <<  dim << "D "<< "Matrix:\n";

  os.setf(std::ios_base::left);
  os << "           ";
  switch(dim){
    case 1:
      os << std::setw(outWidth) << "Col 0" << std::endl;
      break;
    case 2:
      os << std::setw(outWidth) << "Col 0" << std::setw(outWidth) << "Col 1" << std::endl;
      break;
    case 3:
      os << std::setw(outWidth) << "Col 0" << std::setw(outWidth) << "Col 1" << std::setw(outWidth) << "Col 2" << std::endl;
      break;
    default:
      TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  }

  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.setf(std::ios_base::right);

  for(int i = 0; i < dim; ++i){
    os << " Row "<< i<< " ";
    for(int j = 0; j < dim; ++j) {
      os << std::setw(outWidth) << matrix.getElement(i,j);
    }
    os << std::endl;
  }

  // reset format state of os
  os.copyfmt(oldFormatState);

  return os;
}

// End member, friend, and related function definitions of class Matrix.

} // namespace Intrepid
