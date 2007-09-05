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


/** \file
    \brief Definition file for auxiliary classes providing linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


#include <cfloat>

namespace Intrepid {

// Member function definitions of the class Point.

template<class Scalar>
Point<Scalar>::Point(const Point<Scalar>& right) {
    data.resize(right.getDim());
    data = right.data;
    point_type = right.point_type;
}

template<class Scalar>
Point<Scalar>::Point(const int dim, 
                     const PointType point_type_): point_type(point_type_) {
  assert(0<dim && dim<=INTREPID_MAX_DIMENSION);
  data.assign(dim, (Scalar)0);
}

template<class Scalar>
Point<Scalar>::Point(const Scalar x,
                     const PointType point_type_) : point_type(point_type_) {
  data.resize(1);
  data[0] = x;
}

template<class Scalar>
Point<Scalar>::Point(const Scalar x, 
                     const Scalar y, 
                     const PointType point_type_) : point_type(point_type_)  {
  data.resize(2);
  data[0] = x;
  data[1] = y;
}

template<class Scalar>
Point<Scalar>::Point(const int x, 
                     const int y, 
                     const PointType point_type_) : point_type(point_type_)  {
  data.resize(2);
  data[0] = (Scalar)x;
  data[1] = (Scalar)y;
}

template<class Scalar>
Point<Scalar>::Point(const Scalar x, 
                     const Scalar y, 
                     const Scalar z, 
                     const PointType point_type_) : point_type(point_type_)  {
  data.resize(3);
  data[0] = x;
  data[1] = y;
  data[2] = z;
}

template<class Scalar>
Point<Scalar>::Point(const Scalar* dataptr, 
                     const int dim, 
                     const PointType point_type_) : point_type(point_type_)  {
  assert(0<dim && dim<=INTREPID_MAX_DIMENSION);
  data.assign(dataptr, dataptr+dim);
}

template<class Scalar>
inline int Point<Scalar>::getDim() const {
  return data.size();
}

template<class Scalar>
inline PointType Point<Scalar>::getType() const {
  return point_type;
}

template<class Scalar>
inline const char* Point<Scalar>::getName() const {
  return PointNames[point_type];
}

template<class Scalar>
inline const std::vector<Scalar> & Point<Scalar>::getData() const {
  return data;
}

template<class Scalar>
inline void Point<Scalar>::setData(const Scalar* dataptr, 
                                   const int target_dim) {
  assert(this -> getDim() == target_dim);
  for(int dim = 0; dim < target_dim; dim++){
    data[dim] = dataptr[dim];
  }
}

template<class Scalar>
inline void Point<Scalar>::setType(const PointType target_type) {
  point_type = target_type;
}


template<class Scalar>
ErrorCode Point<Scalar>::inRefCell(const CellType cell_type) const{
  ErrorCode error_code = SUCCESS;
  if(point_type != REFERENCE) error_code = FAILURE;
  switch(cell_type) {
    case EDGE:
      if( !(-1.0 <= data[0] && data[0] <= 1.0)) error_code = FAILURE;
      break;
    case TRI:{
      Scalar threshold = 3. * std::abs( DBL_EPSILON );
      Scalar distance = max( max( -data[0], -data[1] ), data[0] + data[1] - 1 );
      if( distance > threshold ) error_code = FAILURE;
      break;
      }
    case QUAD:
      if(!((-1.0 <= data[0] && data[0] <= 1.0) && \
           (-1.0 <= data[1] && data[1] <= 1.0))) error_code = FAILURE;   
      break;
    case TET:{
      Scalar threshold = 4. * std::abs( DBL_EPSILON );
      Scalar distance = max(  max(-data[0],-data[1]), max(-data[2],data[0] + data[1] + data[2] - 1)  );

      if( distance > threshold ) error_code = FAILURE;
      break;
      }
    case HEX:
      if(!((-1.0 <= data[0] && data[0] <= 1.0) && \
           (-1.0 <= data[1] && data[1] <= 1.0) && \
           (-1.0 <= data[2] && data[2] <= 1.0))) error_code = FAILURE;
        break;
    case PRISM:
      break;
    case PYRAMID:
      break;
    default:
      std::cerr<< "Point::inRefCell error: invalid cell type \n";
      exit(1);
  }
  return error_code;
}


template<class Scalar>
const Scalar & Point<Scalar>::operator [] (const int coordinate_id_) const {
  return data[coordinate_id_];
}


template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator = (const Point<Scalar>& right)
{
    assert(this != &right);            // check for self-assignment
    assert(getDim()==right.getDim());  // dimensions must match    
    data = right.data;
    point_type = right.point_type;      // changes point type of left!
    return *this;
}

template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator += (const Point<Scalar>& right)
{
    assert(this != &right);  // check for self-assignment
    switch (getDim()) {      // does not change point type of left!
    case 3:
      data[0] += right.data[0];
      data[1] += right.data[1];
      data[2] += right.data[2];
      return *this;
    case 2:
      data[0] += right.data[0];
      data[1] += right.data[1];
      return *this;
    default: // case 1
      data[0] += right.data[0];
      return *this;
    }
}

template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator ^= (const Point<Scalar>& right)
{                                       // Does not change point_type!
    assert(this != &right);             // check for self-assignment
    assert(getDim()==3);                // cross product only in 3D
    std::vector<Scalar> tmp(3);         
    tmp[0] = data[1]*right.data[2]-data[2]*right.data[1];
    tmp[1] = data[2]*right.data[0]-data[0]*right.data[2];
    tmp[2] = data[0]*right.data[1]-data[1]*right.data[0];
    data   = tmp;
    return *this;
}

template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator -= (const Point<Scalar>& right)
{                                         // Does not change point_type!
    assert(this != &right);               // check for self-assignment
    switch (getDim()) {
    case 3:
      data[0] -= right.data[0];
      data[1] -= right.data[1];
      data[2] -= right.data[2];
      return *this;
    case 2:
      data[0] -= right.data[0];
      data[1] -= right.data[1];
      return *this;
    default: // case 1
      data[0] -= right.data[0];
      return *this;
    }
}

template<class Scalar>
inline Point<Scalar> & Point<Scalar>::operator *= (const Scalar left)
{                                         // Does not change point_type!
    switch (getDim()) {
    case 3:
      data[0] *= left;
      data[1] *= left;
      data[2] *= left;
      return *this;
    case 2:
      data[0] *= left;
      data[1] *= left;
      return *this;
    default: // case 1
      data[0] *= left;
      return *this;
    }
}

// END of member definitions; START friends and related

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
  assert(left.getDim()==right.getDim());
  Scalar result;
  switch (left.getDim()) {
  case 3:
    result = left.getData()[0]*right.getData()[0] +
             left.getData()[1]*right.getData()[1] +
             left.getData()[2]*right.getData()[2];
    return result;
  case 2:
    result = left.getData()[0]*right.getData()[0] +
             left.getData()[1]*right.getData()[1];
    return result;
  default: // case 1
    result = left.getData()[0]*right.getData()[0];
    return result;
  }
}

template<class Scalar>
std::ostream& operator << (std::ostream& os, const Point<Scalar>& point) {
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.precision(6);

  int dim = point.getDim();
  os << "  " <<  dim << "D "<< point.getName() <<" Point (";
  for (int j=0; j < dim; j++) {
    os << std::setw(14) << point.getData()[j];
  }
  std::cout << ")";
  return os;
}

// End member, friend, and related function definitions of class Point.



// Member function definitions of the class LinearMap.

template<class Scalar>
LinearMap<Scalar>::LinearMap(const LinearMap<Scalar>& right){
    elements_.assign(right.elements_.begin(), right.elements_.end());
}

template<class Scalar>
LinearMap<Scalar>::LinearMap(const int dim){
    assert(0<dim && dim<=INTREPID_MAX_DIMENSION);
    std::vector<Scalar> tmp(dim);
    elements_.assign(dim,tmp);
}

template<class Scalar>
LinearMap<Scalar>::LinearMap(const Scalar* dataptr, const int dim){
    assert(0<dim && dim<=INTREPID_MAX_DIMENSION);
    std::vector<Scalar> tmp(dim);
    elements_.assign(dim,tmp);
    for(int i=0; i < dim; ++i){
	for(int j=0; j < dim; ++j){
	    elements_[i][j]=*(dataptr+i*dim+j);
	}
    }
}

template<class Scalar>
void LinearMap<Scalar>::setData(const Scalar* dataptr, const int dim){
  assert(this -> getDim() == dim);
  for(int i=0; i < dim; ++i){
	for(int j=0; j < dim; ++j){
      elements_[i][j]=*(dataptr+i*dim+j);
	}
  }
}


template<class Scalar>
inline LinearMap<Scalar> & LinearMap<Scalar>::operator = (const LinearMap<Scalar>& right)
{
    assert(this != &right);                            // check for self-assignment
    assert(getDim()==right.getDim());  // dimensions must match
    elements_ = right.elements_;
    return *this;
}

template<class Scalar>
int LinearMap<Scalar>::getDim() const {
  return elements_.size();
}

template<class Scalar>
Scalar LinearMap<Scalar>::getElement(int rowID, int colID) const {
    assert(0<=rowID && rowID<(int)elements_.size());
    assert(0<=colID && colID<(int)elements_.size());
    return( elements_[rowID][colID]);
}

template<class Scalar>
Point<Scalar> LinearMap<Scalar>::getColumn(int colID) const {
    assert(0<=colID && colID<elements_.size());
    switch (getDim()) {
    case 3:
      return Point<Scalar>(elements_[0][colID],
                           elements_[1][colID],
		           elements_[2][colID]);
    case 2:
      return Point<Scalar>(elements_[0][colID],
                           elements_[1][colID]);
    default: // case 1
      return Point<Scalar>(elements_[0][colID]);
    }
}

template<class Scalar>
Point<Scalar> LinearMap<Scalar>::getRow(int rowID) const {
    assert(0<=rowID && rowID<elements_.size());
    switch (getDim()) {
    case 3:
      return Point<Scalar>(elements_[rowID][0],
                           elements_[rowID][1],
                           elements_[rowID][2]);
    case 2:
      return Point<Scalar>(elements_[rowID][0],
                           elements_[rowID][1]);
    case 1: // case 1
      return Point<Scalar>(elements_[rowID][0]);
    }
}

template<class Scalar>
Scalar LinearMap<Scalar>::Det() const {
    int i,j,rowID = 0, colID = 0, rowperm[3]={0,1,2}, colperm[3]={0,1,2}; // Complete pivoting
    Scalar emax(0), determinant(0);
    switch (getDim()) {
    case 3:
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( elements_[i][j] ) >  emax){
             rowID = i;  colID = j; emax = std::abs( elements_[i][j] );
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
            B[i][j] = elements_[rowperm[i]][colperm[j]];
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
      return(determinant); // vulnerable to underflow and overflow
    case 2:
      return(elements_[0][0]*elements_[1][1]-
             elements_[0][1]*elements_[1][0]);
    default: // case 1
      return(elements_[0][0]);
    }
}

template<class Scalar>
LinearMap<Scalar> LinearMap<Scalar>::getTranspose() const {
    LinearMap<Scalar> transpose(getDim());
    for(int i=0; i < getDim(); ++i){
        transpose.elements_[i][i]=elements_[i][i];
	for(int j=i+1; j < getDim(); ++j){
            transpose.elements_[i][j]=elements_[j][i];
            transpose.elements_[j][i]=elements_[i][j];
	}
    }
    return transpose;
}

template<class Scalar>
void LinearMap<Scalar>::Transpose() {
    Scalar temp(0);
    for(int i=0; i < getDim(); ++i){
	for(int j=i+1; j < getDim(); ++j){
	    temp=elements_[i][j];
	    elements_[i][j]=elements_[j][i];
	    elements_[j][i]=temp;
	}
    }
}

template<class Scalar>
LinearMap<Scalar> LinearMap<Scalar>::getInverse() const {
    int dim = getDim(), i,j,rowID = 0, colID = 0, rowperm[3]={0,1,2}, colperm[3]={0,1,2}; // Complete pivoting
    Scalar emax(0), determinant(0);

    if (dim==3) {
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( elements_[i][j] ) >  emax){
             rowID = i;  colID = j; emax = std::abs( elements_[i][j] );
          }
        }
      }
      if( emax == 0 ) std::cerr <<" LinearMap getInverse: Zero matrix\n";
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
          B[i][j] = elements_[rowperm[i]][colperm[j]];
        }
      }
      B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
      for(i=0; i < 2; ++i){
        for(j=0; j < 2; ++j){
          S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
        }
      }
      Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
      if( detS == 0 ) std::cerr <<" LinearMap getInverse : Singular matrix\n";

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
      return LinearMap<Scalar>(&Inverse[0][0],dim);
    }
    else if (dim==2) {
      Scalar Inverse[2][2];
      determinant= Det();
      Inverse[0][0] =   elements_[1][1]/determinant;
      Inverse[0][1] = - elements_[0][1]/determinant;
      //
      Inverse[1][0] = - elements_[1][0]/determinant;
      Inverse[1][1] =   elements_[0][0]/determinant;
      return LinearMap<Scalar>(&Inverse[0][0],dim);
    }
    else { // dim==1
      Scalar Inverse[1][1];
      determinant= Det();
      Inverse[0][0] = (Scalar)1 / elements_[0][0];
      return LinearMap<Scalar>(&Inverse[0][0],dim);
    }
}

template<class Scalar>
void LinearMap<Scalar>::Invert() {
    int dim = getDim();
    Scalar determinant = Det();

    if (dim==3) {
      Scalar Inverse[3][3];
      Inverse[0][0] = (-elements_[1][2]*elements_[2][1] + elements_[1][1]*elements_[2][2])/determinant;
      Inverse[0][1] = ( elements_[0][2]*elements_[2][1] - elements_[0][1]*elements_[2][2])/determinant;
      Inverse[0][2] = (-elements_[0][2]*elements_[1][1] + elements_[0][1]*elements_[1][2])/determinant;
      //
      Inverse[1][0] = ( elements_[1][2]*elements_[2][0] - elements_[1][0]*elements_[2][2])/determinant;
      Inverse[1][1] = (-elements_[0][2]*elements_[2][0] + elements_[0][0]*elements_[2][2])/determinant;
      Inverse[1][2] = ( elements_[0][2]*elements_[1][0] - elements_[0][0]*elements_[1][2])/determinant;
      //
      Inverse[2][0] = (-elements_[1][1]*elements_[2][0] + elements_[1][0]*elements_[2][1])/determinant;
      Inverse[2][1] = ( elements_[0][1]*elements_[2][0] - elements_[0][0]*elements_[2][1])/determinant;
      Inverse[2][2] = (-elements_[0][1]*elements_[1][0] + elements_[0][0]*elements_[1][1])/determinant;
      for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
          elements_[i][j] = Inverse[i][j];
    }
    else if (dim==2) {
      Scalar Inverse[2][2];
      Inverse[0][0] =   elements_[1][1]/determinant;
      Inverse[0][1] = - elements_[0][1]/determinant;
      //
      Inverse[1][0] = - elements_[1][0]/determinant;
      Inverse[1][1] =   elements_[0][0]/determinant;
      for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
          elements_[i][j] = Inverse[i][j];
    }
    else { // dim==1
      elements_[0][0] = (Scalar)1 / elements_[0][0];
    }
}

// END of member definitions; START friends and related

template<class Scalar>
const Point<Scalar> operator * (const LinearMap<Scalar>& mat, const Point<Scalar>& vec) {
    Scalar Product[3];
    assert(mat.getDim()==vec.getDim());
    int dim = vec.getDim();
    for (int i = 0; i < dim; ++i)
    {
        Product[i] = 0.0;
	for (int j = 0; j < dim; ++j)
		Product[i] += mat.getElement(i,j)*vec.getData()[j];
    }
    return Point<Scalar>(Product, dim);
}

template<class Scalar>
const LinearMap<Scalar> operator * (const LinearMap<Scalar>& lmat, const LinearMap<Scalar>& rmat) {
    assert(lmat.getDim()==rmat.getDim());
    int dim = lmat.getDim();
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
    return LinearMap<Scalar>(Product, dim);
}

template<class Scalar>
std::ostream& operator << (std::ostream& os, const LinearMap<Scalar>& matrix)
{
  short dim = matrix.getDim();
  os  << "\nLinearMap info: ambient dimension = " << dim << "D\n";
  switch(dim){
    case 1:
      os << "               Col 0" << std::endl;
      break;
    case 2:
      os << "               Col 0          Col 1" << std::endl;
      break;
    case 3:
      os << "               Col 0          Col 1          Col 2" << std::endl;
      break;
    default:
      os << "LinearMap error: invalid ambient dimension\n";
      exit(1);
  }
  os 	<< std::setprecision(6) << std::setiosflags(std::ios::scientific);
  for(int i = 0; i < dim; ++i){
    os << " Row "<< i<< " ";
    for(int j = 0; j < dim; ++j) {
      os << std::setiosflags(std::ios::right) << std::setw(16) << matrix.getElement(i,j);
    }
    os << std::endl;
  }
  os << std::resetiosflags(std::ios::floatfield);
  return os;
}

// End member, friend, and related function definitions of class LinearMap.


} // namespace Intrepid
