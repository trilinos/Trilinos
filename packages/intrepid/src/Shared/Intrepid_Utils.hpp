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

/** \file   Intrepid_Utils.hpp
    \brief  Intrepid utilities.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_UTILS_HPP
#define INTREPID_UTILS_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Rank.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
namespace Intrepid {

/***************************************************************************************************
 ***************************************************************************************************
 **                                                                                               **
 **  Declarations of non-templated utility functions for order and cardinality of operators       **
 **                                                                                               **
 ***************************************************************************************************
 ***************************************************************************************************/


  /** \brief  Returns the rank of fields in a function space of the specified type.

  Field rank is defined as the number of indices needed to specify function value and
  equals 0, 1,or 2 for scalars, vectors and tensors, respectively. The scalar field
  spaces in Intrepid are FUNCTION_SPACE_HGRAD and FUNCTION_SPACE_HVOL. The vector field
  spaces are FUNCTION_SPACE_HCURL, FUNCTION_SPACE_HDIV and FUNCTION_SPACE_VECTOR_HGRAD.
  FUNCTION_SPACE_TENSOR_HGRAD contains rank-2 tensors.

  \param  spaceType         [in]     -  function space type
  \return rank of the fields in the specified function space
  */
  int getFieldRank(const EFunctionSpace spaceType);



  /** \brief  Returns rank of an operator.

  When an operator acts on a field of a certain rank, the result can be a field with the
  same or a different rank. Operator rank is defined the difference between the ranks of
  the output field and the input field:
  \verbatim
  Rank(OPERATOR) = Rank(OPERATOR(FIELD)) - Rank(FIELD)
  \endverbatim
  Therefore, operator rank allows us to figure out the rank of the result:
  \verbatim
  Rank(OPERATOR(FIELD)) = Rank(FIELD) + Rank(OPERATOR)
  \endverbatim
  and provides means to size properly arrays for output results. The following table
  summarizes operator ranks (~ denotes undefined, below slash means 3D).
  By default, in 1D any operator other than VALUE has rank 1, i.e., GRAD, CURL and DIV
  reduce to d/dx and Dk are the higher-order derivatives d^k/dx^k. Only scalar functions
  are allowed in 1D.

  \verbatim
  |========|======|============================|=========|==========|==========|==========|
  | field  | rank |  FUNCTION_SPACE_[type]     |  VALUE  | GRAD, Dk |   CURL   |    DIV   |
  |--------|------|----------------------------|---------|----------|----------|----------|
  | scalar |   0  |  HGRAD, HVOL               |    0    |     1    | 3-dim/~  |     ~    |
  | vector |   1  |  HCURL, HDIV, VECTOR_HGRAD |    0    |     1    | dim - 3  |    -1    |
  | tensor |   2  |  TENSOR_HGRAD              |    0    |     1    | dim - 3  |    -1    |
  |--------|------|----------------------------|---------|----------|----------|----------|
  |   1D   |   0  |  HGRAD, HVOL only          |    0    |     1    |     1    |     1    |
  |=======================================================================================|
  \endverbatim

  \param  spaceType        [in]    - function space type
  \param  operatorType     [in]    - the operator acting on the specified function space
  \param  spaceDim         [in]    - spatial dimension
  \return rank of the operator as defined in the table
  */
  int getOperatorRank(const EFunctionSpace spaceType,
                      const EOperator      operatorType,
                      const int            spaceDim);



  /** \brief  Returns order of an operator.

    \param  operatorType       [in]    - type of the operator whose order we want to know
    \return result ranges from 0 (for OPERATOR_VALUE) to 10 (OPERATOR_D10)
    */
  int getOperatorOrder(const EOperator operatorType);



  /** \brief  Returns the ordinal of a partial derivative of order k based on the multiplicities of
    the partials dx, dy, and dz.

    By default, any implementation of Intrepid::Basis method returns partials of order k
  (specified by OPERATOR_Dk) as a multiset ordered by the lexicographical order of the
    partial derivatives multiplicities. For example, the 10 derivatives of order 3 in 3D
    are enumerated as:
    \verbatim
    D3={(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),(0,3,0),(0,2,1),(0,1,2),(0,0,3)}
  \endverbatim
    The enumeration formula for this lexicographical order is
    <table>
    <tr> <td>\f$i(xMult)            = 0\f$                            </td> <td>in 1D (only 1 derivative)</td> </tr>
    <tr> <td>\f$i(xMult,yMult)      =yMult\f$                         </td> <td>in 2D</td> </tr>
    <tr> <td>\f$i(xMult,yMult,zMult)=zMult+\sum_{r = 0}^{k-xMult} r\f$</td> <td>in 3D</td> </tr>
    </table>
    where the order k of Dk is implicitly defined by xMult + yMult + zMult. Space dimension is
    implicitly defined by the default values of the multiplicities of y and z derivatives.

  \param  xMult            [in]    - multiplicity of dx
    \param  yMult            [in]    - multiplicity of dy (default = -1)
  \param  zMult            [in]    - multiplicity of dz (default = -1)
  \return the ordinal of partial derivative of order k as function of dx, dy, dz multiplicities
    */
  int getDkEnumeration(const int xMult,
                       const int yMult = -1,
                       const int zMult = -1);



  /** \brief  Returns multiplicities of dx, dy, and dz based on the enumeration of the partial
    derivative, its order and the space dimension. Inverse of the getDkEnumeration() method.

    \param  partialMult      [out]    - array with the multiplicities f dx, dy and dz
    \param  derivativeEnum   [in]     - enumeration of the partial derivative
    \param  operatorType     [in]     - k-th partial derivative Dk
    \param  spaceDim         [in]     - space dimension
    */
  void getDkMultiplicities(Teuchos::Array<int>&  partialMult,
                           const int             derivativeEnum,
                           const EOperator       operatorType,
                           const int             spaceDim);



  /** \brief  Returns cardinality of Dk, i.e., the number of all derivatives of order k.

    The set of all partial derivatives of order k is isomorphic to the set of all multisets
    of cardinality k with elements taken from the sets {x}, {x,y}, and {x,y,z} in 1D, 2D,
    and 3D respectively. For example, the partial derivative
    \f$\displaystyle D\{1,2,1\}f = \frac{d^4 f}{dx dy^2 dz}\f$  maps to the multiset
    \f$\{x, y, z\}\f$ with multiplicities \f$\{1,2,1\}\f$. The number of all such multisets
    is given by the binomial coefficient
    \f[       \begin{pmatrix} spaceDim + k - 1 \\ spaceDim - 1 \end{pmatrix}              \f]
Therefore:
    \li     in 1D: cardinality = 1\n
    \li     in 2D: cardinality = k + 1\n
    \li     in 3D: cardinality = (k + 1)*(k + 2)/2

    \param  operatorType     [in]     - k-th derivative operator Dk
    \param  spaceDim         [in]     - space dimension
    \return the number of all partial derivatives of order k
    */
  int getDkCardinality(const EOperator operatorType,
                       const int       spaceDim);



/***************************************************************************************************
 ***************************************************************************************************
 **                                                                                               **
 **                      Declarations of helper functions for the basis class                     **
 **                                                                                               **
 ***************************************************************************************************
 ***************************************************************************************************/

  /** \brief  Fills <var>ordinalToTag_</var> and <var>tagToOrdinal_</var> by basis-specific tag data

    \param  tagToOrdinal     [out]  - Lookup table for the DoF's ordinal by its tag
    \param  ordinalToTag     [out]  - Lookup table for the DoF's tag by its ordinal
    \param  tags             [in]   - a set of basis-dependent tags in flat (rank-1) array format.
    \param  basisCard        [in]   - cardinality of the basis
    \param  tagSize          [in]   - number of fields in a DoF tag
    \param  posScDim         [in]   - position in the tag, counting from 0, of the subcell dim
    \param  posScOrd         [in]   - position in the tag, counting from 0, of the subcell ordinal
    \param  posDfOrd         [in]   - position in the tag, counting from 0, of DoF ordinal relative to the subcell
    */
  void setOrdinalTagData(std::vector<std::vector<std::vector<int> > >   &tagToOrdinal,
                         std::vector<std::vector<int> >                 &ordinalToTag,
                         const int                                      *tags,
                         const int                                      basisCard,
                         const int                                      tagSize,
                         const int                                      posScDim,
                         const int                                      posScOrd,
                         const int                                      posDfOrd);



/***************************************************************************************************
 ***************************************************************************************************
 **                                                                                               **
 **                      Declarations of templated utility functions                              **
 **                                                                                               **
 ***************************************************************************************************
 ***************************************************************************************************/

  enum TypeOfExactData{
    INTREPID_UTILS_FRACTION=0,
    INTREPID_UTILS_SCALAR
  };

/***************************************************************************************************
 *                                                                                                 *
 *               Utility functions for handling external data in tests                             *
 *                                                                                                 *
 ***************************************************************************************************/

/** \brief  Compares the values in the test matrix <var><b>testMat</b></var> to precomputed
            analytic values stored in a file, where the input matrix is an array of arrays.

    \param  testMat          [in]     -  test matrix
    \param  inputFile        [in]     -  input file
    \param  reltol           [in]     -  relative tolerance for equality comparisons
    \param  iprint           [in]     -  if 0, no output; if 1, details are printed
    \param  analyticDataType [in]     -  type of analytic data for comparison:
                                         \li if INTREPID_UTILS_FRACTION, analytic fractions are parsed and computed
                                         \li if INTREPID_UTILS_SCALAR, high-precision scalar data is read
    \return 0 if pass; error code otherwise
 */
template<class Scalar>
int compareToAnalytic(const Teuchos::Array< Teuchos::Array<Scalar> > testMat,
                      std::ifstream & inputFile,
                      Scalar reltol,
                      int iprint,
                      TypeOfExactData analyticDataType = INTREPID_UTILS_FRACTION);

/** \brief  Compares the values in the test matrix <var><b>testMat</b></var> to precomputed
            analytic values stored in a file, where the input matrix is a single contiguous
            array.

    \param  testMat          [in]     -  test matrix
    \param  inputFile        [in]     -  input file
    \param  reltol           [in]     -  relative tolerance for equality comparisons
    \param  iprint           [in]     -  if 0, no output; if 1, details are printed
    \param  analyticDataType [in]     -  type of analytic data for comparison:
                                         \li if INTREPID_UTILS_FRACTION, analytic fractions are parsed and computed
                                         \li if INTREPID_UTILS_SCALAR, high-precision scalar data is read
    \return 0 if pass; error code otherwise
 */
template<class Scalar>
int compareToAnalytic(const Scalar * testMat,
                      std::ifstream & inputFile,
                      Scalar reltol,
                      int iprint,
                      TypeOfExactData analyticDataType = INTREPID_UTILS_FRACTION);



/** \brief  Loads analytic values stored in a file into the matrix <var><b>testMat</b></var>,
            where the output matrix is an array of arrays.

    \param  testMat          [in]     -  test matrix
    \param  inputFile        [in]     -  input file
    \param  analyticDataType [in]     -  type of analytic data for comparison:
                                         \li if INTREPID_UTILS_FRACTION, analytic fractions are parsed and computed
                                         \li if INTREPID_UTILS_SCALAR, high-precision scalar data is read
 */
template<class Scalar>
void getAnalytic(Teuchos::Array< Teuchos::Array<Scalar> > & testMat,
                 std::ifstream & inputFile,
                 TypeOfExactData analyticDataType = INTREPID_UTILS_FRACTION);

/** \brief  Loads analytic values stored in a file into the matrix <var><b>testMat</b></var>,
            where the output matrix is a single contiguous array.

    \param  testMat          [in]     -  test matrix
    \param  inputFile        [in]     -  input file
    \param  analyticDataType [in]     -  type of analytic data for comparison:
                                         \li if INTREPID_UTILS_FRACTION, analytic fractions are parsed and computed
                                         \li if INTREPID_UTILS_SCALAR, high-precision scalar data is read
 */
template<class Scalar>
void getAnalytic(Scalar * testMat,
                 std::ifstream & inputFile,
                 TypeOfExactData analyticDataType = INTREPID_UTILS_FRACTION);



/***************************************************************************************************
 *                                                                                                 *
 *     Utility functions for checking requirements on ranks and dimensions of array arguments      *
 *                                                                                                 *
 ***************************************************************************************************/


/** \brief  Checks if the rank of the array argument is in the required range.

    \param  errmsg          [out] - error message
    \param  array           [in]  - array argument
    \param  lowerBound      [in]  - lower bound for the rank of the array
    \param  upperBound      [in]  - upper bound for the rank of the array

    \return true if lowerBound <= array.rank() <= rankR, false otherwise
  */
  template<class Array>
  bool requireRankRange(std::string&   errmsg,
                        const Array&   array,
                        const int      lowerBound,
                        const int      upperBound);



  /** \brief  Checks if two arrays have matching ranks.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  array2          [in]  - second array argument

      \return true if array.rank1() == array2.rank(), false otherwise
    */
  template<class Array1, class Array2>
  bool requireRankMatch(std::string&   errmsg,
                        const Array1&  array1,
                        const Array2&  array2);



  /** \brief  Checks if the specified array dimension is in the required range.

      \param  errmsg          [out] - error message
      \param  array           [in]  - array argument
      \param  dim             [in]  - dimension ordinal, 0 <= dim < array
      \param  lowerBound      [in]  - lower bound for dimension <var>dim</var>
      \param  upperBound      [in]  - upper bound for dimension <var>dim</var>

      \return true if lowerBound <= array.dimension(dim) <= upperBound, false otherwise
    */
  template<class Array>
  bool requireDimensionRange(std::string&  errmsg,
                             const Array&  array,
                             const int     dim,
                             const int     lowerBound,
                             const int     upperBound);



  /** \brief  Checks arrays for a single matching dimension.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  a1_dim0         [in]  - dimension ordinal for first array
      \param  array2          [in]  - second array argument
      \param  a2_dim0         [in]  - dimension ordinal for second array

      \return true if array1.dimension(a1_dim0) == array2.dimension(a2_dim0), false otherwise
    */
  template<class Array1, class Array2>
  bool requireDimensionMatch(std::string&   errmsg,
                             const Array1&  array1,
                             const int      a1_dim0,
                             const Array2&  array2,
                             const int      a2_dim0);



  /** \brief  Checks arrays for two matching dimensions.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  a1_dim0         [in]  - 1st dimension ordinal for first array
      \param  a1_dim1         [in]  - 2nd dimension ordinal for first array
      \param  array2          [in]  - second array argument
      \param  a2_dim0         [in]  - 1st dimension ordinal for second array
      \param  a2_dim1         [in]  - 2nd dimension ordinal for second array

      \return true if array1.dimension(a1_dim*) == array2.dimension(a2_dim*) for *={0,1}, false otherwise
    */
  template<class Array1, class Array2>
  bool requireDimensionMatch(std::string&   errmsg,
                             const Array1&  array1,
                             const int      a1_dim0, const int a1_dim1,
                             const Array2&  array2,
                             const int      a2_dim0, const int a2_dim1);



  /** \brief  Checks arrays for three matching dimensions.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  a1_dim0         [in]  - 1st dimension ordinal for first array
      \param  a1_dim1         [in]  - 2nd dimension ordinal for first array
      \param  a1_dim2         [in]  - 3rd dimension ordinal for first array
      \param  array2          [in]  - second array argument
      \param  a2_dim0         [in]  - 1st dimension ordinal for second array
      \param  a2_dim1         [in]  - 2nd dimension ordinal for second array
      \param  a2_dim2         [in]  - 3rd dimension ordinal for second array

      \return true if array1.dimension(a1_dim*) == array2.dimension(a2_dim*) for *={0,1,2}, false otherwise
    */
  template<class Array1, class Array2>
  bool requireDimensionMatch(std::string&   errmsg,
                             const Array1&  array1,
                             const int      a1_dim0, const int a1_dim1, const int a1_dim2,
                             const Array2&  array2,
                             const int      a2_dim0, const int a2_dim1, const int a2_dim2);



  /** \brief  Checks arrays for four matching dimensions.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  a1_dim0         [in]  - 1st dimension ordinal for first array
      \param  a1_dim1         [in]  - 2nd dimension ordinal for first array
      \param  a1_dim2         [in]  - 3rd dimension ordinal for first array
      \param  a1_dim3         [in]  - 4th dimension ordinal for first array
      \param  array2          [in]  - second array argument
      \param  a2_dim0         [in]  - 1st dimension ordinal for second array
      \param  a2_dim1         [in]  - 2nd dimension ordinal for second array
      \param  a2_dim2         [in]  - 3rd dimension ordinal for second array
      \param  a2_dim3         [in]  - 4th dimension ordinal for second array

      \return true if array1.dimension(a1_dim*) == array2.dimension(a2_dim*) for *={0,1,2,3}, false otherwise
    */
  template<class Array1, class Array2>
  bool requireDimensionMatch(std::string&   errmsg,
                             const Array1&  array1,
                             const int      a1_dim0, const int a1_dim1, const int a1_dim2, const int a1_dim3,
                             const Array2&  array2,
                             const int      a2_dim0, const int a2_dim1, const int a2_dim2, const int a2_dim3);



  /** \brief  Checks arrays for five matching dimensions.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  a1_dim0         [in]  - 1st dimension ordinal for first array
      \param  a1_dim1         [in]  - 2nd dimension ordinal for first array
      \param  a1_dim2         [in]  - 3rd dimension ordinal for first array
      \param  a1_dim3         [in]  - 4th dimension ordinal for first array
      \param  a1_dim4         [in]  - 5th dimension ordinal for first array
      \param  array2          [in]  - second array argument
      \param  a2_dim0         [in]  - 1st dimension ordinal for second array
      \param  a2_dim1         [in]  - 2nd dimension ordinal for second array
      \param  a2_dim2         [in]  - 3rd dimension ordinal for second array
      \param  a2_dim3         [in]  - 4th dimension ordinal for second array
      \param  a2_dim4         [in]  - 5th dimension ordinal for second array

      \return true if array1.dimension(a1_dim*) == array2.dimension(a2_dim*) for *={0,1,2,3,4}, false otherwise
    */
  template<class Array1, class Array2>
  bool requireDimensionMatch(std::string&   errmsg,
                             const Array1&  array1,
                             const int      a1_dim0, const int a1_dim1,
                             const int      a1_dim2, const int a1_dim3, const int a1_dim4,
                             const Array2&  array2,
                             const int      a2_dim0, const int a2_dim1,
                             const int      a2_dim2, const int a2_dim3, const int a2_dim4);



  /** \brief  Checks arrays for all their dimensions match. Arrays with equal ranks required.

      \param  errmsg          [out] - error message
      \param  array1          [in]  - first array argument
      \param  array2          [in]  - second array argument

      \return true if array1.dimension(i) == array2.dimension(i), 0 <= i < rank, false otherwise.
   */
  template<class Array1, class Array2>
  bool requireDimensionMatch(std::string&   errmsg,
                             const Array1&  array1,
                             const Array2&  array2);



/***************************************************************************************************
 ***************************************************************************************************
 **                                                                                               **
 **                           Definitions of templated functions                                  **
 **                                                                                               **
 ***************************************************************************************************
 ***************************************************************************************************/


/***************************************************************************************************
 *                                                                                                 *
 *               Utility functions for handling external data in tests                             *
 *                                                                                                 *
 ***************************************************************************************************/

template<class Scalar>
int compareToAnalytic(const Teuchos::Array< Teuchos::Array<Scalar> > testMat,
                      std::ifstream & inputFile,
                      Scalar reltol,
                      int iprint,
                      TypeOfExactData analyticDataType ) {

  // This little trick lets us print to std::cout only if
  // iprint > 0.
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  std::string line;
  Scalar testentry;
  Scalar abstol;
  Scalar absdiff;
  int i=0, j=0;
  int err = 0;

  while (! inputFile.eof() )
    {
    std::getline (inputFile,line);
    std::istringstream linestream(line);
    std::string chunk;
    j = 0;
    while( linestream >> chunk ) {
      int num1;
      int num2;
      std::string::size_type loc = chunk.find( "/", 0);
      if( loc != std::string::npos ) {
        chunk.replace( loc, 1, " ");
        std::istringstream chunkstream(chunk);
        chunkstream >> num1;
        chunkstream >> num2;
        testentry = (Scalar)(num1)/(Scalar)(num2);
        abstol = ( std::fabs(testentry) < reltol ? reltol : std::fabs(reltol*testentry) );
        absdiff = std::fabs(testentry - testMat[i][j]);
        if (absdiff > abstol) {
          err++;
          *outStream << "FAILURE --> ";
        }
        *outStream << "entry[" << i << "," << j << "]:" << "   "
          << testMat[i][j] << "   " << num1 << "/" << num2 << "   "
          << absdiff << "   " << "<?" << "   " << abstol << "\n";
      }
      else {
        std::istringstream chunkstream(chunk);
        if (analyticDataType == INTREPID_UTILS_FRACTION) {
          chunkstream >> num1;
          testentry = (Scalar)(num1);
        }
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
          chunkstream >> testentry;
        abstol = ( std::fabs(testentry) < reltol ?reltol : std::fabs(reltol*testentry) );
        absdiff = std::fabs(testentry - testMat[i][j]);
        if (absdiff > abstol) {
          err++;
          *outStream << "FAILURE --> ";
        }
        *outStream << "entry[" << i << "," << j << "]:" << "   "
          << testMat[i][j] << "   " << testentry << "   "
          << absdiff << "   " << "<?" << "   " << abstol << "\n";
      }
      j++;
    }
    i++;
    }

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return err;
} // end compareToAnalytic



template<class Scalar>
int compareToAnalytic(const Scalar * testMat,
                      std::ifstream & inputFile,
                      Scalar reltol,
                      int iprint,
                      TypeOfExactData analyticDataType ) {

  // This little trick lets us print to std::cout only if
  // iprint > 0.
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  std::string line;
  Scalar testentry;
  Scalar abstol;
  Scalar absdiff;
  int i=0, j=0, offset=0;
  int err = 0;

  while (! inputFile.eof() )
    {
    std::getline (inputFile,line);
    std::istringstream linestream(line);
    std::string chunk;
    j = 0;
    while( linestream >> chunk ) {
      int num1;
      int num2;
      std::string::size_type loc = chunk.find( "/", 0);
      if( loc != std::string::npos ) {
        chunk.replace( loc, 1, " ");
        std::istringstream chunkstream(chunk);
        chunkstream >> num1;
        chunkstream >> num2;
        testentry = (Scalar)(num1)/(Scalar)(num2);
        abstol = ( std::fabs(testentry) < reltol ? reltol : std::fabs(reltol*testentry) );
        absdiff = std::fabs(testentry - testMat[i*offset+j]);
        if (absdiff > abstol) {
          err++;
          *outStream << "FAILURE --> ";
        }
        *outStream << "entry[" << i << "," << j << "]:" << "   "
          << testMat[i*offset+j] << "   " << num1 << "/" << num2 << "   "
          << absdiff << "   " << "<?" << "   " << abstol << "\n";
      }
      else {
        std::istringstream chunkstream(chunk);
        if (analyticDataType == INTREPID_UTILS_FRACTION) {
          chunkstream >> num1;
          testentry = (Scalar)(num1);
        }
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
          chunkstream >> testentry;
        abstol = ( std::fabs(testentry) < reltol ?reltol : std::fabs(reltol*testentry) );
        absdiff = std::fabs(testentry - testMat[i*offset+j]);
        if (absdiff > abstol) {
          err++;
          *outStream << "FAILURE --> ";
        }
        *outStream << "entry[" << i << "," << j << "]:" << "   "
          << testMat[i*offset+j] << "   " << testentry << "   "
          << absdiff << "   " << "<?" << "   " << abstol << "\n";
      }
      j++;
    }
    i++;
    offset = j;
    }

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return err;
} // end compareToAnalytic



template<class Scalar>
void getAnalytic(Teuchos::Array< Teuchos::Array<Scalar> > & testMat,
                 std::ifstream & inputFile,
                 TypeOfExactData analyticDataType ) {

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  std::string line;
  Scalar testentry;
  int i=0, j=0;

  while (! inputFile.eof() )
    {
    std::getline (inputFile,line);
    std::istringstream linestream(line);
    std::string chunk;
    j = 0;
    while( linestream >> chunk ) {
      int num1;
      int num2;
      std::string::size_type loc = chunk.find( "/", 0);
      if( loc != std::string::npos ) {
        chunk.replace( loc, 1, " ");
        std::istringstream chunkstream(chunk);
        chunkstream >> num1;
        chunkstream >> num2;
        testentry = (Scalar)(num1)/(Scalar)(num2);
        testMat[i][j] = testentry;
      }
      else {
        std::istringstream chunkstream(chunk);
        if (analyticDataType == INTREPID_UTILS_FRACTION) {
          chunkstream >> num1;
          testentry = (Scalar)(num1);
        }
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
          chunkstream >> testentry;
        testMat[i][j] = testentry;
      }
      j++;
    }
    i++;
    }

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
} // end getAnalytic



template<class Scalar>
void getAnalytic(Scalar * testMat,
                 std::ifstream & inputFile,
                 TypeOfExactData analyticDataType) {

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  std::string line;
  Scalar testentry;
  int i=0, j=0, offset=0;

  while (! inputFile.eof() )
    {
    std::getline (inputFile,line);
    std::istringstream linestream(line);
    std::string chunk;
    j = 0;
    while( linestream >> chunk ) {
      int num1;
      int num2;
      std::string::size_type loc = chunk.find( "/", 0);
      if( loc != std::string::npos ) {
        chunk.replace( loc, 1, " ");
        std::istringstream chunkstream(chunk);
        chunkstream >> num1;
        chunkstream >> num2;
        testentry = (Scalar)(num1)/(Scalar)(num2);
        testMat[i*offset+j] = testentry;
      }
      else {
        std::istringstream chunkstream(chunk);
        if (analyticDataType == INTREPID_UTILS_FRACTION) {
          chunkstream >> num1;
          testentry = (Scalar)(num1);
        }
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
          chunkstream >> testentry;
        testMat[i*offset+j] = testentry;
      }
      j++;
    }
    i++;
    offset = j;
    }

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
} // end getAnalytic


/***************************************************************************************************
 *                                                                                                 *
 *     Utility functions for checking requirements on ranks and dimensions of array arguments      *
 *                                                                                                 *
 ***************************************************************************************************/


template<class Array>
bool requireRankRange(std::string&   errmsg,
                      const Array&   array,
                      const int      lowerBound,
                      const int      upperBound){

  TEUCHOS_TEST_FOR_EXCEPTION( (lowerBound > upperBound) , std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireRankRange): lowerBound <= upperBound required!");

  bool OK = true;
  if( (lowerBound == upperBound) && !(getrank(array) == (size_t)lowerBound) ) {
    errmsg += "\n>>> Array rank = ";
    errmsg += (char)(48 + getrank(array) );
    errmsg += " while rank-";
    errmsg += (char) (48 + lowerBound);
    errmsg += " array required.";
    OK = false;
  }
  else if ( (lowerBound < upperBound) &&  !( ((size_t)lowerBound <= getrank(array) ) && (getrank(array) <= (size_t)upperBound)  ) ){
    errmsg += "\n>>> Array rank = ";
    errmsg += (char)(48 + getrank(array) );
    errmsg += " while a rank between ";
    errmsg += (char) (48 + lowerBound);
    errmsg += " and ";
    errmsg += (char) (48 + upperBound);
    errmsg += " is required.";
    OK = false;
  }
  return OK;
}


template<class Array1, class Array2>
bool requireRankMatch(std::string&   errmsg,
                      const Array1&  array1,
                      const Array2&  array2){
  bool OK = true;
  if(getrank(array1) != getrank(array2) ) {
    errmsg += "\n>>> Array ranks are required to match.";
    OK = false;
  }
  return OK;
}


template<class Array>
bool requireDimensionRange(std::string&  errmsg,
                           const Array&  array,
                           const int     dim,
                           const int     lowerBound,
                           const int     upperBound){

  TEUCHOS_TEST_FOR_EXCEPTION( (lowerBound > upperBound) , std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionRange): lowerBound <= upperBound required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= dim) && ((size_t)dim < getrank(array) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionRange): 0 <= dim < array.rank() required!");

  bool OK = true;
  if( (lowerBound > upperBound) || ( (size_t)dim >= getrank(array) ) ) {
    errmsg += "\n>>> Unexpected error: ";
    OK = false;
  }
  if( (lowerBound == upperBound) && !(static_cast<int>(array.dimension(dim)) == lowerBound) ) {
    errmsg += "\n>>> dimension(";
    errmsg += (char)(48 + dim);
    errmsg += ") = ";
    errmsg += (char)(48 + array.dimension(dim) );
    errmsg += " while dimension(";
    errmsg += (char)(48 + dim);
    errmsg += ") = ";
    errmsg += (char)(48 + lowerBound);
    errmsg += " required.";
    OK = false;
  }
  else if( (lowerBound < upperBound) &&
           !( ((size_t)lowerBound <= (size_t)array.dimension(dim) ) && (static_cast<size_t>(array.dimension(dim)) <= (size_t)upperBound) ) ){
    errmsg += "\n>>> dimension(";
    errmsg += (char)(48 + dim);
    errmsg += ") = ";
    errmsg += (char)(48 + array.dimension(dim) );
    errmsg += " while ";
    errmsg += (char)(48 + lowerBound);
    errmsg += " <= dimension(";
    errmsg += (char)(48 + dim);
    errmsg += ") <= ";
    errmsg +=(char)(48 + upperBound);
    errmsg +=" required.";
    OK = false;
  }
  return OK;
}



template<class Array1, class Array2>
bool requireDimensionMatch(std::string&   errmsg,
                           const Array1&  array1,
                           const int      a1_dim0,
                           const Array2&  array2,
                           const int      a2_dim0){

  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim0) && ((size_t)a1_dim0 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim0 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim0) && ((size_t)a2_dim0 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim0 < array2.rank() required!");

  bool OK = true;
  if(static_cast<int>(array1.dimension(a1_dim0)) != static_cast<int>(array2.dimension(a2_dim0)) ){
    errmsg += "\n>>> dimension(";
    errmsg += (char)(48 + a1_dim0);
    errmsg += ") of 1st array and dimension(";
    errmsg += (char)(48 + a2_dim0);
    errmsg += ") of 2nd array are required to match.";
    OK = false;
  }
  return OK;
}



template<class Array1, class Array2>
bool requireDimensionMatch(std::string&   errmsg,
                           const Array1&  array1,
                           const int      a1_dim0, const int a1_dim1,
                           const Array2&  array2,
                           const int      a2_dim0, const int a2_dim1){

  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim0) && ((size_t)a1_dim0 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim0 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim1) && ((size_t)a1_dim1 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim1 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim0) && ((size_t)a2_dim0 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim0 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim1) && ((size_t)a2_dim1 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim1 < array2.rank() required!");

  bool OK = true;
  if( !requireDimensionMatch(errmsg, array1, a1_dim0, array2, a2_dim0) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim1, array2, a2_dim1) ){
    OK = false;
  }
  return OK;
}



template<class Array1, class Array2>
bool requireDimensionMatch(std::string&   errmsg,
                           const Array1&  array1,
                           const int      a1_dim0, const int a1_dim1, const int a1_dim2,
                           const Array2&  array2,
                           const int      a2_dim0, const int a2_dim1, const int a2_dim2){

  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim0) && ((size_t)a1_dim0 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim0 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim1) && ((size_t)a1_dim1 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim1 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim2) && ((size_t)a1_dim2 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim2 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim0) && ((size_t)a2_dim0 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim0 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim1) && ((size_t)a2_dim1 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim1 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim2) && ((size_t)a2_dim2 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim2 < array2.rank() required!");


  bool OK = true;
  if( !requireDimensionMatch(errmsg, array1, a1_dim0, array2, a2_dim0) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim1, array2, a2_dim1) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim2, array2, a2_dim2) ){
    OK = false;
  }
  return OK;
}



template<class Array1, class Array2>
bool requireDimensionMatch(std::string&   errmsg,
                           const Array1&  array1,
                           const int      a1_dim0, const int a1_dim1, const int a1_dim2, const int a1_dim3,
                           const Array2&  array2,
                           const int      a2_dim0, const int a2_dim1, const int a2_dim2, const int a2_dim3){

  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim0) && ((size_t)a1_dim0 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim0 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim1) && ((size_t)a1_dim1 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim1 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim2) && ((size_t)a1_dim2 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim2 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim3) && ((size_t)a1_dim3 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim3 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim0) && ((size_t)a2_dim0 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim0 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim1) && ((size_t)a2_dim1 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim1 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim2) && ((size_t)a2_dim2 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim2 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim3) && ((size_t)a2_dim3 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim3 < array2.rank() required!");
  bool OK = true;
  if( !requireDimensionMatch(errmsg, array1, static_cast<int>(a1_dim0), array2, static_cast<int>(a2_dim0)) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim1, array2, a2_dim1) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim2, array2, a2_dim2) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim3, array2, a2_dim3) ){
    OK = false;
  }
  return OK;
}



template<class Array1, class Array2>
bool requireDimensionMatch(std::string&   errmsg,
                           const Array1&  array1,
                           const int      a1_dim0, const int a1_dim1, const int a1_dim2,
                           const int      a1_dim3, const int a1_dim4,
                           const Array2&  array2,
                           const int      a2_dim0, const int a2_dim1, const int a2_dim2,
                           const int      a2_dim3, const int a2_dim4){

  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim0) && ((size_t)a1_dim0 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim0 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim1) && ((size_t)a1_dim1 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim1 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim2) && ((size_t)a1_dim2 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim2 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim3) && ((size_t)a1_dim3 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim3 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a1_dim4) && ((size_t)a1_dim4 < getrank(array1) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a1_dim4 < array1.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim0) && ((size_t)a2_dim0 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim0 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim1) && ((size_t)a2_dim1 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim1 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim2) && ((size_t)a2_dim2 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim2 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim3) && ((size_t)a2_dim3 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim3 < array2.rank() required!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (0 <= a2_dim4) && ((size_t)a2_dim4 < getrank(array2) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): 0 <= a2_dim4 < array2.rank() required!");

  bool OK = true;
  if( !requireDimensionMatch(errmsg, array1, a1_dim0, array2, a2_dim0) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim1, array2, a2_dim1) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim2, array2, a2_dim2) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim3, array2, a2_dim3) ){
    OK = false;
  }
  if( !requireDimensionMatch(errmsg, array1, a1_dim4, array2, a2_dim4) ){
    OK = false;
  }
  return OK;
}



template<class Array1, class Array2>
bool requireDimensionMatch(std::string&   errmsg,
                           const Array1&  array1,
                           const Array2&  array2){

  TEUCHOS_TEST_FOR_EXCEPTION( !requireRankMatch(errmsg, array1, array2 ), std::invalid_argument,
                      ">>> ERROR (Intrepid_Utils::requireDimensionMatch): Arrays with equal ranks are required to test for all dimensions match." )

  bool OK = true;
  for(size_t dim = 0; dim < getrank(array1); dim++){
    if( !requireDimensionMatch(errmsg, array1, dim, array2, dim) ){
      OK = false;
      break;
    }
  }
  return OK;
}


} // end namespace Intrepid

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

