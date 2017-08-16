#ifndef INTREPID2_ORIENTATIONTOOLS_HPP
#define INTREPID2_ORIENTATIONTOOLS_HPP

// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_OrientationTools.hpp
    \brief  Header file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_FieldContainer_Kokkos.hpp"
#include "Intrepid2_RealSpaceTools.hpp"

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_BasisSet.hpp"

#include <Intrepid2_KokkosRank.hpp>
#include "Kokkos_Core.hpp"

namespace Intrepid2 {

  // Orientation encoding and decoding
  class Orientation {
  private:
    unsigned int _edgeOrt, _faceOrt;

  public:
    template<typename NodeType>
    static void getElementNodeMap(NodeType *subCellVerts,
                                  int & numVertex,
                                  const shards::CellTopology & cellTopo,
                                  const NodeType *elemNodes,
                                  const int subCellDim,
                                  const int subCellOrd);
    
    // orientation is always computed on the right-handed coordinates
    template<typename NodeType>
    static int getOrientation(const NodeType *subCellVerts,
                              const int numVertex);

    template<typename NodeType>
    static Orientation getOrientation(const shards::CellTopology & cellTopo,
                                      const NodeType *elemNodes);
    
    Orientation();
    Orientation(const Orientation &b) = default;

    bool isAlignedToReference() const;

    void setEdgeOrientation(const int numEdge, const int edgeOrt[]);
    void getEdgeOrientation(int *edgeOrt, const int numEdge) const;
    void setFaceOrientation(const int numFace, const int faceOrt[]);
    void getFaceOrientation(int *faceOrt, const int numFace) const;
  };

  template<class Scalar>
  class OrientationTools {
  private:

    class DenseMatrix {
    private:
      int _offm, _offn, _m, _n, _cs, _rs;
      Kokkos::View<Scalar*> _a;

    public:
      DenseMatrix();
      DenseMatrix(const DenseMatrix &b) = default;

      DenseMatrix(const int m, const int n);
      void setView(const DenseMatrix &b,
                   const int offm, const int m,
                   const int offn, const int n);

      int NumRows() const;
      int NumCols() const;
      int RowStride() const;
      int ColStride() const;
      Scalar* ValuePtr() const;

      Scalar& Value(const int i,
                    const int j);

      Scalar Value(const int i,
                   const int j) const;

      size_t countNumNonZeros(const Scalar epsilon) const;

      std::ostream& showMe(std::ostream &os) const;
    };

    class CoeffMatrix {
    private:
      int _m, _n;

      Kokkos::View<size_t*>  _ap;      //!< pointers to column index and values
      Kokkos::View<int*>     _aj;      //!< column index compressed format
      Kokkos::View<Scalar*>  _ax;      //!< values

      void createInternalArrays(const int m,
                                const int n,
                                const size_t nnz);

    public:
      CoeffMatrix();
      CoeffMatrix(const CoeffMatrix &b) = default;

      void import(const DenseMatrix &b,
                  const bool transpose);

      int NumRows() const;
      int NumCols() const;

      size_t RowPtr(const int i) const;
      int* ColsInRow(const int i) const;
      Scalar* ValuesInRow(const int i) const;
      int NumNonZerosInRow(const int i) const;

      std::ostream& showMe(std::ostream &os) const;
    };

    /** \brief  Computes modified point for line segment.

        \param  ot       [out] - modified point value
        \param  pt       [in]  - input point in [-1.0 , 1.0]
        \param  ort      [in]  - orientation number between 0 and 1
    */
    static void getModifiedLinePoint(double &ot,
                                     const double pt,
                                     const int ort);

    /** \brief  Computes modified point for triangle.

        \param  ot0      [out] - modified coordinate 0
        \param  ot1      [out] - modified coordinate 1
        \param  pt0      [out] - input coordinate 0
        \param  pt1      [out] - input coordinate 1
        \param  ort      [in]  - orientation number between 0 and 5
    */
    static void getModifiedTrianglePoint(double &ot0,
                                         double &ot1,
                                         const double pt0,
                                         const double pt1,
                                         const int ort);

    /** \brief  Computes modified point for quadrilateral.

        \param  ot0      [out] - modified coordinate 0
        \param  ot1      [out] - modified coordinate 1
        \param  pt0      [out] - input coordinate 0
        \param  pt1      [out] - input coordinate 1
        \param  ort      [in]  - orientation number between 0 and 7
    */
    static void getModifiedQuadrilateralPoint(double &ot0,
                                              double &ot1,
                                              const double pt0,
                                              const double pt1,
                                              const int ort);

    template<class ArrayType>
    static void getEdgeCoeffMatrix_HGRAD(CoeffMatrix &                   C,
                                         const Basis<Scalar,ArrayType> & lineBasis,
                                         const Basis<Scalar,ArrayType> & cellBasis,
                                         const int                       edgeId,
                                         const int                       edgeOrt);
    
    template<class ArrayType>
    static void getEdgeCoeffMatrix_HCURL(CoeffMatrix &                   C,
                                         const Basis<Scalar,ArrayType> & cellBasis,
                                         const Basis<Scalar,ArrayType> & lineBasis,
                                         const int                       edgeId,
                                         const int                       edgeOrt) {}

    template<class ArrayType>
    static void getEdgeCoeffMatrix_HDIV(CoeffMatrix &                   C,
                                        const Basis<Scalar,ArrayType> & lineBasis,
                                        const Basis<Scalar,ArrayType> & cellBasis,
                                        const int                       edgeId,
                                        const int                       edgeOrt) {}

    template<class ArrayType>
    static void getTriangleCoeffMatrix_HGRAD(CoeffMatrix &                   C,
                                             const Basis<Scalar,ArrayType> & faceBasis,
                                             const Basis<Scalar,ArrayType> & cellBasis,
                                             const int                       faceId,
                                             const int                       faceOrt);

    template<class ArrayType>
    static void getTriangleCoeffMatrix_HCURL(CoeffMatrix &                   C,
                                             const Basis<Scalar,ArrayType> & faceBasis,
                                             const Basis<Scalar,ArrayType> & cellBasis,
                                             const int                       faceId,
                                             const int                       faceOrt) {}

    template<class ArrayType>
    static void getTriangleCoeffMatrix_HDIV(CoeffMatrix &                   C,
                                            const Basis<Scalar,ArrayType> & faceBasis,
                                            const Basis<Scalar,ArrayType> & cellBasis,
                                            const int                       faceId,
                                            const int                       faceOrt) {}

    template<class ArrayType>
    static void getQuadrilateralCoeffMatrix_HGRAD(CoeffMatrix &                   C,
                                                  const Basis<Scalar,ArrayType> & faceBasis,
                                                  const Basis<Scalar,ArrayType> & cellBasis,
                                                  const int                       faceId,
                                                  const int                       faceOrt) {}

    template<class ArrayType>
    static void getQuadrilateralCoeffMatrix_HCURL(CoeffMatrix &                   C,
                                                  const Basis<Scalar,ArrayType> & faceBasis,
                                                  const Basis<Scalar,ArrayType> & cellBasis,
                                                  const int                       faceId,
                                                  const int                       faceOrt) {}

    template<class ArrayType>
    static void getQuadrilateralCoeffMatrix_HDIV(CoeffMatrix &                   C,
                                                 const Basis<Scalar,ArrayType> & faceBasis,
                                                 const Basis<Scalar,ArrayType> & cellBasis,
                                                 const int                       faceId,
                                                 const int                       faceOrt) {}

    // all edges and faces should apply this coeff matrix
    template<class ArrayType>
    static void applyCoeffMatrix(ArrayType &         outValues,
                                 const ArrayType &   refValues,
                                 const CoeffMatrix & C,
                                 const unsigned int  offset,
                                 const unsigned int  numDofs);

    // vertex and interior DOFs only
    template<class ArrayType>
    static void copyBasisValues(ArrayType &         outValues,
                                const ArrayType &   refValues,
                                const unsigned int  offRows, const unsigned int  numRows,
                                const unsigned int  offCols, const unsigned int  numCols);

  public:
    
    template<class ArrayType>
    static bool isLeftHandedCell(const ArrayType & pts);
    
    /** \brief  Computes modified parameterization maps of 1- and 2-subcells with orientation.

        \param  outPoints       [out] - rank-2 (P,D2) array with points in 1D or 2D modified domain with orientation
        \param  refPoints       [in]  - rank-2 (P,D2) array with points in 1D or 2D parameter domain
        \param  cellTopo        [in]  - cell topology of the parameterized domain (1- and 2-subcells)
        \param  cellOrt         [in]  - cell orientation number (zero is aligned with shards default configuration
    */
    template<class ArrayType>
    static void mapToModifiedReference(ArrayType &                  outPoints,
                                       const ArrayType &            refPoints,
                                       const shards::CellTopology & cellTopo,
                                       const int                    cellOrt = 0);


    // basis object indlues std::vector which needs copied by reerence.
    // this should be removed.
    template<class ArrayType>
    static void getBasisFunctionsByTopology(ArrayType &                     outValues,
                                            const ArrayType &               refValues,
                                            const Basis<Scalar,ArrayType> & basis);


    template<class ArrayType>
    static void getModifiedBasisFunctions(ArrayType &                        outValues,
                                          const ArrayType &                  refValues,
                                          const BasisSet<Scalar,ArrayType> & basis,
                                          const Orientation                  ort);


    static bool verbose;
    static bool reverse;
    static std::ostream* verboseStreamPtr;
  };

}

// include templated function definitions
#include "Intrepid2_OrientationToolsDef.hpp"

#endif

#endif
