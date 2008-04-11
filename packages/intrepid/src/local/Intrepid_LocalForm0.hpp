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

/** \file   Intrepid_LocalForm0.hpp
    \brief  Header file for the Intrepid::LocalForm0 class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_LOCAL_FORM_0_HPP
#define INTREPID_LOCAL_FORM_0_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_LocalField.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_Cubature.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"


namespace Intrepid {

/** \class Intrepid::LocalForm0
    \brief Defines the base class for the representation of differential 0 forms in Intrepid.

    The local field interface relies on two member functions, <var>getOperator</var>
    and <var>getFunctional</var>, with varying signatures. Detailed info here ...
*/
template<class Scalar>
class LocalForm0 : public LocalField<Scalar> {
  private:

  const static EField                                                      fieldType_ = FIELD_FORM_0;

  Teuchos::RCP<Basis<Scalar> >                                             basis_;
  ECell                                                                    basisCell_;
  ECoordinates                                                             basisCoordSystem_;
  EBasis                                                                   basisType_;
  int                                                                      basisDegree_;
  Teuchos::Array<LocalDofTag>                                              basisDofTags_;
  int                                                                      basisNumDofs_;
  Teuchos::Array<Teuchos::Array<Teuchos::Array<VarContainer<Scalar> > > >  basisVals_;

  Teuchos::Array<Teuchos::Array<Teuchos::RCP<Cubature<Scalar> > > >        cubature_;
  Teuchos::Array<Teuchos::Array<Teuchos::Array<Point<Scalar> > > >         cubPoints_; 
  Teuchos::Array<Teuchos::Array<Teuchos::Array<Scalar> > >                 cubWeights_;
  Teuchos::Array<Teuchos::Array<int> >                                     numCubPoints_;

  ECompEngine                                                              compEngine_;

  public:

  //LocalForm0() {}

  /** brief  Constructor.

      \param basis         [in]     - Refcount pointer to an Intrepid::Basis object.
      \param cubature      [in]     - Refcount pointer to an Intrepid::Cubature object.
      \param compEngine    [in]     - Computational engine (manual, BLAS, etc.).
  */
  LocalForm0(const Teuchos::RCP<Basis<Scalar> > basis,
             const Teuchos::Array<Teuchos::Array<Teuchos::RCP<Cubature<Scalar> > > > cubature, 
             const ECompEngine compEngine = COMP_CPP);


  ~LocalForm0() {}


  /** \brief Returns values of a (primitive) operator <var>primOp</var> applied to
             FEM basis functions, evaluated at the array of (preset) cubature points.

      \param primOp           [in]     - Input operator (primitive).
      \param subDim           [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subCellId        [in]     - Subcell id.

      \return
             - Output array.
  */
  const VarContainer<Scalar> & getOperator(const EOperator  primOp,
                                           const int        subDim,
                                           const int        subCellId);


  /** \brief Returns values of a (primitive) operator <var>primOp</var> applied to
             FEM basis functions, evaluated at the array <var>inputPoints</var> of
             <strong>reference</strong> points.

      The context in which this interface function can be called is limited to FEM
      reconstructions (i.e. via a reference cell). The function returns
      operator values <strong>decoupled from the geometry of a particular cell</strong>
      (untransformed values). 

      \param outputValues    [out]     - Output array.
      \param inputPoints      [in]     - Array of input (reference) points.
      \param primOp           [in]     - Input operator (primitive).
  */
  void getOperator(VarContainer<Scalar> &                  outputValues,
                   const Teuchos::Array<Point<Scalar> > &  inputPoints,
                   const EOperator                         primOp);


  /** \brief Returns values of a (primitive) operator <var>primOp</var> applied to
             FEM or FVD basis functions, evaluated at the array <var>inputPoints</var>
             of <strong>physical</strong> points.

      The context in which this interface function can be called is two-fold.
      For FEM reconstructions (i.e. those based on a reference cell), the function returns
      the values of the operator applied to basis functions in <strong>physical</strong> space,
      i.e. relevant geometric transformations will be incorporated into the evaluation.
      For FV/D reconstructions, the function computes operator values <strong>directly</strong>
      on physical cells.

      \param outputValues    [out]     - Output array.
      \param inputPoints      [in]     - Array of input (physical) points.
      \param primOp           [in]     - Input operator (primitive).
      \param cell             [in]     - Physical cell.
  */
  void getOperator(LexContainer<Scalar> &                  outputValues,
                   const Teuchos::Array<Point<Scalar> > &  inputPoints,
                   const EOperator                         primOp,
                   const Cell<Scalar> &                    cell);


  /** \brief Returns discrete representation (matrix) of integral quantities (one for every
             cell in the multicell <var>mCell</var>) involving
             a left differential operator <var>leftOp</var> applied to basis
             functions, external input data <var>inputData</var>, and a right
             differential operator <var>rightOp</var> applied to basis functions,
             where left and right basis functions belong to the same basis.

      The integral to be computed is

      \f$ \displaystyle \int_{\Omega} (leftOp\;\phi) (inputData) (rightOp\;\widehat{\phi}) d\Omega\f$

      where \f$\Omega\f$ is any valid integration domain (line, surface, cell),
      \f$\phi\f$ and \f$\widehat{\phi}\f$ are basis functions from the same basis,
      and <var>inputData</var> is a data array that should be formatted as below.\n
      Note:\n
      \f$s\f$ denotes a scalar quantity,\n
      \f$v_i\f$ denotes the i-th component of a vector \f$v\f$,\n
      \f$T_{ij}\f$ denotes the ij-th component of a tensor \f$T\f$,\n
      \f$q_k\f$ denotes the k-th integration point,\n
      \f$N\f$ denotes the number of integration points,\n
      \f$d\f$ denotes the space dimension

      <table>
        <tr> <td>\f$ value(s, cell\_id, q_k) = inputData[cell\_id \cdot N + k] \f$</td>
             <td>if <var>inputFormat==DATA_SCALAR</var></td> </tr>
        <tr> <td>\f$ value(v_i, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d + k \cdot d + i] \f$</td>
             <td>if <var>inputFormat==DATA_VECTOR</var></td> </tr>
        <tr> <td>\f$ value(T_{ij}, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d^2 + k \cdot d^2 + i \cdot d + j] \f$</td>
             <td>if <var>inputFormat==DATA_TENSOR</var></td> </tr>
      </table>

      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param mCell            [in]     - Multicell.
      \param inputData        [in]     - Input data.
      \param inputFormat      [in]     - Format of input data (scalar, vector, etc.).
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  void getOperator(LexContainer<Scalar> &          outputValues,
                   const EOperator                 leftOp,
                   const EOperator                 rightOp,
                   MultiCell<Scalar> &             mCell,
                   const Teuchos::Array<Scalar> &  inputData,
                   const EDataFormat               inputFormat,
                   const EIntegrationDomain        intDomain = INTEGRATION_DOMAIN_CELL);


  /** \brief Returns discrete representation (matrix) of integral quantities (one for every
             cell in the multicell <var>mCell</var>) involving
             a left differential operator <var>leftOp</var> and a right
             differential operator <var>rightOp</var> applied to basis functions,
             where left and right basis functions belong to the same basis.

      The integral to be computed is

      \f$ \displaystyle \int_{\Omega} (leftOp\;\phi) (rightOp\;\widehat{\phi}) d\Omega\f$

      where \f$\Omega\f$ is any valid integration domain (line, surface, cell) and
      \f$\phi\f$ and \f$\widehat{\phi}\f$ are basis functions from the same basis.

      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param mCell            [in]     - Multicell.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  void getOperator(LexContainer<Scalar> &      outputValues,
                   const EOperator             leftOp,
                   const EOperator             rightOp,
                   MultiCell <Scalar> &        mCell,
                   const EIntegrationDomain    intDomain = INTEGRATION_DOMAIN_CELL);


  /** \brief Returns discrete representation (matrix) of integral quantities (one for every
             cell in the multicell <var>mCell</var>) involving
             a left differential operator <var>leftOp</var> applied to basis
             functions from the local field \f$\Lambda^m\f$, external input data <var>inputData</var>,
             and a right differential operator <var>rightOp</var> applied to basis functions
             from the local field \f$\Lambda^n\f$.

      The integral to be computed is

      \f$ \displaystyle \int_{\Omega} (leftOp\;\phi) (inputData) (rightOp\; \mbox{\boldmath$\varphi$}) d\Omega\f$

      where \f$\Omega\f$ is any valid integration domain (line, surface, cell),
      \f$\phi\f$ and \f$\mbox{\boldmath$\varphi$}\f$ are basis functions that span
      potentially different local fields,
      and <var>inputData</var> is a data array that should be formatted as below.
      The domain of <var>leftOp</var> is the parent local field. The domain, i.e. the local field for
      <var>rightOp</var> must be specified via <var>rightOpField</var>.\n
      Note:\n
      \f$s\f$ denotes a scalar quantity,\n
      \f$v_i\f$ denotes the i-th component of a vector \f$v\f$,\n
      \f$T_{ij}\f$ denotes the ij-th component of a tensor \f$T\f$,\n
      \f$q_k\f$ denotes the k-th integration point,\n
      \f$N\f$ denotes the number of integration points,\n
      \f$d\f$ denotes the space dimension

      <table>
        <tr> <td>\f$ value(s, cell\_id, q_k) = inputData[cell\_id \cdot N + k] \f$</td>
             <td>if <var>inputFormat==DATA_SCALAR</var></td> </tr>
        <tr> <td>\f$ value(v_i, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d + k \cdot d + i] \f$</td>
             <td>if <var>inputFormat==DATA_VECTOR</var></td> </tr>
        <tr> <td>\f$ value(T_{ij}, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d^2 + k \cdot d^2 + i \cdot d + j] \f$</td>
             <td>if <var>inputFormat==DATA_TENSOR</var></td> </tr>
      </table>

      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param rightOpField     [in]     - Local field of the right operator.
      \param mCell            [in]     - Multicell.
      \param inputData        [in]     - Input data.
      \param inputFormat      [in]     - Format of input data (scalar, vector, etc.).
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  void getOperator(LexContainer<Scalar> &           outputValues,
                   const EOperator                  leftOp,
                   const EOperator                  rightOp,
                   const LocalField<Scalar> &       rightOpField,
                   MultiCell<Scalar> &              mCell,
                   const Teuchos::Array<Scalar> &   inputData,
                   const EDataFormat                inputFormat,
                   const EIntegrationDomain         intDomain = INTEGRATION_DOMAIN_CELL);


  /** \brief Fills the Intrepid::LexContainer <var>leftValues</var> with a discrete
             representation of <var>leftOp</var> that is convenient for fast integration.

      The discrete representation depends on <var>leftOp</var>, <var>rightOp</var>, 
      <var>rightOpField</var>, and various other "convenience" factors.

      \param leftValues      [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param rightOpField     [in]     - Local field of the right operator.
      \param mCell            [in]     - Multicell.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  void fillLeft(LexContainer<Scalar> &      leftValues,
                const EOperator             leftOp,
                const EOperator             rightOp,
                const LocalField<Scalar> &  rightOpField,
                MultiCell<Scalar> &         mCell,
                const EIntegrationDomain    intDomain);


  /** \brief Fills the Intrepid::LexContainer <var>rightValues</var> with a discrete
             representation of <var>rightOp</var> that is convenient for fast integration.

      The discrete representation depends on <var>leftOp</var>, <var>rightOp</var>, 
      <var>rightOpField</var>, and various other "convenience" factors.

      \param rightValues     [out]     - Output array.
      \param rightOp          [in]     - Right operator.
      \param leftOp           [in]     - Left operator.
      \param rightOpField     [in]     - Local field of the right operator.
      \param mCell            [in]     - Multicell.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  void fillRight(LexContainer<Scalar> &      leftValues,
                 const EOperator             leftOp,
                 const EOperator             rightOp,
                 const LocalField<Scalar> &  rightOpField,
                 MultiCell<Scalar> &         mCell,
                 const EIntegrationDomain    intDomain);


  /** \brief Performs numerical integration (note that geometric transformations
             and integration weights have already been included in <var>leftValues</var>
             and <var>rightValues</var>).

      \param outputValues     [out]     - Output array.
      \param leftValues       [out]     - Left array.
      \param rightValues      [out]     - Right array.
  */
  void integrate(LexContainer<Scalar> &        outputValues,
                 const LexContainer<Scalar> &  leftValues,
                 const LexContainer<Scalar> &  rightValues) const;


  /** \brief Returns field type.
  */
  EField getFieldType() const {return fieldType_;}

}; // class LocalForm0

}// end namespace Intrepid

#include "Intrepid_LocalForm0Def.hpp"

#endif
