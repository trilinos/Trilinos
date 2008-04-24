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

/** \file   Intrepid_LocalField.hpp
    \brief  Header file for the Intrepid::LocalField class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_LOCAL_FIELD_HPP
#define INTREPID_LOCAL_FIELD_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Intrepid_Cell.hpp"
#include "Intrepid_LexContainer.hpp"
#include "Intrepid_VarContainer.hpp"


namespace Intrepid {

/** \class Intrepid::LocalField
    \brief Defines the base class for the representation of local fields in Intrepid. 

    The local field interface relies on two member functions, <var>getOperator</var>
    and <var>getFunctional</var>, with varying signatures. Detailed info here ...
*/
template<class Scalar>
class LocalField {
  private:
    
  public:

  virtual ~LocalField() {}


  /** \brief Returns VarContainer with a multi-indexed quantity representing the values of a 
  (primitive) operator <var>primOp</var> applied to FEM basis functions, evaluated at the array 
  <var>inputPoints</var> of <strong>reference</strong> points. The rank of the output VarContainer
  and its index range vary depending on the type of <var>primOp</var> and the concrete field; see
  implementation for details.

  The context in which this interface function can be called is limited to FEM reconstructions 
  (i.e. via a reference cell). The function returns operator values <strong>decoupled from the 
  geometry of a particular cell</strong> (untransformed values). 

      \param outputValues    [out]     - Output container.
      \param inputPoints      [in]     - Array of input (reference) points.
      \param primOp           [in]     - Input operator (primitive).
  */
  virtual void getOperator(VarContainer<Scalar> &                  outputValues,
                           const Teuchos::Array<Point<Scalar> > &  inputPoints,
                           const EOperator                         primOp) = 0;


  /** \brief Returns VarContainer with a multi-indexed quantity representing the values of a 
  (primitive) operator <var>primOp</var> applied to FEM or FVD basis functions, evaluated at the 
  array <var>inputPoints</var> of <strong>physical</strong> points. The rank of the output VarContainer
  and its index range vary depending on the type of <var>primOp</var> and the concrete field; see
  implementation for details.  

  The context in which this interface function can be called is two-fold. For FEM reconstructions 
  (i.e. those based on a reference cell), the function returns the values of the operator applied to 
  basis functions in <strong>physical</strong> space, i.e. relevant geometric transformations will 
  be incorporated into the evaluation. For FV/D reconstructions, the function computes operator 
  values <strong>directly</strong> on physical cells.

      \param outputValues    [out]     - Output container.
      \param inputPoints      [in]     - Array of input (reference) points.
      \param primOp           [in]     - Input operator (primitive).
      \param cell             [in]     - Physical cell.
  */
  virtual void getOperator(VarContainer<Scalar> &                  outputValues,
                           const Teuchos::Array<Point<Scalar> > &  inputPoints,
                           const EOperator                         primOp,
                           const Cell<Scalar> &                    cell) = 0;


  /** \brief Returns a rank-3 LexContainer such that <var>outputValues{c,l,r}</var> = \f$\mathbf{A}^c_{lr}\f$
    where:
    
    \li \f$ \displaystyle \mathbf{A}^c_{lr} = \int_{\Omega_c} (leftOp\;\phi_l) (inputData) (rightOp\;{\phi_r}) d\Omega \f$    
    \li  \f$\phi_l\f$ and \f${\phi_r}\f$ are basis functions from the <strong>native</strong> local field \f$\mathcal{F}\f$;
    \li  <var>leftOp</var> and <var>rightOp</var> are admissible operator types for \f$\mathcal{F}\f$;
    \li  <var>inputData</var> is a user-supplied data array;
    \li  \f$\Omega_c\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The multi-index <var>{c,l,r}</var> has the following upper bounds:
    \li \f$ 0 \le c < \f$ number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le l,r < \dim \mathcal{F} \f$, i.e., the number of basis functions in the basis.
    
    The user-supplied data array should be formatted as follows:\n

      <table>
        <tr> <td>\f$ value(s, cell\_id, q_k) = inputData[cell\_id \cdot N + k] \f$</td>
             <td>if <var>inputFormat==DATA_SCALAR</var></td> </tr>
        <tr> <td>\f$ value(v_i, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d + k \cdot d + i] \f$</td>
             <td>if <var>inputFormat==DATA_VECTOR</var></td> </tr>
        <tr> <td>\f$ value(T_{ij}, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d^2 + k \cdot d^2 + i \cdot d + j] \f$</td>
             <td>if <var>inputFormat==DATA_TENSOR</var></td> </tr>
      </table>
    
    where \n
    \li   \f$s\f$ denotes a scalar quantity,\n
    \li   \f$v_i\f$ denotes the i-th component of a vector \f$v\f$,\n
    \li   \f$T_{ij}\f$ denotes the ij-th component of a tensor \f$T\f$,\n
    \li   \f$q_k\f$ denotes the k-th integration point,\n
    \li   \f$N\f$ denotes the number of integration points,\n
    \li   \f$d\f$ denotes the space dimension
    
    
      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param mCell            [in]     - Multicell.
      \param inputData        [in]     - Input data.
      \param inputFormat      [in]     - Format of input data (scalar, vector, etc.).
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(LexContainer<Scalar> &          outputValues,
                           const EOperator                 leftOp,
                           const EOperator                 rightOp,
                           MultiCell<Scalar> &             mCell,
                           const Teuchos::Array<Scalar> &  inputData,
                           const EDataFormat               inputFormat,
                           const EIntegrationDomain        intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief  Returns a rank-3 LexContainer such that <var>outputValues{c,l,r}</var> = \f$\mathbf{A}^c_{lr}\f$
  where:
    
    \li \f$ \displaystyle \mathbf{A}^c_{lr} = \int_{\Omega_c} (leftOp\;\phi_l) (rightOp\;{\phi_r}) d\Omega \f$    
    \li  \f$\phi_l\f$ and \f${\phi_r}\f$ are basis functions from the <strong>native</strong> local field \f$\mathcal{F}\f$;
    \li  <var>leftOp</var> and <var>rightOp</var> are admissible operator types for \f$\mathcal{F}\f$;
    \li  \f$\Omega_c\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The multi-index <var>{c,l,r}</var> has the following upper bounds:
    \li \f$ 0 \le c < \f$ number of integration domains (number of cells in the multicell is the default)
    \li \f$ 0 \le l,r < \dim \mathcal{F} \f$, i.e., the number of basis functions in the basis.
    
      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param mCell            [in]     - Multicell.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(LexContainer<Scalar> &      outputValues,
                           const EOperator             leftOp,
                           const EOperator             rightOp,
                           MultiCell <Scalar> &        mCell,
                           const EIntegrationDomain    intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief Returns a rank-3 LexContainer such that <var>outputValues{c,l,r}</var> = \f$\mathbf{A}^c_{lr}\f$
  where:
    
    \li \f$ \displaystyle \mathbf{A}^c_{lr} = \int_{\Omega_c} (leftOp\;\phi_l) (inputData) (rightOp\;{\varphi_r}) d\Omega \f$    
    \li  \f$\phi_l\f$ is basis function from the <strong>native</strong> local field \f$\mathcal{F}_l\f$
    \li  \f$\varphi_r\f$ is basis function from an <strong>auxiliary</strong> local field \f$\mathcal{F}_r\f$
    \li  <var>leftOp</var> is an admissible operator for functions from \f$\mathcal{F}_l\f$
    \li  <var>rightOp</var> is an admissible operator for functions from \f$\mathcal{F}_r\f$
    \li  <var>inputData</var> is a user-supplied data array;
    \li  \f$\Omega_c\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    Note that the domain of <var>leftOp</var> is always the parent local field. The domain, i.e. the 
    local field for <var>rightOp</var> must be specified via <var>rightOpField</var>.\n

    The multi-index <var>{c,l,r}</var> has the following upper bounds:
    
    \li \f$ 0 \le c < \f$ number of integration domains (number of cells in the multicell is the default)
    \li \f$ 0 \le l < \dim \mathcal{F}_l \f$, i.e., the number of basis functions in the parent (left) local field.
    \li \f$ 0 \le r < \dim \mathcal{F}_r \f$, i.e., the number of basis functions in the auxiliary (right) local field.
    
    The user-supplied data array should be formatted as follows:\n
    
    <table>
    <tr> <td>\f$ value(s, cell\_id, q_k) = inputData[cell\_id \cdot N + k] \f$</td>
    <td>if <var>inputFormat==DATA_SCALAR</var></td> </tr>
    <tr> <td>\f$ value(v_i, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d + k \cdot d + i] \f$</td>
    <td>if <var>inputFormat==DATA_VECTOR</var></td> </tr>
    <tr> <td>\f$ value(T_{ij}, cell\_id, q_k) = inputData[cell\_id \cdot N \cdot d^2 + k \cdot d^2 + i \cdot d + j] \f$</td>
    <td>if <var>inputFormat==DATA_TENSOR</var></td> </tr>
    </table>
    
    where \n
    \li   \f$s\f$ denotes a scalar quantity,\n
    \li   \f$v_i\f$ denotes the i-th component of a vector \f$v\f$,\n
    \li   \f$T_{ij}\f$ denotes the ij-th component of a tensor \f$T\f$,\n
    \li   \f$q_k\f$ denotes the k-th integration point,\n
    \li   \f$N\f$ denotes the number of integration points,\n
    \li   \f$d\f$ denotes the space dimension
    
      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param rightOpField     [in]     - Local field of the right operator.
      \param mCell            [in]     - Multicell.
      \param inputData        [in]     - Input data.
      \param inputFormat      [in]     - Format of input data (scalar, vector, etc.).
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(LexContainer<Scalar> &           outputValues,
                           const EOperator                  leftOp,
                           const EOperator                  rightOp,
                           const LocalField<Scalar> &       rightOpField,
                           MultiCell<Scalar> &              mCell,
                           const Teuchos::Array<Scalar> &   inputData,
                           const EDataFormat                inputFormat,
                           const EIntegrationDomain         intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief Returns field type.
  */
  virtual EField getFieldType() const = 0;

}; // class LocalField 

}// end namespace Intrepid

#endif
