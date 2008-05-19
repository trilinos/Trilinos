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
#include "Intrepid_FieldContainer.hpp"


namespace Intrepid {

/** \class Intrepid::LocalField
    \brief Defines the base class for the representation of local fields in Intrepid. 

    The local field interface relies on two member functions, <var>getOperator</var>
    and <var>getFunctional</var>, with varying signatures. Detailed info here ...
*/
template<class Scalar, class ArrayType = FieldContainer<Scalar> >
class LocalField {
  private:
    
  public:

  virtual ~LocalField() {}

  /** \brief Returns <var>ArrayType</var> with a multi-indexed quantity representing the values of a 
  (primitive) operator <var>primOp</var> applied to FEM basis functions, evaluated at the array 
  <var>inputPoints</var> of <strong>reference</strong> points. The rank of the output <var>ArrayType</var>
  and its dimensions vary depending on the type of <var>primOp</var> and the concrete field; see
  implementation for details.

  The context in which this interface function can be called is limited to FEM reconstructions 
  (i.e. via a reference cell). The function returns operator values <strong>decoupled from the 
  geometry of a particular cell</strong> (untransformed values). 

      \param outputValues    [out]     - Output container.
      \param inputPoints      [in]     - Array of input (reference) points.
      \param primOp           [in]     - Input operator (primitive).
  */
  virtual void getOperator(ArrayType &                             outputValues,
                           const Teuchos::Array<Point<Scalar> > &  inputPoints,
                           const EOperator                         primOp) = 0;


  /** \brief Returns <var>ArrayType</var> with a multi-indexed quantity representing the values of a 
  (primitive) operator <var>primOp</var> applied to FEM or FVD basis functions, evaluated at the 
  array <var>inputPoints</var> of <strong>physical</strong> points. The rank of the output <var>ArrayType</var>
  and its dimensions vary depending on the type of <var>primOp</var> and the concrete field; see
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
  virtual void getOperator(ArrayType &                             outputValues,
                           const Teuchos::Array<Point<Scalar> > &  inputPoints,
                           const EOperator                         primOp,
                           const Cell<Scalar> &                    cell) = 0;


  /** \brief Returns a rank-3 <var>ArrayType</var> such that <var>outputValues(C,L,R)</var> = \f$\mathbf{A}^C_{LR}\f$ where:
    
    \li \f$ \displaystyle \mathbf{A}^C_{LR} = \int_{\Omega_C} (leftOp\;\phi_L) (inputData) (rightOp\;{\phi_R}) d\Omega \f$    
    \li \f$\phi_L\f$ and \f${\phi_R}\f$ are basis functions from the <strong><strong>native</strong></strong> local field \f$\mathcal{F}\f$;
    \li <var>leftOp</var> and <var>rightOp</var> are admissible operator types for \f$\mathcal{F}\f$;
    \li <var>inputData</var> is a user-supplied data container;
    \li \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    Dimensions of the <var>outputValues</var> container are as follows:
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le L,R < \dim \mathcal{F} \f$, i.e., the number of basis functions in the <strong>native</strong> basis.
    
    The user-supplied data container should be formatted as follows:
    
    \li \f$ inputData(C,P)     = s(q_P)     |_{\Omega_C} \f$ if <var>s</var> is a rank-0 field (scalar function)
    \li \f$ inputData(C,P,I)   = v_I(q_P)   |_{\Omega_C} \f$ if <var>v</var> is a rank-1 field (vector function)
    \li \f$ inputData(C,P,I,J) = T_{IJ}(q_P)|_{\Omega_C} \f$ if <var>T</var> is a rank-2 field (tensor function)
        
    The rank of the user-supplied field is implicitly determined from the rank of <var>inputData</var>.
    Dimensions of the <var>inputData</var> are as follows:

    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le P < \f$ the number of cubature points on the reference integration domain
    \li \f$ 0 \le I,J < d \f$ the dimension of the admissible cell type of the LocalField (the space dimensions)    
    
      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param mCell            [in]     - Multicell.
      \param inputData        [in]     - Input data.
      \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(ArrayType &              outputValues,
                           const EOperator          leftOp,
                           const EOperator          rightOp,
                           MultiCell<Scalar> &      mCell,
                           const ArrayType &        inputData,
                           const bool               reuseJacobians = false,
                           const EIntegrationDomain intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief A simpler version of the previous method which does not take user-specified data. Returns 
    a rank-3 <var>ArrayType</var> such that <var>outputValues(C,L,R)</var> = \f$\mathbf{A}^C_{LR}\f$ where:
    
    \li \f$ \displaystyle \mathbf{A}^C_{LR} = \int_{\Omega_C} (leftOp\;\phi_L) (rightOp\;{\phi_R}) d\Omega \f$    
    \li  \f$\phi_L\f$ and \f${\phi_R}\f$ are basis functions from the <strong><strong>native</strong></strong> local field \f$\mathcal{F}\f$;
    \li  <var>leftOp</var> and <var>rightOp</var> are admissible operator types for \f$\mathcal{F}\f$;
    \li  \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    Dimensions of the <var>outputValues</var> container are as follows:
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le L,R < \dim \mathcal{F} \f$, i.e., the number of basis functions in the <strong>native</strong> basis.
    
      \param outputValues    [out]     - Output array.
      \param leftOp           [in]     - Left operator.
      \param rightOp          [in]     - Right operator.
      \param mCell            [in]     - Multicell.
      \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(ArrayType &                 outputValues,
                           const EOperator             leftOp,
                           const EOperator             rightOp,
                           MultiCell <Scalar> &        mCell,
                           const bool                  reuseJacobians = false,
                           const EIntegrationDomain    intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief Returns a rank-3 <var>ArrayType</var> such that <var>outputValues(C,L,R)</var> = \f$\mathbf{A}^C_{LR}\f$ where:
    
    \li \f$ \displaystyle \mathbf{A}^C_{LR} = \int_{\Omega_C} (leftOp\;\phi_L) (inputData) (rightOp\;{\varphi_R}) d\Omega \f$    
    \li  \f$\phi_L\f$ is basis function from the <strong><strong>native</strong></strong> local field \f$\mathcal{F}_L\f$
    \li  \f$\varphi_R\f$ is basis function from an <strong>auxiliary</strong> local field \f$\mathcal{F}_R\f$
    \li  <var>leftOp</var> is an admissible operator for functions from \f$\mathcal{F}_L\f$
    \li  <var>rightOp</var> is an admissible operator for functions from \f$\mathcal{F}_R\f$
    \li  <var>inputData</var> is a user-supplied data array;
    \li  \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    Note that the domain of <var>leftOp</var> is always the <strong>native</strong> local field. The domain, i.e. the 
    local field for <var>rightOp</var> must be specified via <var>rightOpField</var>.\n
    
    Dimensions of the <var>outputValues</var> container are as follows:
    
    \li \f$ 0 \le C < \f$ number of integration domains (number of cells in the multicell is the default)
    \li \f$ 0 \le L < \dim \mathcal{F}_L \f$, i.e., the number of basis functions in the <strong>native</strong> (left) local field.
    \li \f$ 0 \le R < \dim \mathcal{F}_R \f$, i.e., the number of basis functions in the <strong>auxiliary</strong> (right) local field.
    
    \warning The <strong><strong>native</strong></strong> and <strong>auxiliary</strong> LocalField objects 
    <strong>must be instantiated with identical cell types and cubature sets!</strong> 
    
    The user-supplied data container should be formatted as follows:
    
    \li \f$ inputData(C,P)     = s(q_P)     |_{\Omega_C} \f$ if <var>s</var> is a rank-0 field (scalar function)
    \li \f$ inputData(C,P,I)   = v_I(q_P)   |_{\Omega_C} \f$ if <var>v</var> is a rank-1 field (vector function)
    \li \f$ inputData(C,P,I,J) = T_{IJ}(q_P)|_{\Omega_C} \f$ if <var>T</var> is a rank-2 field (tensor function)
    
    The rank of the user-supplied field is implicitly determined from the rank of <var>inputData</var>.
    Dimensions of the <var>inputData</var> are as follows:
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le P < \f$ the number of cubature points on the reference integration domain
    \li \f$ 0 \le I,J < d \f$ the dimension of the admissible cell type of the LocalField (the space dimensions)    
    
    \param outputValues    [out]     - Output array.
    \param leftOp           [in]     - Left operator.
    \param rightOp          [in]     - Right operator.
    \param rightOpField     [in]     - Local field of the right operator.
    \param mCell            [in]     - Multicell.
    \param inputData        [in]     - Input data.
    \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    */
  virtual void getOperator(ArrayType &                outputValues,
                           const EOperator            leftOp,
                           const EOperator            rightOp,
                           const LocalField<Scalar> & rightOpField,
                           MultiCell<Scalar> &        mCell,
                           const ArrayType &          inputData,
                           const bool                 reuseJacobians = false,
                           const EIntegrationDomain   intDomain = INTEGRATION_DOMAIN_CELL) = 0;
  
  
  /** \brief A simpler version of the previous method which does not take user specified data. Returns 
    a rank-3 <var>ArrayType</var> such that <var>outputValues(C,L,R)</var> = \f$\mathbf{A}^C_{LR}\f$ where:
    
    \li \f$ \displaystyle \mathbf{A}^C_{LR} = \int_{\Omega_C} (leftOp\;\phi_L) (inputData) (rightOp\;{\varphi_R}) d\Omega \f$    
    \li  \f$\phi_L\f$ is basis function from the <strong><strong>native</strong></strong> local field \f$\mathcal{F}_L\f$
    \li  \f$\varphi_R\f$ is basis function from an <strong>auxiliary</strong> local field \f$\mathcal{F}_R\f$
    \li  <var>leftOp</var> is an admissible operator for functions from \f$\mathcal{F}_L\f$
    \li  <var>rightOp</var> is an admissible operator for functions from \f$\mathcal{F}_R\f$
    \li  \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    Note that the domain of <var>leftOp</var> is always the <strong>native</strong> local field. The domain, i.e. the 
    local field for <var>rightOp</var> must be specified via <var>rightOpField</var>.\n
    
    Dimensions of the <var>outputValues</var> container are as follows:
    
    \li \f$ 0 \le C < \f$ number of integration domains (number of cells in the multicell is the default)
    \li \f$ 0 \le L < \dim \mathcal{F}_L \f$, i.e., the number of basis functions in the <strong>native</strong> (left) local field.
    \li \f$ 0 \le R < \dim \mathcal{F}_R \f$, i.e., the number of basis functions in the <strong>auxiliary</strong> (right) local field.
    
    \warning The <strong><strong>native</strong></strong> and <strong>auxiliary</strong> LocalField objects 
    <strong>must be instantiated with identical cell types and cubature sets!</strong> 

    \param outputValues    [out]     - Output array.
    \param leftOp           [in]     - Left operator.
    \param rightOp          [in]     - Right operator.
    \param rightOpField     [in]     - Local field of the right operator.
    \param mCell            [in]     - Multicell.
    \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    */
  virtual void getOperator(ArrayType &                  outputValues,
                           const EOperator              leftOp,
                           const EOperator              rightOp,
                           const LocalField<Scalar> &   rightOpField,
                           MultiCell<Scalar> &          mCell,
                           const bool                   reuseJacobians = false,
                           const EIntegrationDomain     intDomain = INTEGRATION_DOMAIN_CELL) = 0;
  

  /** \brief Returns field type.
  */
  virtual EField getFieldType() const = 0;
  
  
  /** \brief Returns cell type for whch the field was instantiated
  */
  virtual ECell getCellType() const = 0;
  
  
  /** \brief Returns number of cubature points associated with a subcell.
    
    \param subcellDim       [in]     - dimension of the subcell
    \param subcellId        [in]     - order of the subcell relative to the cell
  */
  virtual int getNumCubPoints(const int subcellDim,
                              const int subcellId) const = 0;
  
}; // class LocalField 

}// end namespace Intrepid

#endif
