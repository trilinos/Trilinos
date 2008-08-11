// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copytest (2007) Sandia Corporation
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
    and <var>getFunctional</var>, with varying signatures. The second template
    parameter specifies the type of the output container to be used by these methods. The default
    is Intrepid's FieldContainer class. Detailed info here ...
*/
template<class Scalar, class ArrayType = FieldContainer<Scalar> >
class LocalField {
  private:
    
  public:

  virtual ~LocalField() {}

  /** \brief Fills <var>outputValues</var> with values of <var>primOp</var>, acting on a set of FEM 
  basis functions, at a set of <strong>reference cell</strong> points <var>inputPoints</var>. Rank of the
  output ArrayType container depends on the ranks of the basis fields and the specified operator and
  ranges from 2 to 5. The first two dimensions are always the number of points in <var>inputPoints</var> 
  and dimension of the basis set.

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


  /** \brief Fills <var>outputValues</var> with values of a (primitive) operator <var>primOp</var> 
  acting on FEM or FVD basis functions, at a set of <strong>physical cell</strong> points 
  <var>inputPoints</var>. Rank of the output ArrayType container depends on the ranks of the basis 
  fields and the specified operator and ranges from 2 to 5. The first two dimensions are always the 
  number of points in <var>inputPoints</var> and dimension of the basis set.
 
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


  /** \brief Computes a local discrete operator acting on the </strong>native</strong> basis set and 
    whsoe range has the same dimension as </strong>native</strong> basis set, i.e., a square \f$ N\times N\f$
    matrix where <var>N</var> is the number of </strong>native</strong> basis functions.  
    The elements of the \f$ N\times N\f$ matrix are defined as follows:
    
    \f[\displaystyle\mathbf{A}^C_{IJ}=\int_{\Omega_C}(trialOp\;\phi_J)(inputData)(testOp\;{\phi_I})d\Omega\f]
    
    where  
    \li \f$\phi_J\f$ and \f${\phi_I}\f$ are basis functions in the <strong>native</strong> local field \f$\mathcal{F}_n\f$;
    \li <var>trialOp</var> and <var>testOp</var> are admissible operator types for \f$\mathcal{F}_n\f$;
    \li <var>inputData</var> is a user-specified data container;
    \li \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The set of local operators is returned as a rank-3 <var>ArrayType</var> container such that 
    <var>outputValues(C,I,J)</var> = \f$\mathbf{A}^C_{IJ}\f$. Dimensions of the container are as follows

    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le I,J < N =\dim \mathcal{F}_n\f$, i.e., dimension of the <strong>native</strong> local field basis.
    
    \remarks 
    The user-specified field data in <var>inputData</var> must be compatible with <var>trialOp</var> and 
    <var>testOp</var> in the sense that the operation under the integral must contract to a scalar value.
    The user-specified data container should be formatted as follows:
    
    \li \f$ inputData(C,P)       = s(q_P)        |_{\Omega_C} \f$ if <var>s</var> is a rank-0 field (scalar function)
    \li \f$ inputData(C,P,D1)    = v_{D1}(q_P)   |_{\Omega_C} \f$ if <var>v</var> is a rank-1 field (vector function)
    \li \f$ inputData(C,P,D1,D2) = T_{D1,D2}(q_P)|_{\Omega_C} \f$ if <var>T</var> is a rank-2 field (tensor function)
      
    Dimensions of the <var>inputData</var> are as follows:
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le P < \f$ the number of cubature points on the reference integration domain
    \li \f$ 0 \le D1,D2 < D \f$ the dimension of the admissible cell type of the LocalField (the space dimensions)    
    
    The rank of the user-specified field is implicitly determined from the rank of <var>inputData</var>.
    
      \param outputValues    [out]     - Output container.
      \param trialOp          [in]     - Trial operator.
      \param testOp           [in]     - Test operator.
      \param mCell            [in]     - Multicell.
      \param inputData        [in]     - Input data container.
      \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(ArrayType &              outputValues,
                           const EOperator          trialOp,
                           const EOperator          testOp,
                           MultiCell<Scalar> &      mCell,
                           const ArrayType &        inputData,
                           const bool               reuseJacobians = false,
                           const EIntegrationDomain intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief A simpler version of the previous method in which definition of the local discrete
    operator does not require user-specified data. The elements of the \f$ N\times N\f$ matrix are defined as follows:
    
    \f[\displaystyle\mathbf{A}^C_{IJ}=\int_{\Omega_C}(trialOp\;\phi_J)(testOp\;{\phi_I}) d\Omega \f]
    
    where
    \li  \f$\phi_J\f$ and \f${\phi_I}\f$ are basis functions in the <strong>native</strong> local field \f$\mathcal{F}_n\f$;
    \li  <var>trialOp</var> and <var>testOp</var> are admissible operator types for \f$\mathcal{F}_n\f$; 
    \li  \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The set of local operators is returned as a rank-3 <var>ArrayType</var> container such that 
    <var>outputValues(C,I,J)</var> = \f$\mathbf{A}^C_{IJ}\f$. Dimensions of the container are as follows
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le I,J < N =\dim \mathcal{F}_n \f$, i.e., dimension of the <strong>native</strong> local field basis.
    
    \remarks
    <var>trialOp</var> and <var>testOp</var> must be compatible with each other in the sense that the 
    operation under the integral must contract to a scalar value.
    
      \param outputValues    [out]     - Output container.
      \param trialOp          [in]     - Trial operator.
      \param testOp           [in]     - Test operator.
      \param mCell            [in]     - Multicell.
      \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
      \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getOperator(ArrayType &                 outputValues,
                           const EOperator             trialOp,
                           const EOperator             testOp,
                           MultiCell <Scalar> &        mCell,
                           const bool                  reuseJacobians = false,
                           const EIntegrationDomain    intDomain = INTEGRATION_DOMAIN_CELL) = 0;


  /** \brief  Computes a local discrete operator acting on the <strong>native</strong> basis set and 
    whose range has the dimension of the <strong>auxiliary</strong> basis set, i.e., a rectangular
    \f$ A\times N\f$ matrix where <var>N</var> and <var>A</var> are the numbers of <strong>native</strong>
    and <strong>auxiliary</strong> basis functions, respectively.  The elements of the \f$ A\times N\f$ 
    matrix are defined as follows:
    
    \f[\displaystyle\mathbf{A}^C_{IJ}=\int_{\Omega_C}(trialOp\;\phi_J)(inputData)(testOp\;{\varphi_I})d\Omega\f]
    
    where
    \li  \f$\phi_J\f$ is basis function from the <strong>native</strong> local field \f$\mathcal{F}_n\f$
    \li  \f$\varphi_I\f$ is basis function from an <strong>auxiliary</strong> local field \f$\mathcal{F}_a\f$
    \li  <var>trialOp</var> is an admissible operator for functions in \f$\mathcal{F}_n\f$
    \li  <var>testOp</var> is an admissible operator for functions in \f$\mathcal{F}_a\f$
    \li  <var>inputData</var> is a user-specified data container;
    \li  \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The set of local operators is returned as a rank-3 <var>ArrayType</var> container such that 
    <var>outputValues(C,I,J)</var> = \f$\mathbf{A}^C_{IJ}\f$. Dimensions of the container are as follows
    
    \li \f$ 0 \le C < \f$ number of integration domains (number of cells in the multicell is the default)
    \li \f$ 0 \le J < N = \dim \mathcal{F}_n \f$, i.e., dimension of the <strong>native</strong> (trial)  local field basis.
    \li \f$ 0 \le I < A = \dim \mathcal{F}_a \f$, i.e., dimension of the <strong>auxiliary</strong> (test)  local field basis.
    
    \remarks 
    The user-specified field data in <var>inputData</var> must be compatible with <var>trialOp</var> and 
    <var>testOp</var> in the sense that the operation under the integral must contract to a scalar value.
    The user-specified data container should be formatted as follows:
    
    \li \f$ inputData(C,P)       = s(q_P)        |_{\Omega_C} \f$ if <var>s</var> is a rank-0 field (scalar function)
    \li \f$ inputData(C,P,D1)    = v_{D1}(q_P)   |_{\Omega_C} \f$ if <var>v</var> is a rank-1 field (vector function)
    \li \f$ inputData(C,P,D1,D2) = T_{D1,D2}(q_P)|_{\Omega_C} \f$ if <var>T</var> is a rank-2 field (tensor function)
    
    Dimensions of the <var>inputData</var> are as follows:
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le P < \f$ the number of cubature points on the reference integration domain
    \li \f$ 0 \le D1,D2 < D \f$ the dimension of the admissible cell type of the LocalField (the space dimensions)    
    
    The rank of the user-specified field is implicitly determined from the rank of <var>inputData</var>.
    
    The domain of <var>trialOp</var> is always the <strong>native</strong> LocalField, whereas the 
    domain of <var>testOp</var> is an <strong>auxiliary</strong> LocalField that must be specified 
    in <var>testOpField</var>. 
    
    \warning The <strong><strong>native</strong></strong> and <strong>auxiliary</strong> LocalField objects 
    <strong>must be instantiated with identical cell types and cubature sets!</strong> 
    
    \param outputValues    [out]     - Output container.
    \param trialOp          [in]     - Trial operator.
    \param testOp           [in]     - Test operator.
    \param testOpField      [in]     - Local field of the test operator.
    \param mCell            [in]     - Multicell.
    \param inputData        [in]     - Input data container.
    \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    */
  virtual void getOperator(ArrayType &                outputValues,
                           const EOperator            trialOp,
                           const EOperator            testOp,
                           const LocalField<Scalar> & testOpField,
                           MultiCell<Scalar> &        mCell,
                           const ArrayType &          inputData,
                           const bool                 reuseJacobians = false,
                           const EIntegrationDomain   intDomain = INTEGRATION_DOMAIN_CELL) = 0;
  
  
  /** \brief A simpler version of the previous method in which definition of the local discrete
    operator does not require user-specified data. The elements of the \f$ A\times N\f$ matrix are defined as follows:
    
    \f[\displaystyle\mathbf{A}^C_{IJ}=\int_{\Omega_C}(trialOp\;\phi_J)(inputData)(testOp\;{\varphi_I})d\Omega\f]
    
    where
    \li  \f$\phi_J\f$ is basis function from the <strong>native</strong> local field \f$\mathcal{F}_n\f$
    \li  \f$\varphi_I\f$ is basis function from an <strong>auxiliary</strong> local field \f$\mathcal{F}_a\f$
    \li  <var>trialOp</var> is an admissible operator for functions in \f$\mathcal{F}_n\f$
    \li  <var>testOp</var> is an admissible operator for functions in \f$\mathcal{F}_a\f$
    \li  \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The set of local operators is returned as a rank-3 <var>ArrayType</var> container such that 
    <var>outputValues(C,I,J)</var> = \f$\mathbf{A}^C_{IJ}\f$. Dimensions of the container are as follows
    
    \li \f$ 0 \le C < \f$ number of integration domains (number of cells in the multicell is the default)
    \li \f$ 0 \le J < N = \dim \mathcal{F}_n \f$, i.e., dimension of the <strong>native</strong> (trial)  local field basis.
    \li \f$ 0 \le I < A = \dim \mathcal{F}_a \f$, i.e., dimension of the <strong>auxiliary</strong> (test)  local field basis.
    
    \remarks
    <var>trialOp</var> and <var>testOp</var> must be compatible with each other in the sense that the 
    operation under the integral must contract to a scalar value.
    
    The domain of <var>trialOp</var> is always the <strong>native</strong> LocalField, whereas the 
    domain of <var>testOp</var> is an <strong>auxiliary</strong> LocalField that must be specified 
    in <var>testOpField</var>. 
    
    \warning The <strong><strong>native</strong></strong> and <strong>auxiliary</strong> LocalField objects 
    <strong>must be instantiated with identical cell types and cubature sets!</strong> 

    \param outputValues    [out]     - Output container.
    \param trialOp          [in]     - Trial operator.
    \param testOp           [in]     - Test operator.
    \param testOpField      [in]     - Local field of the test operator.
    \param mCell            [in]     - Multicell.
    \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    */
  virtual void getOperator(ArrayType &                  outputValues,
                           const EOperator              trialOp,
                           const EOperator              testOp,
                           const LocalField<Scalar> &   testOpField,
                           MultiCell<Scalar> &          mCell,
                           const bool                   reuseJacobians = false,
                           const EIntegrationDomain     intDomain = INTEGRATION_DOMAIN_CELL) = 0;
  
  
  /** \brief Computes a local discrete functional acting on the <strong>native</strong> basis set, i.e.,
    a vector with \f$N\f$ components where <var>N</var> is the number of </strong>native</strong> basis 
    functions. The elements of this <var>N</var>-dimensional vector are defined as follows:
    
    \f[ \mathbf{f}^C_{I} = \int_{\Omega_C} (trialData) (testOp\;{\phi_I}) d\Omega \f]
    
    where
    \li <var>trialData</var> is a user-specified data container;
    \li \f${\phi_I}\f$ is basis function from the <strong><strong>native</strong></strong> local field \f$\mathcal{F}_n\f$;
    \li <var>testOp</var> is an admissible primitive operator for \f$\mathcal{F}_n\f$;
    \li \f$\Omega_C\f$ is a subcell that represents a valid integration domain for the resulting expression.
    
    The set of local functionals is returned as a rank-2 <var>ArrayType</var> container such that 
    <var>outputValues(C,I)</var> = \f$\mathbf{f}^C_{I}\f$. Dimensions of the container are as follows:
  
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le I < N = \dim \mathcal{F}_n\f$, i.e., the dimension of the <strong>native</strong> local field basis.
    
    \remarks
    The user-specified data in <var>trialData</var> must be compatible with <var>testOp</var> in the sense
    that the operation under the integral must contract to a scalar value. 
    The user-specified data container should be formatted as follows:
    
    \li \f$ inputData(C,P)       = s(q_P)        |_{\Omega_C} \f$ if <var>s</var> is a rank-0 field (scalar function)
    \li \f$ inputData(C,P,D1)    = v_{D1}(q_P)   |_{\Omega_C} \f$ if <var>v</var> is a rank-1 field (vector function)
    \li \f$ inputData(C,P,D1,D2) = T_{D1,D2}(q_P)|_{\Omega_C} \f$ if <var>T</var> is a rank-2 field (tensor function)
    
    Dimensions of the <var>inputData</var> are as follows:
    
    \li \f$ 0 \le C < \f$ the number of integration domains (default is number of cells in the multicell)
    \li \f$ 0 \le P < \f$ the number of cubature points on the reference integration domain
    \li \f$ 0 \le D1,D2 < D \f$ the dimension of the admissible cell type of the LocalField (the space dimensions)  
    
    The rank of the user-specified field is implicitly determined from the rank of <var>inputData</var>.
    
    \param outputValues    [out]     - Output container.
    \param trialData        [in]     - User provided data for the "trial" argument.
    \param testOp           [in]     - primitive operator for the "test" argument.
    \param mCell            [in]     - Multicell.
    \param reuseJacobians   [in]     - Forces reuse of Jacobian and subcell measure values at cub. pts.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  virtual void getFunctional(ArrayType &              outputValues,
                             const ArrayType &        inputData,
                             const EOperator          testOp,
                             MultiCell<Scalar> &      mCell,
                             const bool               reuseJacobians = false,
                             const EIntegrationDomain intDomain = INTEGRATION_DOMAIN_CELL) = 0;
  
  
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
