// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __INTREPID_HGRAD_C0_FEM_HPP__
#define __INTREPID_HGRAD_C0_FEM_HPP__

#include "Intrepid_Basis.hpp"

namespace Intrepid{
  template<class Scalar, class ArrayScalar>
  class Basis_HGRAD_C0_FEM : public Basis<Scalar, ArrayScalar>, public DofCoordsInterface<ArrayScalar> {
  private:

    /** \brief Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
     */
    void initializeTags();

  public:

    /** \brief Constructor.
    */
    Basis_HGRAD_C0_FEM();
    /** \brief  FEM basis evaluation on a <strong>reference Quadrilateral</strong> cell.

                Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
                points in the <strong>reference Quadrilateral</strong> cell. For rank and dimensions of
                I/O array arguments see Section \ref basis_md_array_sec .

        \param  outputValues      [out] - rank-2 or 3 array with the computed basis values
        \param  inputPoints       [in]  - rank-2 array with dimensions (P,D) containing reference points
        \param  operatorType      [in]  - operator applied to basis functions
     */
    void getValues(ArrayScalar &          outputValues,
                   const ArrayScalar &    inputPoints,
                   const EOperator        operatorType) const;


    /**  \brief  FVD basis evaluation: invocation of this method throws an exception.
     */
    void getValues(ArrayScalar &          outputValues,
                   const ArrayScalar &    inputPoints,
                   const ArrayScalar &    cellVertices,
                   const EOperator        operatorType = OPERATOR_VALUE) const;

    /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
                <strong>reference Quadrilateral</strong>.

        \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
                                       dimensioned (F,D)
    */
    void getDofCoords(ArrayScalar & DofCoords) const;

  };


}


#include "Intrepid_HGRAD_C0_FEMDef.hpp"
#endif
