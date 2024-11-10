// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_LagrangianInterpolation.hpp
    \brief  Header file for the Intrepid2::LagrangianInterpolation class.
    \author Created by Mauro Perego
 */
#ifndef __INTREPID2_LAGRANGIANINTERPOLATION_HPP__
#define __INTREPID2_LAGRANGIANINTERPOLATION_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_PointTools.hpp"

#include "Intrepid2_Basis.hpp"

// -- HGRAD family
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

// -- HCURL family
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"

#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_LINE_Cn_FEM.hpp"

// -- HDIV family
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"

#include "Intrepid2_HDIV_TRI_In_FEM.hpp"
#include "Intrepid2_HDIV_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"

// -- Lower order family
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"

#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HCURL_WEDGE_I1_FEM.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HDIV_WEDGE_I1_FEM.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Intrepid2_OrientationTools.hpp"


namespace Intrepid2 {

/** \class  Intrepid2::LagrangianInterpolation and LagrangianTools classes
    \brief  A class providing static members to perform Lagrangian interpolation on a finite element.


    The Lagrangian interpolant \f$\mathcal I_h\f$ is defined as
    \f[
    \mathcal I_h(f) = \sum_i \alpha^f_i \phi_i, \quad \alpha_i^f := L_i(f),
    \f]
    where \f$\{\phi_i\}\f$ is the basis of the finite element, \f$\alpha_i^f\f$ are the
    <var><b>basisCoeffs</b></var> and \f$L_i\f$ are the DOFs of the basis.
    This class works only for Lagrangian elements that provide a set of orthonormal DOFs defined as
    \f[
    L_i(f) := f(\mathbf x_i) \cdot \beta_i, \quad L_i(\phi_j) = \delta_{ij},
    \f]
    where \f$\beta_i\f$ are referred to as <var><b>dofCoeffs</b></var>, and \f$\mathbf x_i\f$ are the coordinates of the basis nodes.

    In order to perform the interpolation, evaluate the function \f$f\f$ at the set of points \f$\{\mathbf x_i\}\f$ computed using the basis method <var><b>getDoCoords</b></var>
    and then obtain the basis coefficients \f$\alpha_i^f\f$ by calling the function <var><b>getBasisCoeffs</b></var>.

    \remark The interpolation is performed at the reference element. Therefore, the function \f$f\f$,
            which is contravariant, needs to mapped to the reference space, using the inverse operation of a pullback, before calling <var><b>getBasisCoeffs</b></var>.
 */

template<typename DeviceType>
class LagrangianInterpolation {
public:

  /** \brief  Computes the basis weights of the function interpolation.

      \code
      C  - num. cells
      F  - num. fields
      D  - spatial dimension
      \endcode

      \param  basisCoeffs         [out] - rank-2 view (C,F) that will contain the basis coefficients of the interpolation.
      \param  functionAtDofCoords [in]  - variable rank view that contains the function evaluated at reference DOF coordinates.
      \param  cellBasis           [in]  - pointer to the basis for the interpolation
      \param  cellOrientations    [in]  - rank-1 view (C) containing the Orientation objects at each cell

      \remark The output views need to be pre-allocated. <var><b>dofCoeffs</b></var> and <var><b>functionAtDofCoords</b></var> have
              rank 2, (C,F) for scalar basis and  3, (C,F,D) for vector basis.
              <var><b>functValsAtDofCoords</b></var> contains the function evaluated at the reference <var><b>dofCoords</b></var> and contravariantly transformed
              to the reference element. The reference <var><b>dofCoords</b></var> are obtained with the method cellBasis->getDofCoords(...), NOT with the getOrientedDofCoords() method
   */
  template<typename basisCoeffsViewType,
  typename funcViewType,
  typename BasisType,
  typename ortViewType>
  static void
  getBasisCoeffs(basisCoeffsViewType basisCoeffs,
      const funcViewType functionAtDofCoords,
      const BasisType* cellBasis,
      const ortViewType orts);
  };

/** \class  Intrepid2::LagrangianTools
  \brief  A class providing tools for Lagrangian elements as static members.

  Lagrangian orthonormal DOFs are defined as
  \f[
  L_i(f) := f(\mathbf x_i) \cdot \beta_i, \quad L_i(\phi_j) = \delta_{ij},
  \f]
  where \f$\beta_i\f$ are referred to as <var><b>dofCoeffs</b></var>, and \f$\mathbf x_i\f$ as the <var><b>dofCoords</b></var>.

  This class provides tools to compute <var><b>dofCoords</b></var> and <var><b>dofCoeffs</b></var> for the <b>oriented</b> reference element.
*/

template<typename DeviceType>
class LagrangianTools {
public:

  /** \brief  Computes the coordinates associated with the basis DOFs for the reference oriented element

      \code
      C  - num. cells
      F  - num. fields
      D  - spatial dimension
      \endcode

      \param  dofCoords        [out] - rank-3 view (C,F,D), that will contain coordinates associated with the basis DOFs.
      \param  cellBasis        [in]  - pointer to the basis for the interpolation
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell

      \remark the output views need to be pre-allocated.
   */
  template<typename BasisType,
  class ...coordsProperties,
  typename ortValueType, class ...ortProperties>
  static void
  getOrientedDofCoords(
      Kokkos::DynRankView<typename BasisType::scalarType, coordsProperties...> dofCoords,
      const BasisType* cellBasis,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations
  );

      /** \brief  Computes the coefficients associated with the basis DOFs for the reference oriented element

      \code
      C  - num. cells
      F  - num. fields
      D  - spatial dimension
      \endcode

      \param  dofCoeffs        [out] - variable rank view that will contain coefficients associated with the basis DOFs.
      \param  cellBasis        [in]  - pointer to the basis for the interpolation
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell

      \remark the output views need to be pre-allocated. <var><b>dofCoeffs</b></var> has rank 2, (C,F) for scalar basis and  3,
              (C,F,D) for vector basis.
   */
  template<typename BasisType,
  class ...coeffsProperties,
  typename ortValueType, class ...ortProperties>
  static void
  getOrientedDofCoeffs(
      Kokkos::DynRankView<typename BasisType::scalarType, coeffsProperties...> dofCoeffs,
      const BasisType* cellBasis,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations
  );
};

} // Intrepid2 namespace

// include templated function definitions
#include "Intrepid2_LagrangianInterpolationDef.hpp"

#endif





