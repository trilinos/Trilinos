// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasisFamily.hpp
    \brief  Stateless class representing a family of basis functions, templated on H(vol) and H(grad) on the line.  Only hypercube topologies are supported at the moment, but the intent is ultimately to support all standard topologies.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_DerivedBasisFamily_h
#define Intrepid2_DerivedBasisFamily_h

#include "Intrepid2_Basis.hpp"

#include "Intrepid2_DerivedBasis_HGRAD_QUAD.hpp"
#include "Intrepid2_DerivedBasis_HCURL_QUAD.hpp"
#include "Intrepid2_DerivedBasis_HDIV_QUAD.hpp"
#include "Intrepid2_DerivedBasis_HVOL_QUAD.hpp"

#include "Intrepid2_DerivedBasis_HGRAD_HEX.hpp"
#include "Intrepid2_DerivedBasis_HCURL_HEX.hpp"
#include "Intrepid2_DerivedBasis_HDIV_HEX.hpp"
#include "Intrepid2_DerivedBasis_HVOL_HEX.hpp"

#include "Intrepid2_DerivedBasis_HGRAD_WEDGE.hpp"
#include "Intrepid2_DerivedBasis_HCURL_WEDGE.hpp"
#include "Intrepid2_DerivedBasis_HDIV_WEDGE.hpp"
#include "Intrepid2_DerivedBasis_HVOL_WEDGE.hpp"

#include "Intrepid2_SerendipityBasis.hpp"

namespace Intrepid2
{
  //! \brief EmptyBasisFamily allows us to set a default void family for a given topology
    class EmptyBasisFamily
    {
    public:
      using HGRAD = void;
      using HCURL = void;
      using HDIV  = void;
      using HVOL  = void;
    };
  
  /** \class Intrepid2::DerivedBasisFamily
      \brief A family of basis functions, constructed from H(vol) and H(grad) bases on the line.
   
   At present, only hypercube topologies (line, quadrilateral, hexahedron) are supported, but other topologies will be supported in the future.
  */
  template<class LineBasisHGRAD, class LineBasisHVOL, class TriangleBasisFamily = EmptyBasisFamily, class TetrahedronBasisFamily = EmptyBasisFamily, class PyramidBasisFamily = EmptyBasisFamily>
  class DerivedBasisFamily
  {
  public:
    using ExecutionSpace  = typename LineBasisHGRAD::ExecutionSpace;
    using OutputValueType = typename LineBasisHGRAD::OutputValueType;
    using PointValueType  = typename LineBasisHGRAD::PointValueType;
    
    using Basis      = typename LineBasisHGRAD::BasisBase;
    using BasisPtr   = Teuchos::RCP<Basis>;
    using DeviceType = typename Basis::DeviceType;
    
    // line bases
    using HGRAD_LINE = LineBasisHGRAD;
    using HVOL_LINE  = LineBasisHVOL;
    
    // quadrilateral bases
    using HGRAD_QUAD = Basis_Derived_HGRAD_QUAD<HGRAD_LINE>;
    using HCURL_QUAD = Basis_Derived_HCURL_QUAD<HGRAD_LINE, HVOL_LINE>;
    using HDIV_QUAD  = Basis_Derived_HDIV_QUAD <HGRAD_LINE, HVOL_LINE>;
    using HVOL_QUAD  = Basis_Derived_HVOL_QUAD <HVOL_LINE>;
    
    // hexahedron bases
    using HGRAD_HEX = Basis_Derived_HGRAD_HEX<HGRAD_LINE>;
    using HCURL_HEX = Basis_Derived_HCURL_HEX<HGRAD_LINE, HVOL_LINE>;
    using HDIV_HEX  = Basis_Derived_HDIV_HEX <HGRAD_LINE, HVOL_LINE>;
    using HVOL_HEX  = Basis_Derived_HVOL_HEX <HVOL_LINE>;
    
    // triangle bases
    using HGRAD_TRI = typename TriangleBasisFamily::HGRAD;
    using HCURL_TRI = typename TriangleBasisFamily::HCURL;
    using HDIV_TRI = typename TriangleBasisFamily::HDIV;
    using HVOL_TRI = typename TriangleBasisFamily::HVOL;
    
    // tetrahedron bases
    using HGRAD_TET = typename TetrahedronBasisFamily::HGRAD;
    using HCURL_TET = typename TetrahedronBasisFamily::HCURL;
    using HDIV_TET = typename TetrahedronBasisFamily::HDIV;
    using HVOL_TET = typename TetrahedronBasisFamily::HVOL;
    
    // wedge bases
    using HGRAD_WEDGE = Basis_Derived_HGRAD_WEDGE<HGRAD_TRI, HGRAD_LINE>;
    using HCURL_WEDGE = Basis_Derived_HCURL_WEDGE<HGRAD_TRI, HCURL_TRI, HGRAD_LINE, HVOL_LINE>;
    using HDIV_WEDGE  = Basis_Derived_HDIV_WEDGE < HDIV_TRI, HVOL_TRI,  HGRAD_LINE, HVOL_LINE>;
    using HVOL_WEDGE  = Basis_Derived_HVOL_WEDGE < HVOL_TRI, HVOL_LINE>;
    
    // pyramid bases
    using HGRAD_PYR = typename PyramidBasisFamily::HGRAD;
    using HCURL_PYR = typename PyramidBasisFamily::HCURL;
    using HDIV_PYR  = typename PyramidBasisFamily::HDIV;
    using HVOL_PYR  = typename PyramidBasisFamily::HVOL;
  };
  
  /** \brief  Factory method for line bases in the given family.
      \param [in] fs        - the function space for the basis.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getLineBasis(Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_LINE (polyOrder,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_LINE(polyOrder,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }

  /** \brief  Factory method for isotropic quadrilateral bases in the given family.
      \param [in] fs        - the function space for the basis.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getQuadrilateralBasis(Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_QUAD (polyOrder,pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_QUAD(polyOrder,pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_QUAD (polyOrder,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_QUAD(polyOrder,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for potentially anisotropic quadrilateral bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder_x - the polynomial order of the basis in the x dimension.
      \param [in] polyOrder_y - the polynomial order of the basis in the y dimension.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getQuadrilateralBasis(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_QUAD (polyOrder_x,polyOrder_y,pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_QUAD(polyOrder_x,polyOrder_y,pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_QUAD (polyOrder_x,polyOrder_y,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_QUAD(polyOrder_x,polyOrder_y,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic bases on the hexahedron in the given family.
      \param [in] fs        - the function space for the basis.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getHexahedronBasis(Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_HEX (polyOrder,pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_HEX(polyOrder,pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_HEX (polyOrder,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_HEX(polyOrder,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }

  /** \brief  Factory method for isotropic HGRAD bases on a hypercube for the given family.  Note that this will return a Line<2> for its base cell topology.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] spaceDim - the number of dimensions for the hypercube on which the basis is defined.
      \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getHypercubeBasis_HGRAD(int polyOrder, int spaceDim, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    using BasisPtr = typename BasisFamily::BasisPtr;
    
    BasisPtr lineBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
    BasisPtr tensorBasis = lineBasis;
    
    for (int d=1; d<spaceDim; d++)
    {
      tensorBasis = Teuchos::rcp(new Basis_TensorBasis<BasisBase>(tensorBasis, lineBasis, FUNCTION_SPACE_HGRAD));
    }
    
    return tensorBasis;
  }

  /** \brief  Factory method for isotropic HVOL bases on a hypercube for the given family.  Note that this will return a Line<2> for its base cell topology.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] spaceDim - the number of dimensions for the hypercube on which the basis is defined.
      \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getHypercubeBasis_HVOL(int polyOrder, int spaceDim, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    using BasisPtr = typename BasisFamily::BasisPtr;
    
    BasisPtr lineBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
    BasisPtr tensorBasis = lineBasis;
    
    for (int d=1; d<spaceDim; d++)
    {
      tensorBasis = Teuchos::rcp(new Basis_TensorBasis<BasisBase>(tensorBasis, lineBasis, FUNCTION_SPACE_HVOL));
    }
    
    return tensorBasis;
  }
  
  /** \brief  Factory method for potentially anisotropic hexahedron bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder_x - the polynomial order of the basis in the x dimension.
      \param [in] polyOrder_y - the polynomial order of the basis in the y dimension.
      \param [in] polyOrder_z - the polynomial order of the basis in the z dimension.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getHexahedronBasis(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_HEX (polyOrder_x,polyOrder_y,polyOrder_z,pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_HEX(polyOrder_x,polyOrder_y,polyOrder_z,pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_HEX (polyOrder_x,polyOrder_y,polyOrder_z,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_HEX(polyOrder_x,polyOrder_y,polyOrder_z,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }

  /** \brief  Factory method for isotropic HGRAD Serendipity bases on a hypercube for the given family.  Note that this will return a Line<2> for its base cell topology.  Note also that the family must use hierarchical bases.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] spaceDim - the number of dimensions for the hypercube on which the basis is defined.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getSerendipityBasis_HGRAD(int polyOrder, int spaceDim)
  {
    auto fullBasis = getHypercubeBasis_HGRAD<BasisFamily>(polyOrder, spaceDim);

    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    auto serendipityBasis = Teuchos::rcp( new SerendipityBasis<BasisBase>(fullBasis) );
    return serendipityBasis;
  }

  /** \brief  Factory method for isotropic HGRAD Serendipity bases on a hypercube for the given family.  Note that this will return a Line<2> for its base cell topology.  Note also that the family must use hierarchical bases.
      \param [in] polyOrder - the polynomial order of the basis.
      \param [in] spaceDim - the number of dimensions for the hypercube on which the basis is defined.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getSerendipityBasis_HVOL(int polyOrder, int spaceDim)
  {
    auto fullBasis = getHypercubeBasis_HVOL<BasisFamily>(polyOrder, spaceDim);

    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    auto serendipityBasis = Teuchos::rcp( new SerendipityBasis<BasisBase>(fullBasis) );
    return serendipityBasis;
  }
  
  /** \brief  Factory method for isotropic tetrahedron bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getTetrahedronBasis(Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_TET (polyOrder,pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_TET(polyOrder,pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_TET (polyOrder,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_TET(polyOrder,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic triangle bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getTriangleBasis(Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_TRI (polyOrder,pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_TRI(polyOrder,pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_TRI (polyOrder,pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_TRI(polyOrder,pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }

  /** \brief  Factory method for isotropic wedge bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getWedgeBasis(Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_WEDGE (polyOrder, pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_WEDGE(polyOrder, pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_WEDGE (polyOrder, pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_WEDGE(polyOrder, pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }

  /** \brief  Factory method for anisotropic wedge bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getWedgeBasis(Intrepid2::EFunctionSpace fs, ordinal_type polyOrder_xy, ordinal_type polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_WEDGE (polyOrder_xy, polyOrder_z, pointType));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_WEDGE(polyOrder_xy, polyOrder_z, pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_WEDGE (polyOrder_xy, polyOrder_z, pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_WEDGE(polyOrder_xy, polyOrder_z, pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }

  /** \brief  Factory method for pyramid bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getPyramidBasis(Intrepid2::EFunctionSpace fs, ordinal_type polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_PYR (polyOrder, pointType));
//      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_PYR(polyOrder, pointType));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_PYR (polyOrder, pointType));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_PYR(polyOrder, pointType));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic bases in the given family on the specified cell topology.
      \param [in] cellTopo    - the cell topology on which the basis is defined.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
      \param [in] pointType   - type of lattice used for creating the DoF coordinates.
   
   At present, only hypercube topologies are supported.  Once basis families support other element types, this method can
   be updated so that it also supports other element types.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getBasis(const shards::CellTopology &cellTopo, Intrepid2::EFunctionSpace fs, int polyOrder, const EPointType pointType = POINTTYPE_DEFAULT)
  {
    using Teuchos::rcp;
    switch (cellTopo.getBaseKey())
    {
      case shards::Line<>::key:          return getLineBasis<BasisFamily>(fs,polyOrder, pointType);
      case shards::Quadrilateral<>::key: return getQuadrilateralBasis<BasisFamily>(fs,polyOrder,pointType);
      case shards::Triangle<>::key:      return getTriangleBasis<BasisFamily>(fs,polyOrder,pointType);
      case shards::Hexahedron<>::key:    return getHexahedronBasis<BasisFamily>(fs,polyOrder,pointType);
      case shards::Tetrahedron<>::key:   return getTetrahedronBasis<BasisFamily>(fs,polyOrder,pointType);
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported cell topology");
    }
  }
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasisFamily_h */
