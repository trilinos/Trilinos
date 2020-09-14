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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
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
  template<class LineBasisHGRAD, class LineBasisHVOL, class TriangleBasisFamily = EmptyBasisFamily, class TetrahedronBasisFamily = EmptyBasisFamily>
  class DerivedBasisFamily
  {
  public:
    using ExecutionSpace  = typename LineBasisHGRAD::ExecutionSpace;
    using OutputValueType = typename LineBasisHGRAD::OutputValueType;
    using PointValueType  = typename LineBasisHGRAD::PointValueType;
    
    using Basis    = ::Intrepid2::Basis<ExecutionSpace,OutputValueType,PointValueType>;
    using BasisPtr = Teuchos::RCP<Basis>;
    
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
  };
  
  /** \brief  Factory method for line bases in the given family.
      \param [in] fs        - the function space for the basis.
      \param [in] polyOrder - the polynomial order of the basis.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getLineBasis(Intrepid2::EFunctionSpace fs, int polyOrder)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_LINE (polyOrder));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_LINE(polyOrder));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic quadrilateral bases in the given family.
      \param [in] fs        - the function space for the basis.
      \param [in] polyOrder - the polynomial order of the basis.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getQuadrilateralBasis(Intrepid2::EFunctionSpace fs, int polyOrder)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_QUAD (polyOrder));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_QUAD(polyOrder));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_QUAD (polyOrder));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_QUAD(polyOrder));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for potentially anisotropic quadrilateral bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder_x - the polynomial order of the basis in the x dimension.
      \param [in] polyOrder_y - the polynomial order of the basis in the y dimension.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getQuadrilateralBasis(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_QUAD (polyOrder_x,polyOrder_y));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_QUAD(polyOrder_x,polyOrder_y));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_QUAD (polyOrder_x,polyOrder_y));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_QUAD(polyOrder_x,polyOrder_y));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic bases on the hexahedron in the given family.
      \param [in] fs        - the function space for the basis.
      \param [in] polyOrder - the polynomial order of the basis.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getHexahedronBasis(Intrepid2::EFunctionSpace fs, int polyOrder)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_HEX (polyOrder));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_HEX(polyOrder));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_HEX (polyOrder));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_HEX(polyOrder));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for potentially anisotropic hexahedron bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder_x - the polynomial order of the basis in the x dimension.
      \param [in] polyOrder_y - the polynomial order of the basis in the y dimension.
      \param [in] polyOrder_z - the polynomial order of the basis in the z dimension.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getHexahedronBasis(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, int polyOrder_z)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_HEX (polyOrder_x,polyOrder_y,polyOrder_z));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_HEX(polyOrder_x,polyOrder_y,polyOrder_z));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_HEX (polyOrder_x,polyOrder_y,polyOrder_z));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_HEX(polyOrder_x,polyOrder_y,polyOrder_z));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic tetrahedron bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getTetrahedronBasis(Intrepid2::EFunctionSpace fs, int polyOrder)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      //Note: only HGRAD is available for Hierarchical basis at the moment
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_TET (polyOrder));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_TET(polyOrder));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_TET (polyOrder));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_TET(polyOrder));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic triangle bases in the given family.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getTriangleBasis(Intrepid2::EFunctionSpace fs, int polyOrder)
  {
    using Teuchos::rcp;
    switch (fs)
    {
      //Note: only HGRAD is available for Hierarchical basis at the moment
      case FUNCTION_SPACE_HVOL:  return rcp(new typename BasisFamily::HVOL_TRI (polyOrder));
      case FUNCTION_SPACE_HCURL: return rcp(new typename BasisFamily::HCURL_TRI(polyOrder));
      case FUNCTION_SPACE_HDIV:  return rcp(new typename BasisFamily::HDIV_TRI (polyOrder));
      case FUNCTION_SPACE_HGRAD: return rcp(new typename BasisFamily::HGRAD_TRI(polyOrder));
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported function space");
    }
  }
  
  /** \brief  Factory method for isotropic bases in the given family on the specified cell topology.
      \param [in] cellTopo    - the cell topology on which the basis is defined.
      \param [in] fs          - the function space for the basis.
      \param [in] polyOrder   - the polynomial order of the basis.
   
   At present, only hypercube topologies are supported.  Once basis families support other element types, this method can
   be updated so that it also supports other element types.
     */
  template<class BasisFamily>
  static typename BasisFamily::BasisPtr getBasis(shards::CellTopology &cellTopo, Intrepid2::EFunctionSpace fs, int polyOrder)
  {
    using Teuchos::rcp;
    switch (cellTopo.getBaseKey())
    {
      case shards::Line<>::key:          return getLineBasis<BasisFamily>(fs,polyOrder);
      case shards::Quadrilateral<>::key: return getQuadrilateralBasis<BasisFamily>(fs,polyOrder);
      case shards::Triangle<>::key:      return getTriangleBasis<BasisFamily>(fs,polyOrder);
      case shards::Hexahedron<>::key:    return getHexahedronBasis<BasisFamily>(fs,polyOrder);
      case shards::Tetrahedron<>::key:   return getTetrahedronBasis<BasisFamily>(fs,polyOrder);
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported cell topology");
    }
  }
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasisFamily_h */
