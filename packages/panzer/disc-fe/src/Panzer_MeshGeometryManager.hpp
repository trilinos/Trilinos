// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_MESH_GEOMETRY_MANAGER_HPP
#define PANZER_MESH_GEOMETRY_MANAGER_HPP

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Panzer_BasisDescriptor.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Intrepid2_Basis.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_WEDGE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_PYR_C1_FEM.hpp"

#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TET_COMP12_FEM.hpp"

#include "Intrepid2_HGRAD_WEDGE_C2_FEM.hpp"

namespace panzer {
  
  /** \brief Holds information about the mesh geometry and a "mesh basis"
   * 
   * Recognizing that the mesh topology and geometry are really distinct, this
   * class contains information tied to the mesh geometry required in
   * FE reference-to-physical transformations.
   *  
  */ 

  class MeshGeometryManager {

  public:

    MeshGeometryManager(const Teuchos::RCP<const shards::CellTopology> & ct)
      : cell_topo_(ct)
    {
      // Create a "mesh basis", or the HGrad basis consistent with the mesh geometry
      // (as defined by the extended cell topology)

      mesh_basis_ = createHGradBasis<double,double>(*ct);
      basis_descriptor_ = createBasisDescriptor(*ct);

    }
 
    //! Get (extended) CellTopology
    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return cell_topo_; }

    //! Get HGrad basis associated with the mesh geometry
    Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > getMeshBasis() const
    { return mesh_basis_; }

    //! Get HGrad PureBasis associated with the mesh geometry
    Teuchos::RCP<const panzer::PureBasis> getMeshPureBasis(const int & num_cells) const 
    { return Teuchos::rcp(new panzer::PureBasis(basis_descriptor_,cell_topo_,num_cells)); }

    //! 
    
    // TODO BWR Also will want HCURL and HDIV bases consistent with the mesh

    // TODO BWR be careful about storing information that is housed in worksets
    // TODO BWR there is a lot of memory management happening and dont want lingering
    // TODO BWR RCPs open here

    // TODO BWR get set orts??? NO 
    // TODO BWR get set cell nodes??? NO REMEMBER WE ARE ON AN ELEMENT BLOCK

  private:
  
    //! (Extended) mesh topology
    Teuchos::RCP<const shards::CellTopology> cell_topo_;
    //! HGrad basis for the mesh geometry
    Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > mesh_basis_;
    //! BasisDescriptor for creating a PureBasis when needed 
    panzer::BasisDescriptor basis_descriptor_;

     /** \brief Generates default HGrad basis based on cell topology 
        \param cellTopo            [in] - cell topology 
  
        \note Taken from Intrepid2 for now.
     */

    // TODO BWR remove static? Prolly not necessary if just for this purpose
    // TODO BWR do this as pure basis?
    template<typename outputValueType, 
             typename pointValueType>
    static Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,outputValueType,pointValueType> >
    createHGradBasis( const shards::CellTopology & cellTopo ) {
      Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,outputValueType,pointValueType> > r_val;

      switch (cellTopo.getKey()) {
      case shards::Line<2>::key:          r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_LINE_C1_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Line<3>::key:          r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_LINE_C2_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Triangle<3>::key:      r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<4>::key: r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<4>::key:   r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<8>::key:    r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Wedge<6>::key:         r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_WEDGE_C1_FEM  <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Pyramid<5>::key:       r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_PYR_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;

      case shards::Triangle<6>::key:      r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<9>::key: r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<10>::key:  r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_C2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<11>::key:  r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_COMP12_FEM<PHX::exec_space,outputValueType,pointValueType>()); break;
        //case shards::Hexahedron<20>::key:   r_val = Teuchos::rcp(new Basis_HGRAD_HEX_I2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<27>::key:   r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
        //case shards::Wedge<15>::key:        r_val = Teuchos::rcp(new Basis_HGRAD_WEDGE_I2_FEM  <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Wedge<18>::key:        r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_WEDGE_C2_FEM  <PHX::exec_space,outputValueType,pointValueType>()); break;
        //case shards::Pyramid<13>::key:      r_val = Teuchos::rcp(new Basis_HGRAD_PYR_I2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;

      case shards::Quadrilateral<8>::key: 
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key: 
      default: {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                   ">>> ERROR (MeshGeometryManager::createHGradBasis): Cell topology not supported.");
      }
      }
      return r_val;
    }

    /* @brief Given the mesh cell topology, make a \ref panzer::BasisDescriptor 
    *         required to create a \ref panzer::PureBasis
    */
    panzer::BasisDescriptor createBasisDescriptor( const shards::CellTopology & cellTopo ) {

      int order = 1;

      switch (cellTopo.getKey()) {
      case shards::Line<2>::key:          order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_LINE_C1_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Line<3>::key:          order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_LINE_C2_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Triangle<3>::key:      order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<4>::key: order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<4>::key:   order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<8>::key:    order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Wedge<6>::key:         order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_WEDGE_C1_FEM  <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Pyramid<5>::key:       order = 1; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_PYR_C1_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;

      case shards::Triangle<6>::key:      order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<9>::key: order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM   <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<10>::key:  order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_C2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      //case shards::Tetrahedron<11>::key:  order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_TET_COMP12_FEM<PHX::exec_space,outputValueType,pointValueType>()); break;
        //case shards::Hexahedron<20>::key:   r_val = Teuchos::rcp(new Basis_HGRAD_HEX_I2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<27>::key:   order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;
        //case shards::Wedge<15>::key:        r_val = Teuchos::rcp(new Basis_HGRAD_WEDGE_I2_FEM  <PHX::exec_space,outputValueType,pointValueType>()); break;
      case shards::Wedge<18>::key:        order = 2; r_val = Teuchos::rcp(new Intrepid2::Basis_HGRAD_WEDGE_C2_FEM  <PHX::exec_space,outputValueType,pointValueType>()); break;
        //case shards::Pyramid<13>::key:      r_val = Teuchos::rcp(new Basis_HGRAD_PYR_I2_FEM    <PHX::exec_space,outputValueType,pointValueType>()); break;

      case shards::Quadrilateral<8>::key: 
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key: 
      default: {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                   ">>> ERROR MeshGeometryManager: Cell topology not supported.");
      }
      }

      return panzer::BasisDescriptor(order,"HGrad");
    }
  };

}

#endif
