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

#ifndef __Panzer_STK_SurfaceNormals_hpp__
#define __Panzer_STK_SurfaceNormals_hpp__

#include "Panzer_STK_config.hpp"
#include "Teuchos_RCP.hpp"
#include "Intrepid_FieldContainer.hpp"
#include <boost/unordered_map.hpp>
#include <string>
#include <vector>

namespace panzer_stk_classic {

  class STK_Interface;
  

  /** \brief Computes the normals for all nodes associated with a sideset surface

      Computes the node normals for a given side set.  This
      computation will use ghosted element contributions for all
      surfaces associated with the node, but will only use elements
      associated with the sideset.  So if two sidesets adjoin at the
      seam, the nodes on the seam will NOT have the same node normal.
      A simple future addition would be to allow for a list of
      sidesets to be provided to allow for the union computation.

      \param[out] normals Map of the node normals (including ghosted nodes).  Key is the global id from the stk mesh node entity id and value is vector of the size of the parent element dimension containing the normal vector components.  Vector components will be resized on calling this method.
      \param[in] mesh (Required) Panzer stk mesh 
      \param[in] sidesetName (Required) Name of the sideset that the normals will be computed on
      \param[in] elementBlockName (Required) Name of the element block that the outward facing normals will be computed on
      \param[in] out (Optional) The ostream used for serial debug output on print process only.  If non-null this will print debug info.
      \param[in] pout (Optional) The ostream used for parallel debug output by all processes.  If non-null this will print debug info.
  */
  void computeSidesetNodeNormals(boost::unordered_map<unsigned,std::vector<double> >& normals,
				 const Teuchos::RCP<const panzer_stk_classic::STK_Interface>& mesh,
				 const std::string& sidesetName,
				 const std::string& elementBlockName,
				 std::ostream* out = NULL,
				 std::ostream* pout = NULL);

  /** \brief Computes the normals for all nodes associated with a sideset surface

      Computes the node normals for a given side set.  This
      computation will use ghosted element contributions for all
      surfaces associated with the node, but will only use elements
      associated with the sideset.  So if two sidesets adjoin at the
      seam, the nodes on the seam will NOT have the same node normal.
      A simple future addition would be to allow for a list of
      sidesets to be provided to allow for the union computation.

      \param[out] normals Map of the node normals (including ghosted nodes).  Key is the panzer_stk_classic::STK_Interface local element id and value is a multidimensional array of the size of the number of nodes times the parent element dimension containing the normal vector components.  Vector components will be resized on calling this method.
      \param[in] mesh (Required) Panzer stk mesh 
      \param[in] sidesetName (Required) Name of the sideset that the normals will be computed on
      \param[in] elementBlockName (Required) Name of the element block that the outward facing normals will be computed on
      \param[in] out (Optional) The ostream used for serial debug output on print process only.  If non-null this will print debug info.
      \param[in] pout (Optional) The ostream used for parallel debug output by all processes.  If non-null this will print debug info.
  */
  void computeSidesetNodeNormals(boost::unordered_map<std::size_t,Intrepid::FieldContainer<double> >& elementToNormalMap,
				 const Teuchos::RCP<const panzer_stk_classic::STK_Interface>& mesh,
				 const std::string& sidesetName,
				 const std::string& elementBlockName,
				 std::ostream* out = NULL,
				 std::ostream* pout = NULL);


}

#endif
