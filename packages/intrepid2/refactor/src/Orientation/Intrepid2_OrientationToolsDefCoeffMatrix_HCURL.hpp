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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HCURL_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HCURL_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  namespace Impl {

    template<typename outputViewType,
             typename subcellBasisType,
             typename cellBasisType>
    inline
    void
    OrientationTools::
    getCoeffMatrix_HCURL(outputViewType &output,
                         const subcellBasisType subcellBasis,
                         const cellBasisType cellBasis,
                         const ordinal_type subcellId,
                         const ordinal_type subcellOrt) {
      typedef typename outputViewType::execution_space space_type;
      typedef typename outputViewType::value_type value_type;
      
      // with shards, everything should be computed on host space
      typedef typename
        Kokkos::Impl::is_space<space_type>::host_mirror_space::execution_space host_space_type;

      typedef Kokkos::DynRankView<value_type,host_space_type> DynRankViewHostType;

      ///
      /// Topology
      ///

      // populate points on a subcell and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();

      const ordinal_type cellDim = cellTopo.getDimension();
      const ordinal_type subcellDim = subcellTopo.getDimension();
      
      INTREPID2_TEST_FOR_EXCEPTION( subcellDim >= cellDim,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
                                    "cellDim must be greater than subcellDim.");

      const auto subcellBaseKey = subcellTopo.getBaseKey();
      const auto cellBaseKey = cellTopo.getBaseKey();

      INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key &&
                                    subcellBaseKey != shards::Quadrilateral<>::key &&
                                    subcellBaseKey != shards::Triangle<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
                                    "subcellBasis must have line, quad, or triangle topology.");

      INTREPID2_TEST_FOR_EXCEPTION( cellBaseKey != shards::Quadrilateral<>::key &&
                                    cellBaseKey != shards::Triangle<>::key &&
                                    cellBaseKey != shards::Hexahedron<>::key &&
                                    cellBaseKey != shards::Tetrahedron<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
                                    "cellBasis must have quad, triangle, hexhedron or tetrahedron topology.");

      // if node map has left handed system, orientation should be re-enumerated.
      ordinal_type ort = -1;
      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  2)
          ort = subcellOrt;
        break;
      }
      case shards::Triangle<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  6) {
          // in the basis of tet, it uses map to reference subcell and accounts for the left handed face
          //const ordinal_type leftHanded = cellTopo.getNodeMap(2, subcellId, 1) > cellTopo.getNodeMap(2, subcellId, 2);
          //const ordinal_type leftOrt[] = { 0, 2, 1, 3, 5, 4 };
          ort = subcellOrt; //(leftHanded ? leftOrt[subcellOrt] : subcellOrt);
        }
        break;
      }
      case shards::Quadrilateral<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  8) {
        //some faces require special treatment because the dofs of the face bases are not consistent with the corresponding dofs of the hexahedron
        //TODO: modify reference Hexahedron element so that the dofs of the subcell are consistent with those of the cell.
        switch (subcellId) {
            case 3:
            case 4: {
              const ordinal_type modifiedOrt[] = { 0, 3, 2, 1, 4, 7, 6, 5 }; //left hand orientation
              ort = modifiedOrt[subcellOrt];
              break;
            }
            case 2: {
              const ordinal_type modifiedOrt[] = { 0, 3, 2, 1, 6, 5, 4, 7 };
              ort = modifiedOrt[subcellOrt];
              break;
            }
            default:
              ort = subcellOrt;
          }
        }
        break;
      }
      }
      INTREPID2_TEST_FOR_EXCEPTION( ort == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
                                    "Orientation is not properly setup.");
      ///
      /// Scale is computed from Jacobian between reference subcell coord and oriented subcell coord
      ///
      double scale = -1.0;
      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        if        (cellBaseKey == shards::Quadrilateral<>::key || 
                   cellBaseKey == shards::Hexahedron<>::key) {
          scale = 1.0;
        } else if (cellBaseKey == shards::Triangle<>::key) {
          const double scale_edge[] = { 1.0, sqrt(2), 1.0 };
          scale = scale_edge[subcellId];
        } else if (cellBaseKey == shards::Tetrahedron<>::key) {
          const double scale_edge[] = { 1.0, sqrt(2.0), 1.0, 
                                        1.0, sqrt(2.0), sqrt(2.0) };
          scale = scale_edge[subcellId];
        } 
        break;
      }
      case shards::Triangle<>::key: {
        if (cellBaseKey == shards::Tetrahedron<>::key) {
          const double scale_face[] = { 1.0, sqrt(2.0), 1.0, 1.0 };
          scale = scale_face[subcellId];          
        }         
        break;
      }
      case shards::Quadrilateral<>::key: {
        if (cellBaseKey == shards::Hexahedron<>::key) {
          scale = 1.0;
        }         
        break;
      }
      }
      INTREPID2_TEST_FOR_EXCEPTION( scale < 0,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
                                    "Scale is not properly setup.");

      ///
      /// Function space
      ///      
      {
        const std::string cellBasisName(cellBasis.getName());
        if (cellBasisName.find("HCURL") != std::string::npos) {
          const std::string subcellBasisName(subcellBasis.getName());
          // edge hcurl is hgrad with gauss legendre points
          switch (subcellDim) {
          case 1: {
            //TODO: Hex, QUAD, TET and TRI element should have the same 1d basis
            if ((cellBasisName.find("HEX") != std::string::npos) || (cellBasisName.find("QUAD") != std::string::npos)) {
              INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HGRAD") == std::string::npos,
                                          std::logic_error,
                                          ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " 
                                          "subcellBasis function space (1d) is not consistent to cellBasis.");
            } else if ((cellBasisName.find("TET") != std::string::npos) || (cellBasisName.find("TRI") != std::string::npos)) {
              INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("L2") == std::string::npos,
                                          std::logic_error,
                                          ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " 
                                          "subcellBasis function space (1d) is not consistent to cellBasis.");
            }
            break;
          }
          case 2: {
            INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HCURL") == std::string::npos,
                                          std::logic_error,
                                          ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " 
                                          "subcellBasis function space (2d) is not consistent to cellBasis.");
            break;
          }
          }
        }
      }

      ///
      /// Collocation points
      ///     
      const ordinal_type degree = cellBasis.getDegree();

      const ordinal_type numCellBasis = cellBasis.getCardinality();
      const ordinal_type numSubcellBasis = subcellBasis.getCardinality();

      const ordinal_type ordSubcell = cellBasis.getDofOrdinal(subcellDim, subcellId, 0);
      INTREPID2_TEST_FOR_EXCEPTION( ordSubcell == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
                                    "Invalid subcellId returns -1 ordSubcell.");

      const ordinal_type ndofSubcell = cellBasis.getDofTag(ordSubcell)(3);

      // reference points on a subcell
      DynRankViewHostType refPtsSubcell;

      switch (subcellBaseKey) {
      case shards::Line<>::key:
      case shards::Triangle<>::key: {
        // hcurl ndof: p-1, the interior points of p+1 with offset 1: p-1.
        const ordinal_type ndofLine = PointTools::getLatticeSize(subcellTopo, degree+1, 1);
        refPtsSubcell = DynRankViewHostType("refPtsSubcell", ndofLine, subcellDim);
        PointTools::getLattice(refPtsSubcell,
                               subcellTopo,
                               degree+1,
                               1, // offset by 1 so the points are located inside
                               POINTTYPE_EQUISPACED);
        break;
      }
      case shards::Quadrilateral<>::key: {
        // hcurl ndof: (p-1)*p*2, tensor product of interior points of p+2 with offset 1: (p)*(p).

        // tensor product of lines
        const auto lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
        const ordinal_type ndofLine = PointTools::getLatticeSize(lineTopo, degree+2, 1);
        DynRankViewHostType refPtsLine("refPtsLine", ndofLine, 1);
        PointTools::getLattice(refPtsLine,
                               lineTopo,
                               degree+2,
                               1, // offset by 1 so the points are located inside
                               POINTTYPE_EQUISPACED);

        refPtsSubcell = DynRankViewHostType("refPtsSubcell", ndofLine*ndofLine, subcellDim);
        ordinal_type idx = 0;
        for (ordinal_type j=0;j<ndofLine;++j) { // y
          for (ordinal_type i=0;i<ndofLine;++i,++idx) { // x
            refPtsSubcell(idx, 0) = refPtsLine(i,0);
            refPtsSubcell(idx, 1) = refPtsLine(j,0);
          }
        }
        INTREPID2_TEST_FOR_EXCEPTION( idx < ndofSubcell,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                      "counted points on quad is less than ndofSubcell.");
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key ||
                                      subcellBaseKey != shards::Quadrilateral<>::key ||
                                      subcellBaseKey != shards::Triangle<>::key,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                      "subcellBasis must have line, quad, or triangle topology.");
      }
      }

      const ordinal_type nptsSubcell = refPtsSubcell.dimension(0);

      // modified points with orientation
      DynRankViewHostType ortPtsSubcell("ortPtsSubcell", nptsSubcell, subcellDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsSubcell,
                                                     refPtsSubcell,
                                                     subcellTopo,
                                                     ort);
      
      // map to reference coordinates
      DynRankViewHostType refPtsCell("refPtsCell", nptsSubcell, cellDim);
      CellTools<host_space_type>::mapToReferenceSubcell(refPtsCell,
                                                        refPtsSubcell,
                                                        subcellDim,
                                                        subcellId,
                                                        cellTopo);

      ///
      /// Basis evaluation on the collocation points
      ///

      // evaluate values on the reference cell
      DynRankViewHostType refValues("refValues", numCellBasis, nptsSubcell, cellDim);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);

      // evaluate values on the modified cell
      DynRankViewHostType outValues("outValues", numSubcellBasis, nptsSubcell, subcellDim);

      auto orient_values = [ndofSubcell,&subcellBasis,subcellDim,nptsSubcell,&outValues](const ordinal_type c[][2]) {
        for (ordinal_type i=0;i<ndofSubcell;++i) {
          const ordinal_type ii = subcellBasis.getDofOrdinal(subcellDim, 0, i);
          for (ordinal_type j=0;j<nptsSubcell;++j) {
            value_type tmp[2] = {};
            for (ordinal_type k=0;k<subcellDim;++k) 
              for (ordinal_type l=0;l<subcellDim;++l) 
                tmp[k] += c[k][l] * outValues(ii,j,l);
            for (ordinal_type k=0;k<subcellDim;++k) 
              outValues(ii,j,k) = tmp[k];
          }
        }
      };

      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        auto out = Kokkos::subview(outValues, Kokkos::ALL(), Kokkos::ALL(), 0);
        subcellBasis.getValues(out, ortPtsSubcell);

        // second dimension is dummy
        const ordinal_type c[2][1][2] = { { {  1, 0 } },     // 0
                                          { { -1, 0 } } };   // 1
        orient_values(c[ort]);
        break;
      }
      case shards::Triangle<>::key: {
        //Inverse of the gradinent of the orientation map
        subcellBasis.getValues(outValues, ortPtsSubcell);
        const ordinal_type c[6][2][2] = { { {  1,  0 },
                                            {  0,  1 } }, // 0
                                          { {  0,  1 },
                                            { -1, -1 } }, // 1
                                          { { -1, -1 },
                                            {  1,  0 } }, // 2
                                          { {  0,  1 },
                                            {  1,  0 } }, // 3
                                          { { -1, -1 },
                                            {  0,  1 } }, // 4
                                          { {  1,  0 },
                                            { -1, -1 } } }; // 5
        orient_values(c[ort]);
        break;
      }
      case shards::Quadrilateral<>::key: {
        subcellBasis.getValues(outValues, ortPtsSubcell);
        //Transpose of the gradient of the orientation map
        const ordinal_type c[8][2][2] = { { {  1,  0 },
                                            {  0,  1 } }, // 0
                                          { {  0,  1 },
                                            { -1,  0 } }, // 1
                                          { { -1,  0 },
                                            {  0, -1 } }, // 2
                                          { {  0, -1 },
                                            {  1,  0 } }, // 3
                                          { {  0,  1 },
                                            {  1,  0 } }, // 4
                                          { { -1,  0 },
                                            {  0,  1 } }, // 5
                                          { {  0, -1 },
                                            { -1,  0 } }, // 6
                                          { {  1,  0 },
                                            {  0, -1 } } }; // 7
        orient_values(c[ort]);
        break;
      }
      }
      
      ///
      /// Restrict vector valued basis functions on the subcell dimensions
      ///
      auto normalize = [&](DynRankViewHostType v) { 
        value_type norm = 0.0;
        const ordinal_type iend = v.dimension(0);
        for (ordinal_type i=0;i<iend;++i)
          norm += v(i)*v(i);
        norm = std::sqrt(norm);
        for (ordinal_type i=0;i<iend;++i)
          v(i) /= norm;
      };
      
      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        DynRankViewHostType edgeTan("edgeTan", cellDim);
        CellTools<host_space_type>::getReferenceEdgeTangent(edgeTan, subcellId, cellTopo);

        normalize(edgeTan);
        
        DynRankViewHostType tmpValues("tmpValues", numCellBasis, nptsSubcell, subcellDim);
        for (ordinal_type i=0;i<numCellBasis;++i)
          for (ordinal_type j=0;j<nptsSubcell;++j)
            for (ordinal_type k=0;k<cellDim;++k)
              tmpValues(i,j,0) += refValues(i,j,k)*edgeTan(k);
        refValues = tmpValues;
        break;
      }
      case shards::Quadrilateral<>::key: 
      case shards::Triangle<>::key: {
        DynRankViewHostType faceTanU("faceTanU", cellDim), faceTanV("faceTanV", cellDim);
        CellTools<host_space_type>::getReferenceFaceTangents(faceTanU, faceTanV,subcellId, cellTopo);

        normalize(faceTanU);
        normalize(faceTanV);
        
        DynRankViewHostType tmpValues("tmpValues", numCellBasis, nptsSubcell, subcellDim);
        for (ordinal_type i=0;i<numCellBasis;++i)
          for (ordinal_type j=0;j<nptsSubcell;++j)
            for (ordinal_type k=0;k<cellDim;++k) {
              tmpValues(i,j,0) += refValues(i,j,k)*faceTanU(k);
              tmpValues(i,j,1) += refValues(i,j,k)*faceTanV(k);
            }
        refValues = tmpValues;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, "Should not come here" );        
      }
      }
      
      ///
      /// Construct collocation matrix and solve problems
      ///

      // construct collocation matrix; using lapack, it should be left layout
      const ordinal_type dofDim = outValues.dimension(2);
      Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type>
        refMat("refMat", nptsSubcell*dofDim, ndofSubcell),
        ortMat("ortMat", nptsSubcell*dofDim, ndofSubcell),
        lwork("lwork", std::max(nptsSubcell*dofDim, 10), ndofSubcell);

      for (ordinal_type i=0;i<ndofSubcell;++i) {
        const ordinal_type iref = cellBasis.getDofOrdinal(subcellDim, subcellId, i);
        const ordinal_type iout = subcellBasis.getDofOrdinal(subcellDim, 0, i);

        for (ordinal_type j=0;j<nptsSubcell;++j) {
          const ordinal_type joff = j*dofDim;
          for (ordinal_type k=0;k<dofDim;++k) {
            refMat(joff+k,i) = refValues(iref,j,k)*scale;
            ortMat(joff+k,i) = outValues(iout,j,k);
          }
        }
      }

      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,value_type> lapack;
        ordinal_type info = 0;
        lapack.GELS('N', // trans
                    nptsSubcell*dofDim, // m
                    ndofSubcell, // n
                    ndofSubcell, // nrhs
                    refMat.data(), refMat.stride_1(), // A
                    ortMat.data(), ortMat.stride_1(), // B
                    lwork.data(), lwork.span(), // work space
                    &info);

        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
        }

        // clean up numerical noise            
        {
          const double eps = threshold();
          const ordinal_type 
            iend = ortMat.dimension(0),
            jend = ortMat.dimension(1);
          for (ordinal_type i=0;i<iend;++i)
            for (ordinal_type j=0;j<jend;++j)
              if (std::abs(ortMat(i,j)) < eps) ortMat(i,j) = 0;
        }
      }

      {
        // move the data to original device memory
        const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofSubcell);
        Kokkos::deep_copy(Kokkos::subview(output, range, range),
                          Kokkos::subview(ortMat, range, range));
      }
    }
  }

}
#endif
