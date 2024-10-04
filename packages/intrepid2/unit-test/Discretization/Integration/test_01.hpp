// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Unit test (CubatureDirect,CubatureTensor):
            correctness of volume computations for reference cells.
    \author Created by P. Bochev, D. Ridzal and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_CubatureTensor.hpp"

#include "Intrepid2_CubatureDirectLineGaussJacobi20.hpp"
#include "Intrepid2_CubatureDirectTriDefault.hpp"
#include "Intrepid2_CubatureDirectTriSymmetric.hpp"
#include "Intrepid2_CubatureDirectTetDefault.hpp"
#include "Intrepid2_CubatureDirectTetSymmetric.hpp"
//#include "Intrepid2_CubatureTensorPyr.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "test_util.hpp"

namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    } catch (std::logic_error &err) {                                    \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  \

    template<typename ValueType, typename DeviceType>
    int Integration_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      using DeviceSpaceType = typename DeviceType::execution_space;
      typedef typename
        Kokkos::DefaultHostExecutionSpace HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                  Unit Test (CubatureDirect,CubatureTensor)                  |\n"
        << "|                                                                             |\n"
        << "|     1) Computing volumes of reference cells                                 |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubatureDirectLineGauss        <DeviceType,pointValueType,weightValueType> CubatureLineType;
      typedef CubatureDirectTriDefault       <DeviceType,pointValueType,weightValueType> CubatureTriType;
      typedef CubatureDirectTriSymmetric     <DeviceType,pointValueType,weightValueType> CubatureTriSymType;
      typedef CubatureDirectTetDefault       <DeviceType,pointValueType,weightValueType> CubatureTetType;
      typedef CubatureDirectTetSymmetric     <DeviceType,pointValueType,weightValueType> CubatureTetSymType;
      typedef CubatureTensor                 <DeviceType,pointValueType,weightValueType> CubatureTensorType;
      //typedef CubatureTensorPyr              <DeviceType,pointValueType,weightValueType> CubatureTensorPyrType;

      const auto tol = 100.0 * tolerence();

      int errorFlag  = 0;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: exception                                                           |\n"
        << "===============================================================================\n";

      try {
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        *outStream << "-> Line testing\n\n";
        {
          INTREPID2_TEST_ERROR_EXPECTED( CubatureLineType(-1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureLineType(Parameters::MaxCubatureDegreeEdge+1) );
        }

        *outStream << "-> Triangle testing\n\n";
        {
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTriType triCub(-1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTriType triCub(Parameters::MaxCubatureDegreeTri+1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTriSymType triSymCub(-1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTriSymType triSymCub(Parameters::MaxCubatureDegreeTri+1) );
        }

        *outStream << "-> Tetrahedron testing\n\n";
        {
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTetType tetCub(-1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTetType tetCub(Parameters::MaxCubatureDegreeTet+1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTetSymType tetSymCub(-1) );
          INTREPID2_TEST_ERROR_EXPECTED( CubatureTetSymType tetSymCub(Parameters::MaxCubatureDegreeTet+1) );
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }

      } catch (std::logic_error &err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
      };


      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: basic functionalities                                               |\n"
        << "===============================================================================\n";

      try {

        *outStream << "-> Line testing\n\n";
        {
          CubatureLineType lineCub(4);
          INTREPID2_TEST_FOR_EXCEPTION( lineCub.getDimension() != 1, std::logic_error,
                                        ">>> ERROR (Integration::Test01): line cubature must have 1 dimension.");
          INTREPID2_TEST_FOR_EXCEPTION( lineCub.getAccuracy() != 4, std::logic_error,
                                        ">>> ERROR (Integration::Test01): line cubature reports wrong accuracy.");
        }

        *outStream << "-> Triangle testing\n\n";
        {
          CubatureTriType triCub(17);
          INTREPID2_TEST_FOR_EXCEPTION( triCub.getNumPoints() != 61, std::logic_error,
                                        ">>> ERROR (Integration::Test01): triangle cubature reports a wrong number of points.");
          INTREPID2_TEST_FOR_EXCEPTION( triCub.getDimension() != 2, std::logic_error,
                                        ">>> ERROR (Integration::Test01): triangle cubature reports a wrong dimension.");
          
          CubatureTriSymType triSymCub(17);
          INTREPID2_TEST_FOR_EXCEPTION( triSymCub.getNumPoints() != 60, std::logic_error,
                                        ">>> ERROR (Integration::Test01): triangle symmetric cubature reports a wrong number of points.");
          INTREPID2_TEST_FOR_EXCEPTION( triSymCub.getDimension() != 2, std::logic_error,
                                        ">>> ERROR (Integration::Test01): triangle symmetric cubature reports a wrong dimension.");
        }

        *outStream << "-> Tetrahedron testing\n\n";
        {
          CubatureTetType tetCub(17);
          INTREPID2_TEST_FOR_EXCEPTION( tetCub.getNumPoints() != 495, std::logic_error,
                                        ">>> ERROR (Integration::Test01): tetrahedron cubature reports a wrong number of points.");
          INTREPID2_TEST_FOR_EXCEPTION( tetCub.getDimension() != 3, std::logic_error,
                                        ">>> ERROR (Integration::Test01): tetrahedron cubature reports a wrong dimension.");
        }
        
        *outStream << "-> Tetrahedron testing\n\n";
        {
          CubatureTetSymType tetSymCub(17);
          INTREPID2_TEST_FOR_EXCEPTION( tetSymCub.getNumPoints() != 364, std::logic_error,
                                        ">>> ERROR (Integration::Test01): tetrahedron symmetric cubature reports a wrong number of points.");
          INTREPID2_TEST_FOR_EXCEPTION( tetSymCub.getDimension() != 3, std::logic_error,
                                        ">>> ERROR (Integration::Test01): tetrahedron symmetric cubature reports a wrong dimension.");
        }

        *outStream << "-> Quad testing\n\n";
        {
          CubatureTensorType quadCub( CubatureLineType(3), CubatureLineType(7) );

          INTREPID2_TEST_FOR_EXCEPTION( quadCub.getDimension() != 2, std::logic_error,
                                        ">>> ERROR (Integration::Test01): quad cubature must have 2 dimension.");

          ordinal_type accuracy[Parameters::MaxDimension];
          quadCub.getAccuracy( accuracy );
          INTREPID2_TEST_FOR_EXCEPTION( accuracy[0] != 3 || accuracy[1] != 7, std::logic_error,
                                        ">>> ERROR (Integration::Test01): quad cubature reports wrong accuracy.");

        }

        *outStream << "-> Hex testing\n\n";
        {
          CubatureTensorType hexCub( CubatureLineType(1), CubatureLineType(4), CubatureLineType(2) );

          INTREPID2_TEST_FOR_EXCEPTION( hexCub.getDimension() != 3, std::logic_error,
                                        ">>> ERROR (Integration::Test01): hex cubature must have 3 dimension.");

          ordinal_type accuracy[Parameters::MaxDimension];
          hexCub.getAccuracy( accuracy );
          INTREPID2_TEST_FOR_EXCEPTION( accuracy[0] != 1 || accuracy[1] != 4 || accuracy[2] != 2, std::logic_error,
                                        ">>> ERROR (Integration::Test01): hex cubature reports wrong accuracy.");
        }

        *outStream << "-> Prism testing\n\n";
        {
          CubatureTensorType prismCub( CubatureTriType(4), CubatureLineType(3) );

          INTREPID2_TEST_FOR_EXCEPTION( prismCub.getDimension() != 3, std::logic_error,
                                        ">>> ERROR (Integration::Test01): prism cubature must have 3 dimension.");

          ordinal_type accuracy[Parameters::MaxDimension];
          prismCub.getAccuracy( accuracy );
          INTREPID2_TEST_FOR_EXCEPTION( accuracy[0] != 4 || accuracy[1] != 3, std::logic_error,
                                        ">>> ERROR (Integration::Test01): prism cubature reports wrong accuracy.");


          CubatureTensorType prismSymCub( CubatureTriSymType(4), CubatureLineType(3) );
          prismSymCub.getAccuracy( accuracy );
          INTREPID2_TEST_FOR_EXCEPTION( accuracy[0] != 4 || accuracy[1] != 3, std::logic_error,
                                        ">>> ERROR (Integration::Test01): prism symmetric cubature reports wrong accuracy.");
        }

      } catch (std::logic_error &err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
      };

      *outStream
        << "===============================================================================\n"
        << "| TEST 3: volume computations                                                 |\n"
        << "===============================================================================\n";

      try {

        DynRankView ConstructWithLabel(cubPoints,  Parameters::MaxIntegrationPoints, Parameters::MaxDimension);
        DynRankView ConstructWithLabel(cubWeights, Parameters::MaxIntegrationPoints);

        int maxTotalCubatureDegree = Parameters::MaxCubatureDegreeEdge; // sum in all dimensions
        
        *outStream << "-> Line testing\n\n";
        {
          for (ordinal_type deg=0;deg<=maxTotalCubatureDegree;++deg) {
            CubatureLineType cub(deg);
            cub.getCubature(cubPoints, cubWeights);
            const auto npts = cub.getNumPoints();

            const auto testVol = computeRefVolume(npts, cubWeights);
            const auto refVol  = 2.0;
            if (std::abs(testVol - refVol) > tol) {
              *outStream << std::setw(30) << "Line volume --> " << std::setw(10) << std::scientific << testVol <<
                std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
              ++errorFlag;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            }
          }
        }

        *outStream << "-> Triangle testing\n\n";
        {
          for (auto deg=0;deg<=20;++deg) {
            CubatureTriType cub(deg);
            cub.getCubature(cubPoints, cubWeights);
            const auto npts = cub.getNumPoints();

            const auto testVol = computeRefVolume(npts, cubWeights);
            const auto refVol  = 0.5;
            if (std::abs(testVol - refVol) > tol) {
              *outStream << std::setw(30) << "Triangle volume computed with cubature of degree " << deg << " --> " << std::setw(10) << std::scientific << testVol <<
                std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
              ++errorFlag;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            }
          }

          *outStream << "-> Triangle symmetric cubature testing\n\n";
          for (auto deg=0;deg<=Parameters::MaxCubatureDegreeTri;++deg) {
            CubatureTriSymType cub(deg);
            cub.getCubature(cubPoints, cubWeights);
            const auto npts = cub.getNumPoints();

            const auto testVol = computeRefVolume(npts, cubWeights);
            const auto refVol  = 0.5;
            if (std::abs(testVol - refVol) > tol) {
              *outStream << std::setw(30) << "Triangle volume computed with symmetric cubature of degree " << deg << " --> " << std::setw(10) << std::scientific << testVol <<
                std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
              ++errorFlag;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            }
          }


        }

        *outStream << "-> Quad testing\n\n";
        {
          for (ordinal_type y_deg=0;y_deg<=maxTotalCubatureDegree-1;++y_deg)
            for (ordinal_type x_deg=0;x_deg<=maxTotalCubatureDegree-y_deg;++x_deg) {
              const auto x_line = CubatureLineType(x_deg);
              const auto y_line = CubatureLineType(y_deg);
              CubatureTensorType cub( x_line, y_line );

              cub.getCubature(cubPoints, cubWeights);
              const auto npts = cub.getNumPoints();

              const auto testVol = computeRefVolume(npts, cubWeights);
              const auto refVol  = 4.0;
              if (std::abs(testVol - refVol) > tol) {
                *outStream << std::setw(30) << "Quadrilateral volume computed with tensor-product cubature of degree (" << x_deg << ", " << y_deg << ") --> " << std::setw(10) << std::scientific << testVol <<
                  std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";

                ++errorFlag;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              }
            }
        }
        
        *outStream << "-> Tetrahedron testing\n\n";
        {
          for (auto deg=0;deg<=Parameters::MaxCubatureDegreeTet;++deg) {
            CubatureTetType cub(deg);
            
            cub.getCubature(cubPoints, cubWeights);
            const auto npts = cub.getNumPoints();
            
            const auto testVol = computeRefVolume(npts, cubWeights);
            const auto refVol  = 1.0/6.0;
            if (std::abs(testVol - refVol) > tol) {
              *outStream << std::setw(30) << "Tetrahedron volume computed with cubature of degree " << deg << " --> " << std::setw(10) << std::scientific << testVol <<
                std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
              
              ++errorFlag;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            }
          }

          *outStream << "-> Tetrahedron symmetric cubature testing\n\n";
          for (auto deg=0;deg<=Parameters::MaxCubatureDegreeTet;++deg) {
            CubatureTetSymType cub(deg);
            
            cub.getCubature(cubPoints, cubWeights);
            const auto npts = cub.getNumPoints();
            
            const auto testVol = computeRefVolume(npts, cubWeights);
            const auto refVol  = 1.0/6.0;
            if (std::abs(testVol - refVol) > tol) {
              *outStream << std::setw(30) << "Tetrahedron volume computed with symmetric cubature of degree " << deg << " --> " << std::setw(10) << std::scientific << testVol <<
                std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
              
              ++errorFlag;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            }
          }
        }
        
        *outStream << "-> Hexahedron testing\n\n";
        {
          // test up to 10th order in each dimension with standard (non-tensor) points; we'll test up to max degree with tensor points below.
          for (ordinal_type z_deg=0;z_deg<10;++z_deg)
            for (ordinal_type y_deg=0;y_deg<10;++y_deg)
              for (ordinal_type x_deg=0;x_deg<10;++x_deg) {
                const auto x_line = CubatureLineType(x_deg);
                const auto y_line = CubatureLineType(y_deg);
                const auto z_line = CubatureLineType(z_deg);
                CubatureTensorType cub( x_line, y_line, z_line );

                cub.getCubature(cubPoints, cubWeights);
                const auto npts = cub.getNumPoints();

                const auto testVol = computeRefVolume(npts, cubWeights);
                const auto refVol  = 8.0;
                if (std::abs(testVol - refVol) > tol) {
                  *outStream << std::setw(30) << "Hexahedron volume computed with tensor-product cubature of degree (" << x_deg << ", " << y_deg << ", " << z_deg << ") --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";

                  ++errorFlag;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                }
              }
        }
        
        {
          for (ordinal_type z_deg=0;z_deg<=maxTotalCubatureDegree-2;++z_deg)
            for (ordinal_type y_deg=0;y_deg<=maxTotalCubatureDegree-z_deg;++y_deg)
              for (ordinal_type x_deg=0;x_deg<=maxTotalCubatureDegree-z_deg-y_deg;++x_deg) {
                const auto x_line = CubatureLineType(x_deg);
                const auto y_line = CubatureLineType(y_deg);
                const auto z_line = CubatureLineType(z_deg);
                CubatureTensorType cub( x_line, y_line, z_line );

                auto cubTensorPoints  = cub.allocateCubaturePoints();
                auto cubTensorWeights = cub.allocateCubatureWeights();
                cub.getCubature(cubTensorPoints, cubTensorWeights);
                const auto npts = cub.getNumPoints();

                const auto testVol = computeRefVolume(npts, cubTensorWeights);
                const auto refVol  = 8.0;
                if (std::abs(testVol - refVol) > tol) {
                  *outStream << std::setw(30) << "Hexahedron volume computed with tensor-product cubature of degree (" << x_deg << ", " << y_deg << ", " << z_deg << ") --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";

                  ++errorFlag;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                }
              }
        }

        *outStream << "-> Prism testing\n\n";
        {
          for (auto z_deg=0;z_deg<Parameters::MaxCubatureDegreeEdge;++z_deg)
            for (auto xy_deg=0;xy_deg<20;++xy_deg) {
              const auto xy_tri = CubatureTriType(xy_deg);
              const auto z_line = CubatureLineType(z_deg);
              CubatureTensorType cub( xy_tri, z_line );
              
              cub.getCubature(cubPoints, cubWeights);
              const auto npts = cub.getNumPoints();
              
              const auto testVol = computeRefVolume(npts, cubWeights);
              const auto refVol  = 1.0;
              if (std::abs(testVol - refVol) > tol) {
                *outStream << std::setw(30) << "Wedge volume computed with tensor-product cubature of degrees (" << xy_deg << ", " << z_deg << ") --> " << std::setw(10) << std::scientific << testVol <<
                  std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
                ++errorFlag;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              }
            }

          *outStream << "-> Prism symmetric quadrature testing\n\n";
          // test up to 10th order in the extrusion dimension with standard (non-tensor) points; we'll test up to max degree with tensor points below.
          for (auto z_deg=0;z_deg<10;++z_deg)
            for (auto xy_deg=0;xy_deg<Parameters::MaxCubatureDegreeTri;++xy_deg) {
              const auto xy_tri = CubatureTriSymType(xy_deg);
              const auto z_line = CubatureLineType(z_deg);
              CubatureTensorType cub( xy_tri, z_line );
              
              cub.getCubature(cubPoints, cubWeights);
              const auto npts = cub.getNumPoints();
              
              const auto testVol = computeRefVolume(npts, cubWeights);
              const auto refVol  = 1.0;
              if (std::abs(testVol - refVol) > tol) {
                *outStream << std::setw(30) << "Wedge volume computed with symmetric tensor-product cubature of degree (" << xy_deg << ", " << z_deg << ") --> " << std::setw(10) << std::scientific << testVol <<
                  std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
                ++errorFlag;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              }
            }
          
          for (auto z_deg=0;z_deg<Parameters::MaxCubatureDegreeEdge;++z_deg)
            for (auto xy_deg=0;xy_deg<Parameters::MaxCubatureDegreeTri;++xy_deg) {
              const auto xy_tri = CubatureTriSymType(xy_deg);
              const auto z_line = CubatureLineType(z_deg);
              CubatureTensorType cub( xy_tri, z_line );
              
              auto cubTensorPoints  = cub.allocateCubaturePoints();
              auto cubTensorWeights = cub.allocateCubatureWeights();
              cub.getCubature(cubTensorPoints, cubTensorWeights);
              const auto npts = cub.getNumPoints();
              
              const auto testVol = computeRefVolume(npts, cubTensorWeights);
              const auto refVol  = 1.0;
              if (std::abs(testVol - refVol) > tol) {
                *outStream << std::setw(30) << "Wedge volume computed with symmetric tensor-product cubature of degree (" << xy_deg << ", " << z_deg << ") --> " << std::setw(10) << std::scientific << testVol <<
                  std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
                ++errorFlag;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              }
            }
        }
        
        // *outStream << "-> Pyramid testing: over-integration by 2 (due to duffy transformation) \n\n";
        // {
        //   for (auto deg=0;deg<=Parameters::MaxCubatureDegreePyr;++deg) {
        //     const auto xy_line = CubatureLineType(deg);
        //     const auto z_line  = CubatureLineJacobiType(deg);
        //     CubatureTensorPyrType cub( xy_line, xy_line, z_line );
        //     cub.getCubature(cubPoints, cubWeights);
        //     const auto npts = cub.getNumPoints();
            
        //     const auto testVol = computeRefVolume(npts, cubWeights);
        //     const auto refVol  = 4.0/3.0;
        //     if (std::abs(testVol - refVol) > tol) {              
        //       *outStream << std::setw(30) << "Pyramid volume --> " << std::setw(10) << std::scientific << testVol <<
        //         std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - refVol) << "\n";
              
        //       ++errorFlag;
        //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        //     }
        //   }
        // }

      }  catch (std::logic_error &err) {
        *outStream << err.what() << "\n";
        errorFlag = -1;
      };


      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";

      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      return errorFlag;
    }



  } // end of namespace test
} // end of namespace intrepid2
