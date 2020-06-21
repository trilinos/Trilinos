// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_Typedefs.hpp
 *  \brief keep typedefs that commonly appear in many places localized
 */

#ifndef ZOLTAN2_TYPEDEFS
#define ZOLTAN2_TYPEDEFS

#include "Zoltan2_TestHelpers.hpp"
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Zoltan2_EvaluateBaseClass.hpp>
//#include <Tpetra_Map.hpp>
//#include <Xpetra_Vector_decl.hpp>
//#include <Xpetra_CrsMatrix_decl.hpp>
//#include <Xpetra_CrsGraph_decl.hpp>
//
////#ifdef HAVE_ZOLTAN2_GALERI
////#include <Galeri_XpetraProblemFactory.hpp>
////#include <Galeri_XpetraParameters.hpp>
////#endif
//
////#include <Kokkos_DefaultNode.hpp>
//#include "GeometricGenerator.hpp"

// Forward declaration for classes in the Tpetra namespace
// namespace Tpetra {
// 
//   template<typename T1, typename T2, typename T3, typename T4, const bool T5>
//   class CrsMatrix;  
// 
//   template<typename T1, typename T2, typename T3, const bool T4>
//   class CrsGraph;  
// 
//   template<typename T1, typename T2, typename T3, typename T4, const bool T5>
//   class Vector; 
// 
//   template<typename T1, typename T2, typename T3, typename T4, const bool T5>
//   class MultiVector;
// }
// 
// // Forward declaration for classes in the Xpetra namespace
// namespace Xpetra {
// 
//   template<typename T1, typename T2, typename T3, typename T4>
//   class CrsMatrix;  
// 
//   template<typename T1, typename T2, typename T3>
//   class CrsGraph;  
// 
//   template<typename T1, typename T2, typename T3, typename T4>
//   class Vector; 
// 
//   template<typename T1, typename T2, typename T3, typename T4>
//   class MultiVector;
// }
 
 // Forward declaration for classes in the GeometricGen namespace
 namespace GeometricGen {
   template<typename T1, typename T2, typename T3, typename T4>
   class GeometricGenerator;
 }
 
// Forward declaration of classes in the Zoltan2 namespace
 namespace Zoltan2 {
   template<typename T1, typename T2, typename T3, typename T4>
   class BasicUserTypes;
 
   template<typename T1>
   class BaseAdapter;
 
   template<typename T1>
   class BasicIdentifierAdapter;
 
   template<typename T1>
   class XpetraMultiVectorAdapter;
 
   template<typename T1, typename T2>
   class XpetraCrsGraphAdapter;
 
   template<typename T1, typename T2>
   class XpetraCrsMatrixAdapter;
 
   template<typename T1>
   class BasicVectorAdapter; 
 
 #ifdef HAVE_ZOLTAN2_PAMGEN
   template<typename T1>
   class PamgenMeshAdapter;
 #endif
 
   template<typename T1>
   class Problem;
   
   template<typename T1>
   class PartitioningProblem;
   
   template<typename T1>
   class OrderingProblem;
   
   template<typename T1>
   class ColoringProblem;
 }

namespace Zoltan2_TestingFramework {

  // Data types
  typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>     tcrsMatrix_t;
  typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t>                 tcrsGraph_t;
  typedef Tpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t>        tVector_t;
  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t>   tMVector_t;

  typedef Xpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>     xcrsMatrix_t;
  typedef Xpetra::CrsGraph<zlno_t, zgno_t, znode_t>                 xcrsGraph_t;
  typedef Xpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t>        xVector_t;
  typedef Xpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t>   xMVector_t;

  typedef GeometricGen::GeometricGenerator<zscalar_t, zlno_t, zgno_t, znode_t>
  geometricgen_t;

  // Adapter types 
  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t>    userTypes_t;
  typedef Zoltan2::BaseAdapter<userTypes_t>                     base_adapter_t;
  typedef Zoltan2::BasicIdentifierAdapter<userTypes_t>          basic_id_t;

  typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t>              xMV_tMV_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tMVector_t>    xCG_tCG_t;
  typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t, tMVector_t>  xCM_tCM_t;

  typedef Zoltan2::XpetraMultiVectorAdapter<xMVector_t>              xMV_xMV_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<xcrsGraph_t, tMVector_t>    xCG_xCG_t;
  typedef Zoltan2::XpetraCrsMatrixAdapter<xcrsMatrix_t, tMVector_t>  xCM_xCM_t;

#ifdef HAVE_EPETRA_DATA_TYPES
  typedef Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector>         xMV_eMV_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<Epetra_CrsGraph, tMVector_t>   xCG_eCG_t;
  typedef Zoltan2::XpetraCrsMatrixAdapter<Epetra_CrsMatrix, tMVector_t> xCM_eCM_t;
#else // temp compiler issues - dummy place holders
  typedef Zoltan2::BasicVectorAdapter<tMVector_t> xMV_eMV_t;
  typedef Zoltan2::BasicVectorAdapter<tMVector_t> xCG_eCG_t;
  typedef Zoltan2::BasicVectorAdapter<tMVector_t> xCM_eCM_t;
#endif

  typedef Zoltan2::BasicVectorAdapter<tMVector_t>         basic_vector_adapter;

#ifdef HAVE_ZOLTAN2_PAMGEN
  typedef Zoltan2::PamgenMeshAdapter<userTypes_t>          pamgen_adapter_t;
#else
  // This typedef exists only to satisfy the compiler.
  // PamgenMeshAdapter cannot be used when Trilinos is not built with Pamgen
  typedef Zoltan2::BasicVectorAdapter<userTypes_t>         pamgen_adapter_t;
#endif

// when test objects are created they set an enum to define the adapter type
// then when used they can be cast to that type safely using Z2_TEST_UPCAST
enum EAdapterType {
  AT_basic_id_t,
  AT_xMV_tMV_t,
  AT_xMV_xMV_t,
  AT_xMV_eMV_t,
  AT_xCG_tCG_t,
  AT_xCG_xCG_t,
  AT_xCG_eCG_t,
  AT_xCM_tCM_t,
  AT_xCM_xCM_t,
  AT_xCM_eCM_t,
  AT_basic_vector_adapter,
  AT_pamgen_adapter_t
};

#define Z2_TEST_UPCAST(adptr, TEMPLATE_ACTION)                                 \
switch(adptr) {                                                                \
  case AT_basic_id_t: {TEMPLATE_ACTION(basic_id_t)} break;                     \
  case AT_xMV_tMV_t: {TEMPLATE_ACTION(xMV_tMV_t)} break;                       \
  case AT_xMV_xMV_t: {TEMPLATE_ACTION(xMV_xMV_t)} break;                       \
  case AT_xMV_eMV_t: {TEMPLATE_ACTION(xMV_eMV_t)} break;                       \
  case AT_xCG_tCG_t: {TEMPLATE_ACTION(xCG_tCG_t)} break;                       \
  case AT_xCG_xCG_t: {TEMPLATE_ACTION(xCG_xCG_t)} break;                       \
  case AT_xCG_eCG_t: {TEMPLATE_ACTION(xCG_eCG_t)} break;                       \
  case AT_xCM_tCM_t: {TEMPLATE_ACTION(xCM_tCM_t)} break;                       \
  case AT_xCM_xCM_t: {TEMPLATE_ACTION(xCM_xCM_t)} break;                       \
  case AT_xCM_eCM_t: {TEMPLATE_ACTION(xCM_eCM_t)} break;                       \
  case AT_basic_vector_adapter: {TEMPLATE_ACTION(basic_vector_adapter)} break; \
  case AT_pamgen_adapter_t: {TEMPLATE_ACTION(pamgen_adapter_t)} break;         \
  default: throw std::logic_error( "Bad Z2_TEST_UPCAST" );                     \
};

#define Z2_TEST_UPCAST_COORDS(adptr, TEMPLATE_ACTION)                          \
switch(adptr) {                                                                \
  case AT_xMV_tMV_t: {TEMPLATE_ACTION(xMV_tMV_t)} break;                       \
  default: throw std::logic_error( "Bad Z2_TEST_UPCAST_COORDINATES" );         \
};

} // namespace Zoltan2_TestingFramework

#endif // ZOLTAN2_TYPEDEFS
