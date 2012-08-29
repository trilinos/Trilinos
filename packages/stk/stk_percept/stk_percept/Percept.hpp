#ifndef stk_percept_Percept_hpp
#define stk_percept_Percept_hpp

//#define HAVE_INTREPID_DEBUG 0

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- PGI Compiler bug workaround (switch statements in template code)
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

#define USE_PGI_7_1_COMPILER_BUG_WORKAROUND

#if !defined(PGI_INSTANTIATION_FILE) && defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND)
// ensure bases get included only once
#define INTREPID_HGRAD_HEX_C1_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_FEMDEF_HPP

#define INTREPID_HGRAD_LINE_C1_FEMDEF_HPP

#define INTREPID_HGRAD_TRI_C1_FEMDEF_HPP
#define INTREPID_HGRAD_TRI_C2_FEMDEF_HPP

#define INTREPID_HGRAD_QUAD_C1_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C2_FEMDEF_HPP

#define INTREPID_HGRAD_HEX_C1_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_FEMDEF_HPP

#define INTREPID_HGRAD_TET_C1_FEMDEF_HPP
#define INTREPID_HGRAD_TET_C2_FEMDEF_HPP

#define INTREPID_HGRAD_WEDGE_C1_FEMDEF_HPP

#define INTREPID_HGRAD_QUAD_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_WEDGE_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C2_SERENDIPITY_FEMDEF_HPP

      // Shells
#define INTREPID_HGRAD_TRI_C1_FEMDEF_HPP
#define INTREPID_HGRAD_TRI_C2_FEMDEF_HPP

#define INTREPID_HGRAD_QUAD_C1_FEMDEF_HPP


#define INTREPID_CELLTOOLSDEF_HPP

#endif


//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- PGI Compiler old-ness problem - unsupported on boost::ublas
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

/* The following is to avoid the following error from boost::ublas libraries:

"/scratch/srkenno/code/TPLs_src/boost/boost/numeric/ublas/detail/config.hpp", line 170: catastrophic error: 
          #error directive: Your compiler and/or configuration is unsupported
          by this verions of uBLAS. Define BOOST_UBLAS_UNSUPPORTED_COMPILER=0
          to override this message. Boost 1.32.0 includes uBLAS with support
          for many older compilers.
  #error Your compiler and/or configuration is unsupported by this verions of uBLAS. Define BOOST_UBLAS_UNSUPPORTED_COMPILER=0 to 
  override this message. Boost 1.32.0 includes uBLAS with support for many older compilers.
*/

#if defined(__PGI) 
#define BOOST_UBLAS_UNSUPPORTED_COMPILER 0
#endif


//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- Geometry and Mesquite configuration
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

#if defined(STK_BUILT_IN_SIERRA) && !defined(STK_PERCEPT_HAS_GEOMETRY)
#define STK_PERCEPT_HAS_GEOMETRY
#endif

#if defined(STK_BUILT_IN_SIERRA) && !defined(STK_PERCEPT_HAS_MESQUITE)
#define STK_PERCEPT_HAS_MESQUITE
#endif

#if defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)
#undef STK_PERCEPT_HAS_MESQUITE
#endif

#include <stk_percept/ExceptionWatch.hpp>


#endif
