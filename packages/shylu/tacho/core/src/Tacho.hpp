#ifndef __TACHO_HPP__
#define __TACHO_HPP__

/// \file Tacho.hpp
/// \brief Main header file
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLUTacho_config.h"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_Partition.hpp"

#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_DenseMatrixView.hpp"

#include "TachoExp_MatrixMarket.hpp"           

#include "TachoExp_Graph.hpp"
#include "TachoExp_GraphTools_CAMD.hpp"        
#include "TachoExp_GraphTools_Metis.hpp"       
#include "TachoExp_GraphTools_Scotch.hpp"      

#include "TachoExp_SupernodeInfo.hpp"
#include "TachoExp_SymbolicTools.hpp"

#include "TachoExp_Chol.hpp"
#include "TachoExp_Chol_External.hpp"
#include "TachoExp_Chol_ByBlocks.hpp"

#include "TachoExp_Trsm.hpp"
#include "TachoExp_Trsm_External.hpp"
#include "TachoExp_Trsm_ByBlocks.hpp"

#include "TachoExp_Herk.hpp"
#include "TachoExp_Herk_External.hpp"
#include "TachoExp_Herk_ByBlocks.hpp"          

#include "TachoExp_Gemm.hpp"
#include "TachoExp_Gemm_External.hpp"
#include "TachoExp_Gemm_ByBlocks.hpp"

#include "TachoExp_Trsv.hpp"
#include "TachoExp_Trsv_External.hpp"

#include "TachoExp_Gemv.hpp"
#include "TachoExp_Gemv_External.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

#include "TachoExp_TaskFunctor_FactorizeChol.hpp"
#include "TachoExp_TaskFunctor_FactorizeCholByBlocks.hpp"

#include "TachoExp_TaskFunctor_SolveLowerChol.hpp"
#include "TachoExp_TaskFunctor_SolveUpperChol.hpp"

#include "TachoExp_NumericTools.hpp"

// Do not include this. 
// In a gcc (4.9.x), this causes some multiple definition link error with gcc headers.
// No idea yet why it happens as the code is guarded by Tacho::Experimental namespace.
//#include "TachoExp_CommandLineParser.hpp" 

#endif
