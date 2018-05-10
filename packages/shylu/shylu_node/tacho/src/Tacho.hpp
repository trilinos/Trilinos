#ifndef __TACHO_HPP__
#define __TACHO_HPP__

/// \file Tacho.hpp
/// \brief Main header file
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLU_NodeTacho_config.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_Partition.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_DenseMatrixView.hpp"

#include "Tacho_MatrixMarket.hpp"           

#include "Tacho_Graph.hpp"
#include "Tacho_GraphTools_CAMD.hpp"        
#include "Tacho_GraphTools_Metis.hpp"       
//#include "Tacho_GraphTools_MetisMT.hpp"       
#include "Tacho_GraphTools_Scotch.hpp"      

#include "Tacho_SupernodeInfo.hpp"
#include "Tacho_SymbolicTools.hpp"

#include "Tacho_Lapack_External.hpp"
#include "Tacho_Lapack_Team.hpp"

#include "Tacho_Blas_External.hpp"
#include "Tacho_Blas_Team.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_External.hpp"
#include "Tacho_Chol_Internal.hpp"
#include "Tacho_Chol_ByBlocks.hpp"

#include "Tacho_Trsm.hpp"
#include "Tacho_Trsm_External.hpp"
#include "Tacho_Trsm_Internal.hpp"
#include "Tacho_Trsm_ByBlocks.hpp"

#include "Tacho_Herk.hpp"
#include "Tacho_Herk_External.hpp"
#include "Tacho_Herk_Internal.hpp"
#include "Tacho_Herk_ByBlocks.hpp"          

#include "Tacho_Gemm.hpp"
#include "Tacho_Gemm_External.hpp"
#include "Tacho_Gemm_Internal.hpp"
#include "Tacho_Gemm_ByBlocks.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_External.hpp"
#include "Tacho_Trsv_Internal.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_External.hpp"
#include "Tacho_Gemv_Internal.hpp"

#include "Tacho_CholSupernodes.hpp"
#include "Tacho_CholSupernodes_Serial.hpp"
#include "Tacho_CholSupernodes_SerialPanel.hpp"

#include "Tacho_TaskFunctor_FactorizeChol.hpp"
#include "Tacho_TaskFunctor_FactorizeCholPanel.hpp"
#include "Tacho_TaskFunctor_FactorizeCholByBlocks.hpp"
#include "Tacho_TaskFunctor_FactorizeCholByBlocksPanel.hpp"

#include "Tacho_TaskFunctor_SolveLowerChol.hpp"
#include "Tacho_TaskFunctor_SolveUpperChol.hpp"

#include "Tacho_NumericTools.hpp"

// Do not include this. 
// In a gcc (4.9.x), this causes some multiple definition link error with gcc headers.
// No idea yet why it happens as the code is guarded by Tacho::Experimental namespace.
//#include "Tacho_CommandLineParser.hpp" 

#endif
