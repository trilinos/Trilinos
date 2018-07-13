#pragma once
#ifndef __TEST_SUITE_HPP__
#define __TEST_SUITE_HPP__

/// \file test_suite.hpp
/// \brief Simple test suite
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "test_crs_matrix_base.hpp"
#include "test_crs_matrix_base_io.hpp"
#include "test_crs_matrix_view.hpp"
#include "test_crs_row_view.hpp"

#include "test_chol_unblocked.hpp"
//#include "test_chol_blocked.hpp"

#include "test_dense_matrix_base.hpp"
#include "test_dense_matrix_view.hpp"

#include "test_tri_solve_unblocked.hpp"
//#include "test_tri_solve_blocked.hpp"

#include "test_crs_hier_base.hpp"
#include "test_crs_task_view.hpp"

#include "test_chol_by_blocks.hpp"

#include "test_dense_hier_base.hpp"
#include "test_dense_task_view.hpp"

#include "test_tri_solve_by_blocks.hpp"

#include "test_chol_tri_solve_by_blocks.hpp"

namespace Tacho { 
  
  using namespace std;
  
  template<typename VT, typename OT, typename ST = OT,
           typename SpT = void, typename MeT = void>
  class TestSuite : public Disp {
  public:
    static string label;

    static int doUnitTests() {
      int r_val = 0;

      cout << label << "::doUnitTests::Begin" << endl;
      {
        const OT blk_cnt = 6, blks[6] = { 1, 2, 4, 8, 12, 16 };
        const OT nrhs_cnt = 6, nrhs[6] = { 1, 2, 4, 8, 12, 16 };
        // ============================================================
        r_val += testCrsMatrixBase<VT,OT,ST,SpT,MeT>(0,0);
        r_val += testCrsMatrixBase<VT,OT,ST,SpT,MeT>(3,3);
        // ============================================================
        r_val += testCrsMatrixBaseIO<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx");
        // ============================================================
        r_val += testCrsMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,2,2);
        r_val += testCrsMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,0,2,2);
        r_val += testCrsMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  0,1,2,2);
        r_val += testCrsMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,1,0);
        r_val += testCrsMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,0,1);
        // ============================================================
        r_val += testCrsRowView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,2,2);
        r_val += testCrsRowView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,0,2,2);
        r_val += testCrsRowView<VT,OT,ST,SpT,MeT>(3,3, /**/  0,1,2,2);
        r_val += testCrsRowView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,1,0);
        r_val += testCrsRowView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,0,1);
        // ============================================================
        r_val += testCholUnblocked<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
                                                      "mm_crs_chol.mtx");
        // ============================================================
        // for (OT i=0;i<blk_cnt;++i) 
        //   r_val += testCholBlocked<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
        //                                               blks[i],
        //                                               "mm_crs_chol.mtx");
        // ============================================================
        r_val += testDenseMatrixBase<VT,OT,ST,SpT,MeT>(0,0);
        r_val += testDenseMatrixBase<VT,OT,ST,SpT,MeT>(3,3);
        // ============================================================
        r_val += testDenseMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,2,2);
        r_val += testDenseMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,0,2,2);
        r_val += testDenseMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  0,1,2,2);
        r_val += testDenseMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,1,0);
        r_val += testDenseMatrixView<VT,OT,ST,SpT,MeT>(3,3, /**/  1,1,0,1);
        // ============================================================
        for (OT i=0;i<nrhs_cnt;++i) 
          r_val += testTriSolveUnblocked<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
                                                           nrhs[i]);
        // ============================================================      
        // for (OT j=0;j<blk_cnt;++j) 
        //   for (OT i=0;i<nrhs_cnt;++i) 
        //     r_val += testTriSolveBlocked<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
        //                                                    blks[j], nrhs[i]);
        // ============================================================ 
        r_val += testCrsHierBase<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx");
        r_val += testCrsTaskView<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx");
        // ============================================================ 
        r_val += testCholByBlocks<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx");
        // ============================================================ 
        for (OT j=0;j<blk_cnt;++j) 
          for (OT i=0;i<nrhs_cnt;++i) 
            r_val += testDenseHierBase<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
                                                         blks[j], nrhs[i]);
        for (OT j=0;j<blk_cnt;++j) 
          for (OT i=0;i<nrhs_cnt;++i) 
            r_val += testDenseTaskView<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
                                                         blks[j], nrhs[i]);
        // ============================================================ 
        for (OT j=0;j<blk_cnt;++j) 
          for (OT i=0;i<nrhs_cnt;++i) 
            r_val += testTriSolveByBlocks<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
                                                            blks[j], nrhs[i]);
        // ============================================================ 
        for (OT j=0;j<blk_cnt;++j) 
          for (OT i=0;i<nrhs_cnt;++i) 
            r_val += testCholTriSolveByBlocks<VT,OT,ST,SpT,MeT>("mm_crs_input.mtx",
                                                                 blks[j], nrhs[i]);
        // ============================================================ 
      }
      cout << label << "::doUnitTests::End" << endl;

      string eval;
      __EVAL_STRING__(r_val, eval);

      __DOT_LINE__;
      cout << label << "::doUnitTests::Eval::" << eval << endl;
      return r_val;
    }
  };

  template<typename VT, typename OT, typename ST, 
           typename SpT, typename MeT> 
  string TestSuite<VT, OT, ST,
                   SpT, MeT>::label = "NoName";
}

#endif
