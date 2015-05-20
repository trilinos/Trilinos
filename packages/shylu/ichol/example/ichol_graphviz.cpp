#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include "util.hpp"
#include "graph_helper_scotch.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "sequential_for.hpp"
#include "parallel_for.hpp"

#include "task_policy_graphviz.hpp"
#include "task_factory.hpp"
#include "team_factory.hpp"                                                                                   
#include "task_team_factory.hpp"    

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space 
typedef Kokkos::Serial space_type;
                                                                                                              
using namespace Example;

// exec space
typedef space_type ExecSpace;

// flat matrix 
typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;                        
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;  

// scotch reordering
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;                                                

typedef TaskTeamFactory<TaskPolicy,Future,TeamThreadLoopRegion> TaskFactoryType;                              
typedef SequentialFor ForType;                                                                                

// block representation for CrsMatrix                                                                         
typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;                                       
                                                                                                              
// hier matrix                                                                                                
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierBaseType;                     
typedef CrsMatrixView<CrsHierBaseType> CrsHierViewType;                     

typedef TaskView<CrsHierViewType,TaskFactoryType> CrsHierTaskType;   

// ---------------------------------------------------------------------------------
int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }

  // initialization
  // --------------------------------------
  ExecSpace::initialize();                                                                            
  cout << "default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  // import a matrix
  // --------------------------------------
  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[2] << endl;
    return -1;
  }

  CrsMatrixBaseType AA("AA");
  AA.importMatrixMarket(in);

  // reorder a matrix using Scotch
  // --------------------------------------
  GraphHelperType S(AA);
  S.computeOrdering();

  cout << "# of rows = " << S.NumRows() << ", # of blocks " << S.NumBlocks() 
       << endl;

  // permute a matrix with a given ordering
  // --------------------------------------
  CrsMatrixBaseType PA("Permuted A");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);

  // run ichol by blocks
  // --------------------------------------
  CrsMatrixBaseType RR("RR");

  RR.copy(Uplo::Upper, PA);
  {
    CrsHierBaseType HH("HH");
    CrsMatrixHelper::flat2hier(Uplo::Upper, RR, HH,
                               S.NumBlocks(),
                               S.RangeVector(),
                               S.TreeVector());    

    CrsHierTaskType H(&HH);  
    
    int r_val = 0;                                                                                            
    typedef typename CrsTaskViewType::policy_type policy_type;                                                
                                                                                                              
    IChol<Uplo::Upper,AlgoIChol::ByBlocks>::                                                             
      TaskFunctor<ForType,CrsHierTaskType>(H).apply(policy_type::member_null(), r_val);   
    
    ofstream out;
    out.open("graph.gv");
    if (!out.good()) {
      cout << "Error in open the file: graph.gv" << endl;
      return -1;
    }
    
    policy_type policy;
    policy.graphviz(out);
  }

  // finalization
  // --------------------------------------  
  ExecSpace::finalize();
  
  return 0;
}
