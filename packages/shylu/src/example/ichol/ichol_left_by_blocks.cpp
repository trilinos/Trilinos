#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "symbolic_task.hpp"
#include "crs_task_view.hpp"

#include "ichol_left_by_blocks.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// flat matrix 
typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

// taskfy blocks
typedef Example::SymbolicTask Task;
typedef Example::CrsTaskView<CrsMatrixBase,Task> CrsTaskView;

// hier matrix
typedef Example::CrsMatrixBase<CrsTaskView,ordinal_type,size_type> CrsHierBase;
typedef Example::CrsTaskView<CrsHierBase,Task> CrsHierView;

typedef Example::Uplo Uplo;

int main (int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " filename" << " blksize" << endl;
    return -1;
  }

  CrsMatrixBase A;

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  A.importMatrixMarket(in);

  CrsMatrixBase L(A, Uplo::Lower);

  CrsHierBase H(L, 1, 1);
  CrsHierView HH(H);

  int r_val = Example::ichol_left_by_blocks_lower(HH);
  if (r_val != 0) 
    cout << " Error = " << r_val << endl;
  
  cout << HH << endl;
  
  Task::queue::showMe(cout);
  
  return 0;
}
