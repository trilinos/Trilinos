#include <iostream>

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"

using namespace std;

typedef double value_type;
typedef int ordinal_type;
typedef size_t size_type;

struct coo {
  ordinal_type i, j;
  value_type v;
  coo(const ordinal_type ii,
      const ordinal_type jj,
      const value_type vv) : i(ii), j(jj), v(vv) { }

  coo& operator=(const coo &y) {
    this->i = y.i;
    this->j = y.j;
    this->v = y.v;

    return *this;
  }

  bool operator<(const coo &y) const {
    ordinal_type r_val = (this->i - y.i);
    return (r_val == 0 ? this->j < y.j : r_val < 0);
  }
  
  bool operator==(const coo &y) const {
    return (this->i == y.i) && (this->j == y.j);
  }

  bool operator!=(const coo &y) const {
    return !(*this == y);
  }
};

typedef struct coo ijv_type;

int main(int argc, char *argv[]) {

  // ===================================================================================
  int nrank, irank;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &nrank);
  MPI_Comm_rank(comm, &irank);

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Hypre.\n");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int niter = 10;
  clp.setOption("niter", &niter, "Number of iterations for testing");

  clp.recogniseAllOptions(false);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );
  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  auto t_mm_read         = Teuchos::TimeMonitor::getNewTimer("mm::read");
  auto t_pcg_create      = Teuchos::TimeMonitor::getNewTimer("pcg::create");
  auto t_euclid_create   = Teuchos::TimeMonitor::getNewTimer("euclid::create");
  auto t_pcg_set_precond = Teuchos::TimeMonitor::getNewTimer("pcg::setprecond");
  auto t_pcg_setup       = Teuchos::TimeMonitor::getNewTimer("pcg::setup");
  auto t_pcg_solve       = Teuchos::TimeMonitor::getNewTimer("pcg::solve");
  auto t_euclid_destroy  = Teuchos::TimeMonitor::getNewTimer("euclid::destroy");
  auto t_pcg_destroy     = Teuchos::TimeMonitor::getNewTimer("pcg::destroy");
  
  // ===================================================================================
  ordinal_type m = 0, n = 0;
  size_type nnz = 0;

  HYPRE_IJMatrix ijv;
  HYPRE_IJVector vb, vx;

  HYPRE_ParCSRMatrix A;
  HYPRE_ParVector b, x;

  ordinal_type ilower, iupper, jlower, jupper;

  vector<ordinal_type> ncols_ijv, cols_ijv, rows_ijv, rows_vb, rows_vx;
  vector<value_type> vals_ijv, vals_vb, vals_vx;

  // read matrix market and construct hypre matrix ojbects
  {
    Teuchos::TimeMonitor evalTimerMonitor(*t_mm_read);

    ifstream in;
    in.open(file_input);
    if (!in.good()) {
      cout << "Failed in open the file: " << file_input << endl;
      return -1;
    }

    const ordinal_type base = 1;
    {
      string header;
      if (in.is_open()) {
        getline(in, header);
        while (in.good()) {
          char c = in.peek();
          if (c == '%' || c == '\n') {
            in.ignore(256, '\n');
            continue;
          }
          break;
        }
      }

      // check the header
      bool symmetry = (header.find("symmetric") != string::npos);

      // read matrix specification
      in >> m >> n >> nnz;

      // partition of matrix
      const ordinal_type range = m/nrank;

      ilower = range*irank;
      iupper = (irank == (nrank-1) ? (m-1) : range*(irank+1)-1);

      jlower = 0;
      jupper = (n-1);

      HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &ijv);
      HYPRE_IJMatrixSetObjectType(ijv, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(ijv);

      vector<ijv_type> tmp;
      for (size_type k=0;k<nnz;++k) {
        ordinal_type row, col;
        value_type val;
        in >> row >> col >> val;
        
        row -= base;
        col -= base;

        if ((row >= ilower && row <= iupper) && 
            (col >= jlower && col <= jupper)) 
          tmp.push_back(ijv_type(row, col, val));

        if ((symmetry && row != col) && 
            (col >= ilower && col <= iupper) && 
            (row >= jlower && row <= jupper)) 
          tmp.push_back(ijv_type(col, row, val));
      }
      sort(tmp.begin(), tmp.end(), less<ijv_type>());
      
      {
        const size_type ksize = tmp.size();
        size_type k = 0;

        for (ordinal_type icnt=0,i=ilower;i<=iupper;++i,++icnt) {
          rows_ijv.push_back(i);
          size_type jcnt = 0;
          for ( ;k<ksize;++k) {
            if (tmp[k].i != i) {
              break;
            } else {
              cols_ijv.push_back(tmp[k].j);
              vals_ijv.push_back(tmp[k].v);
              ++jcnt;
            }            
          }
          ncols_ijv.push_back(jcnt);
        }
      }

      HYPRE_IJMatrixSetValues(ijv, (iupper-ilower+1), &ncols_ijv[0], &rows_ijv[0], &cols_ijv[0], &vals_ijv[0]);
      HYPRE_IJMatrixAssemble(ijv);
      HYPRE_IJMatrixGetObject(ijv, (void **) &A);

      HYPRE_IJVectorCreate(comm, jlower, jupper, &vb);
      HYPRE_IJVectorCreate(comm, jlower, jupper, &vx);

      HYPRE_IJVectorSetObjectType(vb, HYPRE_PARCSR);
      HYPRE_IJVectorSetObjectType(vx, HYPRE_PARCSR);

      HYPRE_IJVectorInitialize(vb);
      HYPRE_IJVectorInitialize(vx);

      for (ordinal_type j=jlower;j<=jupper;++j) {
        rows_vb.push_back(j);
        rows_vx.push_back(j);
        vals_vb.push_back(1.0);
        vals_vx.push_back(0.0);
      }

      HYPRE_IJVectorSetValues(vb, rows_vb.size(), &rows_vb[0], &vals_vb[0]);
      HYPRE_IJVectorSetValues(vx, rows_vx.size(), &rows_vx[0], &vals_vx[0]);

      HYPRE_IJVectorAssemble(vb);
      HYPRE_IJVectorAssemble(vx);

      HYPRE_IJVectorGetObject(vb, (void **) &b);
      HYPRE_IJVectorGetObject(vx, (void **) &x);
    }
  }
  // ===================================================================================
  {
    HYPRE_Solver pcg;
    {
      Teuchos::TimeMonitor evalTimerMonitor(*t_pcg_create);
      HYPRE_ParCSRPCGCreate(comm, &pcg);
    }

    HYPRE_Solver eu;
    {
      Teuchos::TimeMonitor evalTimerMonitor(*t_euclid_create);
      HYPRE_EuclidCreate(comm, &eu);
      HYPRE_EuclidSetParams(eu, argc, argv);
    }
    {
      Teuchos::TimeMonitor evalTimerMonitor(*t_pcg_set_precond);
      HYPRE_PCGSetPrecond(pcg,
                          (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve,
                          (HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup,
                          eu);
    }
    {
      Teuchos::TimeMonitor evalTimerMonitor(*t_pcg_setup);
      HYPRE_PCGSetup(pcg, (HYPRE_Matrix)A, (HYPRE_Vector)b, (HYPRE_Vector)x);
    }
    // {
    //   Teuchos::TimeMonitor evalTimerMonitor(*t_pcg_solve);
    //   HYPRE_PCGSolve(pcg, (HYPRE_Matrix)A, (HYPRE_Vector)b, (HYPRE_Vector)x);
    // }

    {
      Teuchos::TimeMonitor evalTimerMonitor(*t_euclid_destroy);
      HYPRE_EuclidDestroy(eu);
    }
    {
      Teuchos::TimeMonitor evalTimerMonitor(*t_pcg_destroy);
      HYPRE_ParCSRPCGDestroy(pcg);
    }

    Teuchos::TimeMonitor::summarize(cout);
  }

  // Finalize MPI
  MPI_Finalize();

  return 0;
}

