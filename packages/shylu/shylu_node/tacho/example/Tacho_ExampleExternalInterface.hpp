#ifndef __TACHO_EXAMPLE_EXTERNALINTERFACE_HPP__
#define __TACHO_EXAMPLE_EXTERNALINTERFACE_HPP__
#include "Tacho_config.h"

#if defined(TACHO_USE_INT_INT) 

#include "Tacho.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"
#include "Tacho_Solver.hpp"
#include "Tacho_CommandLineParser.hpp"

namespace tacho {

  enum TACHO_PARAM_INDICES 
    {USEDEFAULTSOLVERPARAMETERS,
     VERBOSITY,
     SMALLPROBLEMTHRESHOLDSIZE,
     TASKING_OPTION_BLOCKSIZE,
     TASKING_OPTION_PANELSIZE,
     TASKING_OPTION_MAXNUMSUPERBLOCKS, 
     LEVELSET_OPTION_SCHEDULING,
     LEVELSET_OPTION_DEVICE_LEVEL_CUT,
     LEVELSET_OPTION_DEVICE_FACTOR_THRES,
     LEVELSET_OPTION_DEVICE_SOLVE_THRES,
     LEVELSET_OPTION_NSTREAMS,
     LEVELSET_OPTION_VARIANT,
     INDEX_LENGTH
    };

  using device_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::device_type;
  using exec_space = typename device_type::execution_space;
  
  using host_device_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type;
  using host_space = typename host_device_type::execution_space;

  using ViewVectorType = Kokkos::View<double*,device_type>; 
  using ViewVectorTypeInt = Kokkos::View<int*,device_type>; 

  template <class SX> class tachoSolver
  {
  public:
    using sched_type = Kokkos::TaskSchedulerMultiple<exec_space>;
    typedef Tacho::Solver<SX,sched_type> solver_type;

    typedef Tacho::ordinal_type ordinal_type;
    typedef Tacho::size_type size_type;

    typedef typename solver_type::value_type value_type;
    typedef typename solver_type::ordinal_type_array ordinal_type_array;
    typedef typename solver_type::size_type_array size_type_array;
    typedef typename solver_type::value_type_array value_type_array;
    typedef typename solver_type::value_type_matrix value_type_matrix;

    // host typedefs
    typedef typename solver_type::ordinal_type_array_host ordinal_type_array_host;
    typedef typename solver_type::size_type_array_host size_type_array_host;
    typedef typename solver_type::value_type_array_host value_type_array_host;
    typedef typename solver_type::value_type_matrix_host value_type_matrix_host;
    
    tachoSolver(const int* solverParams) :
      m_numRows(0),
      m_Solver() 
    {
      setSolverParameters(solverParams);
    }

    ~tachoSolver() 
    {
      m_Solver.release();
    }

    int Initialize(int numRows,
                   /// with TACHO_ENABLE_INT_INT, size_type is "int"
		   int* rowBegin,
		   int* columns,
		   SX* values)
    {
      m_numRows = numRows;
      if (m_numRows == 0) return 0;
      
      const int numTerms = rowBegin[numRows];
      size_type_array_host    ap_host((size_type*)   rowBegin, numRows+1);
      ordinal_type_array_host aj_host((ordinal_type*)columns,  numTerms);

#if defined (KOKKOS_ENABLE_CUDA)
      /// transfer A into device
      value_type_array ax(Kokkos::ViewAllocateWithoutInitializing("ax"), 
			  numTerms);
      value_type_array_host ax_host(values, numTerms);
      Kokkos::deep_copy(ax, ax_host);
#else
      /// wrap pointer on host 
      value_type_array ax(values, numTerms);
#endif

      Kokkos::Impl::Timer timer; 
      {
        timer.reset();
	m_Solver.analyze(numRows, ap_host, aj_host);
        const double t = timer.seconds();
        std::cout << "ExternalInterface:: analyze time " << t << std::endl;
      }
      {
        timer.reset();
        m_Solver.initialize();
        const double t = timer.seconds();
        std::cout << "ExternalInterface:: initialize time " << t << std::endl;
      }
      /// I recommend to separate factorization from the solver initialization
      {
        timer.reset();
	m_Solver.factorize(ax);
        const double t = timer.seconds();
        std::cout << "ExternalInterface:: factorize time " << t << std::endl;
      }

      return 0;
    }

    void exportSupernodes(std::vector<int> &supernodes) {
      const auto supernodes_device = m_Solver.getSupernodes();
      auto supernodes_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), supernodes_device); 
      {
        const int n = supernodes_host.extent(0);
        supernodes.resize(n);
        std::copy(supernodes_host.data(), supernodes_host.data()+n, supernodes.data());
      }
    }

    void exportUpperTriangularFactorsToCrsMatrix(std::vector<int> &rowBeginU,
                                                 std::vector<int> &columnsU,
                                                 std::vector<SX> &valuesU,
                                                 std::vector<int> &perm) {
      {
        const auto perm_device = m_Solver.getPermutationVector();
        auto perm_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), perm_device);
        {
          const int n = perm_host.extent(0);
          perm.resize(n);
          std::copy(perm_host.data(), perm_host.data()+n, perm.data());
        }
      }
      
      {
        typename solver_type::crs_matrix_type A;
        m_Solver.exportFactorsToCrsMatrix(A);

        typename solver_type::crs_matrix_type_host A_host;
        A_host.createMirror(A);
        A_host.copy(A); 
        {
          const auto ap = A_host.RowPtr();
          const int n = ap.extent(0);
          rowBeginU.resize(n);
          std::copy(ap.data(), ap.data()+n, rowBeginU.data());
        }
        {
          const auto aj = A_host.Cols();
          const int n = aj.extent(0);
          columnsU.resize(n);
          std::copy(aj.data(), aj.data()+n, columnsU.data());
        }
        {
          const auto ax = A_host.Values();
          const int n = ax.extent(0);
          valuesU.resize(n);
          std::copy(ax.data(), ax.data()+n, valuesU.data());
        }
      }
    }

    void MySolve(int NRHS,
                 value_type_matrix &b, 
                 value_type_matrix &x)
    {
      if (m_numRows == 0) return;

      const int m = m_numRows;      
      if (static_cast<int>(m_TempRhs.extent(0)) < m || 
          static_cast<int>(m_TempRhs.extent(1)) < NRHS)
        m_TempRhs = value_type_matrix("temp rhs", m, NRHS);      
      
      m_Solver.solve(x, b, m_TempRhs);
    }

  private:
    int m_numRows;
    solver_type m_Solver;
    value_type_matrix m_TempRhs;

    void setSolverParameters(const int* solverParams)
    {
      if (solverParams[USEDEFAULTSOLVERPARAMETERS]) return;
      // common options
      m_Solver.setVerbose                       (solverParams[VERBOSITY]);
      m_Solver.setSmallProblemThresholdsize     (solverParams[SMALLPROBLEMTHRESHOLDSIZE]);

      // tasking options
      m_Solver.setBlocksize                     (solverParams[TASKING_OPTION_BLOCKSIZE]);
      m_Solver.setPanelsize                     (solverParams[TASKING_OPTION_PANELSIZE]);
      m_Solver.setMaxNumberOfSuperblocks        (solverParams[TASKING_OPTION_MAXNUMSUPERBLOCKS]);

      // levelset options
      m_Solver.setLevelSetScheduling            (solverParams[LEVELSET_OPTION_SCHEDULING]);      
      m_Solver.setLevelSetOptionDeviceLevelCut  (solverParams[LEVELSET_OPTION_DEVICE_LEVEL_CUT]);
      m_Solver.setLevelSetOptionDeviceFunctionThreshold
        (solverParams[LEVELSET_OPTION_DEVICE_FACTOR_THRES], 
         solverParams[LEVELSET_OPTION_DEVICE_SOLVE_THRES]);
      m_Solver.setLevelSetOptionNumStreams      (solverParams[LEVELSET_OPTION_NSTREAMS]);
      m_Solver.setLevelSetOptionAlgorithmVariant(solverParams[LEVELSET_OPTION_VARIANT]);
    }

  };

} // namespace tacho

#endif
#endif // INTERFACE_TACHO_SOLVER_H
