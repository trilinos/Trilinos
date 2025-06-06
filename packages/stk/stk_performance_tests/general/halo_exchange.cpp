#include <math.h>
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternBlocking.hpp"
#include "stk_unit_test_utils/timer.hpp"
#include <vector>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <string>
#include <iostream>
#include <memory>
#include <fstream>
#include "stk_unit_test_utils/getOption.h"
#include "gtest/gtest.h"

//#define WRITE_TIMING_FILES

namespace {

std::vector<double> get_stencil(int npts);

class TimerWrapper
{
  public:
    explicit TimerWrapper(int nruns) :
      m_timer(MPI_COMM_WORLD),
      m_nruns(nruns)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      m_timer.start_timing(); 
    }

    ~TimerWrapper()
    {
      MPI_Barrier(MPI_COMM_WORLD);
      m_timer.update_timing();
      m_timer.print_timing(m_nruns);
    }

  private:
    stk::unit_test_util::Timer m_timer;
    int m_nruns = 1;
};

//-----------------------------------------------------------------------------
// Exact solution

class ExactSolution
{
  public:
    virtual ~ExactSolution() {}

    virtual double get_exact_solution(double x, double y) const = 0;

    virtual double get_exact_solution_xx(double x, double y) const = 0;

    virtual double get_exact_solution_yy(double x, double y) const = 0;
};

class ExactSolutionPoly : public ExactSolution
{
  public:
    explicit ExactSolutionPoly(int degree) :
      m_degree(degree)
    {}

    double get_exact_solution(double x, double y) const override
    {
      if (m_degree == 0)
        return 1;
      else
      {
        return std::pow(x, m_degree) + std::pow(y, m_degree);
      }
    }

    double get_exact_solution_xx(double x, double /*y*/) const override
    {
      return get_Uxx(x);
    }

    double get_exact_solution_yy(double /*x*/, double y) const override
    {
      return get_Uxx(y);
    }
    
  private:
    double get_Uxx(double coord) const
    {
      if (m_degree <= 1)
        return 0;
      else if (m_degree == 2)
        return 2;
      else
        return m_degree * (m_degree - 1) * std::pow(coord, m_degree - 2);
    }

    int m_degree;
};


class ExactSolutionPeriodic : public ExactSolution
{
  public:
    ExactSolutionPeriodic(double xmin, double xmax, double ymin, double ymax)
    {
      double pi = std::atan(1)*4;
      m_betaX = 2*pi/(xmax - xmin);
      m_betaY = 2*pi/(ymax - ymin);
    }


    double get_exact_solution(double x, double y) const override
    {
      return std::sin(m_betaX * x) + std::sin(m_betaY * y);
    }

    double get_exact_solution_xx(double /*x*/, double /*y*/) const override
    {
      throw std::runtime_error("unsupported function");
    }

    double get_exact_solution_yy(double /*x*/, double /*y*/) const override
    {
      throw std::runtime_error("unsupported function");
    }

  private:
    double m_betaX;
    double m_betaY;
};

//-----------------------------------------------------------------------------
// Input configuration

struct RunConfig
{
  MPI_Comm comm;
  int nptsPerDirection;
  int nghost;
  int stencilSize;
  std::shared_ptr<ExactSolution> exSol;
  double xMin;
  double xMax;
  double yMin;
  double yMax;
  double loadImbalance;
};

RunConfig get_performance_run_config()
{
    int stencilSize = 5;
    int nptsPerDirection = stk::unit_test_util::get_command_line_option("-npts_per_direction", 250);
    int nghost = stk::unit_test_util::get_command_line_option("-nghost", 10);
    double xMin = 0;
    double xMax = 1;
    double yMin = 0;
    double yMax = 1;

    // deal with the off-by-one issue when using periodic BCs
    int nprocs_per_direction = std::round(std::sqrt(stk::parallel_machine_size(MPI_COMM_WORLD)));
    double dx = (xMax - xMin)/(nprocs_per_direction * nptsPerDirection - 1);
    double dy = (yMax - yMin)/(nprocs_per_direction * nptsPerDirection - 1);
    auto exSol = std::make_shared<ExactSolutionPeriodic>(xMin, xMax + dx, yMin, yMax + dy);
    double loadImbalanceFac = stk::unit_test_util::get_command_line_option("-load_imbalance_fac", 0.0);
    double loadImbalance = 1.0 + loadImbalanceFac * stk::parallel_machine_rank(MPI_COMM_WORLD);

    return {MPI_COMM_WORLD, nptsPerDirection, nghost, stencilSize, exSol, xMin, xMax, yMin, yMax, loadImbalance};
}


//-----------------------------------------------------------------------------
// Timing

struct TimingStats
{
  double tInteriorCompute = 0;
  double tBoundaryCompute = 0;
  double tStartComm = 0;
  double tPack = 0;

  double tFinishComm = 0;
  double tTest = 0;
  double tSetup = 0;
  double tTotal = 0;
};

std::vector<std::string> TimingStatsNames{"t_interior_compute", "t_boundary_compute", "t_start_comm", "t_pack",
                                          "t_finish_comm", "t_test", "t_setup", "t_total"};

std::string format_percent(double val)
{
  const int bufSize = 64;
  char buf[bufSize];
  std::snprintf(buf, bufSize, "%3.2f", val);

  return std::string(buf) + "%";
}

std::string print_line(std::string name, double val, double total, int nameMaxLength=20, int valueMaxLength=10)
{
  std::string val_s = std::to_string(val);
  int namePaddingLen = nameMaxLength - name.size();
  int valuePaddingLen = valueMaxLength - val_s.size();

  for (int i=0; i < namePaddingLen; ++i)
    name += " ";

  for (int i=0; i < valuePaddingLen; ++i)
    val_s += " ";

  return std::string("  ") + name + ": " + val_s + " (" + format_percent(100 * val / total) + ")";
}

void print_timing_stats(std::ostream& os, const TimingStats& timer)
{

  os << "Total time: " << timer.tTotal << " s" << std::endl;
  os << print_line("t_setup", timer.tSetup, timer.tTotal) << std::endl;
  os << print_line("t_interior_compute", timer.tInteriorCompute, timer.tTotal) << std::endl;
  os << print_line("t_boundary_compute", timer.tBoundaryCompute, timer.tTotal) << std::endl;
  os << print_line("t_start_comm", timer.tStartComm, timer.tTotal) << std::endl;
  os << print_line("t_pack", timer.tPack, timer.tTotal) << std::endl;
  os << print_line("t_finish_comm", timer.tFinishComm, timer.tTotal) << std::endl;
  os << print_line("t_test", timer.tTest, timer.tTotal) << std::endl;

  // Note: tBoundaryCompute occurs inside tFinishComm, so it is intentionally omitted here
  double tOther = timer.tTotal - (timer.tInteriorCompute + timer.tStartComm +
                                    timer.tFinishComm + timer.tTest +
                                    timer.tSetup);
                
  os << print_line("t_other", tOther, timer.tTotal) << std::endl;
}

#ifdef WRITE_TIMING_FILES
void write_timing_files(const TimingStats& timer, MPI_Comm comm, const std::string& prefix)
{
  // collect all data on root process
  int commRank = stk::parallel_machine_rank(comm);
  int commSize = stk::parallel_machine_size(comm);

  int nvals = 8;
  double sendData[] = {timer.tInteriorCompute, timer.tBoundaryCompute, timer.tStartComm, timer.tPack,
                        timer.tFinishComm, timer.tTest, timer.tSetup, timer.tTotal};
  int root = 0;
  if (commRank != root)
  {
    MPI_Gather(sendData, nvals, MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, root, comm);
  } else
  {
    std::vector<double> recvData(nvals * commSize);
    MPI_Gather(sendData, nvals, MPI_DOUBLE, recvData.data(), nvals, MPI_DOUBLE, root, comm);

    // write names file
    std::ofstream name_file(prefix + "_field_names.txt");
    for (auto name : TimingStatsNames)
      name_file << name << "\n";
    name_file.close();

    // write data file
    std::ofstream dataFile(prefix + "_field_values.txt");
    dataFile << std::scientific << std::setprecision(16);

    for (int field=0; field < nvals; ++field)
    {
      for (int i=0; i < commSize; ++i)
      {
        dataFile << recvData[i * nvals + field];
        if (i < commSize - 1)
          dataFile << ", ";
      }
      dataFile << "\n";
    }
    dataFile.close();
  }


  MPI_Barrier(comm);  // not really needed, but helps to keep processes aligned when running
                      // several performance tests in a row

}
#endif

//-----------------------------------------------------------------------------
// Heat equation solver

// mock up of 2D heat equation solved with finite difference, to test different methods
// doing parallel data exchange
// Currently, the number of processors must be a square (ie 2^2, 3^2, etc.)
class HeatEquationHalo
{
  public: 
    explicit HeatEquationHalo(const RunConfig& config) :
        m_comm(config.comm),
        m_commSize(stk::parallel_machine_size(config.comm)),
        m_commRank(stk::parallel_machine_rank(config.comm)),
        m_nprocsPerDirection(std::round(std::sqrt(m_commSize))),
        mProcX(m_commRank % m_nprocsPerDirection),
        mProcY((m_commRank - mProcX) / m_nprocsPerDirection),
        m_nptsPerDirection(config.nptsPerDirection),
        m_nghost(config.nghost),
        m_stencilSize(config.stencilSize),
        m_exSol(config.exSol),
        m_xmin(config.xMin),
        m_xmax(config.xMax),
        m_ymin(config.yMin),
        m_ymax(config.yMax),
        m_loadImbalance(config.loadImbalance),
        m_stencil(get_stencil(config.stencilSize)),
        m_sol(std::pow( config.nptsPerDirection + 2*config.nghost, 2)),
        m_res(config.nptsPerDirection * config.nptsPerDirection),
        m_sendBufs(4),
        m_recvBufsExact(4)
    {
        if (!(m_nprocsPerDirection * m_nprocsPerDirection == m_commSize))
          throw std::runtime_error("number of MPI processes must be a square");

        if (m_commSize == 4)
          throw std::runtime_error("number of MPI processes cannot be 4");

        if (!(get_proc_idx(mProcX, mProcY) == m_commRank))
          throw std::runtime_error("error in computing cartesian decomposition");

        if (!(config.stencilSize % 2 == 1))
          throw std::runtime_error("stencil size must be odd");

        if (!(config.nghost >= (config.stencilSize - 1)/2))
          throw std::runtime_error("nghost must be >= half of stencil size");

        
        for (int i=0; i < 4; ++i)
        {
          m_sendBufs[i].resize(m_nghost * m_nptsPerDirection);
          m_recvBufsExact[i].resize(m_nghost * m_nptsPerDirection);
        }

        set_recv_buffer_exact_values();
    }

    virtual ~HeatEquationHalo() {}

    const TimingStats& run(int nsteps)
    {
      MPI_Barrier(m_comm);
      double t_start = MPI_Wtime();

      initialize_solution();
      set_send_buffer_values();
      int nrepeats = std::floor(m_loadImbalance);
      double load_fraction = m_loadImbalance - nrepeats;

      m_timer.tSetup += MPI_Wtime() - t_start;

      for (int i=0; i < nsteps; ++i)
      {
        reset_iteration();

        set_communication_timed();

        for (int j=0; j < nrepeats; ++j)
          compute_interior();

        compute_interior_fraction(load_fraction);

        finish_communication_timed();

        test_ghost_values();
        //test_residual();
      }

      
      m_timer.tTotal += MPI_Wtime() - t_start;

      return m_timer;
    }

    void test_residual()
    {
      double t_start = MPI_Wtime();
      for (int j=0; j < m_nptsPerDirection; ++j)
        for (int i=0; i < m_nptsPerDirection; ++i)
        {          
          double x = get_x_coord(i);
          double y = get_y_coord(j);
          double expectedVal = m_exSol->get_exact_solution_xx(x, y) + m_exSol->get_exact_solution_yy(x, y);
          double tol = 1e-10;
          if (std::abs(expectedVal) > 1) 
            tol *= std::abs(expectedVal);

          EXPECT_NEAR(m_res[get_res_idx(i, j)], expectedVal, tol);
        }

      m_timer.tTest += MPI_Wtime() - t_start;
    }

  protected:
    // indexing function for solution
    int get_sol_idx(int i, int j)
    {
      assert(i >= -m_nghost);
      assert(i <=  m_nptsPerDirection + m_nghost - 1);
      assert(j >= -m_nghost);
      assert(j <=  m_nptsPerDirection + m_nghost - 1);

      int i2 = i + m_nghost;
      int j2 = j + m_nghost;
      return i2 + (m_nptsPerDirection + 2*m_nghost) * j2;
    }

    // indexing function for residual
    int get_res_idx(int i, int j)
    {
      assert(i >= 0);
      assert(i <  m_nptsPerDirection);
      assert(j >= 0);
      assert(j <  m_nptsPerDirection);

      return i + m_nptsPerDirection * j;
    }

    // indexing function for MPI ranks
    // allows negative indexing for periodic wrap-around
    int get_proc_idx(int i, int j)
    {
      int i2 = (i + m_nprocsPerDirection) % m_nprocsPerDirection;
      int j2 = (j + m_nprocsPerDirection) % m_nprocsPerDirection;

      return i2 + m_nprocsPerDirection * j2;
    }

    int get_neighbor_proc_rank(int bndry)
    {
      switch (bndry)
      {
        case 0: return get_proc_idx(mProcX,     mProcY - 1);
        case 1: return get_proc_idx(mProcX + 1, mProcY);
        case 2: return get_proc_idx(mProcX,     mProcY + 1);
        case 3: return get_proc_idx(mProcX - 1, mProcY);
        default:
          throw std::invalid_argument("There can be only 4 boundaries");
      }
    }

    double get_dx()
    {
      return (m_xmax - m_xmin) / (m_nprocsPerDirection * m_nptsPerDirection - 1);
    }

    double get_dy()
    {
      return (m_ymax - m_ymin) / (m_nprocsPerDirection * m_nptsPerDirection - 1);
    }

    double get_x_coord(int i)
    {
      int iGlobal = i + mProcX * m_nptsPerDirection;
      int nptsGlobal = m_nprocsPerDirection * m_nptsPerDirection;
      int i2 = (iGlobal + nptsGlobal) % nptsGlobal;
      return i2 * get_dx();
    }

    double get_y_coord(int j)
    {
      int jGlobal = j + mProcY * m_nptsPerDirection;
      int nptsGlobal = m_nprocsPerDirection * m_nptsPerDirection;
      int j2 = (jGlobal + nptsGlobal) % nptsGlobal;
      return j2 * get_dy();
    }

    void initialize_solution()
    {
      for (int j=0; j < m_nptsPerDirection; ++j)
        for (int i=0; i < m_nptsPerDirection; ++i)
        {
          double x = get_x_coord(i);
          double y = get_y_coord(j);
          m_sol[get_sol_idx(i, j)] = m_exSol->get_exact_solution(x, y);
        }
    }

    void set_communication_timed()
    {
      double tStartComm = MPI_Wtime();
      if (m_commSize == 1) 
      {
        set_exact_ghost_values();
      } else
      {
        start_communication();
      }
      m_timer.tStartComm += MPI_Wtime() - tStartComm;
    }


    void finish_communication_timed()
    {
      double tFinishComm = MPI_Wtime();
      if (m_commSize == 1)
      {
        for (int i=0; i < 4; ++i)
          compute_boundary(i);
      }
      else {
        finish_communication(); 
      }
      m_timer.tFinishComm += MPI_Wtime() - tFinishComm;
    }

    void reset_iteration()
    {
      for (auto& v : m_res)
        v = 0;

      sabotage_ghost_values();
    }

    void compute_interior();

    void compute_interior_fraction(double fraction);


    // compute the points that need ghost data
    // bndry identifies which boundary: South, East, North, West
    void compute_boundary(int bndry);

    double apply_stencil_x(int i, int j)
    {
      int nghostUsed = (m_stencilSize - 1) / 2;
      double val = 0;
      int idx = 0;
      for (int i2 = i - nghostUsed; i2 <= i + nghostUsed; ++i2)
        val += m_stencil[idx++] * m_sol[get_sol_idx(i2, j)];

      return val;
    }

    double apply_stencil_y(int i, int j)
    {
      int nghostUsed = (m_stencilSize - 1) / 2;
      double val = 0;
      int idx = 0;
      for (int j2 = j - nghostUsed; j2 <= j + nghostUsed; ++j2)
      {
        val += m_stencil[idx++] * m_sol[get_sol_idx(i, j2)];
      }

      return val;
    }

    void local_periodic_communication()
    {
      for (int i=0; i < m_nptsPerDirection; ++i)
        for (int g=1; g <= m_nghost; ++g)
        {
          m_sol[get_sol_idx(i, -g)] = m_sol[get_sol_idx(i, m_nptsPerDirection - g)];
          m_sol[get_sol_idx(i, m_nptsPerDirection + g - 1)] = m_sol[get_sol_idx(i, g - 1)];
        }

      for (int j=0; j < m_nptsPerDirection; ++j)
        for (int g=1; g <= m_nghost; ++g)
        {
          m_sol[get_sol_idx(-g, j)] = m_sol[get_sol_idx(m_nptsPerDirection - g, j)];
          m_sol[get_sol_idx(m_nptsPerDirection + g - 1, j)] = m_sol[get_sol_idx(g - 1, j)];
        }
    }

    void set_exact_ghost_values()
    {
      for (int i=0; i < 4; ++i)
        unpack_receive_buffer(i, m_recvBufsExact[i], false);
    }

    void set_send_buffer_values()
    {
      for (int i=0; i < m_nptsPerDirection; ++i)
        for (int g=1; g <= m_nghost; ++g)
        {
          m_sendBufs[0][i * m_nghost + g - 1] = m_sol[get_sol_idx(i, g - 1)];
          m_sendBufs[2][i * m_nghost + g - 1] = m_sol[get_sol_idx(i, m_nptsPerDirection - g)];
        }

      for (int j=0; j < m_nptsPerDirection; ++j)
        for (int g=1; g <= m_nghost; ++g)
        {
          m_sendBufs[3][j * m_nghost + g - 1] = m_sol[get_sol_idx(g - 1, j)];
          m_sendBufs[1][j * m_nghost + g - 1] = m_sol[get_sol_idx(m_nptsPerDirection - g, j)];
        }
    }


    void set_recv_buffer_exact_values()
    {
      for (int i=0; i < m_nptsPerDirection; ++i)
        for (int g=1; g <= m_nghost; ++g)
        {
          double x = get_x_coord(i);
          double yLower = get_y_coord(0) - get_dy() * g;
          double yUpper = get_y_coord(m_nptsPerDirection - 1) + get_dy() * g;

          m_recvBufsExact[0][i * m_nghost + g - 1] = m_exSol->get_exact_solution(x, yLower);
          m_recvBufsExact[2][i * m_nghost + g - 1] = m_exSol->get_exact_solution(x, yUpper);
        }

      for (int j=0; j < m_nptsPerDirection; ++j)
        for (int g=1; g <= m_nghost; ++g)
        {
          double y = get_y_coord(j);
          double xLower = get_x_coord(0) - get_dx() * g;
          double xUpper = get_x_coord(m_nptsPerDirection - 1) + get_dx() * g;
          m_recvBufsExact[3][j * m_nghost + g - 1] = m_exSol->get_exact_solution(xLower, y);
          m_recvBufsExact[1][j * m_nghost + g - 1] = m_exSol->get_exact_solution(xUpper, y);
        }
    }

    void sabotage_ghost_values()
    {
      for (int i=0; i < m_nptsPerDirection; ++i)
        for (int g=1; g <= m_nghost; ++g)
        {
          m_sol[get_sol_idx(i, -g)] = 666;
          m_sol[get_sol_idx(i, m_nptsPerDirection + g - 1)] = 667;
        }

      for (int j=0; j < m_nptsPerDirection; ++j)
        for (int g=1; g <= m_nghost; ++g)
        {
          m_sol[get_sol_idx(-g, j)] = 668;
          m_sol[get_sol_idx(m_nptsPerDirection + g - 1, j)] = 669;
        }
    }

    //-----------------------------------------------------
    // Interface for derived classes to implement communication

    double get_send_buffer_value(int bndry, int idx, int ghost)
    {
      assert(bndry >= 0 && bndry <= 4);
      assert(idx >= 0 && idx < m_nptsPerDirection);
      assert(ghost >= 0 && ghost < m_nghost);

      return m_sendBufs[bndry][idx * m_nghost + ghost];
    }

    virtual void unpack_receive_buffer(int bndry, const std::vector<double>& vals, bool compute=true)
    {
      assert(vals.size() == static_cast<size_t>(m_nghost * m_nptsPerDirection));
  
      if (bndry == 0 || bndry == 2)
      {
        int jStart = bndry == 0 ? -1 : m_nptsPerDirection;
        int jEnd   = bndry == 0 ? -m_nghost - 1: m_nptsPerDirection + m_nghost;
        int inc    = bndry == 0 ? -1 : 1;
        
        int idx = 0;
        for (int i=0; i < m_nptsPerDirection; ++i)
          for (int j=jStart; j != jEnd; j = j + inc)
            m_sol[get_sol_idx(i, j)] = vals[idx++];

      } else if (bndry == 1 || bndry == 3)
      {
        int iStart = bndry == 3 ? -1 : m_nptsPerDirection;
        int iEnd   = bndry == 3 ? -m_nghost - 1: m_nptsPerDirection + m_nghost;
        int inc    = bndry == 3 ? -1 : 1;
        
        int idx = 0;
        for (int j=0; j < m_nptsPerDirection; ++j)
          for (int i=iStart; i != iEnd; i = i + inc)
            m_sol[get_sol_idx(i, j)] = vals[idx++]; 
      } else
      {
        throw std::invalid_argument("There can be only 4 boundaries");
      }

      if (compute)
        compute_boundary(bndry);
    }

    void test_ghost_values()
    {
      for (int i=0; i < 4; ++i)
        test_ghost_values(i);
    }


    void test_ghost_values(int bndry)
    {
      double t_start = MPI_Wtime();
      if (bndry == 0 || bndry == 2)
      {
        int jStart = bndry == 0 ? -1 : m_nptsPerDirection;
        int jEnd   = bndry == 0 ? -m_nghost - 1: m_nptsPerDirection + m_nghost;
        int inc    = bndry == 0 ? -1 : 1;
        
        int idx = 0;
        for (int i=0; i < m_nptsPerDirection; ++i)
          for (int j=jStart; j != jEnd; j = j + inc)
          {
            EXPECT_NEAR(m_sol[get_sol_idx(i, j)], m_recvBufsExact[bndry][idx++], 1e-13);
          }

      } else if (bndry == 1 || bndry == 3)
      {
        int iStart = bndry == 3 ? -1 : m_nptsPerDirection;
        int iEnd   = bndry == 3 ? -m_nghost - 1: m_nptsPerDirection + m_nghost;
        int inc    = bndry == 3 ? -1 : 1;
        
        int idx = 0;
        for (int j=0; j < m_nptsPerDirection; ++j)
          for (int i=iStart; i != iEnd; i = i + inc)
            EXPECT_NEAR(m_sol[get_sol_idx(i, j)], m_recvBufsExact[bndry][idx++], 1e-13); 
      } else
      {
        throw std::invalid_argument("There can be only 4 boundaries");
      }

      m_timer.tTest += MPI_Wtime() - t_start;
    }

    // the implementations of this function should use get_send_buffer_value() to get
    // the data to send to other processors and start the MPI communication
    virtual void start_communication() = 0;

    // the implementation of this function should finish the MPI communication and
    // call unpack_receive_buffer().  If unpack_receive_buffer() is called with
    // compute=false, then this function should also call compute_boundary() afterwards
    virtual void finish_communication() = 0;


    MPI_Comm m_comm;
    int m_commSize;
    int m_commRank;
    int m_nprocsPerDirection;
    int mProcX;
    int mProcY;
    int m_nptsPerDirection;
    int m_nghost;
    int m_stencilSize;
    std::shared_ptr<ExactSolution> m_exSol;
    double m_xmin = 0.0;
    double m_xmax = 1.0;
    double m_ymin = 0.0;
    double m_ymax = 1.0;
    double m_loadImbalance;
    std::vector<double> m_stencil;

    std::vector<double> m_sol;
    std::vector<double> m_res;
    TimingStats m_timer;

  private:
    std::vector<std::vector<double>> m_sendBufs;
    std::vector<std::vector<double>> m_recvBufsExact;
};

void HeatEquationHalo::compute_interior()
{
  double tStartInterior = MPI_Wtime();

  int nghostUsed = (m_stencilSize - 1) / 2;
  double fac_x = 1.0/(get_dx() * get_dx());
  double fac_y = 1.0/(get_dy() * get_dy());

  for (int i=0; i < m_nptsPerDirection; ++i)
    for (int j=nghostUsed; j < m_nptsPerDirection - nghostUsed; ++j)
      m_res[get_res_idx(i, j)] += apply_stencil_y(i, j) * fac_y;


  for (int j=0; j < m_nptsPerDirection; ++j)
    for (int i=nghostUsed; i < m_nptsPerDirection - nghostUsed; ++i)
      m_res[get_res_idx(i, j)] += apply_stencil_x(i, j) * fac_x;

  m_timer.tInteriorCompute += MPI_Wtime() - tStartInterior;
}


void HeatEquationHalo::compute_interior_fraction(double fraction)
{
  double tStartInterior = MPI_Wtime();

  int nghostUsed = (m_stencilSize - 1) / 2;
  double facX = 1.0/(get_dx() * get_dx());
  double facY = 1.0/(get_dy() * get_dy());

  int nptsPerDirection = std::round(std::sqrt(fraction) * m_nptsPerDirection);

  for (int i=0; i < nptsPerDirection; ++i)
    for (int j=nghostUsed; j < nptsPerDirection - nghostUsed; ++j)
      m_res[get_res_idx(i, j)] += apply_stencil_y(i, j) * facY;


  for (int j=0; j < nptsPerDirection; ++j)
    for (int i=nghostUsed; i < nptsPerDirection - nghostUsed; ++i)
      m_res[get_res_idx(i, j)] += apply_stencil_x(i, j) * facX;

  m_timer.tInteriorCompute += MPI_Wtime() - tStartInterior;
}

void HeatEquationHalo::compute_boundary(int bndry)
{
  double tStart = MPI_Wtime();

  int nghostUsed = (m_stencilSize - 1) / 2;

  if (bndry == 0 || bndry == 2) // South or North
  {
    double facY = 1.0/(get_dy() * get_dy());
    int jStart = bndry == 0 ? 0           : m_nptsPerDirection - nghostUsed;
    int jEnd   = bndry == 0 ? nghostUsed  : m_nptsPerDirection;

    for (int i=0; i < m_nptsPerDirection; ++i)
      for (int j=jStart; j < jEnd; ++j)
        m_res[get_res_idx(i, j)] += apply_stencil_y(i, j) * facY;

  } else if (bndry == 1 || bndry == 3) // East or West
  {
    double facX = 1.0/(get_dx() * get_dx());
    int iStart = bndry == 3 ? 0           : m_nptsPerDirection - nghostUsed;
    int iEnd   = bndry == 3 ? nghostUsed  : m_nptsPerDirection;

    for (int j=0; j < m_nptsPerDirection; ++j)
      for (int i=iStart; i < iEnd; ++i)
        m_res[get_res_idx(i, j)] += apply_stencil_x(i, j) * facX;

  } else
  {
    throw std::invalid_argument("There can be only 4 boundaries");
  }

  m_timer.tBoundaryCompute += MPI_Wtime() - tStart;
}

std::vector<double> get_stencil(int npts)
{
  switch (npts)
  {
    case 3:
      return {1, -2, 1};
    case 5:
      return {-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0};
    case 7:
      return {1.0/90.0, -3.0/20.0, 3.0/2.0, -49.0/18.0, 3.0/2.0, -3.0/20.0, 1.0/90.0};
    case 9:
      return {-1.0/560.0, 8.0/315.0, -1.0/5.0, 8.0/5.0, -205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0};
    default:
      throw std::runtime_error("unsupported stencil size: " + std::to_string(npts));

  };
}

//-----------------------------------------------------------------------------
// Implementations

class HeatEquationHaloSerial : public HeatEquationHalo
{
  public:
    explicit HeatEquationHaloSerial(const RunConfig& config) :
      HeatEquationHalo(config)
    {}

  protected:
    void start_communication() override
    {
      throw std::runtime_error("This class is for serial only!");
    }

    void finish_communication() override
    {
      throw std::runtime_error("This class is for serial only!");
    }
};


// the most efficient possible implementation
class HeatEquationHaloEfficient : public HeatEquationHalo
{
  public:
    explicit HeatEquationHaloEfficient(const RunConfig& config) :
      HeatEquationHalo(config),
      m_sendBufs(4),
      m_recvBufs(4),
      m_sendReqs(4),
      m_recvReqs(4),
      m_neighborRanks({get_neighbor_proc_rank(0), get_neighbor_proc_rank(1), get_neighbor_proc_rank(2), get_neighbor_proc_rank(3)}),
      m_tag(stk::get_mpi_tag_manager().get_tag(m_comm))
    {
      int npts = m_nptsPerDirection * m_nghost;
      for (int bndry=0; bndry < 4; ++bndry)
      {
        m_sendBufs[bndry].resize(npts);
        m_recvBufs[bndry].resize(npts);
      }
    }

    ~HeatEquationHaloEfficient()
    {
      if (!m_is_first)
        MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
    }

  protected:
    void start_communication() override
    {
      for (int bndry=0; bndry < 4; ++bndry)
      {
        auto& buf = m_recvBufs[bndry];
        MPI_Irecv(buf.data(), buf.size(), MPI_DOUBLE, m_neighborRanks[bndry], m_tag, m_comm, &(m_recvReqs[bndry]));
      }

      if (!m_is_first)
        MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);

      for (int bndry=0; bndry < 4; ++bndry)
      {
        int idx = 0;
        auto& buf = m_sendBufs[bndry];
        for (int i=0; i < m_nptsPerDirection; ++i)
          for (int g=0; g < m_nghost; ++g)
            buf[idx++] = get_send_buffer_value(bndry, i, g);

        MPI_Isend(buf.data(), buf.size(), MPI_DOUBLE, m_neighborRanks[bndry], m_tag, m_comm, &(m_sendReqs[bndry]));
      }

      m_is_first = false;
    }

    void finish_communication() override
    {
      for (int i=0; i < 4; ++i)
      {
        int idx=0;
        MPI_Waitany(m_recvReqs.size(), m_recvReqs.data(), &idx, MPI_STATUS_IGNORE);
        unpack_receive_buffer(idx, m_recvBufs[idx]);
      }
    }

  private:
    std::vector<std::vector<double>> m_sendBufs;
    std::vector<std::vector<double>> m_recvBufs;
    std::vector<MPI_Request> m_sendReqs;
    std::vector<MPI_Request> m_recvReqs;
    std::array<int, 4> m_neighborRanks;
    stk::MPITag m_tag;
    bool m_is_first = true;

};

//-----------------------------------------------------------------------------
// classes that use ParallelDataExchange

class HeatEquationHaloParallelDataExchangeBase : public HeatEquationHalo
{
  public:
    explicit HeatEquationHaloParallelDataExchangeBase(const RunConfig& config) :
      HeatEquationHalo(config),
      m_sendBufs(m_commSize),
      m_recvBufs(m_commSize),
      m_neighborRanks({get_neighbor_proc_rank(0), get_neighbor_proc_rank(1), get_neighbor_proc_rank(2), get_neighbor_proc_rank(3)})
    {
      for (int bndry=0; bndry < 4; ++bndry)
      {
        int rank = get_neighbor_proc_rank(bndry);
        m_sendBufs[rank].resize(m_nptsPerDirection * m_nghost);
        m_recvBufs[rank].resize(m_nptsPerDirection * m_nghost);
      }
    }

  protected:

    void pack_send_buffers()
    {
      double tStart = MPI_Wtime();
      for (int bndry=0; bndry < 4; ++bndry)
      {
        int rank = get_neighbor_proc_rank(bndry);
        int idx = 0;
        for (int i=0; i < m_nptsPerDirection; ++i)
          for (int g=0; g < m_nghost; ++g)
            m_sendBufs[rank][idx++] = get_send_buffer_value(bndry, i, g);
      }
      m_timer.tPack += MPI_Wtime() - tStart;
    }


    int get_rank_boundary(int rank)
    {
      for (int i=0; i < 4; ++i)
      {
        if (m_neighborRanks[i] == rank)
          return i;
      }

      return -1;
    }

    std::vector<std::vector<double>> m_sendBufs;
    std::vector<std::vector<double>> m_recvBufs;
    std::array<int, 4> m_neighborRanks;
};


class HeatEquationHaloParallelDataExchangeBlocking : public HeatEquationHaloParallelDataExchangeBase
{
  public:
    explicit HeatEquationHaloParallelDataExchangeBlocking(const RunConfig& config) :
      HeatEquationHaloParallelDataExchangeBase(config),
      m_exchanger(m_comm)
    {}

  protected:
    void start_communication() override
    {
      pack_send_buffers();
      m_exchanger.execute(m_sendBufs, m_recvBufs);
    }

    void finish_communication() override
    {
      for (int bndry=0; bndry < 4; ++bndry)
      {
        int rank = get_neighbor_proc_rank(bndry);
        unpack_receive_buffer(bndry, m_recvBufs[rank]);
      }
    }
  
  private:
    stk::DataExchangeUnknownPatternBlocking m_exchanger;
};



class HeatEquationHaloParallelDataExchangeNonBlocking : public HeatEquationHaloParallelDataExchangeBase
{
  public:
    explicit HeatEquationHaloParallelDataExchangeNonBlocking(const RunConfig& config) :
      HeatEquationHaloParallelDataExchangeBase(config),
      m_exchanger(m_comm)
    {}

  protected:
    void start_communication() override
    { 
      if (m_exchanger.are_sends_in_progress())
        m_exchanger.complete_sends();

      pack_send_buffers();
      m_exchanger.start_nonblocking(m_sendBufs, m_recvBufs);
    }

    void finish_communication() override
    {
      m_exchanger.post_nonblocking_receives(m_recvBufs);
      
      auto f = [&](int rank, std::vector<double>& buf)
      {
        unpack_receive_buffer(get_rank_boundary(rank), buf);
      };
      m_exchanger.complete_receives(m_recvBufs, f);
    }

  private:
    stk::DataExchangeUnknownPatternNonBlocking m_exchanger;
};



//-----------------------------------------------------------------------------
// classes that use CommSparse

class HeatEquationHaloCommSparseBase: public HeatEquationHalo
{
  public:
    explicit HeatEquationHaloCommSparseBase(const RunConfig& config) :
      HeatEquationHalo(config),
      m_commSparse(config.comm),
      m_sendBufs(m_commSize),
      m_recvBufs(m_commSize),
      m_neighborRanks({get_neighbor_proc_rank(0), get_neighbor_proc_rank(1), get_neighbor_proc_rank(2), get_neighbor_proc_rank(3)})
{}

  protected:
    void pack_send_buffers()
    {
      double t_start = MPI_Wtime();
      for (int bndry=0; bndry < 4; ++bndry)
      {
        stk::CommBuffer& buf = m_commSparse.send_buffer(get_neighbor_proc_rank(bndry));
        for (int i=0; i < m_nptsPerDirection; ++i)
          for (int g=0; g < m_nghost; ++g)
            buf.pack<double>(get_send_buffer_value(bndry, i, g));
      }
      m_timer.tPack += MPI_Wtime() - t_start;
    }

    stk::CommSparse m_commSparse;
    std::vector<std::vector<double>> m_sendBufs;
    std::vector<std::vector<double>> m_recvBufs;
    std::array<int, 4> m_neighborRanks;
};

// unknown comm plan, no unpacking function
class HeatEquationHaloCommSparseUnknown : public HeatEquationHaloCommSparseBase
{
  public:
    explicit HeatEquationHaloCommSparseUnknown(const RunConfig& config) :
      HeatEquationHaloCommSparseBase(config)
    {}

  protected:
    void start_communication() override
    {
      const bool deallocateSendBuffers = false;
      m_commSparse.reset_buffers();
      for (int phase=0; phase < 2; ++phase)
      {
        pack_send_buffers();

        if (phase == 0)
          m_commSparse.allocate_buffers();
        else
          m_commSparse.communicate(deallocateSendBuffers);
      }
    }

    void finish_communication() override
    {
      for (int i=0; i < 4; ++i)
      {
        m_recvBufs[i].resize(0);
        stk::CommBuffer& buf = m_commSparse.recv_buffer(get_neighbor_proc_rank(i));
        while (buf.remaining())
        {
          double val;
          buf.unpack<double>(val);
          m_recvBufs[i].push_back(val);
        }

        unpack_receive_buffer(i, m_recvBufs[i]);
      }

    }
};


// known comm plan, no unpacking function
class HeatEquationHaloCommSparseKnown : public HeatEquationHaloCommSparseBase
{
  public:
    explicit HeatEquationHaloCommSparseKnown(const RunConfig& config) :
      HeatEquationHaloCommSparseBase(config)
    {}

  protected:
    void start_communication() override
    {
      m_commSparse.reset_buffers();
      std::vector<int> procs_vec(m_neighborRanks.begin(), m_neighborRanks.end());
      for (int phase=0; phase < 2; ++phase)
      {
        pack_send_buffers();

        if (phase == 0)
          m_commSparse.allocate_buffers(procs_vec, procs_vec);
        else
          m_commSparse.communicate();
      }
    }

    void finish_communication() override
    {
      for (int i=0; i < 4; ++i)
      {
        m_recvBufs[i].resize(0);
        stk::CommBuffer& buf = m_commSparse.recv_buffer(get_neighbor_proc_rank(i));
        while (buf.remaining())
        {
          double val;
          buf.unpack<double>(val);
          m_recvBufs[i].push_back(val);
        }

        unpack_receive_buffer(i, m_recvBufs[i]);
      }

    }
};

}  // namespace


TEST(HeatEquationHaloTest, Exactness)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    //ExactSolutionPoly exSol(1);
    auto exSol = std::make_shared<ExactSolutionPoly>(1);
    int stencilSize = 5;
    int nptsPerDirection = 10;
    int nghost = 4;
    double xMin = 0;
    double xMax = 1;
    double yMin = 0;
    double yMax = 1;

    RunConfig config{MPI_COMM_WORLD, nptsPerDirection, nghost, stencilSize, exSol, xMin, xMax, yMin, yMax, 0.0};

    HeatEquationHaloSerial tester(config);

    auto timer = tester.run(1);
    tester.test_residual();

    print_timing_stats(std::cout, timer);
  }
}


TEST(HeatEquationHaloTest, PerformanceEfficient)
{
    HeatEquationHaloEfficient tester(get_performance_run_config());
    TimingStats timer;

    {
      TimerWrapper wrapper(1);
      timer = tester.run(1000);
    }

    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
      print_timing_stats(std::cout, timer);

#ifdef WRITE_TIMING_FILES
    write_timing_files(timer, MPI_COMM_WORLD, "efficient");
#endif
}


TEST(HeatEquationHaloTest, PerformanceParallelDataExchangeBlocking)
{
    auto config = get_performance_run_config();
    HeatEquationHaloParallelDataExchangeBlocking tester(config);
    TimingStats timer;

    {
      TimerWrapper wrapper(1);
      timer = tester.run(1000);
    }

    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
    {
      std::cout << "buffer size = " << config.nghost * config.nptsPerDirection * sizeof(double) / 1024 << "Kib" << std::endl;
      print_timing_stats(std::cout, timer);
    }

#ifdef WRITE_TIMING_FILES
    write_timing_files(timer, MPI_COMM_WORLD, "parallel_data_exchange_blocking");
#endif
}


TEST(HeatEquationHaloTest, PerformanceParallelDataExchangeNonBlocking)
{
    HeatEquationHaloParallelDataExchangeNonBlocking tester(get_performance_run_config());
    TimingStats timer;

    {
      TimerWrapper wrapper(1);
      timer = tester.run(1000);
    }

    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
      print_timing_stats(std::cout, timer);

#ifdef WRITE_TIMING_FILES
    write_timing_files(timer, MPI_COMM_WORLD, "parallel_data_exchange_non_blocking");
#endif
}

TEST(HeatEquationHaloTest, PerformanceCommSparseUnknown)
{
    HeatEquationHaloCommSparseUnknown tester(get_performance_run_config());
    TimingStats timer;

    {
      TimerWrapper wrapper(1);
      timer = tester.run(1000);
    }

    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
      print_timing_stats(std::cout, timer);

#ifdef WRITE_TIMING_FILES
    write_timing_files(timer, MPI_COMM_WORLD, "comm_sparse_unknown");
#endif
}

TEST(HeatEquationHaloTest, PerformanceCommSparseKnown)
{
    HeatEquationHaloCommSparseUnknown tester(get_performance_run_config());
    TimingStats timer;

    {
      TimerWrapper wrapper(1);
      timer = tester.run(1000);
    }

    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
      print_timing_stats(std::cout, timer);

#ifdef WRITE_TIMING_FILES
    write_timing_files(timer, MPI_COMM_WORLD, "comm_sparse_known");
#endif
}
