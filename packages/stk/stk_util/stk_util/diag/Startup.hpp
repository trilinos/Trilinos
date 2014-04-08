#ifndef STK_UTIL_DIAG_Startup_h
#define STK_UTIL_DIAG_Startup_h

#include <mpi.h>                        // for MPI_Comm
#include <stddef.h>                     // for NULL
#include <stk_util/environment/EnvData.hpp>
#include <vector>                       // for vector



namespace sierra {
namespace Env {
  /**
   * @brief Initialize MPI related operations for sierra, outputs banner, etc.
   *        returns 1 if MPI was initialized, 0 otherwise
   *
   * @param argc		an <b>int</b> pointer to the main argc.
   *
   * @param argv		a <b>char</b> pointer to the main argv.
   *
   * @param product_name	a <b>char</b> const pointer to the name of the
   *				product.
   *
   * @param build_date_time	a <b>char</b> const pointer to
   *				__DATE__ " " __TIME__
   *
   * @param mpi_key	        an optional <b>ExecType</b> enumeration
   *                            specifying the type of executable.
   *                            Default is a single executable using
   *                            MPI_COMM_WORLD. Other options result
   *                            in a split communicator.
   *
   * @param peer_sizes          an optional <b>std::vector<int></b>
   *                            const pointer containing the number of
   *                            processors that each peer communicator
   *                            (\sa parrallel_peer_comm()) will
   *                            support. The number of peers is
   *                            determined by the size of the vector.
   *                            Only used if mpi_key = EXEC_TYPE_PEER
   */
 bool StartupSierra(int *argc, char ***argv, const char *product_name, const char *build_date_time,
	           ExecType mpi_key = EXEC_TYPE_WORLD, const std::vector<int> *peer_sizes = NULL,
	           const ExecInfo * preSplitCommunicatorInfo = NULL);

 //
 //  Cleanup any MPI stuff that was created from the startup sierra call, pass in the MPI initialization flag
 //  that StartupSierra generated
 //
 void ShutDownSierra(bool mpiInitFlag);


/**
 * @ingroup EnvMPIDetail EnvOutputDetail EnvCommandLineDetail
 * @brief Class <b>Startup</b> is a sentry class for starting the application.  It
 * ensures that the command line arguments, platform and MPI are ready to go at the start
 * of the application.
 *
 */
class Startup
{
public:

  /**
   * @brief Creates a new <b>Startup</b> instance.
   *
   * @param argc		an <b>int</b> pointer to the main argc.
   *
   * @param argv		a <b>char</b> pointer to the main argv.
   *
   * @param product_name	a <b>char</b> const pointer to the name of the
   *				product.
   *
   * @param build_date_time	a <b>char</b> const pointer to
   *				__DATE__ " " __TIME__
   *
   * @param mpi_key	        an optional <b>ExecType</b> enumeration
   *                            specifying the type of executable.
   *                            Default is a single executable using
   *                            MPI_COMM_WORLD. Other options result
   *                            in a split communicator.
   *
   * @param peer_sizes          an optional <b>std::vector<int></b>
   *                            const pointer containing the number of
   *                            processors that each peer communicator
   *                            (\sa parrallel_peer_comm()) will
   *                            support. The number of peers is
   *                            determined by the size of the vector.
   *                            Only used if mpi_key = EXEC_TYPE_PEER
   */
  Startup(int *argc, char ***argv, const char *product_name, const char *build_date_time,
	  ExecType mpi_key = EXEC_TYPE_WORLD, const std::vector<int> *peer_sizes = NULL,
	  const ExecInfo * preSplitCommunicatorInfo = 0);

  /**
   * @brief Destroys a <b>Startup</b> instance.  IT closes all logging output streams and
   * will finalize MPI only if the Startup::Startup initialized the MPI.
   *
   */
  ~Startup();

private:
  /**
   * @brief Member function <b>startup</b>
   *  @li initializes MPI if not already initialized,
   *  @li queries the environment and command line for environment options,
   *  @li broadcasts the environment options from MPI_COMM_WORLD processor #0 to all other
   *      processors,
   *  @li opens the default parallel output stream.
   *
   * @param argc		an <b>int</b> pointer to the main argc.
   *
   * @param argv		a <b>char</b> pointer to the main argv.
   *
   * @param product_name	a <b>char</b> const pointer to the name of the
   *				product.
   *
   * @param build_time		a <b>char</b> const pointer to __DATE__ " " __TIME__
   *
   * @param mpi_key	        an <b>ExecType</b> 
   *
   * @param peer_sizes          a <b> std::vector<int> </b>
   *
   */
  void startup(int *argc, char ***argv, const char *product_name, const char *build_time,
	       ExecType mpi_key, const std::vector<int> *peer_sizes,
	       const ExecInfo * preSplitCommunicatorInfo);

private:
  bool m_mpiInitFlag;			///< True if Startup initialized MPI
};


/**
 * @ingroup EnvMPIDetail
 * @brief Function <b>reset</b> determines new parallel_size and parallel_rank.
 * Flushes, closes, and reopens log files.
 *
 * This is a very dangerous method, if the previous communicator had been queried and
 * given to an external library then Sierra and that external library will be out-of-sync.
 * It is strongly recommended that parser-interface singletons do not query and save the
 * communicator.
 *
 * This capability is provided to support dakota.
 *
 * @param mpi_comm		a <b>MPI_Comm</b> value for the new mpi
 *				communicator.
 */
void reset(MPI_Comm mpi_comm);

}}
#endif
