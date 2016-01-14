#ifndef KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
#define KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP

#include <TeuchosKokkosCompat_config.h>

#include <KokkosCompat_View.hpp>
#include <Kokkos_Core.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef KOKKOS_HAVE_CUDA
  #ifndef KERNEL_PREFIX
    #ifdef __CUDACC__
    #define KERNEL_PREFIX __host__ __device__
    #endif
  #endif
#endif

namespace Kokkos {
namespace Compat {

namespace Details {

/// \brief Get the value of the "Verbose" parameter as a \c bool.
///
/// This method lets the "Verbose" parameter have type either \c int
/// or \c bool, and returns its value as \c bool.  If the "Verbose"
/// parameter does not exist in the given list, return the default
/// value, which is \c false.
bool
getVerboseParameter (const Teuchos::ParameterList& params);

Teuchos::ParameterList getDefaultNodeParameters ();

} // namespace Details

/// \brief Node that wraps a new Kokkos execution space.
/// \tparam ExecutionSpace The type of the Kokkos execution space to wrap.
/// \tparam MemorySpace The Kokkos memory space in which to work.
///   Defaults to the default memory space of ExecutionSpace.
/// \ingroup kokkos_node_api
template<class ExecutionSpace,
         class MemorySpace = typename ExecutionSpace::memory_space>
class KokkosDeviceWrapperNode {
public:
  typedef ExecutionSpace execution_space;
  typedef MemorySpace memory_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  /// \brief This is NOT a "classic" Node type.
  ///
  /// We will deprecate the "classic" Node types with the 11.14
  /// release of Trilinos, and remove them entirely with the 12.0
  /// release.  This Node type is safe to use.
  static const bool classic = false;

  static int count;

  KokkosDeviceWrapperNode (Teuchos::ParameterList& params) {
    if (count == 0) {
      int curNumThreads = -1; // -1 says "let Kokkos pick"
      if (params.isType<int> ("Num Threads")) {
        curNumThreads = params.get<int> ("Num Threads");
      }
      int curNumNUMA = -1; // -1 says "let Kokkos pick"
      if (params.isType<int> ("Num NUMA")) {
        curNumNUMA = params.get<int> ("Num NUMA");
      }
      int curNumCoresPerNUMA = -1; // -1 says "let Kokkos pick"
      if (params.isType<int> ("Num CoresPerNUMA")) {
        curNumCoresPerNUMA = params.get<int> ("Num CoresPerNUMA");
      }
      int curDevice = -1; // -1 says "let Kokkos pick"
      if (params.isType<int> ("Device")) {
        curDevice = params.get<int> ("Device");
      }
      const bool verbose = Details::getVerboseParameter (params);

      if (verbose) {
        std::ostream& out = std::cout;
        out << "DeviceWrapperNode with ExecutionSpace = "
            << typeid (ExecutionSpace).name () << " initializing with "
            << "\"Num Threads\" = " << curNumThreads
            << ", \"Num NUMA\" = " << curNumNUMA
            << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
            << " \"Device\" = " << curDevice << std::endl;
      }
      init (curNumThreads, curNumNUMA, curNumCoresPerNUMA, curDevice);
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! ExecutionSpace::is_initialized (), std::logic_error,
        "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace = "
        << typeid (ExecutionSpace).name () << "> constructor: Kokkos execution "
        << "space initialization failed.");
    }
    count++;
  }

  KokkosDeviceWrapperNode () {
    if (count == 0) {
      const int curNumThreads = -1; // -1 means "let Kokkos pick"
      const int curNumNUMA = -1;
      const int curNumCoresPerNUMA = -1;
      const int curDevice = 0;

      init (curNumThreads, curNumNUMA, curNumCoresPerNUMA, curDevice);
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! ExecutionSpace::is_initialized (), std::logic_error,
        "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace = "
        << typeid (ExecutionSpace).name () << "> constructor: Kokkos execution "
        << "space initialization failed.");
    }
    count++;
  }

  ~KokkosDeviceWrapperNode ();

  static Teuchos::ParameterList getDefaultParameters ()
  {
    return Details::getDefaultNodeParameters ();
  }

  void init (int numthreads, int numnuma, int numcorespernuma, int device);

  void sync () const { ExecutionSpace::fence (); }

  /// \brief Return the human-readable name of this Node.
  ///
  /// See \ref kokkos_node_api "Kokkos Node API"
  static std::string name();
};

#ifdef KOKKOS_HAVE_CUDA
  typedef KokkosDeviceWrapperNode<Kokkos::Cuda> KokkosCudaWrapperNode;
#endif

#ifdef KOKKOS_HAVE_OPENMP
  typedef KokkosDeviceWrapperNode<Kokkos::OpenMP> KokkosOpenMPWrapperNode;
#endif

#ifdef KOKKOS_HAVE_PTHREAD
  typedef KokkosDeviceWrapperNode<Kokkos::Threads> KokkosThreadsWrapperNode;
#endif

#ifdef KOKKOS_HAVE_SERIAL
  typedef KokkosDeviceWrapperNode<Kokkos::Serial> KokkosSerialWrapperNode;
#endif // KOKKOS_HAVE_SERIAL

} // namespace Compat
} // namespace Kokkos
#endif
