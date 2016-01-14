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

  KokkosDeviceWrapperNode (Teuchos::ParameterList &pl) {
    Teuchos::ParameterList params = getDefaultParameters ();
    params.setParameters (pl);
    const int curNumThreads = params.get<int> ("Num Threads");
    const int curNumNUMA = params.get<int> ("Num NUMA");
    const int curNumCoresPerNUMA = params.get<int> ("Num CoresPerNUMA");
    const int curDevice = params.get<int> ("Device");
    const bool verbose = getVerboseParameter (params);

    if (verbose) {
      std::ostream& out = std::cout;
      out << "DeviceWrapperNode with ExecutionSpace = "
          << typeid (ExecutionSpace).name () << " initializing with "
          << "\"Num Threads\" = " << curNumThreads
          << ", \"Num NUMA\" = " << curNumNUMA
          << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
          << " \"Device\" = " << curDevice << std::endl;
    }
    if (count == 0) {
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
    Teuchos::ParameterList params = getDefaultParameters ();
    const int curNumThreads = params.get<int> ("Num Threads");
    const int curNumNUMA = params.get<int> ("Num NUMA");
    const int curNumCoresPerNUMA = params.get<int> ("Num CoresPerNUMA");
    const int curDevice = params.get<int> ("Device");
    const bool verbose = getVerboseParameter (params);

    if (verbose) {
      std::ostream& out = std::cout;
      out << "DeviceWrapperNode with ExecutionSpace = "
          << typeid (ExecutionSpace).name () << " initializing with "
          << "\"Num Threads\" = " << curNumThreads
          << ", \"Num NUMA\" = " << curNumNUMA
          << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
          << " \"Device\" = " << curDevice << std::endl;
    }
    if (count == 0) {
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

private:
  /// \brief Get the value of the "Verbose" parameter as a \c bool.
  ///
  /// This method lets the "Verbose" parameter have type either \c int
  /// or \c bool, and returns its value as \c bool.  If the "Verbose"
  /// parameter does not exist in the given list, return the default
  /// value, which is \c false.
  static bool
  getVerboseParameter (const Teuchos::ParameterList& params)
  {
    const bool defaultValue = false; // default value of the parameter

    try { // Is it a bool?
      return params.get<bool> ("Verbose");
    }
    catch (...) {}

    try { // Is it an int?
      return params.get<int> ("Verbose");
    }
    catch (...) {}

    return defaultValue;
  }

public:

  static Teuchos::ParameterList getDefaultParameters ()
  {
    Teuchos::ParameterList params;
    params.set ("Verbose", 0);
    // -1 says "Let Kokkos pick"
    params.set ("Num Threads", -1);
    params.set ("Num NUMA", -1);
    params.set ("Num CoresPerNUMA", -1);
    params.set ("Device", 0);
    return params;
  }

  void init (int numthreads, int numnuma, int numcorespernuma, int device);

  void sync () const { ExecutionSpace::fence (); }

  //@{
  //! \name Memory management

  //! \brief Return a const view of a buffer for use on the host.
  template <class T> inline
  Teuchos::ArrayRCP<const T> viewBuffer(size_t size, Teuchos::ArrayRCP<const T> buff) {
    return buff.persistingView(0,size);
  }

  /// \brief Return a non-const view of a buffer for use on the host.
  ///
  /// \param rw [in] 0 if read-only, 1 if read-write.  This is an
  ///   <tt>int</tt> and not a KokkosClassic::ReadWriteOption, in
  ///   order to avoid a circular dependency between KokkosCompat
  ///   and TpetraClassic.
  template <class T> inline
  Teuchos::ArrayRCP<T>
  viewBufferNonConst (const int rw, size_t size, const Teuchos::ArrayRCP<T> &buff) {
    (void) rw; // Silence "unused parameter" compiler warning
    return buff.persistingView(0,size);
  }

  /// \brief Return the human-readable name of this Node.
  ///
  /// See \ref kokkos_node_api "Kokkos Node API"
  static std::string name();

  //@}
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
