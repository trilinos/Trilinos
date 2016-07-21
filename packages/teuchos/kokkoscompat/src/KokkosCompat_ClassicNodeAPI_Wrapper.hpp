#ifndef KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
#define KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP

#include "TeuchosKokkosCompat_config.h"
#include "KokkosCompat_View.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_ParameterList.hpp"

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
  //! The Node's Kokkos execution space.
  typedef ExecutionSpace execution_space;
  //! The Node's Kokkos memory space.
  typedef MemorySpace memory_space;
  /// \brief The Node's Kokkos::Device specialization.
  ///
  /// This is just an (execution space, memory space) pair.
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  /// \brief This is NOT a "classic" Node type.
  ///
  /// We will deprecate the "classic" Node types with the 11.14
  /// release of Trilinos, and remove them entirely with the 12.0
  /// release.  This Node type is safe to use.
  static const bool classic = false;

  /// \brief Reference count of Node instances.
  ///
  /// \warning Do NOT read or modify this.  This is an implementation
  ///   detail of the class.
  ///
  /// If Node originally called ExecutionSpace::initialize() (see
  /// Github Issue #510), and if the count reaches zero in the
  /// destructor, we call ExecutionSpace::finalize().
  static int count;

  /// \brief Is the Node responsible for finalizing its execution space?
  ///
  /// In Node's constructor, was it ever true that count == 0 and
  /// execution_space::is_initialized()?  If so, this means that
  /// somebody else (the user, or some other library) initialized the
  /// execution space before Node got to it.  In that case, the Node's
  /// destructor should NOT call (is not responsible for calling)
  /// execution_space::finalize().
  ///
  /// This fixes Github Issue #510.
  static bool nodeResponsibleForFinalizingExecutionSpace_;

  /// \brief Constructor (that takes a Teuchos::ParameterList).
  ///
  /// \param [in/out] params List of Node configuration parameters.
  ///   If empty, we use defaults.
  KokkosDeviceWrapperNode (Teuchos::ParameterList& params) {
    // Fix for Github Issue #510.  If count is zero, yet the execution
    // space is already initialized, then Node is not responsible for
    // finalizing the execution space.
    if (count == 0 && ExecutionSpace::is_initialized ()) {
      nodeResponsibleForFinalizingExecutionSpace_ = false;
    }
#ifdef KOKKOS_HAVE_CUDA
    // The Cuda Node also handles its host execution space.
    if (count == 0 &&
        std::is_same<ExecutionSpace, ::Kokkos::Cuda>::value &&
        ::Kokkos::HostSpace::execution_space::is_initialized ()) {
      KokkosDeviceWrapperNode< ::Kokkos::HostSpace::execution_space>::nodeResponsibleForFinalizingExecutionSpace_ = false;
    }
#endif // KOKKOS_HAVE_CUDA

    // Kokkos insists that if Kokkos::Cuda is initialized, then its
    // host execution space must also be initialized.  Thus, it
    // suffices to check whether Kokkos::Cuda has been initialized; we
    // don't also have to check the host execution space.
    if (count == 0 && nodeResponsibleForFinalizingExecutionSpace_) {
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
      int curDevice = 0; // -1 does NOT say "let Kokkos pick" for Cuda Devices
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
    }
    count++;
    checkConstructorEnd ();
  }

  //! Default constructor (sets default parameters).
  KokkosDeviceWrapperNode () {
    if (count == 0 && nodeResponsibleForFinalizingExecutionSpace_) {
      const int curNumThreads = -1; // -1 means "let Kokkos pick"
      const int curNumNUMA = -1;
      const int curNumCoresPerNUMA = -1;
      const int curDevice = 0;

      init (curNumThreads, curNumNUMA, curNumCoresPerNUMA, curDevice);
    }
    count++;
    checkConstructorEnd ();
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

private:
  //! Make sure that the constructor initialized everything.
  void checkConstructorEnd () const {
    TEUCHOS_TEST_FOR_EXCEPTION
      (count <= 0, std::logic_error,
       "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace="
       << typeid (ExecutionSpace).name () << ", MemorySpace="
       << typeid (MemorySpace).name () << " > constructor: "
       "count = " << count << " <= 0 at end of constructor.  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! ExecutionSpace::is_initialized (), std::logic_error,
       "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace="
       << typeid (ExecutionSpace).name () << ", MemorySpace="
       << typeid (MemorySpace).name () << " > constructor: "
       "Failed to initialize ExecutionSpace.");
#ifdef KOKKOS_HAVE_CUDA
    // This must be outside the macro, since it has a comma.
    constexpr bool isCuda =
      std::is_same<ExecutionSpace, Kokkos::Cuda>::value;
    TEUCHOS_TEST_FOR_EXCEPTION
      (isCuda && ! Kokkos::HostSpace::execution_space::is_initialized (),
       std::logic_error,
       "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace=Kokkos::Cuda, "
       "MemorySpace=" << typeid (MemorySpace).name () << " > constructor: "
       "Failed to initialize Kokkos::HostSpace::execution_space.");
#endif // KOKKOS_HAVE_CUDA
  }

  //! Make sure that the destructor did its job.
  void checkDestructorEnd () const {
    TEUCHOS_TEST_FOR_EXCEPTION
      (count < 0, std::logic_error,
       "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace="
       << typeid (ExecutionSpace).name () << ", MemorySpace="
       << typeid (MemorySpace).name () << " > destructor: "
       "count = " << count << " < 0 at end of destructor.  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (count == 0 &&
       nodeResponsibleForFinalizingExecutionSpace_ &&
       ExecutionSpace::is_initialized (),
       std::logic_error,
       "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace="
       << typeid (ExecutionSpace).name () << ", MemorySpace="
       << typeid (MemorySpace).name () << " > destructor: "
       "count == 0 and the Node is responsible for finalizing ExecutionSpace, "
       "but its destructor did not.  "
       "Please report this bug to the Tpetra developers.");
#ifdef KOKKOS_HAVE_CUDA
    // This must be outside the macro, since it has a comma.
    constexpr bool isCuda =
      std::is_same<ExecutionSpace, Kokkos::Cuda>::value;
    TEUCHOS_TEST_FOR_EXCEPTION
      (isCuda &&
       count == 0 &&
       KokkosDeviceWrapperNode<Kokkos::HostSpace::execution_space>::nodeResponsibleForFinalizingExecutionSpace_ &&
       Kokkos::HostSpace::execution_space::is_initialized (),
       std::logic_error,
       "Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace=Kokkos::Cuda, "
       "MemorySpace=" << typeid (MemorySpace).name () << " > destructor: "
       "count == 0 and the Node is responsible for finalizing "
       "Kokkos::HostSpace::execution_space, but its destructor did not.  "
       "Please report this bug to the Tpetra developers.");
#endif // KOKKOS_HAVE_CUDA
  }
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
