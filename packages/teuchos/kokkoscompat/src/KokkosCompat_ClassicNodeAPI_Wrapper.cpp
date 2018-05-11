#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "KokkosCompat_Details_KokkosInit.hpp"
#include <Kokkos_Core.hpp>
#include <iostream>

namespace Kokkos {
  namespace Compat {

#ifdef KOKKOS_ENABLE_THREADS
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    KokkosDeviceWrapperNode ()
    {
      KokkosCompat::Details::initializeKokkos ();
    }

    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    KokkosDeviceWrapperNode (Teuchos::ParameterList& params)
    {
      using ::KokkosCompat::Details::getNodeParameters;

      int curNumThreads = -1; // -1 says "let Kokkos pick"
      int curNumNUMA = -1; // -1 says "let Kokkos pick"
      int curNumCoresPerNUMA = -1; // -1 says "let Kokkos pick"
      int curDevice = 0; // -1 does NOT say "let Kokkos pick"
      bool verbose = false;
      getNodeParameters (curNumThreads, curNumNUMA, curNumCoresPerNUMA,
                         curDevice, verbose, params);
      if (verbose) {
        std::ostream& out = std::cout;
        out << "DeviceWrapperNode with ExecutionSpace = Kokkos::Threads "
            << " initializing with "
            << "\"Num Threads\" = " << curNumThreads
            << ", \"Num NUMA\" = " << curNumNUMA
            << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
            << " \"Device\" = " << curDevice << std::endl;
      }

      Kokkos::InitArguments args;
      args.num_threads = curNumThreads;
      args.num_numa = curNumNUMA;
      (void) curNumCoresPerNUMA;
      (void) curDevice;
      // initializes at most once
      KokkosCompat::Details::initializeKokkos (args);
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    getDefaultParameters ()
    {
      using ::KokkosCompat::Details::getDefaultNodeParameters;
      return getDefaultNodeParameters ();
    }

    template<>
    void
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    sync () const
    {
      Kokkos::Threads::fence ();
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    KokkosDeviceWrapperNode ()
    {
      KokkosCompat::Details::initializeKokkos ();
    }

    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    KokkosDeviceWrapperNode (Teuchos::ParameterList& params)
    {
      using ::KokkosCompat::Details::getNodeParameters;

      int curNumThreads = -1; // -1 says "let Kokkos pick"
      int curNumNUMA = -1; // -1 says "let Kokkos pick"
      int curNumCoresPerNUMA = -1; // -1 says "let Kokkos pick"
      int curDevice = 0; // -1 does NOT say "let Kokkos pick" for Cuda Devices
      bool verbose = false;

      getNodeParameters (curNumThreads, curNumNUMA, curNumCoresPerNUMA,
                         curDevice, verbose, params);
      if (verbose) {
        std::ostream& out = std::cout;
        out << "DeviceWrapperNode with ExecutionSpace = Kokkos::OpenMP "
            << " initializing with "
            << "\"Num Threads\" = " << curNumThreads
            << ", \"Num NUMA\" = " << curNumNUMA
            << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
            << " \"Device\" = " << curDevice << std::endl;
      }

      Kokkos::InitArguments args;
      args.num_threads = curNumThreads;
      args.num_numa = curNumNUMA;
      (void) curNumCoresPerNUMA;
      (void) curDevice;
      // initializes at most once
      KokkosCompat::Details::initializeKokkos (args);
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    getDefaultParameters ()
    {
      using ::KokkosCompat::Details::getDefaultNodeParameters;
      return getDefaultNodeParameters ();
    }

    template<>
    void
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    sync () const
    {
      Kokkos::OpenMP::fence ();
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }
#endif

#ifdef KOKKOS_ENABLE_SERIAL
    template<>
    KokkosDeviceWrapperNode<Kokkos::Serial>::
    KokkosDeviceWrapperNode ()
    {
      KokkosCompat::Details::initializeKokkos ();
    }

    template<>
    KokkosDeviceWrapperNode<Kokkos::Serial>::
    KokkosDeviceWrapperNode (Teuchos::ParameterList& params)
    {
      using ::KokkosCompat::Details::getNodeParameters;

      int curNumThreads = -1; // -1 says "let Kokkos pick"
      int curNumNUMA = -1; // -1 says "let Kokkos pick"
      int curNumCoresPerNUMA = -1; // -1 says "let Kokkos pick"
      int curDevice = 0; // -1 does NOT say "let Kokkos pick" for Cuda Devices
      bool verbose = false;

      getNodeParameters (curNumThreads, curNumNUMA, curNumCoresPerNUMA,
                         curDevice, verbose, params);
      if (verbose) {
        std::ostream& out = std::cout;
        out << "DeviceWrapperNode with ExecutionSpace = Kokkos::Serial "
            << " initializing with "
            << "\"Num Threads\" = " << curNumThreads
            << ", \"Num NUMA\" = " << curNumNUMA
            << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
            << " \"Device\" = " << curDevice << std::endl;
      }
      (void) curNumThreads;
      (void) curNumNUMA;
      (void) curNumCoresPerNUMA;
      (void) curDevice;
      KokkosCompat::Details::initializeKokkos ();
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Serial>::
    getDefaultParameters ()
    {
      using ::KokkosCompat::Details::getDefaultNodeParameters;
      return getDefaultNodeParameters ();
    }

    template<>
    void
    KokkosDeviceWrapperNode<Kokkos::Serial>::
    sync () const
    {
      Kokkos::Serial::fence ();
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_CUDA
    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::
    KokkosDeviceWrapperNode ()
    {
      KokkosCompat::Details::initializeKokkos ();
    }

    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::
    KokkosDeviceWrapperNode (Teuchos::ParameterList& params)
    {
      using ::KokkosCompat::Details::getNodeParameters;

      int curNumThreads = -1; // -1 says "let Kokkos pick"
      int curNumNUMA = -1; // -1 says "let Kokkos pick"
      int curNumCoresPerNUMA = -1; // -1 says "let Kokkos pick"
      int curDevice = 0; // -1 does NOT say "let Kokkos pick" for Cuda Devices
      bool verbose = false;

      getNodeParameters (curNumThreads, curNumNUMA, curNumCoresPerNUMA,
                         curDevice, verbose, params);
      if (verbose) {
        std::ostream& out = std::cout;
        out << "DeviceWrapperNode with ExecutionSpace = Kokkos::Cuda "
            << " initializing with "
            << "\"Num Threads\" = " << curNumThreads
            << ", \"Num NUMA\" = " << curNumNUMA
            << ", \"Num CoresPerNUMA\" = " << curNumCoresPerNUMA
            << " \"Device\" = " << curDevice << std::endl;
      }

      ::KokkosCompat::Details::setUpEnvironmentForCuda ();
      Kokkos::InitArguments args;
      args.num_threads = curNumThreads;
      args.num_numa = curNumNUMA;
      (void) curNumCoresPerNUMA;
      args.device_id = curDevice;
      KokkosCompat::Details::initializeKokkos (args);
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Cuda>::
    getDefaultParameters ()
    {
      using ::KokkosCompat::Details::getDefaultNodeParameters;
      return getDefaultNodeParameters ();
    }

    template<>
    void
    KokkosDeviceWrapperNode<Kokkos::Cuda>::
    sync () const
    {
      Kokkos::Cuda::fence ();
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }
#endif // KOKKOS_ENABLE_CUDA

  } // namespace Compat
} // namespace Kokkos



