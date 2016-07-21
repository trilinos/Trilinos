#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_Core.hpp>
#include <cstdlib> // getenv
#include <sstream>

namespace { // (anonymous)

  struct CmdLineArgs {
    int numThreads;
    int numNuma;
    int deviceId;
    int numDevices;

    CmdLineArgs () : // Default values let Kokkos decide
      numThreads (-1), numNuma (-1), deviceId (-1), numDevices (-1)
    {}
  };

  CmdLineArgs
  getCommandLineArgs ()
  {
    using std::cerr;
    using std::endl;

    const bool debug = false;
    if (debug) {
      cerr << "Tpetra: get command-line arguments for Kokkos wrapper Nodes"
           << endl;
    }

    CmdLineArgs argsOut;
    // Some flag value that Kokkos::initialize uses; I just imitated
    // its logic.
    int skip_device = 9999;

    // Attempt to read command-line arguments that were stored in
    // Teuchos::GlobalMPISession. User settings override these.

    const auto argv = ::Teuchos::GlobalMPISession::getArgv ();
    bool gotNumDevices = false;

    // Loop over all command-line arguments.  If the same argument
    // name occurs multiple times, we read all of them and use the
    // last provided value.
    const int narg = static_cast<int> (argv.size ());
    for (int iarg = 0; iarg < narg; ++iarg) {
      // getArgv() only promises to return an array whose entries are
      // convertible to std::string.  Thus, it's not legit to assign
      // to const std::string&; we actually must create an std::string
      // from the array entry.
      const std::string curArg = argv[iarg];

      if (debug) {
        cerr << "  Current command-line argument: " << curArg << endl;
      }

      // Find "--" and "=".  These two together, "--" before "=", with
      // something after "=", indicate a valid command-line argument.
      // It might be a Kokkos one, or it might be something else.
      const size_t posDash = curArg.find ("--");
      if (posDash == std::string::npos) {
        continue; // argument name must start with "--"
      }
      const size_t posEqual = curArg.find ("=");
      if (posEqual == std::string::npos ||
          ! (posDash + static_cast<size_t> (2) < posEqual) ||
          posEqual + static_cast<size_t> (1) == curArg.size ()) {
        // Command-line argument must contain "=",
        // "--" must come before "=" with >= 1 characters in between,
        // and >= 1 characters must follow "=".
        continue;
      }

      // Does the command-line argument start with "--kokkos-"?  If
      // so, it's a Kokkos-specific command-line argument.  Otherwise,
      // Kokkos still reads it.
      const size_t posKokkos = curArg.find ("kokkos-", posDash+2);
      const size_t posArgStart =
        (posKokkos == std::string::npos) ? (posDash+2) : (posKokkos+7);

      // Figure out _which_ Kokkos argument we got.  Put "ndevices" in
      // front of "device", since the latter is a substring of the
      // former.
      const char* validArgNames[4] = {"threads", "numa", "ndevices", "device"};
      const int numValidArgNames = 4;
      const int numDevicesInd = 2;

      int argNameInd = 0; // need to save this for below
      for ( ; argNameInd < numValidArgNames; ++argNameInd) {
        const size_t posArgName =
          curArg.find (validArgNames[argNameInd], posArgStart);
        if (posArgName != std::string::npos) {
          if (argNameInd == numDevicesInd) { // "ndevices" is a special case
            gotNumDevices = true;
          }
          break;
        }
      }

      // Get the one or two numbers that follow the "=".  "ndevices"
      // is the only argument take takes two numbers (separated by
      // ",").  All other cases take just one number.  For simplicity,
      // we always look for both numbers, and ignore the second one if
      // the command-line option doesn't use it.

      int firstNum = 0;
      int secondNum = 0;
      bool gotFirstNum = false;
      bool gotSecondNum = false;

      const size_t posComma = curArg.find (",", posEqual+1);
      if (posComma == std::string::npos) { // no comma after "="
        // "--kokkos-$NAME=$NUMBER1" -- just read first number $NUMBER1
        const std::string arg1 = curArg.substr (posEqual + 1, std::string::npos);
        if (arg1.size () != 0) {
          std::istringstream istr1 (arg1);
          bool threw = false;
          try {
            istr1 >> firstNum;
          }
          catch (...) {
            threw = true;
            gotFirstNum = false;
          }
          if (! threw) {
            // Yes, I could write ! (! istr1), but that looks ridiculous.
            gotFirstNum = ! istr1.bad () && ! istr1.fail ();
          }
        }
      }
      else if (posComma == posEqual + 1) { // comma right after "="
        // "--kokkos-$NAME=,$NUMBER2" -- just read second number $NUMBER2
        const std::string arg2 = curArg.substr (posComma + 1, std::string::npos);
        if (arg2.size () != 0) {
          std::istringstream istr2 (arg2);
          bool threw = false;
          try {
            istr2 >> secondNum;
          }
          catch (...) {
            threw = true;
            gotSecondNum = false;
          }
          if (! threw) {
            // Yes, I could write ! (! istr2), but that looks ridiculous.
            gotSecondNum = ! istr2.bad () && ! istr2.fail ();
          }
        }
      }
      else { // comma follows "=", with >= 1 intervening characters
        // "--kokkos-$NAME=$NUMBER1,$NUMBER2" -- read both numbers
        const std::string arg1 = curArg.substr (posEqual + 1, std::string::npos);
        if (arg1.size () != 0) {
          std::istringstream istr1 (arg1);
          bool threw = false;
          try {
            istr1 >> firstNum;
          }
          catch (...) {
            threw = true;
            gotFirstNum = false;
          }
          if (! threw) {
            // Yes, I could write ! (! istr1), but that looks ridiculous.
            gotFirstNum = ! istr1.bad () && ! istr1.fail ();
          }
        }
        const std::string arg2 = curArg.substr (posComma + 1, std::string::npos);
        if (arg2.size () != 0) {
          std::istringstream istr2 (arg2);
          bool threw = false;
          try {
            istr2 >> secondNum;
          }
          catch (...) {
            threw = true;
            gotSecondNum = false;
          }
          if (! threw) {
            // Yes, I could write ! (! istr2), but that looks ridiculous.
            gotSecondNum = ! istr2.bad () && ! istr2.fail ();
          }
        }
      }

      if (argNameInd == 0) {
        argsOut.numThreads = firstNum;
      }
      else if (argNameInd == 1) {
        argsOut.numNuma = firstNum;
      }
      else if (argNameInd == 2) {
        if (gotFirstNum) {
          argsOut.numDevices = firstNum;
        }
        else {
          continue; // TODO REPORT ERROR
        }
        if (gotSecondNum) {
          skip_device = secondNum;
        }

        argsOut.numDevices = secondNum;
        gotNumDevices = true;
      }
      else if (argNameInd == 3) {
        argsOut.deviceId = firstNum;
      }
      else if (argNameInd == numValidArgNames) {
        continue; // name wasn't in the list of valid names
      }
    } // for each command-line argument

    // If Kokkos::initialize doesn't get the --kokkos-ndevices option,
    // it tries to use environment variables to figure this out.
    if (! gotNumDevices) {
      char* str = getenv ("SLURM_LOCALID");
      if (str == NULL) {
        str = getenv ("MV2_COMM_WORLD_LOCAL_RANK");
      }
      if (str == NULL) {
        str = getenv ("OMPI_COMM_WORLD_LOCAL_RANK");
      }

      if (str != NULL) {
        const std::string sstr (str);
        std::istringstream istr (sstr);

        int localRank = -1;
        bool gotLocalRank = false;
        try {
          istr >> localRank;
          gotLocalRank = true;
        }
        catch (...) {
          gotLocalRank = false;
        }
        if (! istr.bad () && ! istr.fail ()) {
          gotLocalRank = true;
        }

        if (gotLocalRank) {
          argsOut.deviceId = localRank % argsOut.numDevices;
          if (argsOut.deviceId >= skip_device) {
            argsOut.deviceId++;
          }
        }
      }

      if (argsOut.deviceId == -1) {
        argsOut.deviceId = 0;
        // mfh 24 Apr 2016: This looks a little funny (why isn't it
        // outside the argsOut.deviceId == -1 test?), but it imitates
        // the implementation of Kokkos::initialize.  See
        // kokkos/core/src/impl/Kokkos_Core.cpp, lines 368-371.
        if (argsOut.deviceId >= skip_device) {
          argsOut.deviceId++;
        }
      }
    }

    if (debug) {
      cerr << "Tpetra: got command-line arguments: {"
           << "numThreads: " << argsOut.numThreads << ','
           << "numNuma: " << argsOut.numNuma << ','
           << "deviceId: " << argsOut.deviceId << ','
           << "numDevices: " << argsOut.numDevices << '}' << endl;
    }
    return argsOut;
  }
} // namespace (anonymous)

namespace Kokkos {
  namespace Compat {
    namespace Details {
      bool
      getVerboseParameter (const Teuchos::ParameterList& params)
      {
        const bool defaultValue = false; // default value of the parameter

        if (params.isParameter ("Verbose")) {
          if (params.isType<bool> ("Verbose")) { // is it a bool?
            return params.get<bool> ("Verbose");
          }
          else if (params.isType<int> ("Verbose")) { // is it an int?
            return params.get<int> ("Verbose");
          }
          // It might be polite to throw at this point with a helpful
          // message explaining that the parameter has the wrong type,
          // but that would change current behavior, so I'll just
          // leave it.
        }
        return defaultValue;
      }

      Teuchos::ParameterList getDefaultNodeParameters ()
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
    } // namespace Details

    // mfh 01 Jan 2014: These definitions of the class variable count
    // need to be inside the namespace.  Declaring them as "template<>
    // int Kokkos::Compat::KokkosCudaWrapperNode::count = 0" in the
    // global namespace is a C++11 extension and results in compiler
    // warnings with Clang 3.2 on MacOS X.
#ifdef KOKKOS_HAVE_CUDA
    template<> int KokkosCudaWrapperNode::count = 0;
    template<> bool KokkosCudaWrapperNode::nodeResponsibleForFinalizingExecutionSpace_ = true;
#endif
#ifdef KOKKOS_HAVE_OPENMP
    template<> int KokkosOpenMPWrapperNode::count = 0;
    template<> bool KokkosOpenMPWrapperNode::nodeResponsibleForFinalizingExecutionSpace_ = true;
#endif
#ifdef KOKKOS_HAVE_PTHREAD
    template<> int KokkosThreadsWrapperNode::count = 0;
    template<> bool KokkosThreadsWrapperNode::nodeResponsibleForFinalizingExecutionSpace_ = true;
#endif
#ifdef KOKKOS_HAVE_SERIAL
    template<> int KokkosSerialWrapperNode::count = 0;
    template<> bool KokkosSerialWrapperNode::nodeResponsibleForFinalizingExecutionSpace_ = true;
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    ~KokkosDeviceWrapperNode<Kokkos::Threads> ()
    {
      count--;
      // Only call Kokkos::Threads::finalize if the Node reference
      // count is zero, if Kokkos::Threads is already initialized, and
      // if Node's constructor was responsible for initializing it.
      // (See Github Issue #510.)
      if (count == 0 && Threads::is_initialized () && nodeResponsibleForFinalizingExecutionSpace_) {
#ifdef KOKKOS_HAVE_CUDA
        // If any instances of KokkosDeviceWrapperNode<Kokkos::Cuda>
        // exist, they are responsible for calling finalize on Cuda's
        // host execution space.
        if (! Impl::is_same<Kokkos::Threads,HostSpace::execution_space>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif // KOKKOS_HAVE_CUDA
          Threads::finalize ();
      }
      checkDestructorEnd ();
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Threads>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int /* Device */) {
      using std::cerr;
      using std::endl;
      const bool debug = false;

      if (! Kokkos::Threads::is_initialized ()) {

        // Attempt to read command-line arguments that were stored in
        // Teuchos::GlobalMPISession. User settings override these.
        CmdLineArgs args = getCommandLineArgs ();

        if (args.numThreads != -1) {
          NumThreads = args.numThreads;
        }
        if (args.numNuma != -1) {
          NumNUMA = args.numNuma;
        }
        if (args.deviceId != -1) {
          NumCoresPerNUMA = args.deviceId;
        }
        // if (args.numDevice != -1) {
        //   Device = args.numDevice; // Threads doesn't need this one
        // }

        if (NumNUMA > 0 && NumCoresPerNUMA > 0) {
          Kokkos::Threads::initialize (NumThreads, NumNUMA, NumCoresPerNUMA);
        }
        else if (NumNUMA > 0) {
          Kokkos::Threads::initialize (NumThreads, NumNUMA);
        }
        else if (NumThreads > 0) {
          if (debug) {
            cerr << "Tpetra Threads wrapper Node: NumThreads = "
                 << NumThreads << endl;
          }
          Kokkos::Threads::initialize (NumThreads);
          if (debug) {
            cerr << "Tpetra Threads wrapper Node: Concurrency = "
                 << Kokkos::Threads::concurrency () << endl;
          }
        }
        else {
          Kokkos::Threads::initialize ();
        }
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::~KokkosDeviceWrapperNode<Kokkos::OpenMP>() {
      count--;
      // Only call Kokkos::OpenMP::finalize if the Node reference
      // count is zero, if Kokkos::OpenMP is already initialized, and
      // if Node's constructor was responsible for initializing it.
      // (See Github Issue #510.)
      if (count == 0 && OpenMP::is_initialized () && nodeResponsibleForFinalizingExecutionSpace_) {
#ifdef KOKKOS_HAVE_CUDA
        // If any instances of KokkosDeviceWrapperNode<Kokkos::Cuda>
        // exist, they are responsible for calling finalize on Cuda's
        // host execution space.
        if (! Impl::is_same<Kokkos::OpenMP, HostSpace::execution_space>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif // KOKKOS_HAVE_CUDA
        OpenMP::finalize ();
      }
      checkDestructorEnd ();
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::OpenMP::is_initialized ()) {
        if (NumNUMA > 0 && NumCoresPerNUMA > 0) {
          Kokkos::OpenMP::initialize (NumThreads, NumNUMA, NumCoresPerNUMA);
        }
        else if (NumNUMA > 0) {
          Kokkos::OpenMP::initialize (NumThreads, NumNUMA);
        }
        else if (NumThreads > 0) {
          Kokkos::OpenMP::initialize (NumThreads);
        }
        else {
          Kokkos::OpenMP::initialize ();
        }
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }
#endif

#ifdef KOKKOS_HAVE_SERIAL
    template<>
    KokkosDeviceWrapperNode<Kokkos::Serial>::~KokkosDeviceWrapperNode<Kokkos::Serial>() {
      count--;
      // Only call Kokkos::Serial::finalize if the Node reference
      // count is zero, if Kokkos::Serial is already initialized, and
      // if Node's constructor was responsible for initializing it.
      // (See Github Issue #510.)
      if (count == 0 && Serial::is_initialized () && nodeResponsibleForFinalizingExecutionSpace_) {
#ifdef KOKKOS_HAVE_CUDA
        // If any instances of KokkosDeviceWrapperNode<Kokkos::Cuda>
        // exist, they are responsible for calling finalize on Cuda's
        // host execution space.
        if (! Impl::is_same<Kokkos::Serial, HostSpace::execution_space>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif // KOKKOS_HAVE_CUDA
          Serial::finalize ();
      }
      checkDestructorEnd ();
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Serial>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::Serial::is_initialized ()) {
        Kokkos::Serial::initialize ();
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_CUDA
    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::~KokkosDeviceWrapperNode<Kokkos::Cuda>() {
      count--;
      if (count == 0) {
        // Don't call finalize on the host execution space until all
        // Node instances corresponding to the host execution space
        // are gone.  Also, don't call finalize on it if Node wasn't
        // responsible for calling initialize on it (see Github Issue
        // #510).
        if (HostSpace::execution_space::is_initialized () &&
            KokkosDeviceWrapperNode<HostSpace::execution_space>::count == 0 &&
            KokkosDeviceWrapperNode<HostSpace::execution_space>::nodeResponsibleForFinalizingExecutionSpace_) {
          HostSpace::execution_space::finalize ();
        }
        // Only call Kokkos::Cuda::finalize if the Node reference
        // count is zero, if Kokkos::Cuda is already initialized, and
        // if Node's constructor was responsible for initializing it.
        // (See Github Issue #510.)
        if (Cuda::is_initialized () && nodeResponsibleForFinalizingExecutionSpace_) {
          Cuda::finalize ();
        }
      }
      checkDestructorEnd ();
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::
    init(int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {

      // Setting (currently) necessary environment variables for NVIDIA UVM
      #ifdef KOKKOS_USE_CUDA_UVM
        putenv("CUDA_LAUNCH_BLOCKING=1");
      #else
        throw std::runtime_error("Using CudaWrapperNode without UVM is not allowed.");
      #endif

      if (! Kokkos::HostSpace::execution_space::is_initialized ()) {
        if (NumNUMA > 0 && NumCoresPerNUMA > 0) {
          Kokkos::HostSpace::execution_space::initialize (NumThreads, NumNUMA, NumCoresPerNUMA);
        }
        else if (NumNUMA > 0) {
          Kokkos::HostSpace::execution_space::initialize (NumThreads, NumNUMA);
        }
        else if (NumThreads > 0) {
          Kokkos::HostSpace::execution_space::initialize (NumThreads);
        }
        else {
          Kokkos::HostSpace::execution_space::initialize ();
        }
      }
      Kokkos::Cuda::SelectDevice select_device (Device);
      if (! Kokkos::Cuda::is_initialized ()) {
        Kokkos::Cuda::initialize (select_device);
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }
#endif // KOKKOS_HAVE_CUDA

  } // namespace Compat
} // namespace Kokkos



