#include "Compadre_KokkosParser.hpp"

using namespace Compadre;

// for command line arguments, pass them directly in to Kokkos
KokkosParser::KokkosParser(int narg, char **arg, bool print_status) :
        _num_threads(-1), _numa(-1), _device(-1), _ngpu(-1) {

    // determine if Kokkos is already initialized
    // if it has been, get the parameters needed from it
    bool preinitialized = Kokkos::is_initialized();

    if (!preinitialized) {
        Kokkos::initialize(narg, arg);
        _called_initialize = true;
    }

    // get parameters
    // set parameters
    retrievePreviouslyInstantiatedKokkosInitArguments();

    if (print_status && preinitialized) {
        printf("Kokkos already initialized.");
    }

    if (print_status) {
        // call status
        this->status();
    }

    // MPI 
    //char *str;
    //if ((str = getenv("SLURM_LOCALID"))) {
    //  int local_rank = atoi(str);
    //  device = local_rank % ngpu;
    //  if (device >= skip_gpu) device++;
    //}
    //if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
    //  int local_rank = atoi(str);
    //  device = local_rank % ngpu;
    //  if (device >= skip_gpu) device++;
    //}
    //if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
    //  int local_rank = atoi(str);
    //  device = local_rank % ngpu;
    //  if (device >= skip_gpu) device++;
    //}
}

KokkosParser::KokkosParser(int num_threads, int numa, int device, int ngpu, bool print_status) :
        _num_threads(num_threads), _numa(numa), _device(device), _ngpu(ngpu) {

//#ifdef KOKKOS_HAVE_CUDA
//  // only written for handling one gpu
//  compadre_assert_release((ngpu == 1) && "Only one GPU supported at this time.");
//#else
//  ngpu = 0;
//#endif
//#ifdef KOKKOS_HAVE_CUDA
//  compadre_assert_release((ngpu > 0) && "Kokkos has been compiled for CUDA but no GPUs are requested"); 
//#endif

    // returns 1 if we initialized
    _called_initialize = this->initialize();
  
    if (print_status && _called_initialize==1) {
        printf("KOKKOS attempted initialization using settings:\n");
        this->status();
        retrievePreviouslyInstantiatedKokkosInitArguments();
        printf("KOKKOS recalculated and initiatized using settings:\n");
        this->status();
    } else if (_called_initialize==0) {
        // could be improved by retrieving the parameters Kokkos was initialized with
        // get parameters
        // set parameters
        retrievePreviouslyInstantiatedKokkosInitArguments();

        if (print_status) {
            printf("Kokkos already initialized.");
            this->status();
        }
    } else if (print_status) {
        printf("Kokkos failed to initialize.");
    }

}

Kokkos::InitArguments KokkosParser::createInitArguments() const {
  Kokkos::InitArguments args;
  args.num_threads = _num_threads;
  args.num_numa = _numa;
  args.device_id = _device;
  return args;
}

void KokkosParser::retrievePreviouslyInstantiatedKokkosInitArguments() {
// NUMA parts are not tested, and only work for 1 numa region
#ifdef COMPADRE_USE_CUDA
    //auto cuda_space = Kokkos::DefaultExecutionSpace;
    _device = 0;//cuda_space.cuda_device();
    _ngpu = 1;
    _numa = 0;
    _num_threads = 0;
#else
    _device = 0;
    _ngpu = 0;
#endif

#ifdef KOKKOS_ENABLE_THREADS
    //_num_threads = Kokkos::HostSpace::execution_space::get_current_max_threads();
    _numa = 1;
    _num_threads = Kokkos::HostSpace::execution_space::concurrency();//impl_get_current_max_threads();
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    //_num_threads = Kokkos::HostSpace::execution_space::get_current_max_threads();
    _numa = 1;
    _num_threads = Kokkos::HostSpace::execution_space::concurrency();//impl_get_current_max_threads();
#endif
}

int KokkosParser::initialize() {
    // return codes:
    // 1  - success
    // 0  - already initialized
    // -1 - failed for some other reason

    // determine if Kokkos is already initialized
    // if it has been, get the parameters needed from it
    bool preinitialized = Kokkos::is_initialized();
    
    // if already initialized, return
    if (preinitialized) {
        return 0;
    } else {
        try {
            auto our_args = this->createInitArguments();
            Kokkos::initialize(our_args);
            return 1;
        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            throw e;
        } catch (...) {
            return -1;
        }
    }
}

int KokkosParser::finalize(bool hard_finalize) {
    if (hard_finalize || _called_initialize==1) {
        try {
            Kokkos::finalize();
            _called_initialize = 0; // reset since we finalized
            return 1;
        } catch (...) {
            return 0;
        }
    } else {
        return 1;
    }
}

void KokkosParser::status() const {
#ifdef COMPADRE_USE_CUDA
  printf("KOKKOS mode is enabled on GPU with nthreads: %d,  numa: %d, device_id: %d\n", getNumberOfThreads(), getNuma(), getDeviceID());
#else
  printf("KOKKOS mode is enabled on CPU with nthreads: %d,  numa: %d\n", getNumberOfThreads(), getNuma());
#endif
}
