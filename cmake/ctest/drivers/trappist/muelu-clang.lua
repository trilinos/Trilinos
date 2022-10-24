load("sems-clang/11.0.1")
load("sems-openmpi")
load("sems-cmake")
load("sems-ninja")
load("sems-git")

load("sems-superlu")

load("sems-yaml-cpp")
load("sems-hdf5")
load("sems-parallel-netcdf")
load("sems-zlib")

load("sems-boost")
load("sems-python")

load("sems-metis")
load("sems-parmetis")

load("sems-cuda")

pushenv("OMP_NUM_THREADS","2")
pushenv("CUDA_LAUNCH_BLOCKING","1")
pushenv("CUDA_MANAGED_FORCE_DEVICE_ALLOC","1")

-- Only run on the Tesla K40, not the Quadro --
pushenv("CUDA_VISIBLE_DEVICES","0")