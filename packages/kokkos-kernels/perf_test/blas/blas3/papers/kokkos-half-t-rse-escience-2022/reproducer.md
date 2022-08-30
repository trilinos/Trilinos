## To reproduce the half precision results for batched-GEMM:
```bash
git clone https://github.com/kokkos/kokkos.git
git clone https://github.com/kokkos/kokkos-kernels.git
cd kokkos-kernels
git checkout tags/papers/us-rse-escience-2022
cd perf_test/blas/blas3
export KOKKOS_SRC_DIR=/path/to/kokkos
export KOKKOSKERNELS_SRC_DIR=/path/to/kokkos-kernels
```

### On V100
```bash
./KokkosBatched_BatchedGemm_benchmark.sh double SNB VOLTA70
./KokkosBatched_BatchedGemm_benchmark.sh float SNB VOLTA70
./KokkosBatched_BatchedGemm_benchmark.sh half SNB VOLTA70
./KokkosBatched_BatchedGemm_benchmark.sh bhalf SNB VOLTA70
```

### On A100
```bash
./KokkosBatched_BatchedGemm_benchmark.sh double DEFAULT AMPERE80
./KokkosBatched_BatchedGemm_benchmark.sh float DEFAULT AMPERE80
./KokkosBatched_BatchedGemm_benchmark.sh half DEFAULT AMPERE80
./KokkosBatched_BatchedGemm_benchmark.sh bhalf DEFAULT AMPERE80
```