pipeline {
    agent none

    options {
        timeout(time: 3, unit: 'HOURS')
    }

    stages {
        stage('Build & Run') {
            parallel {
                stage('SYCL-OneAPI') {
                    agent {
                        dockerfile {
                            filename 'Dockerfile.sycl'
                            dir 'scripts/docker'
                            label 'nvidia-docker && volta'
                            args '-v /tmp/ccache.kokkos:/tmp/ccache'
                        }
                    }
                    steps {
                        sh '''. /opt/intel/oneapi/setvars.sh --include-intel-llvm && \
                              rm -rf kokkos &&
                              git clone -b develop https://github.com/kokkos/kokkos.git && cd kokkos && \
                              mkdir build && cd build && \
                              cmake \
                                -DCMAKE_BUILD_TYPE=Release \
                                -DCMAKE_CXX_COMPILER=/opt/intel/oneapi/compiler/2023.0.0/linux/bin-llvm/clang++ \
                                -DKokkos_ARCH_VOLTA70=ON \
                                -DKokkos_ENABLE_DEPRECATED_CODE_3=OFF \
                                -DKokkos_ENABLE_SYCL=ON \
                                -DKokkos_ENABLE_UNSUPPORTED_ARCHS=ON \
                                -DCMAKE_CXX_STANDARD=17 \
                              .. && \
                              make -j8 && make install && \
                              cd ../.. && rm -rf kokkos'''
                        sh '''. /opt/intel/oneapi/setvars.sh --include-intel-llvm && \
                              rm -rf build && mkdir -p build && cd build && \
                              cmake \
                                -DCMAKE_BUILD_TYPE=Release \
                                -DCMAKE_CXX_COMPILER=/opt/intel/oneapi/compiler/2023.0.0/linux/bin-llvm/clang++ \
                                -DKokkosKernels_ENABLE_TESTS=ON \
                                -DKokkosKernels_ENABLE_EXAMPLES=ON \
                                -DKokkosKernels_INST_DOUBLE=ON \
                                -DKokkosKernels_INST_ORDINAL_INT=ON \
                                -DKokkosKernels_INST_OFFSET_INT=ON \
                              .. && \
                              make -j8'''
                    }
                }

                stage('HIP-ROCm-5.2') {
                    agent {
                        dockerfile {
                            filename 'Dockerfile.hip'
                            dir 'scripts/docker'
                            additionalBuildArgs '--build-arg BASE=rocm/dev-ubuntu-20.04:5.2'
                            label 'rocm-docker && vega'
                            args '-v /tmp/ccache.kokkos:/tmp/ccache --device=/dev/kfd --device=/dev/dri --security-opt seccomp=unconfined --group-add video --env HIP_VISIBLE_DEVICES=$HIP_VISIBLE_DEVICES'
                        }
                    }
                    steps {
                        sh '''rm -rf kokkos &&
                              git clone -b develop https://github.com/kokkos/kokkos.git && cd kokkos && \
                              mkdir build && cd build && \
                              cmake \
                                -DCMAKE_CXX_COMPILER=hipcc \
                                -DCMAKE_CXX_STANDARD=17 \
                                -DCMAKE_CXX_EXTENSIONS=OFF \
                                -DKokkos_ENABLE_HIP=ON \
                              .. && \
                              make -j8 && make install && \
                              cd ../.. && rm -rf kokkos'''
                        sh '''rm -rf build && mkdir -p build && cd build && \
                              cmake \
                                -DCMAKE_BUILD_TYPE=RelWithDebInfo \
                                -DCMAKE_CXX_COMPILER=hipcc \
                                -DCMAKE_CXX_STANDARD=17 \
                                -DCMAKE_CXX_EXTENSIONS=OFF \
                                -DKokkosKernels_ENABLE_TESTS=ON \
                                -DKokkosKernels_ENABLE_EXAMPLES=ON \
                                -DKokkosKernels_INST_DOUBLE=ON \
                                -DKokkosKernels_INST_ORDINAL_INT=ON \
                                -DKokkosKernels_INST_OFFSET_INT=ON \
                              .. && \
                              make -j8 && ctest --verbose'''
                    }
                }
            }
        }
    }
}
