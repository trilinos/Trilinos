pipeline {
    agent none

    stages {
        stage('HIP-ROCm-4.2-C++14') {
            agent {
                dockerfile {
                    filename 'Dockerfile.hip'
                    dir 'scripts/docker'
                    additionalBuildArgs '--build-arg BASE=rocm/dev-ubuntu-20.04:4.2'
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
                        -DCMAKE_CXX_EXTENSIONS=OFF \
                        -DKokkos_ENABLE_HIP=ON \
                        .. && \
                       make -j8 && make install && \
                       cd ../.. && rm -rf kokkos'''
                sh '''rm -rf build && mkdir -p build && cd build && \
                      cmake \
                        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
                        -DCMAKE_CXX_COMPILER=hipcc \
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
