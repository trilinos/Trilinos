#/bin/bash
export GCC_VERSION=10.3.0
export OPENMPI_VERSION=4.1.6

nohup podman build \
        --build-arg=compiler_version="@${GCC_VERSION}" \
        --build-arg=mpi_version="@${OPENMPI_VERSION}" \
        -t gcc-${GCC_VERSION}-openmpi-${OPENMPI_VERSION}-trilinos-env:${USER}-test . \
  &> gcc-${GCC_VERSION}-openmpi-${OPENMPI_VERSION}-trilinos-env.output &

