FROM docker.io/library/debian:stable

RUN apt-get update \
    && apt-get install -y \
    gcc g++ gfortran \
    make cmake ninja-build \
    libssl-dev ca-certificates \
    ccache \
    git less \
    libopenblas-serial-dev \
    mpi-default-bin mpi-default-dev \
    libhdf5-mpi-dev hdf5-tools \
    libblosc-dev libbz2-dev libxml2-dev \
    libnetcdf-dev \
    libnetcdf-mpi-dev \
    libpnetcdf-dev \
    libnetcdf-pnetcdf-dev \
    libsuperlu-dev \
    libsuperlu-dist-dev \
    libboost-dev \
    libboost-all-dev \
    llvm-dev clang libclang-dev libzstd-dev \
    python3 python3-dev python-is-python3 python3-pip python3-numpy python3-scipy python3-matplotlib python3-mpi4py python3-yaml python3-h5py python3-tk \
    python3-pybind11 \
    python3-torch libtorch-dev \
    --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

# Install binder
RUN git clone https://github.com/RosettaCommons/binder.git /binder \
    && cd /binder && mkdir binder-build && cd binder-build \
    && cmake -DCMAKE_CXX_COMPILER="$(which clang++)" .. \
    && make install \
    && cd / && rm -rf /binder

# allow running MPI as root in the container
# bind MPI ranks to hwthreads
ENV OMPI_MCA_hwloc_base_binding_policy=hwthread \
    MPIEXEC_FLAGS=--allow-run-as-root \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    PIP_BREAK_SYSTEM_PACKAGES=1 \
    PIP_ROOT_USER_ACTION=ignore

# Install Jupyterlab
RUN pip install notebook jupyterlab wurlitzer ipyparallel && \
    jupyter labextension disable "@jupyterlab/apputils-extension:announcements"

COPY . /scripts

ARG REPO
ARG BRANCH

ENV TRILINOS_DIR=/workspace/trilinos/source \
    TRILINOS_BUILD_DIR=/workspace/trilinos/build \
    TRILINOS_INSTALL_DIR=/workspace/trilinos/install \
    KOKKOS_TOOLS_DIR=/workspace/kokkos-tools/source \
    KOKKOS_TOOLS_BUILD_DIR=/workspace/kokkos-tools/build \
    KOKKOS_TOOLS_INSTALL_DIR=/workspace/kokkos-tools/install \
    PYTHONPATH=/workspace/trilinos/build/packages/PyTrilinos2

RUN echo "Cloning ${REPO}, branch ${BRANCH}" && git clone -b ${BRANCH} ${REPO} ${TRILINOS_DIR} && \
    git clone -b develop https://github.com/kokkos/kokkos-tools.git ${KOKKOS_TOOLS_DIR} && \
    mkdir -p $TRILINOS_BUILD_DIR && \
    mkdir -p $KOKKOS_TOOLS_BUILD_DIR

SHELL ["/bin/bash", "-c"]

ARG BUILD_PARALLELISM=2

RUN source /scripts/commandsDevContainer.sh && \
    configure_trilinos && \
    build_trilinos
