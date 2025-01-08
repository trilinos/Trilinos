FROM registry.access.redhat.com/ubi8:latest

# Image basename (useful for describing image characteristics, e.g. toolchain)
ARG AT2_image
ENV AT2_IMAGE=${AT2_image}
# Image path as pushed (useful for referencing where to pull the image)
ARG AT2_image_fullpath
ENV AT2_IMAGE_FULLPATH=${AT2_image_fullpath}

# Disable all host repositories so that we're only using things from UBI
# (in the event that you're building a UBI container on a RHEL host, RedHat
# tries to 'help' and let you use your subscription inside the container)
RUN rm -rf /etc/rhsm-host
RUN yum -y --setopt=tsflags=nodocs update && yum clean all

RUN yum -y install gcc gcc-c++ xz bzip2 patch diffutils file make git python38 procps libX11-devel
RUN git clone https://github.com/spack/spack.git /home/runner/spack && cd /home/runner/spack && git checkout develop-2024-10-13
RUN echo -e "export SPACK_ROOT=/spack\nsource /spack/share/spack/setup-env.sh\n" > /etc/profile.d/spack.sh
RUN bash -l -c "spack config add config:build_jobs:64"

ARG compiler_version=
RUN bash -l -c "spack compiler find"
RUN bash -l -c "spack install intel-oneapi-compilers${compiler_version} target=x86_64 > /dev/null & spack install intel-oneapi-compilers${compiler_version} target=x86_64"
RUN echo -e "spack load intel-oneapi-compilers${compiler_version}\n" >> /etc/profile.d/spack.sh
RUN bash -l -c "spack compiler find"

ARG mpi_version=
# Consider using here-documents once we get to Buildah 1.33.0+
# https://github.com/containers/buildah/issues/3474#issuecomment-1821752538
RUN echo -e "\
spack:\n\
  definitions:\n\
  - compiler:\n\
    - oneapi\n\
  - packages:\n\
    - cmake\n\
    - ninja\n\
    - boost +program_options+system\n\
    - cgns\n\
    - hdf5 +hl\n\
    - metis\n\
    - netcdf-c +parallel-netcdf\n\
    - parallel-netcdf\n\
    - parmetis\n\
    - superlu-dist\n\
    - zlib\n\
    - intel-oneapi-mpi${mpi_version}\n\
    - matio\n\
    - intel-oneapi-mkl\n\
  concretizer:\n\
    unify: when_possible\n\
  specs:\n\
  - matrix:\n\
    - [\$packages]\n\
    - [\$%compiler]\n\
  config:\n\
    install_missing_compilers: false\n\
  view: false\n\
  packages:\n\
    all:\n\
      target: [x86_64]\n" > /root/spack.yaml

RUN echo -e "spack env activate /root\n" >> /etc/profile.d/spack.sh

RUN bash -l -c "spack concretize"
RUN bash -l -c "cd /home/runner && spack env depfile -o Makefile"
# Just in case a mirror fails to pull, retry the spack install so that we can preserve
# all of the packages that did build successfully
RUN bash -l -c "cd /home/runner && make -j6 || true"
RUN bash -l -c "cd /home/runner && make -j6"
RUN rm /home/runner/Makefile

CMD [ "/bin/bash" ]

##########################################################################################
# Everything below this line can be removed long-term (legacy environment modules support)
##########################################################################################

RUN yum install -y environment-modules
RUN echo -e "\
  modules:\n\
    default:\n\
      enable:\n\
      - tcl\n\
      tcl:\n\
        hash_length: 0\n\
        all:\n\
          environment:\n\
            set:\n\
              '{name}_ROOT': '{prefix}'\n\
              '{name}_VERSION': '{version}'\n\
              '{name}_BIN': '{prefix.bin}'\n\
              '{name}_INC': '{prefix.include}'\n\
              '{name}_LIB': '{prefix.lib}'\n" >> /root/spack.yaml
RUN bash -l -c "spack module tcl refresh -y"
RUN bash -l -c "echo \"export MODULEPATH=/spack/share/spack/modules/linux-\\\$(spack arch --operating-system)-x86_64\" >> /etc/profile.d/spack.sh"
RUN echo -e "\
module load intel-oneapi-mpi\n\
module load cmake\n\
module load ninja\n\
module load boost\n\
module load cgns\n\
module load hdf5\n\
module load metis\n\
module load netcdf-c\n\
module load parallel-netcdf\n\
module load parmetis\n\
module load superlu-dist\n\
module load zlib\n\
module load matio\n\
module load binder\n\
module load py-mpi4py\n\
module load py-numpy\n\
module load py-pybind11\n\
module load intel-oneapi-mkl\n" >> /etc/profile.d/spack.sh
