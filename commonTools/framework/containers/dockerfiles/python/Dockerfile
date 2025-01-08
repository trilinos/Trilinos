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

RUN yum -y install gcc gcc-c++ gcc-gfortran xz bzip2 patch diffutils file make git python39 procps environment-modules gettext unzip

# FIXME: Remove this someday
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN useradd runner
USER runner

RUN git clone https://github.com/spack/spack.git /home/runner/spack && cd /home/runner/spack && git checkout develop-2024-10-13
RUN echo -e "export SPACK_ROOT=/home/runner/spack\nsource /home/runner/spack/share/spack/setup-env.sh\n" >> /home/runner/.bashrc
RUN bash -l -c "spack config add config:build_jobs:64"
RUN bash -l -c "spack bootstrap now"

RUN bash -l -c "spack compiler find"

RUN mkdir /home/runner/environment
# Consider using here-documents once we get to Buildah 1.33.0+
# https://github.com/containers/buildah/issues/3474#issuecomment-1821752538
RUN echo -e "\
spack:\n\
  include:\n\
  - modules.yaml\n\
  definitions:\n\
  - compiler:\n\
    - gcc\n\
  - packages:\n\
    - cmake\n\
    - ninja\n\
    - py-pytest-cov\n\
    - zlib\n\
    - emacs\n\
    - gh\n\
    - unzip\n\
  concretizer:\n\
    unify: when_possible\n\
  specs:\n\
  - matrix:\n\
    - [\$packages]\n\
    - [\$%compiler]\n\
  view: false\n\
  packages:\n\
    all:\n\
      target: [x86_64]\n" > /home/runner/environment/spack.yaml

RUN echo -e "\
modules:\n\
  prefix_inspections:\n\
    ./bin:\n\
      - PATH\n\
    ./lib:\n\
    - LIBRARY_PATH\n\
    - LD_LIBRARY_PATH\n\
    ./lib64:\n\
    - LIBRARY_PATH\n\
    - LD_LIBRARY_PATH\n\
    ./include:\n\
    - INCLUDE\n\
    ./man:\n\
      - MANPATH\n\
    ./share/man:\n\
      - MANPATH\n\
    ./share/aclocal:\n\
      - ACLOCAL_PATH\n\
    ./lib/pkgconfig:\n\
      - PKG_CONFIG_PATH\n\
    ./lib64/pkgconfig:\n\
      - PKG_CONFIG_PATH\n\
    ./share/pkgconfig:\n\
      - PKG_CONFIG_PATH\n\
    ./:\n\
      - CMAKE_PREFIX_PATH\n\
  default:\n\
    enable:\n\
      - tcl\n\
    tcl:\n\
      hash_length: 0\n\
      projections:\n\
        all: '{name}/{version}'\n\
      all:\n\
        conflict:\n\
        - '{name}'\n\
        environment:\n\
          set:\n\
            '{name}_ROOT': '{prefix}'\n\
            '{name}_VERSION': '{version}'\n\
            '{name}_BIN': '{prefix.bin}'\n\
            '{name}_INC': '{prefix.include}'\n\
            '{name}_LIB': '{prefix.lib}'\n\
        autoload: none\n\
      blas:\n\
        environment:\n\
          set:\n\
            'BLAS_ROOT': '{prefix}'\n\
            'LAPACK_ROOT': '{prefix}'\n\
      gcc:\n\
        environment:\n\
          set:\n\
            'COMPILER_ROOT': '{prefix}'\n\
      llvm:\n\
        environment:\n\
          set:\n\
            'COMPILER_ROOT': '{prefix}'" >> /home/runner/environment/modules.yaml

RUN bash -l -c "spack module tcl refresh -y"

RUN rm -rf /home/runner/.spack
RUN bash -l -c "spack bootstrap now"

RUN echo -e "spack env activate /home/runner/environment\n" >> /home/runner/.bashrc
RUN bash -l -c "echo \"export MODULEPATH=/home/runner/spack/share/spack/modules/linux-\\\$(spack arch --operating-system)-x86_64\" >> /home/runner/.bashrc"
RUN bash -l -c "spack compiler find $COMPILER_ROOT"


RUN bash -l -c "spack concretize"
RUN bash -l -c "cd /home/runner && spack env depfile -o Makefile"
# Just in case a mirror fails to pull, retry the spack install so that we can preserve
# all of the packages that did build successfully
RUN bash -l -c "cd /home/runner && make -j6 || true"
RUN bash -l -c "cd /home/runner && make -j6"
RUN rm /home/runner/Makefile

# ENV OMPI_ALLOW_RUN_AS_ROOT=1
# ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

CMD [ "/bin/bash" ]

RUN bash -l -c "spack module tcl refresh -y"
RUN echo -e "\
module load cmake\n\
module load ninja\n\
module load zlib\n\
module load emacs\n\
module load py-pytest-cov\n\
module load gh\n\
">> /home/runner/.bashrc

USER root
ADD GenConfig /GenConfig
RUN chmod -R ugo+rwX /GenConfig

RUN echo -e "source /home/runner/.bashrc" >> /root/.bashrc
