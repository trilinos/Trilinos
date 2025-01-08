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

RUN yum -y install gcc gcc-c++ xz bzip2 patch diffutils file make git python38 procps
RUN git clone https://github.com/spack/spack.git /home/runner/spack && cd /home/runner/spack && git checkout develop-2024-10-13
RUN echo -e "export SPACK_ROOT=/spack\nsource /spack/share/spack/setup-env.sh\n" > /etc/profile.d/spack.sh
RUN bash -l -c "spack config add config:build_jobs:64"
RUN bash -l -c "spack bootstrap now"

ARG compiler_version=
RUN bash -l -c "spack compiler find"
RUN bash -l -c "spack install llvm-amdgpu${compiler_version} target=x86_64 > /dev/null & spack install llvm-amdgpu${compiler_version} target=x86_64"
RUN bash -l -c "spack view symlink -i /opt/rocmcc-5.7.1 llvm-amdgpu@5.7.1"
RUN echo -e "export PATH=/opt/rocmcc-5.7.1/bin:\$PATH\nexport LD_LIBRARY_PATH=/opt/rocmcc-5.7.1/lib:\$LD_LIBRARY_PATH\nexport MANPATH=/opt/rocmcc-5.7.1/share/man:\$MANPATH\n" >> /etc/profile.d/spack.sh

ARG mpi_version=
# Consider using here-documents once we get to Buildah 1.33.0+
# https://github.com/containers/buildah/issues/3474#issuecomment-1821752538
RUN echo -e "\
spack:\n\
  compilers:\n\
  - compiler:\n\
      spec: rocmcc@5.7.1\n\
      paths:\n\
        cc: /opt/rocmcc-5.7.1/bin/amdclang\n\
        cxx: /opt/rocmcc-5.7.1/bin/amdclang++\n\
        f77: /opt/rocmcc-5.7.1/bin/amdflang\n\
        fc: /opt/rocmcc-5.7.1/bin/amdflang\n\
      modules: []\n\
      operating_system: rhel8\n\
      environment:\n\
        prepend_path: { \"C_INCLUDE_PATH\": /opt/rocmcc-5.7.1/include, \"CPLUS_INCLUDE_PATH\": /opt/rocmcc-5.7.1/include }\n\
  definitions:\n\
  - compiler:\n\
    - rocmcc@5.7.1\n\
  - packages:\n\
    - hip\n\
  concretizer:\n\
    unify: when_possible\n\
  specs:\n\
  - matrix:\n\
    - [\$packages]\n\
    - [\$%compiler]\n\
  view: false\n\
  packages:\n\
    llvm-amdgpu:\n\
      externals:\n\
      - spec: llvm-amdgpu@5.7.1\n\
        prefix: /opt/rocmcc-5.7.1\n\
      buildable: false\n\
    all:\n\
      target: [x86_64]\n" > /root/spack.yaml

RUN echo -e "spack env activate /root\n" >> /etc/profile.d/spack.sh

RUN bash -l -c "spack external find"
# Just in case a mirror fails to pull, retry the spack install so that we can preserve
# all of the packages that did build successfully
RUN bash -l -c "spack install > /dev/null & spack install > /dev/null & spack install > /dev/null & spack install || true"
RUN bash -l -c "spack install > /dev/null & spack install > /dev/null & spack install > /dev/null & spack install"

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
module load cray-mpich\n\
module load cmake\n\
module load ninja\n\
module load openblas\n" >> /etc/profile.d/spack.sh
