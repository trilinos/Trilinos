# trilinos-containers

The registry for trilinos images can be found here:

[Trilinos Container Registry](https://gitlab-ex.sandia.gov/trilinos-project/trilinos-containers/container_registry)

# Current PR development environment images
**registry-ex.sandia.gov/trilinos-project/trilinos-containers/rhel8/trilinos-pr-env:gcc-8.5.0**

# Getting started
1. You need a machines where you can run containers.  You need podman or docker on the machine.  Rhel8 will have podman by default windows and mac users can install docker and these should work.  Currently all of the containers are built for x86-64 and need to be run on capatable hardware.

1. you may need to login to the registry:
`podman login registry-ex.sandia.gov`
`username: `
`password: <your password for gitlab-ex not you kerberose>`

1. pull down the image:
`podman pull registry-ex.sandia.gov/trilinos-project/trilinos-containers/rhel8/trilinos-pr-env:gcc-8.5.0`

1. To run the container interactively:
`podman run --rm -it registry-ex.sandia.gov/trilinos-project/trilinos-containers/rhel8/trilinos-pr-env:gcc-8.5.0`

once in the container, the environment is already set and it looks like a familiar module environment

```
root@trilinos-container-gcc-8:[/]: module list
Currently Loaded Modulefiles:
 1) gcc/8.5.0      6) cuda/11.8.0     11) netcdf-cxx/4.2        16) parallel-netcdf/1.12.3  21) superlu-dist/5.4.0  
 2) boost/1.82.0   7) git/2.40.0      12) netcdf-fortran/4.6.0  17) parmetis/4.0.3          22) superlu/5.3.0       
 3) ccache/4.8     8) hdf5/1.14.1-2   13) netlib-lapack/3.11.0  18) python/3.10.10          23) yaml-cpp/0.6.2      
 4) cgns/4.3.0     9) metis/5.1.0     14) ninja/1.11.1          19) scotch/7.0.3            24) zlib/1.2.13         
 5) cmake/3.26.3  10) netcdf-c/4.9.2  15) openmpi/4.1.5         20) suite-sparse/5.13.0     
root@trilinos-container-gcc-8:[/]: echo $CC
/spack-installs/gcc/8.5.0/gcc/8.5.0/base/k4bt6ni/bin/gcc
root@trilinos-container-gcc-8:[/]: echo $CXX
/spack-installs/gcc/8.5.0/gcc/8.5.0/base/k4bt6ni/bin/g++
root@trilinos-container-gcc-8:[/]: echo $MPICC
/spack-installs/openmpi/4.1.5/gcc/8.5.0/base/fil2rza/bin/mpicc
root@trilinos-container-gcc-8:[/]: echo $MPICXX
/spack-installs/openmpi/4.1.5/gcc/8.5.0/base/fil2rza/bin/mpic++
root@trilinos-container-gcc-8:[/]: which cmake
/spack-installs/cmake/3.26.3/gcc/8.5.0/base/l4oiivv/bin/cmake

```

# To build another image on top of this image: 

 Create a Dockerfile that uses it

```
from registry-ex.sandia.gov/trilinos-project/trilinos-containers/rhel8/trilinos-pr-env:gcc-8.5.0

RUN <whatever commands you want>
```

# Information about your container

The `AT2_IMAGE` environment variable should be set and contain the name of the image.
This should be helpful when debugging or reproducing problems.
