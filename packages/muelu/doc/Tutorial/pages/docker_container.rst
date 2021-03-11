====================
B. Docker Container
====================

This chapter discusses the basics of the Docker container that comes with this tutorial to allow the user to follow above explanations and do its own experiments with MueLu and Trilinos. A docker container has the advantage that it is rather easy to set up for a user. Even though compiling and installing is easier in recent years by using a cmake based build system, it may be difficult for not so experienced users. The Docker container runs on any machine with Docker installed and brings all the necessary tools for a quick start to MueLu.

Preparations
============

To use the Docker container you must perform the following steps.

#. Install **Docker** on your host machine. You can download it from **www.docker.com**.
#. Download the MueLu **docker-images** repository. It is located at **https://github.com/GrahamBenHarper/docker-images**.
#. In the **docker-images/muelu-tutorial** directory, run `build-container.sh`. This will take some time depending on your machine.

    .. warning::

      Insert screen output of Docker building or pulling

#. Once the build is complete, run **run-container.sh** to run a new tutorial container. The tutorial is installed in the `/opt/trilinos/build/packages/muelu/test/tutorial` directory of the Docker image.

    .. warning::

      Insert screen output of Docker running

.. note::

  The container does not have a graphical interface, so some of the visualization aspects of the tutorial may be unavailable.

Software
========

The Docker container is based on CentOS-Stream 9. The container requires approximately 8GB after Trilinos is built.

The following software is pre-installed:

::

    Terminal: bash
    Text editor: nano
    Version control: git


The following system libraries are installed:

::

    MPI: OpenMI 4.1.1
    Python: Python 3.9.18
    Compiler: gcc 11.4.1
    CMake: 3.26.1
    NetCDF
    HDF5
    BLAS
    LAPACK
    Boost
    Atlas

