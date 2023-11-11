====================
B. Docker Container
====================

This chapter discusses the basics of the Docker container that comes with this tutorial to allow the user to follow above explanations and do its own experiments with MueLu and Trilinos. A virtual machine has the advantage that it is rather easy to set up for a user. Even though compiling and installing got easier the last years by using a cmake based build system it is still a nightmare for not so experienced users. The Docker container runs on any machine with Docker installed and brings all the necessary tools for a quick start to MueLu.

Preparations
============

To use the Docker container you must perform the following steps.

#. Install **Docker** on your host machine. You can download it from **www.docker.com**.
#. Download the MueLu **docker-images** repository. It is located at **https://github.com/GrahamBenHarper/docker-images**.
#. In the **docker-images/muelu-tutorial** directory, run `build-container.sh`. This will take some time depending on your machine.

    .. warning::

      Insert screen output

#. Once the build is complete, run **run-container.sh** to run a new tutorial container. The tutorial is installed in the `/opt/trilinos/build/packages/muelu/test/tutorial` directory of the Docker image.

    .. warning::

      Insert screen output

.. note::

  The container does not have a graphical interface, so some of the visualization aspects of the tutorial may be unavailable.

Software
========

The virtual machine is based on a minimal installation of **Lubuntu 14.04**. The image file has 4 GB with about 250 MB free for the user.

The following software is pre-installed:

::

    Web-browser: midori
    PDF-viewer: evince
    Terminal: LXTerminal
    Visualization: paraview, gnuplot
    File manager: PCManFM
    Analysis: FreeMat v4.0
    GNU octave 3.8.1


The following system libraries are installed:

::

    Trilinos: Trilinos (developer branch: Oct 1, 2014)
    Direct solver: SuperLU 4.3
    VTK: VTK 5.8
    MPI: OpenMPI 1.6.5
    Python: Python 2.7.6
    Compiler: gcc 4.8.2


