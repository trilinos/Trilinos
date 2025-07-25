Quick Start
###########

This document gives a brief introduction to start using Kokkos Kernels in your external project. The first step to use Kokkos Kernels is to install the Kokkos programming model, a quick-start on how to obtain and install Kokkos is available `here <https://kokkos.org/kokkos-core-wiki/get-started/quick-start.html>`_. For the rest of this guide, we will assume that :code:`<kokkos-install-directory>` refers to the directory in which you have installed Kokkos.

Download a release
==================

Our releases are available on `GitHub <https://github.com/kokkos/kokkos-kernels/releases>`_, we recommend selecting the same release as the Kokkos version you have installed to ensure best compatibility.

.. code-block:: shell

  # Uncomment according to the type of file you've downloaded (zip or tar)
  unzip kokkos-kernels-x.y.z.zip
  # tar -xzf kokkos-kernels-x.y.z.tar.gz
  cd kokkos-kernels-x.y.z

Simple Configure, Build and Install
===================================

The simplest way to configure Kokkos Kernels is provided below. With this minimal configuration 

.. code-block:: shell

  cmake -S <kokkos-kernels-source-directory> \
        -B <kokkos-kernels-build-directory> \
	-D CMAKE_INSTALL_PREFIX <kokkos-kernels-install-directory> \
	-D Kokkos_ROOT="<kokkos-install-directory>"
  cmake --build <kokkos-kernels-build-directory>
  cmake --install <kokkos-kernels-build-directory>

On an Apple laptop, after installing Kokkos with Serial and Threads support the first 13 lines of configuration outputs using the above configuration line look as follow

.. code-block:: shell

  -- The CXX compiler identification is AppleClang 16.0.0.16000026
  -- Detecting CXX compiler ABI info
  -- Detecting CXX compiler ABI info - done
  -- Check for working CXX compiler: /Library/Developer/CommandLineTools/usr/bin/c++ - skipped
  -- Detecting CXX compile features
  -- Detecting CXX compile features - done
  -- The project name is: KokkosKernels
  -- Enabled Kokkos devices: THREADS;SERIAL
  -- Found Kokkos version 4.3.1 at /Users/lberge/Research/kokkos_eco/installs/kokkos/lib/cmake/Kokkos
  
  ================================
  Kokkos Kernels version: 4.3.01
  ================================

.. note::

   To avoid any issues we strongly recommend using the same compiler to build Kokkos and Kokkos Kernels. Also please note the version and location of the Kokkos installation that is being used on line 9, as well as the version of Kokkos Kernels used on line 12. Matching versions ensure best compatibility.

For more details on how to build Kokkos Kernels, please refer to `Building <building.html>`_ and the `CMake Options <cmake-keywords.html>`_.

  
Getting Help
============

If you need additional help getting started, please join the `Kokkos Slack Workspace <https://kokkosteam.slack.com/>`_, good channels to ask initial questions are #general, #kokkos-kernels, #build. Here are sign up details. Joining Kokkos Slack is the on ramp for becoming a project contributor.
