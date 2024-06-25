# PyTrilinos2: Automatic Python Interfaces to Trilinos Packages

PyTrilinos2 requires both Pybind11 and Binder 

## Pybind11 installation:

Pybind11 can be installed via
python -m pip install pybind11

## Binder installation:

To build Binder execute the following command sequence in shell with `$HOMEBINDER` the location where you want to build Binder:

```
# clone Binder
cd $HOMEBINDER
git clone https://github.com/RosettaCommons/binder.git && cd binder
git checkout teuchos_rcp

# Create build dir
mkdir $HOMEBINDER/prefix && cd $HOMEBINDER/prefix

# Clone  LLVM
git clone https://github.com/llvm/llvm-project.git llvm && cd llvm
git checkout llvmorg-6.0.1

# git log
#  commit d359f2096850c68b708bc25a7baca4282945949f (HEAD, tag: llvmorg-6.0.1-rc3, tag: llvmorg-6.0.1, origin/release/6.x)
#  Author: Tom Stellard <tstellar@redhat.com>
#  Date:   Thu Jun 14 22:33:33 2018 +0000

# Create symlink pointing to binder/src dir
ln -s $HOMEBINDER/binder/source $HOMEBINDER/prefix/llvm/clang-tools-extra/binder

# Create ``llvm/tools/clang/tools/extra/CMakeLists.txt`` file with content: ``add_subdirectory(binder)``
echo 'add_subdirectory(binder)' > $HOMEBINDER/prefix/llvm/clang-tools-extra/CMakeLists.txt

# Build Binder
mkdir $HOMEBINDER/prefix/build && cd $HOMEBINDER/prefix/build
cmake -S ../llvm/llvm -DLLVM_ENABLE_PROJECTS="clang;libcxx;libcxxabi;clang-tools-extra" \
 -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_EH=1 -DLLVM_ENABLE_RTTI=ON && make -j 10

# At this point, if all above steps is successful, binder should be at
# $HOMEBINDER/prefix/build/bin/binder
```

## PyTrilinos2 configuration:

```
-D Trilinos_ENABLE_PyTrilinos2:BOOL=ON \
-D PYTHON_EXECUTABLE=... \
-D PyTrilinos2_ENABLE_TESTS=ON \
-D PyTrilinos2_BINDER_EXECUTABLE=... \
-D PyTrilinos2_BINDER_GCC_TOOLCHAIN=...\
```
Alternatively, the last option can be replace by the following options to include C++ std headers:
```
-D PyTrilinos2_BINDER_clang_include_dirs=... \
-D PyTrilinos2_BINDER_LibClang_include_dir=... \
```
Depending on the enabled exacution space, Binder might need some extra flags.
For example, for OpenMP, binder needs the `-fopenmp` and an include path of where the OpenMP headers can be included.
```
-D PyTrilinos2_BINDER_FLAGS="-fopenmp;-I..." \
```
As a second example, for CUDA, one might need to specify the used architecture and a path to CUDA related headers:
```
-D PyTrilinos2_BINDER_FLAGS="--cuda-gpu-arch=sm_70;--cuda-path=..." \
```

## Report issues:

PyTrilinos2 relies on an automatic generation of the interface files.
Although this automatic generation has been intensively tested, it is unfeasible to test every combinations of options of Trilinos.
Therefore, if you face a configuration error or a compilation error, please submit an issue on https://github.com/trilinos/Trilinos 
with the PyTrilinos2 tag, the used configuration script, the log of the outpupt of the configuration, and the log of the build process.

## Copyright and License
See PyTrilinos2/COPYRIGHT, PyTrilinos2/LICENSE, https://trilinos.github.io/license.html and individual file headers for additional information.

## Questions? 
Contact lead developers:

* PyTrilinos2 team (GitHub handle: @trilinos/pytrilinos2)
* Kim Liegeois     (GitHub handle: [kliegeois](https://github.com/kliegeois) or knliege@sandia.gov)
* Christian Glusa  (GitHub handle: [cgcgcg](https://github.com/cgcgcg) or caglusa@sandia.gov)
