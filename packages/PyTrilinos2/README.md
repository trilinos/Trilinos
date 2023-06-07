# PyTrilinos2

PyTrilinos2 requires Pybind11 which can be installed via
python -m pip install pybind11

If Trilinos is configured with default options, the provided header
files will work out-of-the-box. Otherwise, Binder needs to be
available as well. Instructions on how to install Binder can be found
below.

## Binder installation:

To build Binder execute the following command sequence in shell with `$HOMEBINDER` the location where you want to build Binder:

```
# clone Binder
cd $HOMEBINDER
git clone https://github.com/kliegeois/binder.git && cd binder
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

Basic options:
```
-D Trilinos_ENABLE_PyTrilinos2:BOOL=ON \
-D PYTHON_EXECUTABLE=... \
-D PyTrilinos2_ENABLE_TESTS=ON \
```
Advanced options to regenerate the source files and potentially update the source files:
```
-D PyTrilinos2_ENABLE_Binder=ON \
-D PyTrilinos2_BINDER_EXECUTABLE=... \
-D PyTrilinos2_BINDER_GCC_TOOLCHAIN=...\
-D PyTrilinos2_ENABLE_Update_Binder=ON \
```
The last option specify if the source files have to be updated (to be committed).
