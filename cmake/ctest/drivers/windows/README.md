# Required Software
The following list of software is required by the Windows build scripts. The build scripts
assume this software will be installed in specific locations, so before installing it would
be worthwhile to browse the scripts for these locations. Alternatively, you could download
and install the software wherever you want - just be sure to copy and update the build scripts
with the correct locations.

**[CMake][1]** - required to setup the build configuration. (Tested using version 3.8.1)

**[Ninja][2]** - required to parallelize the build on Windows. (Tested using version 1.7.2)

**[Visual Studio 2015][3]** - required to build packages.

**[Microsoft MPI][4]** - required to create MPI builds of packages. (Tested using version 8.1.12438.1084)

**Perl** - required by some packages for testing. (Tested using [Strawberry Perl][5] version 5.24.1.1)

**[Git][6]** - required for the scripts to update source code.

**CLAPACK** - required TPL dependency for some packages (Tested using version 3.2.1). Note that if a
pre-built version is not available, you must [build it from source code][7].

[1]: https://cmake.org/download/
[2]: https://ninja-build.org/
[3]: https://www.visualstudio.com/
[4]: https://msdn.microsoft.com/en-us/library/bb524831(v=vs.85).aspx
[5]: http://strawberryperl.com/
[6]: https://git-scm.com/
[7]: http://icl.cs.utk.edu/lapack-for-windows/clapack/


# Script Usage
There are two scripts intended for the typical user - *task_driver_windows.bat* and
*create_windows_package.bat*. For details about additional scripts, please refer to the
[Script File Summary](#Script-File-Summary) section below.

##### task_driver_windows.bat
This script builds and tests Trilinos packages for both Debug and Release configurations,
following the Trilinos package-by-package testing paradigm. The packages are built as static
libraries, and the results are submitted to [CDash][8]. This script may be launched from
the command-line as is, without any additional arguments.

##### create_windows_package.bat
This script creates a ZIP package of static libraries and header files from the set of Trilinos
packages specified in the script. The expected syntax to use the script is:

`create_windows_package.bat <Build Type> <Base Directory>`

where `<Build Type>` specifies the build configuration (Release or Debug) and
`<Base Directory>` specifies a root working directory where the script will clone/update
the Trilinos repository as a sub-directory, create a build sub-directory to do the actual build,
and place resulting output files and ZIP package. For example,

```
create_windows_package.bat Debug C:\path\to\DebugPackageDir
```

would create a Debug build of static libraries, and at the end of the script the base directory
would have the following contents:

```
C:\path\to\DebugPackageDir
   - build                     (Build directory)
   - Trilinos                  (Source code directory)
   - update_output.txt         (Output from the Git update step)
   - configure_output.txt      (Output from the Configure step)
   - build_output.txt          (Output from the build step)
   - package_output.txt        (Output from the package step)
   - trilinos-setup-12.3.zip   (Resulting ZIP package)
```

[8]: https://testing-vm.sandia.gov/cdash


# Script File Summary
##### TrilinosCTestDriverCore.windows.msvc.cmake
CMake script that sets the options common to both the Debug and Release configurations of a 
Trilinos build on Windows using Visual Studio, including the specific Trilinos packages to be
built. This script follows the Trilinos package-by-package testing paradigm. This script 
assumes specific build tools exist (e.g., Ninja, MSVC14, MPI, etc.) and assumes where
they will be located on Windows.

##### create_windows_package.bat
Windows batch script that creates a ZIP package of static libraries and header files. It
updates or clones the latest version of Trilinos as necessary, sets up a build environment
for MSVC14, configures and runs CMake, then builds the specified Trilinos packages using Ninja
and MSVC. After the build is complete, it creates the ZIP package using CPack. The output
for each step of the process is recorded in individual files (e.g., configure_output.txt,
build_output.txt, etc.) within the specified working directory.

##### ctest_windows_mpi_debug.cmake
CMake file that sets the build configuration for a Debug build before calling
*TrilinosCTestDriverCore.windows.msvc.cmake*.

##### ctest_windows_mpi_release.cmake
CMake file that sets the build configuration for a Release build before calling
*TrilinosCTestDriverCore.windows.msvc.cmake*

##### task_driver_windows.bat
Windows batch script that sets important environment variables before running
*ctest_windows_mpi_debug.cmake* and *ctest_windows_mpi_release.cmake* using CTest.


# Notes and Issues
- The scripts assume they are starting from a clean environment. All necessary environment
  variables are set by the scripts.

- When installing Perl on Windows, it likes to insert itself into the system or user PATH
  environment variable. This can have unintended consequences when running CMake, the most
  notable of which is that CMake will (wrongly) assume certain executables in Perl are part
  of the compiler. As a result, CMake won't find the correct MSVC compiler and linker, and
  the configure step will fail. If possible, try to keep Perl out of the PATH when using
  these scripts, or else modify the PATH before running the scripts to remove Perl.

- As noted in Trilinos Pull Request [#1197][i1], the /bigobj compiler flag is necessary for
  debug builds, and is included in the build scripts.

- At the time of writing, the Zoltan tests fail on Windows as noted in Trilinos Issue [#1440][i2].

- The packages are currently built as static libraries. If shared libraries are desired, an
  effort would need to be made to resolve certain build issues. Also, care would need to be
  taken to ensure PATHs are setup correctly for the libraries to find each other during the
  build and testing.

[i1]: https://github.com/trilinos/Trilinos/pull/1197
[i2]: https://github.com/trilinos/Trilinos/issues/1440
