#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

from optparse import OptionParser
import os
import sys
import shutil
import re
import subprocess

# Default list of package enables and disables.  Trikota, STK, and
# Optika are disabled due to explicit dependencies which can't always
# be met.  Mesquite is disabled because it doesn't generate a
# MesquiteConfig.cmake file.
DEFAULT_ENABLE_DISABLE_LIST = [
  ("Amesos", True),
  ("Amesos2", False),
  ("Anasazi", True),
  ("AztecOO", True),
  ("Belos", True),
  ("CTrilinos", False),
  ("Didasko", False),
  ("Epetra", True),
  ("EpetraExt", True),
  ("FEI", True),
  ("ForTrilinos", False),
  ("Galeri", True),
  ("GlobiPack", True),
  ("Ifpack", True),
  ("Ifpack2", True),
  ("Intrepid", True),
  ("Isorropia", True),
  ("Kokkos", True),
  ("Komplex", True),
  ("Mesquite", False),
  ("ML", True),
  ("Moertel", True),
  ("MOOCHO", True),
  ("NOX", True),
  ("Optika", False),
  ("OptiPack", True),
  ("Pamgen", True),
  ("Phalanx", True),
  ("Piro", False),
  ("Pliris", True),
  ("PyTrilinos", False),
  ("RBGen", False),
  ("RTOp", True),
  ("Rythmos", True),
  ("Sacado", True),
  ("SEACAS", False),
  ("Shards", True),
  ("STK", False),
  ("Stokhos", False),
  ("Stratimikos", True),
  ("Sundance", False),
  ("Teko", False),
  ("Teuchos", True),
  ("ThreadPool", True),
  ("Thyra", True),
  ("Tpetra", True),
  ("TriKota", False),
  ("TrilinosCouplings", False),
  ("Triutils", True),
  ("Zoltan", True),
]

#Make reporting errors easier to make consistent between all system calls
def reportError(errorCode, errorMessage, errorFileName="", verboseError=False):
  if errorCode != 0:
    print errorMessage
    if verboseError and errorFileName != "":
      errorFile = open(errorFileName)
      for line in errorFile:
        print line.strip()
      errorFile.close()
    sys.exit(1)

def run_command(command, working_directory, outputfile):
  """
  Run a shell command, redirecting output to the given output file.
  """
  return subprocess.call(command, shell=True,
                         stdout=outputfile, stderr=outputfile,
                         cwd=working_directory)

def configure_cmake_build(build_dir, source_dir, cmake_options, outputfile=None):
  """
  Run the configure step for a CMake build in the given build
  directory using the specified source directory. The cmake_options
  argument should be an associative list (a list of pairs) specifying
  variable name and value. The pairs will be turned into
  '-Dname=value' arguments to the CMake command. Output will be
  directed to the file-like object specified by the outputfile
  argument. This function will perform a clean CMake build, that is
  any existing cache will be deleted. The exit code of the CMake
  process is returned.
  """
  print "config dir exists? " + str(os.path.exists(build_dir))
  print "cwd is " + os.getcwd()
  if not os.path.exists(build_dir):
    print "Creating CMake build directory %s" % build_dir
    os.mkdir(build_dir)

  cachepath = os.path.join(build_dir, "CMakeCache.txt")
  if os.path.exists(cachepath):
    os.remove(cachepath)  

  to_cmake = lambda pair: "-D%s=%s" % pair
  configureCmd = "cmake " + " ".join(map(to_cmake, cmake_options)) + " " + source_dir
  print configureCmd
  return run_command(configureCmd, build_dir, outputfile)

def get_options(workingDir):
  parser = OptionParser()

  repositoryDefault="software.sandia.gov:/space/git/nightly/Trilinos"

  parser.add_option("-r", "--repository", dest="repository", action="store",
    default=repositoryDefault, 
    help="Sets the repository to pull Trilinos from. Default=%default")
  parser.add_option("--cxx-compiler", dest="cxx", action="store", default="g++", 
    help="Sets the C++ compiler to use. Default=%default")
  parser.add_option("--c-compiler", dest="cc", action="store", default="gcc", 
    help="Sets the C compiler to use. Default=%default")
  parser.add_option("--fortran-compiler", dest="fortran", action="store", default="gfortran", 
    help="Sets the Fortran compiler to use. Default=%default")
  parser.add_option("--enable-mpi", dest="enableMpi", action="store_true", default=False, 
    help="Enable an MPI build of Trilinos. Default=%default")
  parser.add_option("--mpi-base-path", dest="mpiBaseDir", action="store", default="", 
    help="Sets the base path where we can find the mpi compiler wrappers, run environmnet and headers libraries. Default=%default")
  parser.add_option("--disable-boost", dest="disableBoost", action="store_true", default=False,
                    help="Disable the Boost TPL. Default=%default")
  parser.add_option("--boost-dir", dest="boostDir", action="store", default="/usr/lib", 
    help="Sets the location of the boost library. Default=%default")
  parser.add_option("--disable-netcdf", dest="disableNetcdf", action="store_true", default=False,
                    help="Disable the NetCDF TPL. Default=%default")
  parser.add_option("--netcdf-dir", dest="netcdfDir", action="store", default="/usr/lib", 
    help="Sets the location of the netcdf library. Default=%default")
  parser.add_option("--no-dashboard", dest="doDashboardBuild", action="store_false", default=True, 
    help="Disables the dashboard submission build. The default is to do the dashboard build")
  parser.add_option("--install-dir", dest="installDir", action="store", default=workingDir+"/install", 
    help="Sets the location of where Trilinos will be installed to. Default=%default")
  parser.add_option("--working-dir", dest="workingDir", action="store", default=workingDir, 
    help="Sets the directory that this script will run from. Default=%default")
  parser.add_option("--branch", dest="branch", action="store", default="master", 
    help="Sets the branch that should be used to build the tarball and installation from. Default=%default")
  parser.add_option("--shared", dest="shared", action="store_true", default=False, 
    help="Build shared libraries instead of static. Default=%default")
  parser.add_option("--enable-secondary-stable", dest="secondaryStable", action="store_true", default=False, 
    help="Enable secondary stable code. Default=master")
  parser.add_option("--verbose-errors", dest="verboseErrors", action="store_true", default=False, 
    help="Enable verbose error reporting. This will cause the output of the command that failed to be sent to stdout. Default=master")
  parser.add_option("-j", dest="buildProcs", action="store", default=1, 
    help="Specify the number of processes to run make with. Default=%default")

  return parser.parse_args()

def main(package_enable_disable_list, options):
  workingDir = os.getcwd()

  options.configurePath = "configure"
  options.tarballBuildPath = "tarball_build"
  #strip any relative bits out of the path
  options.workingDir = os.path.abspath(options.workingDir)

  #create the working dir if it doesn't exist
  if not os.path.exists(options.workingDir):
    os.mkdir(options.workingDir)

  #switch to working dir. This is essentially a noop if the working dir = current dir.
  os.chdir(options.workingDir)

  print "Running nightly_create_tarball.py"
  print "Installing to " + options.installDir

  #
  # Build the CMake configure options
  #
  cmake_configure_options = []
  #choosing which compilers to use, mpi vs standard
  compilerPathString = ""
  if options.enableMpi:
    cmake_configure_options.append(("TPL_ENABLE_MPI", "ON"))
    if options.mpiBaseDir != "":
      cmake_configure_options.append(("MPI_BASE_DIR", options.mpiBaseDir))
  else:
    cmake_configure_options.append(("CMAKE_CXX_COMPILER", options.cxx))
    cmake_configure_options.append(("CMAKE_C_COMPILER", options.cc))
    cmake_configure_options.append(("CMAKE_Fortran_COMPILER", options.fortran))

  if options.shared:
    cmake_configure_options.append(("BUILD_SHARED_LIBS:BOOL", "ON"))

  if options.secondaryStable:
    # have to manually disable the teuchos float and complex options
    # since they are enabled by default when not in development mode.
    cmake_configure_options.append(("Trilinos_ENABLE_SECONDARY_STABLE_CODE", "ON"))

  cmake_configure_options.append(("HAVE_GCC_ABI_DEMANGLE", "ON"))
  cmake_configure_options.append(("Trilinos_WARNINGS_AS_ERRORS_FLAGS", ""))
  cmake_configure_options.append(("CMAKE_VERBOSE_MAKEFILE", "TRUE"))
  for package, is_enabled in package_enable_disable_list:
    if is_enabled:
      on_or_off = "ON"
    else:
      on_or_off = "OFF"
    cmake_configure_options.append(("Trilinos_ENABLE_%s" % package, on_or_off))

  cmake_configure_options.append(("Trilinos_ENABLE_EXPLICIT_INSTANTIATION", "ON"))
  if not options.disableBoost:
    cmake_configure_options.append(("TPL_ENABLE_Boost", "ON"))
    cmake_configure_options.append(("Boost_INCLUDE_DIRS", options.boostDir))
  if not options.disableNetcdf:
    cmake_configure_options.append(("TPL_ENABLE_Netcdf", "ON"))
    cmake_configure_options.append(("Netcdf_LIBRARY_DIRS",
                                   os.path.join(options.netcdfDir, "lib")))
    cmake_configure_options.append(("Netcdf_INCLUDE_DIRS",
                                   os.path.join(options.netcdfDir, "include")))
  cmake_configure_options.append(("CMAKE_INSTALL_PREFIX", options.installDir))

  #removing install directory first so that subsequent installs aren't polluted by old installs
  print "attempting to remove the old install dir"
  try:
    shutil.rmtree(options.installDir, True)
  except:
    print "execption while removing " + options.installDir
    print sys.exc_info()[1]

  #removing configure directory first so that subsequent installs aren't polluted by old tarballs or configuration issues
  print "attempting to remove the old configure dir"
  try:
    shutil.rmtree(options.configurePath, True)
  except:
    print "execption while removing " + options.configurePath
    print sys.exc_info()[1]

  print "done trying to remove old install and configure directories"

  gitBranchError = 0
  #clone repo
  if os.path.exists("Trilinos"):
    print "Repository already exists so updating instead of cloning."
    os.chdir("Trilinos")
    if options.branch != "master":
      gitBranchCmd = "git checkout " + options.branch
      gitBranchError = os.system(gitBranchCmd)
    gitPullCmd = "git pull"
    gitError = os.system(gitPullCmd)
    os.chdir(options.workingDir)
  else:
    print "Repository doesn't exist so cloning."
    gitCloneCmd = "git clone " + options.repository + " Trilinos"
    gitError = os.system(gitCloneCmd)
    os.chdir("Trilinos")
    if options.branch != "master":
      gitBranchCmd = "git checkout --track origin/" + options.branch
      gitBranchError = os.system(gitBranchCmd)
    os.chdir(options.workingDir)


  reportError(gitError, "Could not retrieve repository " + options.repository)
  reportError(gitBranchError, "Could not checkout branch " + options.branch)

  #configure and create tarball
  print "Configuring to create a tarball."
  configure_dir = os.path.join(options.workingDir, options.configurePath)
  source_dir = os.path.join(options.workingDir, "Trilinos")
  reportError(
    configure_cmake_build(
      configure_dir,
      source_dir,
      cmake_configure_options,
      file("configure.out", "w")),
    "Failure while configuring for tarball creation. See configure.out for details.",
    options.verboseErrors)
  
  print "Creating a tarball."
  reportError(run_command("make package_source",
                           configure_dir,
                           file("make_tarball.out", "w")),
              "Failure while creating tarball, see make_tarball.out for details.",
              "make_tarball.out",
              options.verboseErrors)
  
  tarball_build_path = os.path.join(options.workingDir, options.tarballBuildPath)
  if os.path.exists(tarball_build_path):
    print "Deleting directory " + tarball_build_path
    shutil.rmtree(tarball_build_path)

  os.mkdir(tarball_build_path)
  fileList = os.listdir(configure_dir)

  #find the tarball we need to copy
  for entry in fileList:
    if re.match("trilinos-.*-Source\.tar\.gz", entry):
      tarballFileName = entry
      print "found tarball: '" + tarballFileName + "'"
      break

  print "Expanding the tarball."
  shutil.copy(os.path.join(configure_dir, tarballFileName), tarball_build_path)
  tarball_path = os.path.join(tarball_build_path, tarballFileName)
  # Use command line tar because of broken tarfile module in python 2.4.
  reportError(run_command("tar -xzvf " + tarball_path, tarball_build_path, file("tar.out", "w")),
              "Could not expand tarball, see tar.out for details.", "tar.out", options.verboseErrors)

  #Make install
  print "Building and installing from the tarball."
  trilinosDir = tarballFileName[:-7]
  tarball_cmake_path = os.path.join(tarball_build_path, "build")
  cmake_configure_options.append(("Trilinos_ENABLE_DEVELOPMENT_MODE", "OFF"))
  error = configure_cmake_build(tarball_cmake_path,
                                os.path.join(tarball_build_path, trilinosDir),
                                cmake_configure_options, file("tarball_configure.out", "w"))

  reportError(error, "Configuring from tarball failed, see tarball_configure.out for details.",
              "tarball_configure.out", options.verboseErrors)

  tarballMakeInstallCmd = "make install -j " + options.buildProcs
  tarballMakeInstallError = run_command(tarballMakeInstallCmd, tarball_cmake_path,
                                        file("tarball_make_install.out", "w"))

  reportError(tarballMakeInstallError,
              "Make install from the tarball failed, see tarball_make_install.out for details.",
              "tarball_make_install.out", options.verboseErrors)

  if options.doDashboardBuild:
    makeExperimentalBuildCmd = "make dashboard -j " + options.buildProcs
    makeExperimentalBuildError = run_command(makeExperimentalBuildCmd, tarball_cmake_path,
                                             file("tarball_make_experimental.out", "w"))

    reportError(makeExperimentalBuildError,
                "Make experimental from the tarball failed, " +
                "see tarball_make_experimental.out for details.",
                "tarball_make_experimental.out", options.verboseErrors)

  print "Installation from a tarball completed successfully."

if __name__ == '__main__':
  options, args = get_options(os.getcwd())
  main(DEFAULT_ENABLE_DISABLE_LIST, options)
