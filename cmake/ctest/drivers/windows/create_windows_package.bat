@rem Script to generate a Trilinos static-library package on Windows.
@rem Expected syntax to use the script is:

@rem     create_windows_package.bat <Build Type> <Base Directory>

@rem <Build Type> is either Debug or Release.
@rem <Base Directory> is the root directory where the script will clone
@rem    Trilinos, create a build directory, do the build, create a package, etc.

setlocal

rem Set the location of Git, Ninja, etc.
set SEMS_DIR=C:\projects\sems\install\win-x86_64
set NINJA_DIR=%SEMS_DIR%\utility\ninja\1.7.2
set CMAKE_DIR=%SEMS_DIR%\utility\cmake\3.12.0\bin
set GIT_EXE=%SEMS_DIR%\utility\git\2.13.0\cmd\git.exe
set PATH=%NINJA_DIR%;%PATH%

rem Get the script arguments
set BUILD_TYPE=%~1
set BASE_DIR=%~2

rem Setup directories
set BUILD_DIR="%BASE_DIR%\build"
set INSTALL_DIR="%BASE_DIR%\install"
set TRILINOS_DIR="%BASE_DIR%\Trilinos"

rem Cleanup old files in the base directory
if not exist "%BASE_DIR%" mkdir "%BASE_DIR%"
del /q "%BASE_DIR%"

rem Cleanup the old build directory
if exist %BUILD_DIR% rmdir /s /q %BUILD_DIR%
mkdir %BUILD_DIR%

rem Cleanup the old install directory. 
@rem Note: The install directory seems to be a necessary part of the 
@rem   cpack process for Trilinos, which is why we clean it up here and
@rem   specify CMAKE_INSTALL_PREFIX in the configure step below.
if exist %INSTALL_DIR% rmdir /s /q %INSTALL_DIR%

rem Clone or update Trilinos
> "%BASE_DIR%\update_output.txt" 2>&1 (

  if exist %TRILINOS_DIR% (
    cd %TRILINOS_DIR%
    %GIT_EXE% pull
  ) else (
    %GIT_EXE% clone git@github.com:trilinos/Trilinos.git %TRILINOS_DIR%
    cd %TRILINOS_DIR%
    %GIT_EXE% checkout develop
  )

)

rem Add the /bigobj compiler flag to get the debug build to work properly. 
rem As this flag doesn't hurt the release build, we don't bother checking the build type.
set CL=/bigobj

rem Setup the build environment before running CMake (this will help CMake find the right compilers)
cd %BUILD_DIR%
call "%SEMS_DIR%\compiler\VisualStudio\14.0\VC\vcvarsall.bat" amd64

rem Configure using CMake
> "%BASE_DIR%\configure_output.txt" 2>&1 (

  %CMAKE_DIR%\cmake.exe -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING=%BUILD_TYPE% ^
  -D CMAKE_INSTALL_PREFIX:PATH=%INSTALL_DIR% ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D Trilinos_ENABLE_CXX11:BOOL=ON ^
  -D Trilinos_ENABLE_TESTS:BOOL=OFF ^
  -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF ^
  -D Trilinos_ENABLE_FORTRAN:BOOL=OFF ^
  -D Trilinos_ENABLE_DEBUG:BOOL=OFF ^
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON ^
  -D Trilinos_GENERATE_REPO_VERSION_FILE:BOOL=OFF ^
  -D CMAKE_OBJECT_PATH_MAX=500 ^
  -D BLAS_LIBRARY_NAMES:STRING="blas;libf2c" ^
  -D BLAS_LIBRARY_DIRS:STRING="%SEMS_DIR%\tpl\clapack\3.2.1\VisualStudio\14.0\base\lib" ^
  -D LAPACK_LIBRARY_NAMES:STRING="lapack" ^
  -D LAPACK_LIBRARY_DIRS:STRING="%SEMS_DIR%\tpl\clapack\3.2.1\VisualStudio\14.0\base\lib" ^
  -D TPL_ENABLE_MPI:BOOL=ON ^
  -D MPI_BASE_DIR:PATH="%SEMS_DIR%\compiler\Microsoft MPI\8.1.12438.1084" ^
  -D PERL_EXECUTABLE:FILEPATH="%SEMS_DIR%\compiler\strawberry_perl\5.24.1.1\perl\bin\perl.exe" ^
  -D TPL_ENABLE_DLlib:BOOL=OFF ^
  -D TPL_ENABLE_Pthread:BOOL=OFF ^
  -D Trilinos_ENABLE_AztecOO:BOOL=ON ^
  -D Trilinos_ENABLE_Epetra:BOOL=ON ^
  -D Trilinos_ENABLE_ML:BOOL=ON ^
  -D Trilinos_ENABLE_MueLu:BOOL=ON ^
  -D Trilinos_ENABLE_Pamgen:BOOL=ON ^
  -D Trilinos_ENABLE_ROL:BOOL=ON ^
  -D Trilinos_ENABLE_Xpetra:BOOL=ON ^
  -D Trilinos_ENABLE_Zoltan:BOOL=ON ^
  -D Trilinos_ENABLE_Triutils:BOOL=ON ^
  -D Trilinos_ENABLE_Teuchos:BOOL=ON ^
  -D Trilinos_ENABLE_Belos:BOOL=ON ^
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON ^
  -D Trilinos_ENABLE_Amesos:BOOL=ON ^
  -D Trilinos_ENABLE_Ifpack:BOOL=ON ^
  -D Trilinos_ENABLE_TrilinosSS:BOOL=ON ^
  -D Trilinos_ENABLE_Galeri:BOOL=OFF ^
  -D Trilinos_ENABLE_Kokkos:BOOL=OFF ^
  -D Trilinos_ENABLE_Sacado:BOOL=OFF ^
  -D Trilinos_ENABLE_Intrepid:BOOL=OFF ^
  -D Trilinos_ENABLE_Thyra:BOOL=OFF ^
  %TRILINOS_DIR%

)

rem Build using Ninja
> "%BASE_DIR%\build_output.txt" 2>&1 (

  ninja.exe

)

rem Package it all up
> "%BASE_DIR%\package_output.txt" 2>&1 (

  %CMAKE_DIR%\cpack.exe -G ZIP
  
  rem Copy the generated zip file elsewhere
  copy *.zip "%BASE_DIR%"
)

endlocal
