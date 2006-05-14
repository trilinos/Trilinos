<? 

# --------------------------------------------------- #
# File Locations                                      #
#                                                     #
# Specify here the location of the matrix to be       #
# uploaded and the script files.                      #
#                                                     #
# `MatrixDirectory' specifies where to store the      #
# uploaded matrices.                                  #
#                                                     #
# `ImageDirectory' specifies where to store the       #
# images created on-the-fly.                          #
#                                                     #
# `TempDirectory' is the directory used to store any  #
# temporary data.                                     #
#                                                     #
# `PythonDirectory' is the directory containing the   #
# step_process.py; typically, this is the ml/python/  #
# MatrixPortal directory.                             #
#                                                     #
# `HTMLImageDirectory' is the a web-accessible        #
# directory (as local path or as an absolute http     #
# address) that will contain the automatically        #
# generated images.                                   #
#                                                     #
# `PYTHONPATH' will be given to tbe python script     #
# before executing; typically it contains the location#
# of PyTrilinos and of PyChart.                       #
#                                                     #
# `LD_LIBRARY_PATH' can be needed to Trilinos         #
# if compiled as shared libraries.                    #
#                                                     #
# `ENABLE_MPI' can be either TRUE of FALSE. If TRUE,  #
# then PyTrilinos and Trilinos have been configured   #
# and compiled with MPI support. mpirun will be used  #
# to fire up all python scripts using PyTrilinos.     #
#                                                     #
# `MAX_PROCS', `MPI_BOOT'; and `MPI_HALT' specify the #
# maximum number of processors, what to call to ini-  #
# tialize and finalize MPI. For LAM/MPI, these values #
# are lamboot and lamhalt.                            #
# --------------------------------------------------- #
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$

$MACHINE = "givens4";

if ($MACHINE == "givens4")
{
  $HBMatrixDirectory = "/home/chinella/Web/MatrixPortal/HBMatrices/";
  $H5MatrixDirectory = "/home/chinella/Web/MatrixPortal/H5Matrices/";
  $ImageDirectory = "/home/chinella/Web/MatrixPortal/tmp/";
  $TempDirectory = "/tmp/";
  $PythonDirectory = "/home/chinella/Web/MatrixPortal/";
  $HTMLImageDirectory = "http://givens4.ethz.ch/MatrixPortal/tmp";
  $PYTHONPATH = "/home/masala/Trilinos/LINUX_SERIAL/lib/python2.3/site-packages/:/home/masala/lib/python2.3/site-packages/";
  $LD_LIBRARY_PATH = "/home/masala/Trilinos/LINUX_SERIAL/lib";
  $ENABLE_MPI = TRUE;
  $MAX_PROCS = 16;
  $MPI_BOOT = "lamboot";
  $MPI_HALT = "lamhalt";
  $MPI_RUN  = "mpirun";

  $INCLUDES = "-I/home/masala/Trilinos/LINUX_MPI/include";
  $LDFLAGS = "-L/home/masala/Trilinos/LINUX_MPI/lib";
  $LIBS = "-lml -lifpack -laztecoo -lamesos -lgaleri -lepetraext -lepetra -lteuchos -llapack -lblas";
  $CXX = "mpic++";
}
else if ($MACHINE == "kythira")
{
  $HBMatrixDirectory = "/var/www/html/MatrixPortal/HBMatrices/";
  $H5MatrixDirectory = "/var/www/html/MatrixPortal/H5Matrices/";
  $ImageDirectory = "/var/www/html/MatrixPortal/tmp/";
  $TempDirectory = "/tmp/";
  $PythonDirectory = "/home/msala/Trilinos/packages/ml/python/MatrixPortal";
  $HTMLImageDirectory = "http://kythira.ethz.ch/MatrixPortal/tmp";
  $PYTHONPATH = "/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/";
  $LD_LIBRARY_PATH = "/home/msala/Trilinos/LINUX_SERIAL/lib";
  $ENABLE_MPI = FALSE;
  $MAX_PROCS = 1;
  $MPI_BOOT = "";
  $MPI_HALT = "";

  $INCLUDES = "-I/home/msala/Trilinos/LINUX_SERIAL/include";
  $LDFLAGS = "-L/home/msala/lib -L/home/msala/Trilinos/LINUX_SERIAL/lib";
  $LIBS = "-lml -lifpack -laztecoo -lamesos -lgaleri -lepetraext -lepetra -lteuchos -llapack -lblas -lf2c -lm -lMAIN__";
  $CXX = "g++";
}
else if ($MACHINE == "aphrodite-serial")
{
  $HBMatrixDirectory = "/var/www/html/MatrixPortal/HBMatrices/";
  $H5MatrixDirectory = "/var/www/html/MatrixPortal/H5Matrices/";
  $ImageDirectory = "/var/www/html/tmp/";
  $TempDirectory = "/tmp/";
  $PythonDirectory = "/home/jhu/Trilinos/development-branch/Trilinos/packages/ml/python/MatrixPortal/";
  $HTMLImageDirectory = "http://localhost/MatrixPortal/tmp";
  $PYTHONPATH = "/home/jhu/Trilinos/development-branch/sandbox-serial/lib/python2.4/site-packages/";
  $LD_LIBRARY_PATH = "/home/jhu/Trilinos/development-branch/sandbox-serial";
  $ENABLE_MPI = FALSE;
  $MPI_ROOT_DIR = "";
  $MPI_ROOT_DIR = "";
  $MAX_PROCS = 1;
  $MPI_BOOT = "";
  $MPI_HALT = "";

  $INCLUDES = "-I/home/masala/Trilinos/LINUX_MPI/include";
  $LDFLAGS = "-L/home/masala/lib -L/home/masala/Trilinos/LINUX_MPI/lib";
  $LIBS = "-lml -lifpack -laztecoo -lamesos -lgaleri -lepetraext -lepetra -lteuchos -llapack -lblas -lf2c";
  $CXX = "mpic++";
}

# --------------------------------------------------- #
# Web Page Styles and Appearance                      #
#                                                     #
# Functions top_left() and top_right() can be used to #
# personalize the appearance of the web page.         #
# --------------------------------------------------- #

function top_left() { ?>
          <a href=<? global $dir; echo $dir; ?>index.html>
            <img src="<?php global $dir; echo $dir; ?>common/matrix_portal_logo.png" height=60 border=0 alt="The Matrix Portal" /></a>
<?  } ?>

<? function top_right() { ?>
<a href="http://www.ethz.ch/">
      <img src="<?php global $dir; echo $dir; ?>common/eth_logo.gif" border=0 alt="ETHZ" /></a>
<?  } ?>
