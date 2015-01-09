#!/bin/bash -e

SOURCE_DIR=${ROL_SOURCE:-`pwd`}
INSTALL_DIR=${ROL_HOME:-"/usr/local"}

echo ""
echo "***************************************************"
echo "*** Installing Rapid Optimization Library (ROL) ***"
echo "***        Header-only installation             ***"
echo "***************************************************"

echo ""

echo "Main source directory (where all ROL directories live):"
echo "---->" $SOURCE_DIR
echo ""
echo "Source /src directory (where we'll get headers from):"
echo "---->" $SOURCE_DIR/src
echo ""
echo "Install directory (the main installation directory):"
echo "---->" $INSTALL_DIR
echo ""
echo "Include directory (where we'll install headers):"
echo "---->" $INSTALL_DIR/include

echo ""
echo "Let's start ..."
echo ""

if [ ! -d "$INSTALL_DIR" ]; then
  echo "WARNING: Install directory does not exist."
  echo "         Creating install directory ..."
  echo "mkdir $INSTALL_DIR"
  mkdir $INSTALL_DIR
  echo ""
fi

if [ ! -w "$INSTALL_DIR" ]; then
  echo "ERROR:   You don't have write permissions in" $INSTALL_DIR
  echo "         Have you set the ROL_HOME environment variable?"
  echo ""
  echo "***************************************************"
  echo "***         Install of ROL failed!!!            ***"
  echo "***************************************************"
  echo ""
  false
fi

if [ ! -d "$INSTALL_DIR/include" ]; then
  echo "WARNING: Include directory does not exist."
  echo "         Creating include directory ..."
  echo "mkdir $INSTALL_DIR/include"
  mkdir $INSTALL_DIR/include
  echo ""
fi

if [ ! -w "$INSTALL_DIR/include" ]; then
    echo "ERROR:   You don't have write permissions in" $INSTALL_DIR/include
    echo "         Have you set the ROL_HOME environment variable?"
    echo ""
    echo "***************************************************"
    echo "***         Install of ROL failed!!!            ***"
    echo "***************************************************"
    echo ""
    false
fi

if [ ! -f "$SOURCE_DIR/src/function/ROL_Objective.hpp" ]; then
  echo "ERROR:   Can't recognize directory " $SOURCE_DIR/src
  echo "         or the header files contained in it."
  echo "         Have you set the ROL_SOURCE environment variable?"
  echo "         If you haven't, you must install from the"
  echo "         directory that contains the ROL /src directory."
  echo ""
  echo "***************************************************"
  echo "***         Install of ROL failed!!!            ***"
  echo "***************************************************"
  echo ""
  false
fi

if [ -d "$INSTALL_DIR/include" ]; then
  echo "--Grabbing header files from: " $SOURCE_DIR/src
  echo "-----Copying header files to: " $INSTALL_DIR/include
fi

echo ""

find $SOURCE_DIR/src/ -name "*.hpp" -exec cp {} $INSTALL_DIR/include/. \;

cat $SOURCE_DIR/ROL.txt

echo "***************************************************"
echo "***       Install of ROL succeeded!!!           ***"
echo "***************************************************"
echo ""
