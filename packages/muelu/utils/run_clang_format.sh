#!/bin/bash

# assume this is run from MueLu's top level
if [ ! -f "utils/run_clang_format.sh" ] ; then
  echo "Error! This script should be run from the packages/muelu directory! Please use"
  echo "  utils/run_clang_format.sh"
  exit 1
fi

# check if clang-format exists on the user's machine. if not, formatting can't be performed
if [ ! command -v clang-format &> /dev/null ] ; then
  echo "Error! clang-format was not found! Please place clang-format in your PATH or current directory"
  echo
  echo "If you are looking to reformat your code due to the MueLu clang-format Github action, please use clang-format 14.0.0"
  echo 'See https://github.com/llvm/llvm-project/releases/tag/llvmorg-14.0.0 to find the appropriate clang-format version'
  echo "Linux x64 developers can find it here:"
  echo '  https://github.com/llvm/llvm-project/releases/download/llvmorg-14.0.0/clang+llvm-14.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz'
  exit 1
else
  if [[ ! "$@" =~ "-f" ]] ; then
    # get the clang-format version by running clang-format --version
    CLANGFORMAT_VERSION=$(clang-format --version)
    # check if clang-format result contains the string 14.0.0. if not, error out
    if [[ ! "${CLANGFORMAT_VERSION}" =~ "14.0.0" ]] ; then
      echo "Error! clang-format is not version 14.0.0, so it may have conflicts with the MueLu style!"
      echo "Detected clang-format version: ${CLANGFORMAT_VERSION}"
      echo "If you would like to run this anyway, you can force run this with"
      echo "  ./utils/run_clang_format.sh -f"
      echo
      echo "If you are looking to reformat your code due to the MueLu clang-format Github action, please use clang-format 14.0.0"
      echo 'See https://github.com/llvm/llvm-project/releases/tag/llvmorg-14.0.0 to find the appropriate clang-format version'
      echo "Linux x64 developers can find it here:"
      echo '  https://github.com/llvm/llvm-project/releases/download/llvmorg-14.0.0/clang+llvm-14.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz'
      exit 1
    fi
  fi
fi

# grab all *.hpp files and *.cpp files in the MueLu src directory
HEADERS=$(find . -name "*.hpp")
SOURCES=$(find . -name "*.cpp")

echo "Running MueLu clang formatting..."

# apply clang-format inline to each **.hpp, utilizing the style file
for HEADER in ${HEADERS} ; do
  echo ${HEADER}
  clang-format -style=file -i ${HEADER}
done

# apply clang-format inline to each **.cpp, utilizing the style file
for SOURCE in ${SOURCES} ; do
  echo ${SOURCE}
  clang-format -style=file -i ${SOURCE}
done

echo "Done running MueLu clang formatting!"
