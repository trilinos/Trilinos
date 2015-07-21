#!/bin/bash

package_dirs=`ls`
#package_dirs="teuchos nox moocho"

base_dir=$PWD

output="pkg_dir, code_size_mb, pkg_size_mb, ratio code/pkg" 
for package_dir in $package_dirs ; do
  #echo "Processing package_${package_dir}' ..."
  if [ -e ${package_dir}/CMakeLists.txt ] ; then 
    cd $package_dir
    src_size_bytes=$(cat $(find . -name "*.cpp" && find . -name "*.hpp" && find . -name "*.h" && find . -name "*.c" && find . -name "*.C" && find . -name "*.H") | wc -c)
    pkg_size_bytes=$(cat $(find . -type f) | wc -c)
    ratio_src_pkg=$(echo "scale=2; ${src_size_bytes}.0/${pkg_size_bytes}.0" | bc)
    src_size_mb=$(echo "scale=2; ${src_size_bytes} / 1024 / 1000" | bc)
    pkg_size_mb=$(echo "scale=2; ${pkg_size_bytes} / 1024 / 1000" | bc)
    output="${output}\n${package_dir},${src_size_mb},${pkg_size_mb},${ratio_src_pkg}"
  fi
  cd $base_dir
done
 
echo -e $output | column -t -s ', '
