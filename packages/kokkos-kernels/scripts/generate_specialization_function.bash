#!/bin/bash

Function=$1             #e.g. abs: function name
FunctionExtended=$2     #e.g. KokkosBlas1_impl_MV_abs: prefix for files etc.
MasterHeader=$3         #e.g. Kokkos_Blas1_MV_impl_abs.hpp: where the actual function definition and declaration lives 
NameSpace=$4            #e.g. KokkosBlas: namespace it lives in 
KokkosKernelsPath=$5
ScalarList="double float Kokkos::complex<double> Kokkos::complex<float>"
LayoutList="LayoutLeft LayoutRight"
ExecMemSpaceList="Cuda,CudaSpace Cuda,CudaUVMSpace OpenMP,HostSpace Threads,HostSpace Serial,HostSpace"

mkdir generated_specializations_hpp
mkdir generated_specializations_cpp/${Function}
filename_hpp=generated_specializations_hpp/${FunctionExtended}_decl_specializations.hpp
Function_UpperCase=`echo ${FunctionExtended} | awk '{print toupper($0)}'`


echo "#ifndef ${Function_UpperCase}_DECL_SPECIALISATION_HPP_" > ${filename_hpp}
echo "#define ${Function_UpperCase}_DECL_SPECIALISATION_HPP_" >> ${filename_hpp}

cat ${KokkosKernelsPath}/scripts/header >> ${filename_hpp}

echo "namespace ${NameSpace} {" >> ${filename_hpp}
echo "namespace Impl {" >> ${filename_hpp}

for Scalar in ${ScalarList}; do
for Layout in ${LayoutList}; do
for ExecMemSpace in ${ExecMemSpaceList}; do
   ExecMemSpaceArray=(${ExecMemSpace//,/ })
   ExecSpace=${ExecMemSpaceArray[0]}
   MemSpace=${ExecMemSpaceArray[1]}
   echo "Generate: " ${FunctionExtended} " " ${Scalar} " " ${Layout} " " ${ExecSpace} " " ${MemSpace}
   ${KokkosKernelsPath}/scripts/generate_specialization_type.bash ${Function} ${FunctionExtended} ${Scalar} ${Layout} ${ExecSpace} ${MemSpace} ${MasterHeader} ${NameSpace} ${KokkosKernelsPath}
done
done
done

echo "} // Impl" >> ${filename_hpp}
echo "} // ${NameSpace}" >> ${filename_hpp}
echo "#endif // ${Function_UpperCase}_DECL_SPECIALISATION_HPP_" >> ${filename_hpp}

