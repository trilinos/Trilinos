#!/bin/bash

Function=$1             #e.g. abs: function name
FunctionExtended=$2     #e.g. KokkosBlas1_abs: prefix for files etc.
MasterHeader=$3         #e.g. Kokkos_Blas1_abs_spec.hpp: where the specialization layer lives
NameSpace=$4            #e.g. KokkosBlas: namespace it lives in
KokkosKernelsPath=$5
OrdinalList="int int64_t"
OffsetList="int size_t"
ScalarList="double float Kokkos::complex<double> Kokkos::complex<float>"
LayoutList="LayoutLeft LayoutRight"
ExecMemSpaceList="Cuda,CudaSpace Cuda,CudaUVMSpace OpenMP,HostSpace Threads,HostSpace Serial,HostSpace OpenMP,Experimental::HBWSpace Threads,Experimental::HBWSpace Serial,Experimental::HBWSpace"

mkdir generated_specializations_hpp
mkdir generated_specializations_cpp/${Function}
filename_hpp_root=generated_specializations_hpp/${FunctionExtended}_eti_spec
filename_spec_avail_hpp=${filename_hpp_root}_avail.hpp
filename_spec_decl_hpp=${filename_hpp_root}_decl.hpp
Function_UpperCase=`echo ${FunctionExtended} | awk '{print toupper($0)}'`

# Start eti_spec_avail file:
echo "#ifndef ${Function_UpperCase}_ETI_SPEC_AVAIL_HPP_" > ${filename_spec_avail_hpp}
echo "#define ${Function_UpperCase}_ETI_SPEC_AVAIL_HPP_" >> ${filename_spec_avail_hpp}

cat ${KokkosKernelsPath}/scripts/header >> ${filename_spec_avail_hpp}

echo "namespace ${NameSpace} {" >> ${filename_spec_avail_hpp}
echo "namespace Impl {" >> ${filename_spec_avail_hpp}

# Start eti_spec_decl file:
echo "#ifndef ${Function_UpperCase}_ETI_SPEC_DECL_HPP_" > ${filename_spec_decl_hpp}
echo "#define ${Function_UpperCase}_ETI_SPEC_DECL_HPP_" >> ${filename_spec_decl_hpp}

cat ${KokkosKernelsPath}/scripts/header >> ${filename_spec_decl_hpp}

echo "namespace ${NameSpace} {" >> ${filename_spec_decl_hpp}
echo "namespace Impl {" >> ${filename_spec_decl_hpp}

for Scalar in ${ScalarList}; do
for Layout in ${LayoutList}; do
for ExecMemSpace in ${ExecMemSpaceList}; do
for Offset in ${OffsetList}; do
for Ordinal in ${OrdinalList}; do
   ExecMemSpaceArray=(${ExecMemSpace//,/ })
   ExecSpace=${ExecMemSpaceArray[0]}
   MemSpace=${ExecMemSpaceArray[1]}
   echo "Generate: " ${FunctionExtended} " " ${Scalar} " " ${Layout} " " ${ExecSpace} " " ${MemSpace} " "  ${Offset} " " ${Ordinal}
   #echo "${KokkosKernelsPath}/scripts/generate_specialization_type_sparse.bash ${Function} ${FunctionExtended} ${Scalar} ${Layout} ${ExecSpace} ${MemSpace} ${Offset} ${Ordinal} ${MasterHeader} ${NameSpace} ${KokkosKernelsPath}"
   ${KokkosKernelsPath}/scripts/generate_specialization_type_sparse.bash ${Function} ${FunctionExtended} ${Scalar} ${Layout} ${ExecSpace} ${MemSpace} ${Offset} ${Ordinal} ${MasterHeader} ${NameSpace} ${KokkosKernelsPath}
done
done
done
done
done

echo "} // Impl" >> ${filename_spec_avail_hpp}
echo "} // ${NameSpace}" >> ${filename_spec_avail_hpp}
echo "#endif // ${Function_UpperCase}_ETI_SPEC_AVAIL_HPP_" >> ${filename_spec_avail_hpp}

echo "} // Impl" >> ${filename_spec_decl_hpp}
echo "} // ${NameSpace}" >> ${filename_spec_decl_hpp}
echo "#endif // ${Function_UpperCase}_ETI_SPEC_DECL_HPP_" >> ${filename_spec_decl_hpp}
