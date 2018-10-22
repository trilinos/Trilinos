#!/bin/bash
Function=$1                 #e.g. abs: function name
FunctionExtended=$2         #e.g. KokkosBlas1_impl_MV_abs: prefix for files etc.
Scalar=$3                   #e.g. double
Layout=$4                   #e.g. LayoutLeft
ExecSpace=$5                #e.g. OpenMP
MemSpace=$6                 #e.g. HBWSPACE
SlowMemSpace=$7                 #e.g. HostSpace
OffsetType=$8
OrdinalType=$9

filename_master_hpp=${10}      #e.g. Kokkos_Blas1_MV_impl_abs.hpp: where the actual function definition and declaration lives 
NameSpace=${11}                #e.g. KokkosBlas: namespace it lives in
KokkosKernelsPath=${12}        

Macro=`echo ${FunctionExtended} | awk '{print toupper($0)}'`
Scalar_UpperCase=`echo ${Scalar} | awk '{print toupper($0)}' | sed 's|\:\:|\_|g' | sed 's|<|_|g' | sed 's|>|_|g'`

Offset_UpperCase=`echo ${OffsetType} | awk '{print toupper($0)}' | sed 's|\:\:|\_|g' | sed 's|<|_|g' | sed 's|>|_|g'`
Ordinal_UpperCase=`echo ${OrdinalType} | awk '{print toupper($0)}' | sed 's|\:\:|\_|g' | sed 's|<|_|g' | sed 's|>|_|g'`

Scalar_FileName=`echo ${Scalar} | sed 's|\:\:|\_|g' | sed 's|<|_|g' | sed 's|>|_|g'`
Layout_UpperCase=`echo ${Layout} | awk '{print toupper($0)}'`
ExecSpace_UpperCase=`echo ${ExecSpace} | awk '{print toupper($0)}'`

prefix="Experimental::"
echo echo ${MemSpace#$prefix}
echo ${SlowMemSpace#$prefix}
MemSpace_UpperCase=`echo ${MemSpace#$prefix} | awk '{print toupper($0)}'`
SlowMemSpace_UpperCase=`echo ${SlowMemSpace#$prefix} | awk '{print toupper($0)}'`

#MemSpace_UpperCase=`echo ${MemSpace} | awk '{print toupper($0)}'`
#SlowMemSpace_UpperCase=`echo ${SlowMemSpace} | awk '{print toupper($0)}'`

OffsetType_FileName=`echo ${OffsetType} | sed 's|\ |\_|g'`
OrdinalType_FileName=`echo ${OrdinalType} | sed 's|\ |\_|g'`

filename_cpp=generated_specializations_cpp/${Function}/${FunctionExtended}_eti_spec_inst_${Scalar_FileName}_${OffsetType_FileName}_${OrdinalType_FileName}_${Layout}_${ExecSpace}_${MemSpace#$prefix}_${SlowMemSpace#$prefix}.cpp
filename_spec_avail_hpp=generated_specializations_hpp/${FunctionExtended}_eti_spec_avail.hpp
filename_spec_decl_hpp=generated_specializations_hpp/${FunctionExtended}_eti_spec_decl.hpp


cat ${KokkosKernelsPath}/scripts/header > ${filename_cpp}
echo "" >> ${filename_cpp}
echo "#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true" >> ${filename_cpp}
echo "#include \"KokkosKernels_config.h\"" >> ${filename_cpp}
echo "" >> ${filename_cpp}
echo "#if defined (KOKKOSKERNELS_INST_${Scalar_UpperCase}) \\" >> ${filename_cpp} 
echo " && defined (KOKKOSKERNELS_INST_${Layout_UpperCase}) \\" >> ${filename_cpp} 
echo " && defined (KOKKOSKERNELS_INST_EXECSPACE_${ExecSpace_UpperCase}) \\" >> ${filename_cpp} 
echo " && defined (KOKKOSKERNELS_INST_MEMSPACE_${MemSpace_UpperCase}) \\" >> ${filename_cpp} 
echo " && defined (KOKKOSKERNELS_INST_MEMSPACE_${SlowMemSpace_UpperCase}) \\" >> ${filename_cpp} 
echo " && defined (KOKKOSKERNELS_INST_ORDINAL_${Ordinal_UpperCase}) \\" >> ${filename_cpp} 
echo " && defined (KOKKOSKERNELS_INST_OFFSET_${Offset_UpperCase}) " >> ${filename_cpp}
echo "#include \"${filename_master_hpp}\"" >> ${filename_cpp}
echo "namespace ${NameSpace} {" >> ${filename_cpp}
echo "namespace Impl {" >> ${filename_cpp}
echo " ${Macro}_ETI_SPEC_INST(${Scalar}, ${OrdinalType}, ${OffsetType}, Kokkos::${Layout}, Kokkos::${ExecSpace}, Kokkos::${MemSpace}, Kokkos::${SlowMemSpace})" >> ${filename_cpp}
echo "} // Impl" >> ${filename_cpp} 
echo "} // ${NameSpace}" >> ${filename_cpp}
echo "#endif" >> ${filename_cpp}

echo "" >> ${filename_spec_avail_hpp}
echo "#if defined (KOKKOSKERNELS_INST_${Scalar_UpperCase}) \\" >> ${filename_spec_avail_hpp}
echo " && defined (KOKKOSKERNELS_INST_${Layout_UpperCase}) \\" >> ${filename_spec_avail_hpp}
echo " && defined (KOKKOSKERNELS_INST_EXECSPACE_${ExecSpace_UpperCase}) \\" >> ${filename_spec_avail_hpp}
echo " && defined (KOKKOSKERNELS_INST_MEMSPACE_${MemSpace_UpperCase}) \\" >> ${filename_spec_avail_hpp}
echo " && defined (KOKKOSKERNELS_INST_MEMSPACE_${SlowMemSpace_UpperCase}) \\" >> ${filename_spec_avail_hpp}
echo " && defined (KOKKOSKERNELS_INST_ORDINAL_${Ordinal_UpperCase}) \\" >> ${filename_spec_avail_hpp} 
echo " && defined (KOKKOSKERNELS_INST_OFFSET_${Offset_UpperCase}) " >> ${filename_spec_avail_hpp}
echo " ${Macro}_ETI_SPEC_AVAIL(${Scalar}, ${OrdinalType}, ${OffsetType}, Kokkos::${Layout}, Kokkos::${ExecSpace}, Kokkos::${MemSpace}, Kokkos::${SlowMemSpace})" >> ${filename_spec_avail_hpp}
echo "#endif" >> ${filename_spec_avail_hpp}

echo "" >> ${filename_spec_decl_hpp}
echo "#if defined (KOKKOSKERNELS_INST_${Scalar_UpperCase}) \\" >> ${filename_spec_decl_hpp}
echo " && defined (KOKKOSKERNELS_INST_${Layout_UpperCase}) \\" >> ${filename_spec_decl_hpp}
echo " && defined (KOKKOSKERNELS_INST_EXECSPACE_${ExecSpace_UpperCase}) \\" >> ${filename_spec_decl_hpp}
echo " && defined (KOKKOSKERNELS_INST_MEMSPACE_${MemSpace_UpperCase}) \\" >> ${filename_spec_decl_hpp}
echo " && defined (KOKKOSKERNELS_INST_MEMSPACE_${SlowMemSpace_UpperCase}) \\" >> ${filename_spec_decl_hpp}
echo " && defined (KOKKOSKERNELS_INST_ORDINAL_${Ordinal_UpperCase}) \\" >> ${filename_spec_decl_hpp} 
echo " && defined (KOKKOSKERNELS_INST_OFFSET_${Offset_UpperCase}) " >> ${filename_spec_decl_hpp}
echo " ${Macro}_ETI_SPEC_DECL(${Scalar}, ${OrdinalType}, ${OffsetType}, Kokkos::${Layout}, Kokkos::${ExecSpace}, Kokkos::${MemSpace}, Kokkos::${SlowMemSpace})" >> ${filename_spec_decl_hpp}
echo "#endif" >> ${filename_spec_decl_hpp}
