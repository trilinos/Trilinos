## HCurl 1st order space
## 

ADD_SUBDIRECTORY(HCURL_QUAD_I1_FEM)
ADD_SUBDIRECTORY(HCURL_HEX_I1_FEM)
ADD_SUBDIRECTORY(HCURL_TRI_I1_FEM)
ADD_SUBDIRECTORY(HCURL_TET_I1_FEM)
ADD_SUBDIRECTORY(HCURL_WEDGE_I1_FEM)

## Hcurl high order space
##

ADD_SUBDIRECTORY(HCURL_QUAD_In_FEM)
ADD_SUBDIRECTORY(HCURL_HEX_In_FEM)
ADD_SUBDIRECTORY(HCURL_TRI_In_FEM)
ADD_SUBDIRECTORY(HCURL_TET_In_FEM)

## HDiv 1st order space
##

ADD_SUBDIRECTORY(HDIV_QUAD_I1_FEM)
ADD_SUBDIRECTORY(HDIV_HEX_I1_FEM)
ADD_SUBDIRECTORY(HDIV_TRI_I1_FEM)
ADD_SUBDIRECTORY(HDIV_TET_I1_FEM)
ADD_SUBDIRECTORY(HDIV_WEDGE_I1_FEM)

## Hdiv high order space
##

ADD_SUBDIRECTORY(HDIV_QUAD_In_FEM)
ADD_SUBDIRECTORY(HDIV_HEX_In_FEM)
ADD_SUBDIRECTORY(HDIV_TRI_In_FEM)
ADD_SUBDIRECTORY(HDIV_TET_In_FEM)

## HGrad 1st order space
##

ADD_SUBDIRECTORY(HGRAD_LINE_C1_FEM)
ADD_SUBDIRECTORY(HGRAD_QUAD_C1_FEM)
ADD_SUBDIRECTORY(HGRAD_HEX_C1_FEM)
ADD_SUBDIRECTORY(HGRAD_TRI_C1_FEM)
ADD_SUBDIRECTORY(HGRAD_TET_C1_FEM)
ADD_SUBDIRECTORY(HGRAD_WEDGE_C1_FEM)
ADD_SUBDIRECTORY(HGRAD_PYR_C1_FEM)

## HGrad 2nd order space
##

ADD_SUBDIRECTORY(HGRAD_LINE_C2_FEM)
ADD_SUBDIRECTORY(HGRAD_QUAD_C2_FEM)
ADD_SUBDIRECTORY(HGRAD_HEX_C2_FEM)
ADD_SUBDIRECTORY(HGRAD_TRI_C2_FEM)
ADD_SUBDIRECTORY(HGRAD_TET_C2_FEM)
ADD_SUBDIRECTORY(HGRAD_WEDGE_C2_FEM)
ADD_SUBDIRECTORY(HGRAD_PYR_I2_FEM)

## HGrad serendipity elements
## 

#ADD_SUBDIRECTORY(HGRAD_HEX_I2_FEM)
#ADD_SUBDIRECTORY(HGRAD_WEDGE_I2_FEM)
#ADD_SUBDIRECTORY(HGRAD_PYR_I2_FEM)

## HGrad high order space
##

ADD_SUBDIRECTORY(HGRAD_LINE_Cn_FEM_JACOBI)
ADD_SUBDIRECTORY(HGRAD_LINE_Cn_FEM)
ADD_SUBDIRECTORY(HGRAD_QUAD_Cn_FEM)
ADD_SUBDIRECTORY(HGRAD_HEX_Cn_FEM)
ADD_SUBDIRECTORY(HGRAD_TRI_Cn_FEM_ORTH)
ADD_SUBDIRECTORY(HGRAD_TRI_Cn_FEM)
ADD_SUBDIRECTORY(HGRAD_TET_Cn_FEM_ORTH)
ADD_SUBDIRECTORY(HGRAD_TET_Cn_FEM)

## HVOL 0th order space
## 

ADD_SUBDIRECTORY(HVOL_C0_FEM)

## HVOL high order space
##

ADD_SUBDIRECTORY(HVOL_LINE_Cn_FEM)
ADD_SUBDIRECTORY(HVOL_QUAD_Cn_FEM)
ADD_SUBDIRECTORY(HVOL_HEX_Cn_FEM)
ADD_SUBDIRECTORY(HVOL_TRI_Cn_FEM)
ADD_SUBDIRECTORY(HVOL_TET_Cn_FEM)

## HGrad composite element

ADD_SUBDIRECTORY(HGRAD_TET_COMP12_FEM)


## DeRham complex

ADD_SUBDIRECTORY(DeRHAM_TET_FEM)
ADD_SUBDIRECTORY(DeRHAM_HEX_FEM)
ADD_SUBDIRECTORY(DeRHAM_QUAD_FEM)
ADD_SUBDIRECTORY(DeRHAM_TRI_FEM)
