
#define PETRA_FORTRAN

#ifdef PETRA_ADDRESS64BIT

#define PETRA_MATRIX long integer
#define PETRA_VECTOR long integer
#define PETRA_MULTIVECTOR long integer
#define PETRA_COMM long integer
#define PETRA_MAP long integer
#define PETRA_LOCALMAP long integer
#define PETRA_BLOCKMAP long integer
#define PETRA_LOCALBLOCKMAP long integer

#else

#define PETRA_MATRIX integer
#define PETRA_VECTOR integer
#define PETRA_MULTIVECTOR integer
#define PETRA_COMM integer
#define PETRA_MAP integer
#define PETRA_LOCALMAP integer
#define PETRA_BLOCKMAP integer
#define PETRA_LOCALBLOCKMAP integer

#endif

#define PETRA_ERROR_FLAG integer

      PETRA_COMM             petra_comm_create
      PETRA_COMM             petra_comm_create_serial

      PETRA_MAP              petra_map_create
      PETRA_LOCALMAP         petra_localmap_create
      PETRA_BLOCKMAP         petra_blockmap_create1
      PETRA_BLOCKMAP         petra_blockmap_create2
      PETRA_LOCALBLOCKMAP    petra_localblockmap_create1
      PETRA_LOCALBLOCKMAP    petra_localblockmap_create2

      PETRA_MATRIX           petra_rdp_dcrs_matrix_create
      PETRA_MATRIX           petra_rdp_dvbr_matrix_create

      PETRA_VECTOR           petra_rdp_vector_create
      PETRA_MULTIVECTOR     petra_rdp_multivector_create

      PETRA_ERROR_FLAG  petra_comm_getmypid
      PETRA_ERROR_FLAG  petra_comm_getnumproc

      PETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_allocate
      PETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_putblockrow
      PETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_fillcomplete
      PETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_matvec

      PETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_allocate
      PETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_putrow
      PETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_sumintodiagonal
      PETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_fillcomplete
      PETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_matvec

      PETRA_ERROR_FLAG  petra_rdp_vector_putvector
      PETRA_ERROR_FLAG  petra_rdp_vector_lincomb
      PETRA_ERROR_FLAG  petra_rdp_vector_norm2
