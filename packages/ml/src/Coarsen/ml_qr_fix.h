/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#ifndef ML_QR_FIX_H
#define ML_QR_FIX_H

/* If we need more than 16 kernel components, define ML_QR_FIX_TYPE
 * as unsigned int, otherwise use unsigned short int to conserve memory */
#define ML_QR_FIX_TYPE long int

typedef struct ML_qr_fix {

  int                 level;
  int                 numDeadNodDof;
 /* -mb: can later replace the following two with a hash structure */ 
  int                 nDeadNodDof; 
  ML_QR_FIX_TYPE     *xDeadNodDof;

} ML_qr_fix;

#ifdef __cplusplus
extern "C" {
  int ML_qr_fix_Create(const int nCoarseNod);

  int ML_qr_fix_Destroy(void);

  int ML_qr_fix_Print(ML_qr_fix* ptr);

  int ML_qr_fix_NumDeadNodDof(void);

  ML_QR_FIX_TYPE ML_qr_fix_getDeadNod(const int inx);

  void ML_qr_fix_setNumDeadNod(int num);

  void ML_qr_fix_setDeadNod( const int inx, ML_QR_FIX_TYPE val);

  int  ML_fixCoarseMtx(
          ML_Operator *Cmat,          /*-up- coarse operator in MSR format   */
          const int    CoarseMtxType  /*-in- coarse-lev mtx storage type     */
  );
 
  int  ML_qr_fix_Bitsize(void);
}
#else

int ML_qr_fix_Create(const int nCoarseNod);

int ML_qr_fix_Destroy(void);

int ML_qr_fix_Print(ML_qr_fix* ptr);

int ML_qr_fix_NumDeadNodDof(void);

ML_QR_FIX_TYPE ML_qr_fix_getDeadNod(const int inx);

void ML_qr_fix_setNumDeadNod(int num);

void ML_qr_fix_setDeadNod( const int inx, ML_QR_FIX_TYPE val);

int  ML_fixCoarseMtx(
        ML_Operator *Cmat,          /*-up- coarse operator in MSR format     */
        const int    CoarseMtxType  /*-in- coarse-lev mtx storage type       */
     );

int  ML_qr_fix_Bitsize(void);

#endif
#endif
