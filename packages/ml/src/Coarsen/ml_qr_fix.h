#ifndef ML_QR_FIX_H
#define ML_QR_FIX_H

typedef struct ML_qr_fix {

  int level;

} ML_qr_fix;

int ML_qr_fix_Create(ML_qr_fix** ptr);

int ML_qr_fix_Destroy(ML_qr_fix** ptr);

int ML_qr_fix_Print(ML_qr_fix* ptr);

#endif
