#ifndef __HILBERT_CONST_H
#define __HILBERT_CONST_H

/* Bits per unsigned word */

#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  \
  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

extern void Zoltan_BSFC_fhsfc2d(double coord[], unsigned *nkey, unsigned key[]);
extern void Zoltan_BSFC_fhsfc3d(double coord[], unsigned *nkey, unsigned key[]);

#endif
