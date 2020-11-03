/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _PE_STR_UTIL_CONST_H_
#define _PE_STR_UTIL_CONST_H_
/* Function prototypes */
extern int token_compare(char *      token, /* The input character string */
                         const char *key    /* The key to compare with token */
);

extern void strip_string(char        inp_str[], /* The string to strip */
                         const char *tokens     /* The tokens to strip from the beginning and
                                                 * end of the input string */
);

extern void clean_string(char        inp_str[], /* The string to clean */
                         const char *tokens     /* The tokens to strip multiple copies of */
);

extern void string_to_lower(char inp_str[], /* The string to convert to lower case */
                            char cstop      /* Character where to stop */
);

#endif /* _PE_STR_UTIL_CONST_H_ */
