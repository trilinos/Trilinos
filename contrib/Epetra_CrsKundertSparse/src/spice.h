/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *	Portability, global externs, kitchen sink (yuk).
 *	In future releases this will include things from misc.h and util.h,
 *	   which duplicate each other in places
 */

#include <math.h>

#ifndef	M_PI
#  define M_PI		3.14159265358979323846
#endif
#ifndef	M_E
#  define M_E  	   2.7182818284590452354
#endif
#ifndef	M_LOG2E
#  define M_LOG2E		1.4426950408889634074
#endif
#ifndef	M_LOG10E
#  define M_LOG10E        0.43429448190325182765
#endif
#define TS_LIMIT 5

#include "port.h"
#include "hw.h"
#include "config.h"
#include "capabil.h"

#define	NUMELEMS(ARRAY)	(sizeof(ARRAY)/sizeof(*ARRAY))

extern char *Spice_Exec_Dir;
extern char *Spice_Lib_Dir;
extern char *Spice_Help_Dir;
extern char *Spice_Model_Dir;
extern char Spice_OptChar;
extern char *Def_Editor;
extern char *Bug_Addr;
extern int AsciiRawFile;
extern char *Spice_Host;
extern char *Spiced_Log;

extern char Spice_Notice[ ];
extern char Spice_Version[ ];
extern char Spice_Build_Date[ ];

extern char *News_File;
extern char *Default_MFB_Cap;
extern char *Spice_Path;
extern char *Help_Path;
extern char *Lib_Path;
extern int  Patch_Level;

#ifdef MAIN_PROGRAM
    int report_interval, new_raw_head, load_mode, mat_dense;
    int device_error, model_error, *timer_flag, dev_math_error;
    double *simulation_time, *logic_break, *timer_calibration;
    double *accepted_simulation_time;
    int exp_num, exp_d;
    void **exp_list, **exp_values;

    char *xfile;
#else
    extern int report_interval, new_raw_head, load_mode, mat_dense;
    extern int device_error, model_error, *timer_flag, dev_math_error;
    extern double *simulation_time, *logic_break, *timer_calibration;
    extern double *accepted_simulation_time;
    extern char *xfile;
    extern int exp_num, exp_d;
    extern void **exp_list, **exp_values;
#endif

#ifdef SHARED_MEM
/*
void update_SM (double *, double);
*/
#endif /* SHARED_MEM */

#define LOAD_NORMAL 1
#define LOAD_ENERGY 2
