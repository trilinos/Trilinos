#include <time.h>

/* Constants used in timer routines */
#define TIME_WALL 1
#define TIME_CPU 2

/* Function prototypes */
double LB_Time();
double LB_Time_Resolution();
void LB_Print_Time (LB *, double, char *);
int LB_Set_Timer_Param(char *, char *);

