
typedef struct mpics_t {
    char *app_pathname;
    int   debug_level;

    int   gsize;
    int   grank;

    /* pcontrol - default is yes*/
    int enabled;

    MPI_Comm split_comm;
} mpics_t;


extern mpics_t mpics_data;


void mpics_log(char *fmt, ...);
