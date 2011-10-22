typedef struct commsplitter_t {
    char *app_pathname;
    int   debug_level;

    int   gsize;
    int   grank;

    /* pcontrol - default is yes*/
    int enabled;

    MPI_Comm split_comm;
} commsplitter_t;


extern commsplitter_t commsplitter_data;


void commsplitter_log(char *fmt, ...);
