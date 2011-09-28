#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <limits.h>
#include <string.h>

#include <zlib.h>

#include <mpi.h>

#include "commsplitter.h"
#include "fsymbols.h"



#define SUBSTITUTE_COMM \
{\
        orig_comm = *comm;\
        MPI_Comm_compare(*comm, MPI_COMM_WORLD, &compare_result);\
        if (compare_result == MPI_IDENT) {\
            *comm = mpics_data.split_comm;\
        }\
}

#define RESTORE_COMM \
{\
    *comm = orig_comm;\
}




char *get_exe_name(char *pathname)
{
    char *slash;

    slash=rindex(pathname, '/')+1;

    return(slash);
}

char *get_app_pathname_from_proc(void)
{
    int exelen, bufsize=PATH_MAX;
    char *buf = NULL;

    buf = malloc(bufsize);
    if (buf == NULL) {
        mpics_log("unable to allocate space for full executable path.\n");
        PMPI_Abort(MPI_COMM_WORLD, -1);
    }

    exelen = readlink("/proc/self/exe", buf, PATH_MAX);
    if (exelen == -1) {
        free(buf);
    } else {
        buf[exelen] = '\0';
        return(buf);
    }

    return NULL;
}

void get_app_args_from_proc(int *argc, char **argv, int max_args)
{
    int i=0;
    char *buf;
    FILE *f;
    char *arg;

    *argc   = 0;
    *argv = NULL;

    buf = malloc(PATH_MAX);
    f = fopen("/proc/self/cmdline", "r");
    if (f != NULL) {
        if (fread(buf, 1, PATH_MAX, f) > 0) {
            arg = buf;
            while(*arg != '\0') {
                argv[i] = strdup(arg);
                arg += strlen(argv[i]) + 1;
                i++;
                if (i==max_args) {
                    mpics_log("ERROR: too many args.  truncating.");
                    break;
                }
            }
        }
        *argc = i;
    }

    free(buf);
    fclose(f);
}


static int mpics_MPI_Init(int *argc, char ***argv)
{
    int rc = 0;
    int enabled_save;
    uint32_t color=0;

    enabled_save = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Init(argc, argv);

    mpics_data.enabled = enabled_save;

    mpics_data.app_pathname = get_app_pathname_from_proc();
    mpics_init(get_exe_name(mpics_data.app_pathname));
    mpics_log("app_pathname is %s\n", mpics_data.app_pathname);


    if(mpics_data.app_pathname == NULL) {
        mpics_log("Unable to determine application's full pathname.  Cannot split MPI_COMM_WORLD.\n");
        PMPI_Abort(MPI_COMM_WORLD, -1);
    } else {
        color = crc32(0L, Z_NULL, 0);
        color = crc32(color, (Bytef *)mpics_data.app_pathname, strlen(mpics_data.app_pathname));
        mpics_log("crc32 for app_pathname(%s) is %lu\n", mpics_data.app_pathname, color);
    }

    enabled_save = mpics_data.enabled;
    mpics_data.enabled = 0;
    rc=PMPI_Comm_split(MPI_COMM_WORLD, color, mpics_data.grank, &mpics_data.split_comm);
    mpics_data.enabled = enabled_save;

    return(rc);
}

extern int MPI_Init(int *argc, char ***argv)
{
    int rc = 0;

    rc = mpics_MPI_Init(argc, argv);

    return(rc);
}

extern void F77_MPI_INIT(int *ierr)
{
    int   rc = 0;
    int   argc;
    char **argv;

    argv=(char **)malloc(256*sizeof(char*));
    get_app_args_from_proc(&argc, argv, 256);

    rc = mpics_MPI_Init(&argc,(char ***) &argv);
    *ierr = rc;

    return;
}


static int mpics_MPI_Init_thread(int *argc, char ***argv, int required, int *provided)
{
    int rc = 0;
    int enabled_save;
    uint32_t color=0;

    enabled_save = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Init_thread(argc, argv, required, provided);

    mpics_data.enabled = enabled_save;

    mpics_data.app_pathname = get_app_pathname_from_proc();
    mpics_init(get_exe_name(mpics_data.app_pathname));
    mpics_log("app_pathname is %s\n", mpics_data.app_pathname);


    if(mpics_data.app_pathname == NULL) {
        mpics_log("Unable to determine application's full pathname.  Cannot split MPI_COMM_WORLD.\n");
        PMPI_Abort(MPI_COMM_WORLD, -1);
    } else {
        color = crc32(0L, Z_NULL, 0);
        color = crc32(color, (Bytef *)mpics_data.app_pathname, strlen(mpics_data.app_pathname));
        mpics_log("crc32 for app_pathname(%s) is %lu\n", mpics_data.app_pathname, color);
    }

    enabled_save = mpics_data.enabled;
    mpics_data.enabled = 0;
    rc=PMPI_Comm_split(MPI_COMM_WORLD, color, mpics_data.grank, &mpics_data.split_comm);
    mpics_data.enabled = enabled_save;

    return(rc);
}

extern int MPI_Init_thread(int *argc, char ***argv, int required, int *provided)
{
    int rc = 0;

    rc = mpics_MPI_Init_thread(argc, argv, required, provided);

    return(rc);
}

extern void F77_MPI_INIT_THREAD(int *required, int *provided, int *ierr)
{
    int   rc = 0;
    int   argc;
    char **argv;

    argv=(char **)malloc(256*sizeof(char*));
    get_app_args_from_proc(&argc, argv, 256);

    rc = mpics_MPI_Init_thread(&argc, (char ***) &argv, *required, provided);
    *ierr = rc;

    return;
}


static int mpics_MPI_Finalize()
{
    int rc = 0;

    mpics_finalize();
    mpics_data.enabled = 0;
    mpics_log("calling PMPI_Finalize\n");
    rc = PMPI_Finalize();
    mpics_log("returning from PMPI_Finalize\n");

    return(rc);
}

extern int MPI_Finalize(void)
{
    int rc = 0;

    rc = mpics_MPI_Finalize();

    return(rc);
}

extern void F77_MPI_FINALIZE(int *ierr)
{
    int rc = 0;

    rc = mpics_MPI_Finalize();
    *ierr = rc;

    return;
}

static int mpics_MPI_Pcontrol(
        const int flag,
        va_list args)
{
    mpics_log("MPI_Pcontrol enter: flag=%d\n", flag);

    if (flag == 0) {
        if (!mpics_data.enabled) {
            mpics_log("WARN: commsplitter already disabled\n");
        }
        mpics_data.enabled = 0;
    } else {
        if (mpics_data.enabled) {
            mpics_log("WARN: commsplitter already enabled\n");
        }
        mpics_data.enabled = 1;
    }

    return MPI_SUCCESS;
}

int MPI_Pcontrol(
        const int flag,
        ...)
{
    int rc=0;
    va_list args;

    va_start(args, flag);
    rc = mpics_MPI_Pcontrol(flag, args);
    va_end(args);

    return(rc);
}

int F77_MPI_PCONTROL(
        int *flag,
        ...)
{
    int rc=0;
    va_list args;

    va_start(args, flag);
    rc = mpics_MPI_Pcontrol(*flag, args);
    va_end(args);

    return(rc);
}




static int mpics_MPI_Allgather(
        void *sendbuf,
        int *sendcount,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcount,
        MPI_Datatype *recvtype,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Allgather(sendbuf, *sendcount, *sendtype, recvbuf, *recvcount, *recvtype, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Allgather(
        void *sendbuf,
        int sendcount,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcount,
        MPI_Datatype recvtype,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Allgather(sendbuf, &sendcount, &sendtype, recvbuf, &recvcount, &recvtype, &comm);

    return(rc);
}


extern void F77_MPI_ALLGATHER(
        void *sendbuf,
        int *sendcount,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcount,
        MPI_Fint *recvtype,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Allgather(sendbuf, sendcount, &c_sendtype, recvbuf, recvcount, &c_recvtype, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Allgatherv(
        void *sendbuf,
        int *sendcount,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcounts,
        int *displs,
        MPI_Datatype *recvtype,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Allgatherv(sendbuf, *sendcount, *sendtype, recvbuf, recvcounts, displs, *recvtype, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Allgatherv(
        void *sendbuf,
        int sendcount,
        MPI_Datatype sendtype,
        void *recvbuf,
        int *recvcounts,
        int *displs,
        MPI_Datatype recvtype,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Allgatherv(sendbuf, &sendcount, &sendtype, recvbuf, recvcounts, displs, &recvtype, &comm);

    return(rc);
}


extern void F77_MPI_ALLGATHERV(
        void *sendbuf,
        int *sendcount,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcounts,
        int *displs,
        MPI_Fint *recvtype,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Allgatherv(sendbuf, sendcount, &c_sendtype, recvbuf, recvcounts, displs, &c_recvtype, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Allreduce(
        void *sendbuf,
        void *recvbuf,
        int *count,
        MPI_Datatype *datatype,
        MPI_Op *op,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Allreduce(sendbuf, recvbuf, *count, *datatype, *op, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Allreduce(
        void *sendbuf,
        void *recvbuf,
        int count,
        MPI_Datatype datatype,
        MPI_Op op,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Allreduce(sendbuf, recvbuf, &count, &datatype, &op, &comm);

    return(rc);
}


extern void F77_MPI_ALLREDUCE(
        void *sendbuf,
        void *recvbuf,
        int *count,
        MPI_Fint *datatype,
        MPI_Fint *op,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Allreduce(sendbuf, recvbuf, count, &c_datatype, &c_op, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Alltoall(
        void *sendbuf,
        int *sendcount,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcnt,
        MPI_Datatype *recvtype,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Alltoall(sendbuf, *sendcount, *sendtype, recvbuf, *recvcnt, *recvtype, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Alltoall(
        void *sendbuf,
        int sendcount,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcnt,
        MPI_Datatype recvtype,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Alltoall(sendbuf, &sendcount, &sendtype, recvbuf, &recvcnt, &recvtype, &comm);

    return(rc);
}


extern void F77_MPI_ALLTOALL(
        void *sendbuf,
        int *sendcount,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcnt,
        MPI_Fint *recvtype,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Alltoall(sendbuf, sendcount, &c_sendtype, recvbuf, recvcnt, &c_recvtype, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Alltoallv(
        void *sendbuf,
        int *sendcnts,
        int *sdispls,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcnts,
        int *rdispls,
        MPI_Datatype *recvtype,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Alltoallv(sendbuf, sendcnts, sdispls, *sendtype, recvbuf, recvcnts, rdispls, *recvtype, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Alltoallv(
        void *sendbuf,
        int *sendcnts,
        int *sdispls,
        MPI_Datatype sendtype,
        void *recvbuf,
        int *recvcnts,
        int *rdispls,
        MPI_Datatype recvtype,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Alltoallv(sendbuf, sendcnts, sdispls, &sendtype, recvbuf, recvcnts, rdispls, &recvtype, &comm);

    return(rc);
}


extern void F77_MPI_ALLTOALLV(
        void *sendbuf,
        int *sendcnts,
        int *sdispls,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcnts,
        int *rdispls,
        MPI_Fint *recvtype,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Alltoallv(sendbuf, sendcnts, sdispls, &c_sendtype, recvbuf, recvcnts, rdispls, &c_recvtype, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Attr_delete(
        MPI_Comm *comm,
        int *keyval)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Attr_delete(*comm, *keyval);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Attr_delete(
        MPI_Comm comm,
        int keyval)
{
    int rc;


    rc = mpics_MPI_Attr_delete(&comm, &keyval);

    return(rc);
}


extern void F77_MPI_ATTR_DELETE(
        MPI_Fint *comm,
        int *keyval,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Attr_delete(&c_comm, keyval);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Attr_get(
        MPI_Comm *comm,
        int *keyval,
        void *attr_value,
        int *flag)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Attr_get(*comm, *keyval, attr_value, flag);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Attr_get(
        MPI_Comm comm,
        int keyval,
        void *attr_value,
        int *flag)
{
    int rc;


    rc = mpics_MPI_Attr_get(&comm, &keyval, attr_value, flag);

    return(rc);
}


extern void F77_MPI_ATTR_GET(
        MPI_Fint *comm,
        int *keyval,
        void *attr_value,
        int *flag,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Attr_get(&c_comm, keyval, attr_value, flag);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Attr_put(
        MPI_Comm *comm,
        int *keyval,
        void *attr_value)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Attr_put(*comm, *keyval, attr_value);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Attr_put(
        MPI_Comm comm,
        int keyval,
        void *attr_value)
{
    int rc;


    rc = mpics_MPI_Attr_put(&comm, &keyval, attr_value);

    return(rc);
}


extern void F77_MPI_ATTR_PUT(
        MPI_Fint *comm,
        int *keyval,
        void *attr_value,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Attr_put(&c_comm, keyval, attr_value);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Barrier(
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Barrier(*comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Barrier(
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Barrier(&comm);

    return(rc);
}


extern void F77_MPI_BARRIER(
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Barrier(&c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Bcast(
        void *buffer,
        int *count,
        MPI_Datatype *datatype,
        int *root,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Bcast(buffer, *count, *datatype, *root, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Bcast(
        void *buffer,
        int count,
        MPI_Datatype datatype,
        int root,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Bcast(buffer, &count, &datatype, &root, &comm);

    return(rc);
}


extern void F77_MPI_BCAST(
        void *buffer,
        int *count,
        MPI_Fint *datatype,
        int *root,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Bcast(buffer, count, &c_datatype, root, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Bsend(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Bsend(buf, *count, *datatype, *dest, *tag, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Bsend(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Bsend(buf, &count, &datatype, &dest, &tag, &comm);

    return(rc);
}


extern void F77_MPI_BSEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Bsend(buf, count, &c_datatype, dest, tag, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Bsend_init(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Bsend_init(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Bsend_init(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Bsend_init(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_BSEND_INIT(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Bsend_init(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Cart_coords(
        MPI_Comm *comm,
        int *rank,
        int *maxdims,
        int *coords)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_coords(*comm, *rank, *maxdims, coords);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_coords(
        MPI_Comm comm,
        int rank,
        int maxdims,
        int *coords)
{
    int rc;


    rc = mpics_MPI_Cart_coords(&comm, &rank, &maxdims, coords);

    return(rc);
}


extern void F77_MPI_CART_COORDS(
        MPI_Fint *comm,
        int *rank,
        int *maxdims,
        int *coords,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_coords(&c_comm, rank, maxdims, coords);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Cart_create(
        MPI_Comm *comm,
        int *ndims,
        int *dims,
        int *periods,
        int *reorder,
        MPI_Comm *comm_cart)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_create(*comm, *ndims, dims, periods, *reorder, comm_cart);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_create(
        MPI_Comm comm,
        int ndims,
        int *dims,
        int *periods,
        int reorder,
        MPI_Comm *comm_cart)
{
    int rc;


    rc = mpics_MPI_Cart_create(&comm, &ndims, dims, periods, &reorder, comm_cart);

    return(rc);
}


extern void F77_MPI_CART_CREATE(
        MPI_Fint *comm,
        int *ndims,
        int *dims,
        int *periods,
        int *reorder,
        MPI_Fint *comm_cart,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Comm c_comm_cart;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_create(&c_comm, ndims, dims, periods, reorder, &c_comm_cart);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *comm_cart = MPI_Comm_c2f(c_comm_cart);
    }
    return;
}





static int mpics_MPI_Cart_get(
        MPI_Comm *comm,
        int *maxdims,
        int *dims,
        int *periods,
        int *coords)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_get(*comm, *maxdims, dims, periods, coords);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_get(
        MPI_Comm comm,
        int maxdims,
        int *dims,
        int *periods,
        int *coords)
{
    int rc;


    rc = mpics_MPI_Cart_get(&comm, &maxdims, dims, periods, coords);

    return(rc);
}


extern void F77_MPI_CART_GET(
        MPI_Fint *comm,
        int *maxdims,
        int *dims,
        int *periods,
        int *coords,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_get(&c_comm, maxdims, dims, periods, coords);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Cart_map(
        MPI_Comm *comm,
        int *ndims,
        int *dims,
        int *periods,
        int *newrank)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_map(*comm, *ndims, dims, periods, newrank);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_map(
        MPI_Comm comm,
        int ndims,
        int *dims,
        int *periods,
        int *newrank)
{
    int rc;


    rc = mpics_MPI_Cart_map(&comm, &ndims, dims, periods, newrank);

    return(rc);
}


extern void F77_MPI_CART_MAP(
        MPI_Fint *comm,
        int *ndims,
        int *dims,
        int *periods,
        int *newrank,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_map(&c_comm, ndims, dims, periods, newrank);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Cart_rank(
        MPI_Comm *comm,
        int *coords,
        int *rank)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_rank(*comm, coords, rank);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_rank(
        MPI_Comm comm,
        int *coords,
        int *rank)
{
    int rc;


    rc = mpics_MPI_Cart_rank(&comm, coords, rank);

    return(rc);
}


extern void F77_MPI_CART_RANK(
        MPI_Fint *comm,
        int *coords,
        int *rank,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_rank(&c_comm, coords, rank);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Cart_shift(
        MPI_Comm *comm,
        int *direction,
        int *displ,
        int *source,
        int *dest)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_shift(*comm, *direction, *displ, source, dest);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_shift(
        MPI_Comm comm,
        int direction,
        int displ,
        int *source,
        int *dest)
{
    int rc;


    rc = mpics_MPI_Cart_shift(&comm, &direction, &displ, source, dest);

    return(rc);
}


extern void F77_MPI_CART_SHIFT(
        MPI_Fint *comm,
        int *direction,
        int *displ,
        int *source,
        int *dest,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_shift(&c_comm, direction, displ, source, dest);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Cart_sub(
        MPI_Comm *comm,
        int *remain_dims,
        MPI_Comm *comm_new)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cart_sub(*comm, remain_dims, comm_new);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cart_sub(
        MPI_Comm comm,
        int *remain_dims,
        MPI_Comm *comm_new)
{
    int rc;


    rc = mpics_MPI_Cart_sub(&comm, remain_dims, comm_new);

    return(rc);
}


extern void F77_MPI_CART_SUB(
        MPI_Fint *comm,
        int *remain_dims,
        MPI_Fint *comm_new,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Comm c_comm_new;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cart_sub(&c_comm, remain_dims, &c_comm_new);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *comm_new = MPI_Comm_c2f(c_comm_new);
    }
    return;
}





static int mpics_MPI_Cartdim_get(
        MPI_Comm *comm,
        int *ndims)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Cartdim_get(*comm, ndims);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Cartdim_get(
        MPI_Comm comm,
        int *ndims)
{
    int rc;


    rc = mpics_MPI_Cartdim_get(&comm, ndims);

    return(rc);
}


extern void F77_MPI_CARTDIM_GET(
        MPI_Fint *comm,
        int *ndims,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Cartdim_get(&c_comm, ndims);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Comm_create(
        MPI_Comm *comm,
        MPI_Group *group,
        MPI_Comm *comm_out)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_create(*comm, *group, comm_out);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_create(
        MPI_Comm comm,
        MPI_Group group,
        MPI_Comm *comm_out)
{
    int rc;


    rc = mpics_MPI_Comm_create(&comm, &group, comm_out);

    return(rc);
}


extern void F77_MPI_COMM_CREATE(
        MPI_Fint *comm,
        MPI_Fint *group,
        MPI_Fint *comm_out,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Group c_group;
    MPI_Comm c_comm_out;

    c_comm = MPI_Comm_f2c(*comm);
    c_group = MPI_Group_f2c(*group);

    rc = mpics_MPI_Comm_create(&c_comm, &c_group, &c_comm_out);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *comm_out = MPI_Comm_c2f(c_comm_out);
    }
    return;
}





static int mpics_MPI_Comm_dup(
        MPI_Comm *comm,
        MPI_Comm *comm_out)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_dup(*comm, comm_out);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_dup(
        MPI_Comm comm,
        MPI_Comm *comm_out)
{
    int rc;


    rc = mpics_MPI_Comm_dup(&comm, comm_out);

    return(rc);
}


extern void F77_MPI_COMM_DUP(
        MPI_Fint *comm,
        MPI_Fint *comm_out,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Comm c_comm_out;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_dup(&c_comm, &c_comm_out);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *comm_out = MPI_Comm_c2f(c_comm_out);
    }
    return;
}





static int mpics_MPI_Comm_group(
        MPI_Comm *comm,
        MPI_Group *group)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_group(*comm, group);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_group(
        MPI_Comm comm,
        MPI_Group *group)
{
    int rc;


    rc = mpics_MPI_Comm_group(&comm, group);

    return(rc);
}


extern void F77_MPI_COMM_GROUP(
        MPI_Fint *comm,
        MPI_Fint *group,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Group c_group;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_group(&c_comm, &c_group);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *group = MPI_Group_c2f(c_group);
    }
    return;
}





static int mpics_MPI_Comm_rank(
        MPI_Comm *comm,
        int *rank)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_rank(*comm, rank);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_rank(
        MPI_Comm comm,
        int *rank)
{
    int rc;


    rc = mpics_MPI_Comm_rank(&comm, rank);

    return(rc);
}


extern void F77_MPI_COMM_RANK(
        MPI_Fint *comm,
        int *rank,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_rank(&c_comm, rank);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Comm_remote_group(
        MPI_Comm *comm,
        MPI_Group *group)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_remote_group(*comm, group);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_remote_group(
        MPI_Comm comm,
        MPI_Group *group)
{
    int rc;


    rc = mpics_MPI_Comm_remote_group(&comm, group);

    return(rc);
}


extern void F77_MPI_COMM_REMOTE_GROUP(
        MPI_Fint *comm,
        MPI_Fint *group,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Group c_group;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_remote_group(&c_comm, &c_group);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *group = MPI_Group_c2f(c_group);
    }
    return;
}





static int mpics_MPI_Comm_remote_size(
        MPI_Comm *comm,
        int *size)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_remote_size(*comm, size);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_remote_size(
        MPI_Comm comm,
        int *size)
{
    int rc;


    rc = mpics_MPI_Comm_remote_size(&comm, size);

    return(rc);
}


extern void F77_MPI_COMM_REMOTE_SIZE(
        MPI_Fint *comm,
        int *size,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_remote_size(&c_comm, size);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Comm_size(
        MPI_Comm *comm,
        int *size)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_size(*comm, size);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_size(
        MPI_Comm comm,
        int *size)
{
    int rc;


    rc = mpics_MPI_Comm_size(&comm, size);

    return(rc);
}


extern void F77_MPI_COMM_SIZE(
        MPI_Fint *comm,
        int *size,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_size(&c_comm, size);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Comm_split(
        MPI_Comm *comm,
        int *color,
        int *key,
        MPI_Comm *comm_out)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_split(*comm, *color, *key, comm_out);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_split(
        MPI_Comm comm,
        int color,
        int key,
        MPI_Comm *comm_out)
{
    int rc;


    rc = mpics_MPI_Comm_split(&comm, &color, &key, comm_out);

    return(rc);
}


extern void F77_MPI_COMM_SPLIT(
        MPI_Fint *comm,
        int *color,
        int *key,
        MPI_Fint *comm_out,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Comm c_comm_out;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_split(&c_comm, color, key, &c_comm_out);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *comm_out = MPI_Comm_c2f(c_comm_out);
    }
    return;
}





static int mpics_MPI_Comm_test_inter(
        MPI_Comm *comm,
        int *flag)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Comm_test_inter(*comm, flag);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Comm_test_inter(
        MPI_Comm comm,
        int *flag)
{
    int rc;


    rc = mpics_MPI_Comm_test_inter(&comm, flag);

    return(rc);
}


extern void F77_MPI_COMM_TEST_INTER(
        MPI_Fint *comm,
        int *flag,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Comm_test_inter(&c_comm, flag);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Errhandler_get(
        MPI_Comm *comm,
        MPI_Errhandler *errhandler)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Errhandler_get(*comm, errhandler);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Errhandler_get(
        MPI_Comm comm,
        MPI_Errhandler *errhandler)
{
    int rc;


    rc = mpics_MPI_Errhandler_get(&comm, errhandler);

    return(rc);
}


extern void F77_MPI_ERRHANDLER_GET(
        MPI_Fint *comm,
        MPI_Errhandler *errhandler,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Errhandler_get(&c_comm, errhandler);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Errhandler_set(
        MPI_Comm *comm,
        MPI_Errhandler *errhandler)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Errhandler_set(*comm, *errhandler);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Errhandler_set(
        MPI_Comm comm,
        MPI_Errhandler errhandler)
{
    int rc;


    rc = mpics_MPI_Errhandler_set(&comm, &errhandler);

    return(rc);
}


extern void F77_MPI_ERRHANDLER_SET(
        MPI_Fint *comm,
        MPI_Errhandler *errhandler,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Errhandler_set(&c_comm, errhandler);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_File_open(
        MPI_Comm *comm,
        char *filename,
        int *amode,
        MPI_Info *info,
        MPI_File *fh)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_File_open(*comm, filename, *amode, *info, fh);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_File_open(
        MPI_Comm comm,
        char *filename,
        int amode,
        MPI_Info info,
        MPI_File *fh)
{
    int rc;


    rc = mpics_MPI_File_open(&comm, filename, &amode, &info, fh);

    return(rc);
}


extern void F77_MPI_FILE_OPEN(
        MPI_Fint *comm,
        char *filename,
        int *amode,
        MPI_Fint *info,
        MPI_Fint *fh,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Info c_info;
    MPI_File c_fh;

    c_comm = MPI_Comm_f2c(*comm);
    c_info = MPI_Info_f2c(*info);

    rc = mpics_MPI_File_open(&c_comm, filename, amode, &c_info, &c_fh);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *fh = MPI_File_c2f(c_fh);
    }
    return;
}





static int mpics_MPI_Gather(
        void *sendbuf,
        int *sendcnt,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcount,
        MPI_Datatype *recvtype,
        int *root,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Gather(sendbuf, *sendcnt, *sendtype, recvbuf, *recvcount, *recvtype, *root, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Gather(
        void *sendbuf,
        int sendcnt,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcount,
        MPI_Datatype recvtype,
        int root,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Gather(sendbuf, &sendcnt, &sendtype, recvbuf, &recvcount, &recvtype, &root, &comm);

    return(rc);
}


extern void F77_MPI_GATHER(
        void *sendbuf,
        int *sendcnt,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcount,
        MPI_Fint *recvtype,
        int *root,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Gather(sendbuf, sendcnt, &c_sendtype, recvbuf, recvcount, &c_recvtype, root, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Gatherv(
        void *sendbuf,
        int *sendcnt,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcnts,
        int *displs,
        MPI_Datatype *recvtype,
        int *root,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Gatherv(sendbuf, *sendcnt, *sendtype, recvbuf, recvcnts, displs, *recvtype, *root, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Gatherv(
        void *sendbuf,
        int sendcnt,
        MPI_Datatype sendtype,
        void *recvbuf,
        int *recvcnts,
        int *displs,
        MPI_Datatype recvtype,
        int root,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Gatherv(sendbuf, &sendcnt, &sendtype, recvbuf, recvcnts, displs, &recvtype, &root, &comm);

    return(rc);
}


extern void F77_MPI_GATHERV(
        void *sendbuf,
        int *sendcnt,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcnts,
        int *displs,
        MPI_Fint *recvtype,
        int *root,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Gatherv(sendbuf, sendcnt, &c_sendtype, recvbuf, recvcnts, displs, &c_recvtype, root, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Graph_create(
        MPI_Comm *comm,
        int *nnodes,
        int *index,
        int *edges,
        int *reorder,
        MPI_Comm *comm_graph)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Graph_create(*comm, *nnodes, index, edges, *reorder, comm_graph);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Graph_create(
        MPI_Comm comm,
        int nnodes,
        int *index,
        int *edges,
        int reorder,
        MPI_Comm *comm_graph)
{
    int rc;


    rc = mpics_MPI_Graph_create(&comm, &nnodes, index, edges, &reorder, comm_graph);

    return(rc);
}


extern void F77_MPI_GRAPH_CREATE(
        MPI_Fint *comm,
        int *nnodes,
        int *index,
        int *edges,
        int *reorder,
        MPI_Fint *comm_graph,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;
    MPI_Comm c_comm_graph;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Graph_create(&c_comm, nnodes, index, edges, reorder, &c_comm_graph);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *comm_graph = MPI_Comm_c2f(c_comm_graph);
    }
    return;
}





static int mpics_MPI_Graph_get(
        MPI_Comm *comm,
        int *maxindex,
        int *maxedges,
        int *index,
        int *edges)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Graph_get(*comm, *maxindex, *maxedges, index, edges);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Graph_get(
        MPI_Comm comm,
        int maxindex,
        int maxedges,
        int *index,
        int *edges)
{
    int rc;


    rc = mpics_MPI_Graph_get(&comm, &maxindex, &maxedges, index, edges);

    return(rc);
}


extern void F77_MPI_GRAPH_GET(
        MPI_Fint *comm,
        int *maxindex,
        int *maxedges,
        int *index,
        int *edges,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Graph_get(&c_comm, maxindex, maxedges, index, edges);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Graph_map(
        MPI_Comm *comm,
        int *nnodes,
        int *index,
        int *edges,
        int *newrank)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Graph_map(*comm, *nnodes, index, edges, newrank);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Graph_map(
        MPI_Comm comm,
        int nnodes,
        int *index,
        int *edges,
        int *newrank)
{
    int rc;


    rc = mpics_MPI_Graph_map(&comm, &nnodes, index, edges, newrank);

    return(rc);
}


extern void F77_MPI_GRAPH_MAP(
        MPI_Fint *comm,
        int *nnodes,
        int *index,
        int *edges,
        int *newrank,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Graph_map(&c_comm, nnodes, index, edges, newrank);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Graph_neighbors(
        MPI_Comm *comm,
        int *rank,
        int *maxneighbors,
        int *neighbors)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Graph_neighbors(*comm, *rank, *maxneighbors, neighbors);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Graph_neighbors(
        MPI_Comm comm,
        int rank,
        int maxneighbors,
        int *neighbors)
{
    int rc;


    rc = mpics_MPI_Graph_neighbors(&comm, &rank, &maxneighbors, neighbors);

    return(rc);
}


extern void F77_MPI_GRAPH_NEIGHBORS(
        MPI_Fint *comm,
        int *rank,
        int *maxneighbors,
        int *neighbors,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Graph_neighbors(&c_comm, rank, maxneighbors, neighbors);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Graph_neighbors_count(
        MPI_Comm *comm,
        int *rank,
        int *nneighbors)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Graph_neighbors_count(*comm, *rank, nneighbors);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Graph_neighbors_count(
        MPI_Comm comm,
        int rank,
        int *nneighbors)
{
    int rc;


    rc = mpics_MPI_Graph_neighbors_count(&comm, &rank, nneighbors);

    return(rc);
}


extern void F77_MPI_GRAPH_NEIGHBORS_COUNT(
        MPI_Fint *comm,
        int *rank,
        int *nneighbors,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Graph_neighbors_count(&c_comm, rank, nneighbors);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Graphdims_get(
        MPI_Comm *comm,
        int *nnodes,
        int *nedges)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Graphdims_get(*comm, nnodes, nedges);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Graphdims_get(
        MPI_Comm comm,
        int *nnodes,
        int *nedges)
{
    int rc;


    rc = mpics_MPI_Graphdims_get(&comm, nnodes, nedges);

    return(rc);
}


extern void F77_MPI_GRAPHDIMS_GET(
        MPI_Fint *comm,
        int *nnodes,
        int *nedges,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Graphdims_get(&c_comm, nnodes, nedges);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Iprobe(
        int *source,
        int *tag,
        MPI_Comm *comm,
        int *flag,
        MPI_Status *status)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Iprobe(*source, *tag, *comm, flag, status);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Iprobe(
        int source,
        int tag,
        MPI_Comm comm,
        int *flag,
        MPI_Status *status)
{
    int rc;


    rc = mpics_MPI_Iprobe(&source, &tag, &comm, flag, status);

    return(rc);
}


extern void F77_MPI_IPROBE(
        int *source,
        int *tag,
        MPI_Fint *comm,
        int *flag,
        MPI_Status *status,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Iprobe(source, tag, &c_comm, flag, status);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Irecv(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *source,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Irecv(buf, *count, *datatype, *source, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Irecv(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int source,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Irecv(buf, &count, &datatype, &source, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_IRECV(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *source,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Irecv(buf, count, &c_datatype, source, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Irsend(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Irsend(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Irsend(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Irsend(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_IRSEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Irsend(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Isend(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Isend(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Isend(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Isend(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_ISEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Isend(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Issend(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Issend(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Issend(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Issend(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_ISSEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Issend(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Pack(
        void *inbuf,
        int *incount,
        MPI_Datatype *datatype,
        void *outbuf,
        int *count,
        int *position,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Pack(inbuf, *incount, *datatype, outbuf, *count, position, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Pack(
        void *inbuf,
        int incount,
        MPI_Datatype datatype,
        void *outbuf,
        int count,
        int *position,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Pack(inbuf, &incount, &datatype, outbuf, &count, position, &comm);

    return(rc);
}


extern void F77_MPI_PACK(
        void *inbuf,
        int *incount,
        MPI_Fint *datatype,
        void *outbuf,
        int *count,
        int *position,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Pack(inbuf, incount, &c_datatype, outbuf, count, position, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Pack_size(
        int *incount,
        MPI_Datatype *datatype,
        MPI_Comm *comm,
        int *size)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Pack_size(*incount, *datatype, *comm, size);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Pack_size(
        int incount,
        MPI_Datatype datatype,
        MPI_Comm comm,
        int *size)
{
    int rc;


    rc = mpics_MPI_Pack_size(&incount, &datatype, &comm, size);

    return(rc);
}


extern void F77_MPI_PACK_SIZE(
        int *incount,
        MPI_Fint *datatype,
        MPI_Fint *comm,
        int *size,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Pack_size(incount, &c_datatype, &c_comm, size);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Probe(
        int *source,
        int *tag,
        MPI_Comm *comm,
        MPI_Status *status)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Probe(*source, *tag, *comm, status);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Probe(
        int source,
        int tag,
        MPI_Comm comm,
        MPI_Status *status)
{
    int rc;


    rc = mpics_MPI_Probe(&source, &tag, &comm, status);

    return(rc);
}


extern void F77_MPI_PROBE(
        int *source,
        int *tag,
        MPI_Fint *comm,
        MPI_Status *status,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Probe(source, tag, &c_comm, status);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Recv(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *source,
        int *tag,
        MPI_Comm *comm,
        MPI_Status *status)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Recv(buf, *count, *datatype, *source, *tag, *comm, status);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Recv(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int source,
        int tag,
        MPI_Comm comm,
        MPI_Status *status)
{
    int rc;


    rc = mpics_MPI_Recv(buf, &count, &datatype, &source, &tag, &comm, status);

    return(rc);
}


extern void F77_MPI_RECV(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *source,
        int *tag,
        MPI_Fint *comm,
        MPI_Status *status,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Recv(buf, count, &c_datatype, source, tag, &c_comm, status);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Recv_init(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *source,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Recv_init(buf, *count, *datatype, *source, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Recv_init(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int source,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Recv_init(buf, &count, &datatype, &source, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_RECV_INIT(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *source,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Recv_init(buf, count, &c_datatype, source, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Reduce(
        void *sendbuf,
        void *recvbuf,
        int *count,
        MPI_Datatype *datatype,
        MPI_Op *op,
        int *root,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Reduce(sendbuf, recvbuf, *count, *datatype, *op, *root, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Reduce(
        void *sendbuf,
        void *recvbuf,
        int count,
        MPI_Datatype datatype,
        MPI_Op op,
        int root,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Reduce(sendbuf, recvbuf, &count, &datatype, &op, &root, &comm);

    return(rc);
}


extern void F77_MPI_REDUCE(
        void *sendbuf,
        void *recvbuf,
        int *count,
        MPI_Fint *datatype,
        MPI_Fint *op,
        int *root,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Reduce(sendbuf, recvbuf, count, &c_datatype, &c_op, root, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Rsend(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Rsend(buf, *count, *datatype, *dest, *tag, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Rsend(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Rsend(buf, &count, &datatype, &dest, &tag, &comm);

    return(rc);
}


extern void F77_MPI_RSEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Rsend(buf, count, &c_datatype, dest, tag, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Rsend_init(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Rsend_init(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Rsend_init(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Rsend_init(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_RSEND_INIT(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Rsend_init(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Scan(
        void *sendbuf,
        void *recvbuf,
        int *count,
        MPI_Datatype *datatype,
        MPI_Op *op,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Scan(sendbuf, recvbuf, *count, *datatype, *op, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Scan(
        void *sendbuf,
        void *recvbuf,
        int count,
        MPI_Datatype datatype,
        MPI_Op op,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Scan(sendbuf, recvbuf, &count, &datatype, &op, &comm);

    return(rc);
}


extern void F77_MPI_SCAN(
        void *sendbuf,
        void *recvbuf,
        int *count,
        MPI_Fint *datatype,
        MPI_Fint *op,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Scan(sendbuf, recvbuf, count, &c_datatype, &c_op, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Scatter(
        void *sendbuf,
        int *sendcnt,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcnt,
        MPI_Datatype *recvtype,
        int *root,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Scatter(sendbuf, *sendcnt, *sendtype, recvbuf, *recvcnt, *recvtype, *root, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Scatter(
        void *sendbuf,
        int sendcnt,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcnt,
        MPI_Datatype recvtype,
        int root,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Scatter(sendbuf, &sendcnt, &sendtype, recvbuf, &recvcnt, &recvtype, &root, &comm);

    return(rc);
}


extern void F77_MPI_SCATTER(
        void *sendbuf,
        int *sendcnt,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcnt,
        MPI_Fint *recvtype,
        int *root,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Scatter(sendbuf, sendcnt, &c_sendtype, recvbuf, recvcnt, &c_recvtype, root, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Scatterv(
        void *sendbuf,
        int *sendcnts,
        int *displs,
        MPI_Datatype *sendtype,
        void *recvbuf,
        int *recvcnt,
        MPI_Datatype *recvtype,
        int *root,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Scatterv(sendbuf, sendcnts, displs, *sendtype, recvbuf, *recvcnt, *recvtype, *root, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Scatterv(
        void *sendbuf,
        int *sendcnts,
        int *displs,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcnt,
        MPI_Datatype recvtype,
        int root,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Scatterv(sendbuf, sendcnts, displs, &sendtype, recvbuf, &recvcnt, &recvtype, &root, &comm);

    return(rc);
}


extern void F77_MPI_SCATTERV(
        void *sendbuf,
        int *sendcnts,
        int *displs,
        MPI_Fint *sendtype,
        void *recvbuf,
        int *recvcnt,
        MPI_Fint *recvtype,
        int *root,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Scatterv(sendbuf, sendcnts, displs, &c_sendtype, recvbuf, recvcnt, &c_recvtype, root, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Send(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Send(buf, *count, *datatype, *dest, *tag, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Send(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Send(buf, &count, &datatype, &dest, &tag, &comm);

    return(rc);
}


extern void F77_MPI_SEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Send(buf, count, &c_datatype, dest, tag, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Send_init(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Send_init(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Send_init(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Send_init(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_SEND_INIT(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Send_init(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Sendrecv(
        void *sendbuf,
        int *sendcount,
        MPI_Datatype *sendtype,
        int *dest,
        int *sendtag,
        void *recvbuf,
        int *recvcount,
        MPI_Datatype *recvtype,
        int *source,
        int *recvtag,
        MPI_Comm *comm,
        MPI_Status *status)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Sendrecv(sendbuf, *sendcount, *sendtype, *dest, *sendtag, recvbuf, *recvcount, *recvtype, *source, *recvtag, *comm, status);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Sendrecv(
        void *sendbuf,
        int sendcount,
        MPI_Datatype sendtype,
        int dest,
        int sendtag,
        void *recvbuf,
        int recvcount,
        MPI_Datatype recvtype,
        int source,
        int recvtag,
        MPI_Comm comm,
        MPI_Status *status)
{
    int rc;


    rc = mpics_MPI_Sendrecv(sendbuf, &sendcount, &sendtype, &dest, &sendtag, recvbuf, &recvcount, &recvtype, &source, &recvtag, &comm, status);

    return(rc);
}


extern void F77_MPI_SENDRECV(
        void *sendbuf,
        int *sendcount,
        MPI_Fint *sendtype,
        int *dest,
        int *sendtag,
        void *recvbuf,
        int *recvcount,
        MPI_Fint *recvtype,
        int *source,
        int *recvtag,
        MPI_Fint *comm,
        MPI_Status *status,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Sendrecv(sendbuf, sendcount, &c_sendtype, dest, sendtag, recvbuf, recvcount, &c_recvtype, source, recvtag, &c_comm, status);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Sendrecv_replace(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *sendtag,
        int *source,
        int *recvtag,
        MPI_Comm *comm,
        MPI_Status *status)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Sendrecv_replace(buf, *count, *datatype, *dest, *sendtag, *source, *recvtag, *comm, status);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Sendrecv_replace(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int sendtag,
        int source,
        int recvtag,
        MPI_Comm comm,
        MPI_Status *status)
{
    int rc;


    rc = mpics_MPI_Sendrecv_replace(buf, &count, &datatype, &dest, &sendtag, &source, &recvtag, &comm, status);

    return(rc);
}


extern void F77_MPI_SENDRECV_REPLACE(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *sendtag,
        int *source,
        int *recvtag,
        MPI_Fint *comm,
        MPI_Status *status,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Sendrecv_replace(buf, count, &c_datatype, dest, sendtag, source, recvtag, &c_comm, status);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Ssend(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Ssend(buf, *count, *datatype, *dest, *tag, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Ssend(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Ssend(buf, &count, &datatype, &dest, &tag, &comm);

    return(rc);
}


extern void F77_MPI_SSEND(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Ssend(buf, count, &c_datatype, dest, tag, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Ssend_init(
        void *buf,
        int *count,
        MPI_Datatype *datatype,
        int *dest,
        int *tag,
        MPI_Comm *comm,
        MPI_Request *request)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Ssend_init(buf, *count, *datatype, *dest, *tag, *comm, request);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Ssend_init(
        void *buf,
        int count,
        MPI_Datatype datatype,
        int dest,
        int tag,
        MPI_Comm comm,
        MPI_Request *request)
{
    int rc;


    rc = mpics_MPI_Ssend_init(buf, &count, &datatype, &dest, &tag, &comm, request);

    return(rc);
}


extern void F77_MPI_SSEND_INIT(
        void *buf,
        int *count,
        MPI_Fint *datatype,
        int *dest,
        int *tag,
        MPI_Fint *comm,
        MPI_Fint *request,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Ssend_init(buf, count, &c_datatype, dest, tag, &c_comm, &c_request);

    *ierr = (MPI_Fint)rc;
    if (rc == MPI_SUCCESS) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}





static int mpics_MPI_Topo_test(
        MPI_Comm *comm,
        int *top_type)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Topo_test(*comm, top_type);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Topo_test(
        MPI_Comm comm,
        int *top_type)
{
    int rc;


    rc = mpics_MPI_Topo_test(&comm, top_type);

    return(rc);
}


extern void F77_MPI_TOPO_TEST(
        MPI_Fint *comm,
        int *top_type,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Topo_test(&c_comm, top_type);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Unpack(
        void *inbuf,
        int *insize,
        int *position,
        void *outbuf,
        int *count,
        MPI_Datatype *datatype,
        MPI_Comm *comm)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Unpack(inbuf, *insize, position, outbuf, *count, *datatype, *comm);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Unpack(
        void *inbuf,
        int insize,
        int *position,
        void *outbuf,
        int count,
        MPI_Datatype datatype,
        MPI_Comm comm)
{
    int rc;


    rc = mpics_MPI_Unpack(inbuf, &insize, position, outbuf, &count, &datatype, &comm);

    return(rc);
}


extern void F77_MPI_UNPACK(
        void *inbuf,
        int *insize,
        int *position,
        void *outbuf,
        int *count,
        MPI_Fint *datatype,
        MPI_Fint *comm,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Unpack(inbuf, insize, position, outbuf, count, &c_datatype, &c_comm);

    *ierr = (MPI_Fint)rc;
    return;
}





static int mpics_MPI_Win_create(
        void *base,
        MPI_Aint *size,
        int *disp_unit,
        MPI_Info *info,
        MPI_Comm *comm,
        MPI_Win *win)
{
    int rc, compare_result, enabled_state;
    MPI_Comm orig_comm;

    if (mpics_data.enabled) {
        SUBSTITUTE_COMM;
    }

    enabled_state = mpics_data.enabled;
    mpics_data.enabled = 0;

    rc = PMPI_Win_create(base, *size, *disp_unit, *info, *comm, win);

    mpics_data.enabled = enabled_state;
    if (mpics_data.enabled) {
        RESTORE_COMM;
    }

    return(rc);
}



extern int MPI_Win_create(
        void *base,
        MPI_Aint size,
        int disp_unit,
        MPI_Info info,
        MPI_Comm comm,
        MPI_Win *win)
{
    int rc;


    rc = mpics_MPI_Win_create(base, &size, &disp_unit, &info, &comm, win);

    return(rc);
}


extern void F77_MPI_WIN_CREATE(
        void *base,
        MPI_Aint *size,
        int *disp_unit,
        MPI_Fint *info,
        MPI_Fint *comm,
        MPI_Win *win,
        MPI_Fint *ierr)
{
    int rc;

    MPI_Info c_info;
    MPI_Comm c_comm;

    c_info = MPI_Info_f2c(*info);
    c_comm = MPI_Comm_f2c(*comm);

    rc = mpics_MPI_Win_create(base, size, disp_unit, &c_info, &c_comm, win);

    *ierr = (MPI_Fint)rc;
    return;
}
