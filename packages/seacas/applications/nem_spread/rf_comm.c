/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/* System Include files */

#include <stdlib.h>
#include <stdio.h>

/* User include files */

#include "rf_comm.h"
#include "rf_salsa.h"

#include "rf_message.h"
#include "rf_allo.h"
#include "rf_mp_const.h"

/************* ROUTINES IN THIS FILE ******************************************
*
*       Name                 Type                 Called_By
*    -----------           ---------          ----------------
*    nwrite_big              int                 rf_load_lb_info.c
*    nread_big               int                 rf_load_lb_info.c
*    brdcst_maxlen ()        void                el_exoII_io.c
*    brdcst_nosync ()        static void         brdcst_maxlen ()
*    brdcst ()               void                Many places
*    psync ()                void                Many places
*    print_sync_start ()     void                Many places
*    print_sync_end ()       void                Many places
*    gsum_int ()             int                 Many places
*    gsum_double ()          double              Many places
*    gmax_double ()          double              Many places
*    gmin_double ()          double              Many places
*    gmax_int ()             int                 Many places
*    gavg_double ()          double              Many places
*    gmin_int ()             int                 Many places
*
******************************************************************************/

/* static function declarations */

static void brdcst_nosync(int node, int nprocs, char *data, int length,
                          int sender);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int nwrite_big(char *buffer, int nbytes, int dest, int type, int *flag)

/*
 *     nwrite_big:
 *         This duplicates the ncube nwrite call.  However, it may be used
 *   for large messages that may be larger than the communications buffer
 *   size, without causing overflow.  Basically, it breaks the message up
 *   into parts, equal to one third the size of the message buffer.  It
 *   then waits for acknowledgement of each part of the message before
 *   sending the next part.
 *         The algorithm probably needs more work.  It currently runs at
 *   1.5 megabytes per second on the nCUBE 2, out of a possible speed of
 *   2.0 megabytes per second.
 */

{
  if (nbytes > (BRCST_COMM_BUFF_SIZE / 2)) {

    int length = BRCST_COMM_BUFF_SIZE / 3, ack, loc_counter = 0,
      type_ack = MT_BIG_WRITE_ACK;

    do {

      if (md_write((buffer+loc_counter), length, dest, type, flag) != 0) {
        (void) fprintf(stderr, "nwrite_big error on node %d", Proc);
        exit (-1);
      }
      loc_counter += length;

      if (md_read((char *) &ack, sizeof(int), &dest, &type_ack, flag) !=
          sizeof(int) ) {
        (void) fprintf(stderr, "nwrite_big ack error on node %d", Proc);
        exit (-1);
      }

    } while ((loc_counter + length) < nbytes);

    length = nbytes - loc_counter;
    if (md_write(&buffer[loc_counter], length, dest, type, flag) != 0) {
      (void) fprintf(stderr, "nwrite_big error on node %d", Proc);
      exit(1);
    }

    return 0;
  }

  return (md_write(buffer, nbytes, dest, type, flag));

} /* nwrite_big */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int nread_big(char *buffer, int nbytes, int *source, int *type, int *flag)

/*
 *  nread_big:
 *     See the writeup of nwrite_big
 *
 *     NOTE:
 *       Currently, nread_big doesn't handle arbitrary sources or types.
 *       These must be explicitly supplied in the calling arguments.
 */

{
  if (nbytes > (BRCST_COMM_BUFF_SIZE / 2)) {

    int length = BRCST_COMM_BUFF_SIZE / 3, ack = 0, loc_counter = 0,
      type_ack = MT_BIG_WRITE_ACK;

    if (*source < 0) {
      (void) fprintf(stderr, "nread_big error");
      exit (-1);
    }

    if (*type < 0) {
      (void) fprintf(stderr, "nread_big error");
      exit (-1);
    }

    do {

      if (md_write((char *) &ack, sizeof (int), *source, type_ack, flag) != 0) {
        (void) fprintf(stderr, "nread_big ack error on node %d", Proc);
        exit (-1);
      }

      if (md_read(&buffer[loc_counter], length, source, type, flag) != length) {
        (void) fprintf(stderr, "nread_big error on node %d", Proc);
        exit (-1);
      }

      loc_counter += length;

    } while ((loc_counter + length) < nbytes);

    length = nbytes - loc_counter;
    if (md_read(&buffer[loc_counter], length, source, type, flag) != length) {
      (void) fprintf(stderr, "nread_big error on node %d", Proc);
      exit (-1);
    }

    return nbytes;
  }

  return (md_read(buffer, nbytes, source, type, flag));

} /* nread_big */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void brdcst_maxlen(
                   const int node,   /* This is the current processor number */
                   const int nprocs, /* The dimension of the logical hypercube
                                      */
                   char     *data,   /* The address of the data to be broadcast
                                      */
                   int       length, /* Length in bytes of the data to be
                                        broadcast */
                   int       sender  /* The broadcasting processor.           */
)

/*
 *       This function broadcasts a vector stored on processor "sender" to all
 *    processors on the machine.  A maximum message size is enforced in this
 *    call; i.e., if length is greater than the max message size, then more than
 *    one broadcast is carried out.
 *       Note that each processor must know the value of length before the
 *    broadcast takes place.
 */

{

  int num_bcasts; /* number of normal broadcasts */
  int istart_pos; /* index into piecs of data */
  int nleft_over; /* data left at end */
  int j;          /* loop counter */

  num_bcasts = length / MAX_MESG_SIZE;

  /* Send out full sized messages. */

  istart_pos = 0;
  for (j = 0; j < num_bcasts; j++) {
    brdcst_nosync(node, nprocs, &data[istart_pos], MAX_MESG_SIZE, sender);
    istart_pos += MAX_MESG_SIZE;

    /* Synchronize if needed to ensure buffers clear. */

    if ((j % MAX_MESG_IN_BUF) == (MAX_MESG_IN_BUF - 1))
      psync(node, nprocs);
  }

  /* Now handle left overs. */

  nleft_over = length - num_bcasts*MAX_MESG_SIZE;
  if (nleft_over)
    brdcst_nosync(node, nprocs, &data[istart_pos], nleft_over, sender);

  /*
   *   psync(node, nprocs);
   */
}

/******************************************************************************
* Note that the hard sync in the j loop could be replaced with a
* "soft" sync.  However, the payoff is not that great compared with doing
* no sync at all.
*   soft sync code:
*      if (Dim > 0) {
*       np = (1 << Dim) - 1;
*       type = BRDCST + 101 + Dim ;
*       if (node == np) {
*         if ( md_write((char *) &p, sizeof(int),   0,  type, &st) != 0) {
*           (void) fprintf(stderr,"brdcst: ERROR on node %d\n", node);
*           (void) fprintf(stderr,"md_write failed, message type %d/n",type);
*           exit (-1);
*         }
*       } else if (node == 0) {
*         if ( md_read((char *) &p, sizeof(int), &np, &type, &st) != sizeof(int)){
*           (void) fprintf(stderr, "brdcst: ERROR on node %d\n", node);
*           (void) fprintf(stderr," md_read failed, message type %d/n",type);
*           exit (-1);
*         }
*       }
*      }
******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void brdcst_nosync(
     int   node,     /* node is the current processor number */
     int   nprocs,   /* number of processors in the current run */
     char *data,      /* *buf is the address of the data to be broadcast */
     int   length,   /* length is the number of bytes to be broadcast */
     int   sender    /* broadcasting processor */
)

/*
 * This function broadcasts a vector storred on processor "sender" to all
 * processors on the machine. This function works on hypercube as well as the
 * 2D mesh architecture.
 */

{
  int        type;          /* type of next message */
  int        left_sender;   /* beginning of left half of processors */
  int        right_sender;  /* beginning of right half of processors */
  int        left_nprocs;   /* length of one sequence of processors */
  int        right_nprocs;  /* length of other sequence of processors */
  int        new_nprocs;    /* recursive number of processors left */
  int        receiver;      /* list element receiving at this step */
  int        cflag;
  static int offset = 50;

  if (nprocs <= 1) return;

  left_sender = 0;
  new_nprocs  = nprocs;
  type = BRDCST + offset;
  offset++;
  if (offset > 99) offset = 50;

  while (new_nprocs > 1) {
    left_nprocs  = new_nprocs / 2;
    right_nprocs = new_nprocs - left_nprocs;
    right_sender = left_sender + left_nprocs;

    if (sender < right_sender) {
      receiver = sender + right_nprocs;
    }
    else {
      receiver = sender - right_nprocs;
      if (receiver < left_sender)
        receiver = left_sender;
    }

    if (sender == node) {
      if (md_write(data, length, receiver, type, &cflag)) {
        (void) fprintf(stderr, "brdcst_nosync: ERROR on node %d\n", node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
    }
    else if (receiver == node) {
      if (md_read(data, length, &sender, &type, &cflag) != length) {
        (void) fprintf(stderr, "brdcst_nosync: ERROR on node %d\n", node);
        (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
        exit(1);
      }
    }

    if (node >= right_sender) {
      left_sender = right_sender;
      new_nprocs  = right_nprocs;
      if (receiver >= right_sender)
        sender = receiver;
    }
    else {
      new_nprocs = left_nprocs;
      if (receiver < right_sender)
        sender = receiver;
    }
  }

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void brdcst(
             int   node,   /* node is the current processor number           */
             int   nprocs, /* length of processor list                       */
             char *data,   /* *buf is the address of the data to be broadcast*/
             int   length, /* length is the number of bytes to be broadcast  */
             int   sender  /* broadcasting processor                         */
             )

/*
 * This function broadcasts a vector storred on processor "sender" to all
 * processors on the machine.  A synchronization is enforced at the end to
 * ensure buffers are cleared.
 */

{
  int       type;             /* type of next message */
  int       left_sender;      /* beginning of left half of processors */
  int       right_sender;     /* beginning of right half of processors */
  int       left_nprocs;      /* length of one sequence of processors */
  int       right_nprocs;     /* length of other sequence of processors */
  int       new_nprocs;       /* recursive number of processors left */
  int       receiver;         /* list element receiving at this step */
  int       cflag;            /* dummy parameter for compatibility */
  static int offset = 0;      /* offset to avoid type conflicts */

  if (nprocs <= 1) return;

  left_sender = 0;
  new_nprocs  = nprocs;
  type = BRDCST + offset;
  offset++;
  if (offset > 49)
    offset = 0;

  while (new_nprocs > 1) {
    left_nprocs  = new_nprocs / 2;
    right_nprocs = new_nprocs - left_nprocs;
    right_sender = left_sender + left_nprocs;

    if (sender < right_sender) {
      receiver = sender + right_nprocs;
    }
    else {
      receiver = sender - right_nprocs;
      if (receiver < left_sender)
        receiver = left_sender;
    }

    if (sender == node) {
      if (md_write(data, length, receiver, type, &cflag)) {
        (void) fprintf(stderr, "brdcst: ERROR on node %d\n", node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
    }
    else if (receiver == node) {
      if (md_read(data, length, &sender, &type, &cflag) != length) {
        (void) fprintf(stderr, "brdcst: ERROR on node %d\n", node);
        (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
        exit(1);
      }
    }

    if (node >= right_sender) {
      left_sender = right_sender;
      new_nprocs  = right_nprocs;
      if (receiver >= right_sender)
        sender = receiver;
    }
    else {
      new_nprocs = left_nprocs;
      if (receiver < right_sender)
        sender = receiver;
    }
  }

  /* Synchronize to ensure message buffers clear. */

  psync(Proc, Num_Proc);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void psync(int node, int nprocs)

{
  int        type;         /* type of next message */
  int        partner;      /* processor I exchange with */
  int        mask;         /* bit pattern identifying partner */
  int        hbit;         /* largest nonzero bit in nprocs */
  int        nprocs_small; /* largest power of 2 <= nprocs */
  int        cflag;        /* dummy argument for compatability */
  static int offset = 0;   /* offset to avoid type conflict */

  type = SYNC + 12*offset;
  offset++;
  if (offset > 8)
    offset = 0;

  /* Find next lower power of 2. */

  for(hbit=0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *)NULL, 0, partner, type, &cflag) != 0) {
      (void) fprintf(stderr, "psync: ERROR on node %d\n", node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", node);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *)NULL, 0, &partner, &type, &cflag) != 0) {
      (void) fprintf(stderr, "psync: ERROR on node %d\n", node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
  }

  /*  Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask=nprocs_small>>1; mask; mask >>= 1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *)NULL, 0, partner, type, &cflag)) {
        (void) fprintf(stderr, "psync: ERROR on node %d\n", node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *)NULL, 0, &partner, &type, &cflag) != 0) {
        (void) fprintf(stderr, "psync: ERROR on node %d\n", node);
        (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
        exit(1);
      }
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *)NULL, 0, &partner, &type, &cflag) != 0) {
      (void) fprintf(stderr, "psync: ERROR on node %d\n", node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_write((char *)NULL, 0, partner, type, &cflag)) {
      (void) fprintf(stderr, "psync: ERROR on node %d\n", node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_sync_start(int proc, int nprocs, int do_print_line)

/* Routine to allow IO between print_sync_start and print_sync_end to be
   printed by each processor entirely before the next processor begins its IO.
   The printing sequence is from proc = 0 to the last processor,
   number_of_procs = nprocs - 1.

   The last argument is a boolean variable.  If true, a line of # is printed to
   indicate the start of a print_sync I/O block.

   NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

   Author: John Shadid (1421, SNL)

*/

{
  int        flag = 1, from, st, type;
  static int offset = 0;

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if ( proc != 0) {
    from = proc -1;
    if ( md_read( (char *) &flag, sizeof(int), &from, &type, &st)
         != sizeof(int) ) {
      (void) fprintf(stderr, "print_sync_start: ERROR on node %d\n", Proc);
      (void) fprintf(stderr, "md_read failed, message type %d\n", type);
      exit (-1);
    }
  }
  else {
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 36; flag++) (void) printf("#");
      (void) printf(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++) (void) printf("#");
      (void) printf("\n");
    }
  }

#ifdef DEBUG_PSYNC
  (void) printf("\t\tSTART OF PRINT_SYNC SECTION, Proc = %4d, Message_Type = "
                "%6d\n", proc, type);
#endif

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_sync_end(int proc, int nprocs, int do_print_line)

/*

  Routine to allow IO between print_sync_start and print_sync_end to be printed
  by each processor entirely before the next processor begins its IO.  The
  printing sequence is from proc = 0 to the last processor, number_of_procs =
  nprocs - 1.

  The last argument is a boolean variable.  If true, a line of # is printed to
  indicate the start of a print_sync I/O block.

  NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

  Author: John Shadid (1421, SNL)

*/

{
  int         st, flag = 1, from, type, to;
  static int  offset = 0;
  extern void sync(int, int);

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

#ifdef DEBUG_PSYNC
  (void) printf("\t\tEND OF PRINT_SYNC SECTION, Proc = %4d, Message_Type = "
                "%6d\n", proc, type);
#endif

  if (proc < nprocs -1)
    to = proc + 1;
  else {
    to = 0;
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 36; flag++) (void) printf("#");
      (void) printf(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) (void) printf("#");
      (void) printf("\n\n");
    }
  }

  if (md_write((char *) &flag, sizeof(int), to, type, &st) != 0 ) {
    (void) fprintf(stderr, "print_sync_end: ERROR on node %d\n", Proc);
    (void) fprintf(stderr, "md_write failed, message type %d\n", type);
    exit (-1);
  }
  if (proc == 0) {
    from = nprocs -1;
    if (md_read((char *) &flag, sizeof(int), &from, &type, &st)
        != sizeof(int) ) {
      (void) fprintf(stderr, "print_sync_end: ERROR on node %d\n", Proc);
      (void) fprintf(stderr, "md_read failed, message type %d/n", type);
      exit (-1);
    }

#ifdef DEBUG_PSYNC
    (void) printf("\t\t\t Proc 0 received message from %5d, type = %5d, flag = "
                  "%d\n", from, type, flag);
#endif
  }

  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from Proc
   * (Num_Proc-1)
   */

  psync (Proc, Num_Proc);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int gsum_int(int val, int node, int nprocs)

{

  int   type;         /* type of next message */
  int   partner;      /* processor I exchange with */
  int   mask;         /* bit pattern identifying partner */
  int   hbit;         /* largest nonzero bit in nprocs */
  int   nprocs_small; /* largest power of 2 <= nprocs */
  int   val2;         /* arriving value to add */
  int   cflag;        /* dummy argument for compatability */
  char *yo_error = "gsum_int: ERROR on node ";

  /*********************** first executable statment *****************/

  type = GSUM_INT;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *) &val2, sizeof(int), &partner, &type, &cflag) !=
        sizeof(int)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
    val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask=nprocs_small>>1; mask; mask>>=1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *) &val2, sizeof(int), &partner, &type, &cflag) !=
          sizeof(int)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      }
      val += val2;
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *) &val, sizeof(int), &partner, &type, &cflag) !=
        sizeof(int)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

  return val;

} /* gsum_int */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifndef SMOS_GSUMD    /* Only compile in if SMOS_GSUMD is NOT defined. */

double gsum_double(double val, int node, int nprocs)

{
  int    type;         /* type of next message */
  int    partner;      /* processor I exchange with */
  int    mask;         /* bit pattern identifying partner */
  int    hbit;         /* largest nonzero bit in nprocs */
  int    nprocs_small; /* largest power of 2 <= nprocs */
  double val2;         /* arriving value to add */
  int    cflag;        /* dummy argument for compatability */
  char  *yo_error = "gsum_double: ERROR on node ";

  type = GSUM_DOUBLE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *) &val2, sizeof(double), &partner, &type, &cflag) !=
        sizeof(double)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
    val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *) &val2, sizeof(double), &partner, &type, &cflag) !=
          sizeof(double)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      val += val2;
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *) &val, sizeof(double), &partner, &type, &cflag) !=
        sizeof(double)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs ) {
    if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

  return val;

} /* gsum_double */

#endif        /* #ifndef SMOS_GSUMD */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double gmax_double(double val, int node, int nprocs)

{

  int    type;         /* type of next message */
  int    partner;      /* processor I exchange with */
  int    mask;         /* bit pattern identifying partner */
  int    hbit;         /* largest nonzero bit in nprocs */
  int    nprocs_small; /* largest power of 2 <= nprocs */
  double val2;         /* arriving value to add */
  int    cflag;        /* dummy argument for compatability */
  char  *yo_error = "gmax_double: ERROR on node ";
  type = GMAX_DOUBLE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *) &val2, sizeof(double), &partner, &type, &cflag) !=
        sizeof(double)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
    if (val2 > val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *) &val2, sizeof(double), &partner, &type, &cflag) !=
          sizeof(double)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
        exit(1);
      }
      if (val2 > val)
        val = val2;
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *) &val, sizeof(double), &partner, &type, &cflag) !=
        sizeof(double)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

  return val;

} /* gmax_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double gmin_double(double val, int node, int nprocs)

{
  int    type;         /* type of next message */
  int    partner;      /* processor I exchange with */
  int    mask;         /* bit pattern identifying partner */
  int    hbit;         /* largest nonzero bit in nprocs */
  int    nprocs_small; /* largest power of 2 <= nprocs */
  double val2;         /* arriving value to add */
  int    cflag;        /* dummy argument for compatability */
  char  *yo_error = "gmin_double: ERROR on node ";

  type = GMAX_DOUBLE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *) &val2, sizeof(double), &partner, &type, &cflag) !=
        sizeof(double)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
    if (val2 < val)
      val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *) &val2, sizeof(double), &partner, &type, &cflag) !=
          sizeof(double)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (val2 < val)
        val = val2;
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *) &val, sizeof(double), &partner, &type, &cflag) !=
        sizeof(double)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_write((char *) &val, sizeof(double), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

  return val;

} /* gmin_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int gmax_int(int val, int node, int nprocs)

{

  int   type;         /* type of next message */
  int   partner;      /* processor I exchange with */
  int   mask;         /* bit pattern identifying partner */
  int   hbit;         /* largest nonzero bit in nprocs */
  int   nprocs_small; /* largest power of 2 <= nprocs */
  int   val2;         /* arriving value to add */
  int   cflag;        /* dummy argument for compatability */
  char *yo_error = "gmax_int: ERROR on node ";

  type = GMAX_DOUBLE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *) &val2, sizeof(int), &partner, &type, &cflag) !=
        sizeof(int)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
    if (val2 > val)
      val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *) &val2, sizeof(int), &partner, &type, &cflag) !=
          sizeof(int)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
        exit(1);
      }
      if (val2 > val)
        val = val2;
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *) &val, sizeof(int), &partner, &type, &cflag) !=
        sizeof(int)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

  return val;

} /* gmax_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double gavg_double(double val, int node, int nprocs)

{
   return(gsum_double(val, node, nprocs)/nprocs);

} /* gavg_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int gmin_int(int val, int node, int nprocs)

{
  int   type;         /* type of next message */
  int   partner;      /* processor I exchange with */
  int   mask;         /* bit pattern identifying partner */
  int   hbit;         /* largest nonzero bit in nprocs */
  int   nprocs_small; /* largest power of 2 <= nprocs */
  int   val2;         /* arriving value to add */
  int   cflag;        /* dummy argument for compatability */
  char *yo_error = "gmin_int: ERROR on node ";

  type = GMAX_DOUBLE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_read((char *) &val2, sizeof(int), &partner, &type, &cflag) !=
        sizeof(int)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
    if (val2 < val)
      val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      type++;
      partner = node ^ mask;
      if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (md_read((char *) &val2, sizeof(int), &partner, &type, &cflag) !=
          sizeof(int)) {
        (void) fprintf(stderr, "%s%d\n", yo_error, node);
        (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
        exit(1);
      }
      if (val2 < val)
        val = val2;
    }
  }
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (md_read((char *) &val, sizeof(int), &partner, &type, &cflag) !=
        sizeof(int)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_read failed, message type = %d\n", type);
      exit(1);
    }
  }
  else if (node+nprocs_small < nprocs) {
    if (md_write((char *) &val, sizeof(int), partner, type, &cflag)) {
      (void) fprintf(stderr, "%s%d\n", yo_error, node);
      (void) fprintf(stderr, "md_write failed, message type = %d\n", type);
      exit(1);
    }
  }

  return val;

} /* gmin_int */

/******************************************************************************/
/*                       END of rf_comm.c                                     */
/******************************************************************************/
