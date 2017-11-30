/*
* Copyright(c) 2005-2017 National Technology &Engineering Solutions
* of Sandia, LLC(NTESS).Under the terms of Contract DE - NA0003525 with
* NTESS, the U.S.Government retains certain rights in this software.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and / or other materials provided
* with the                                                 distribution.
*
* * Neither the name of NTESS nor the names of its
* contributors may be used to endorse or promote products derived
* from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <exodusII.h>
#include <stdio.h>

int main()
{
  float version = 0.0;

  ex_opts(EX_VERBOSE | EX_ABORT);
  int CPU_word_size = 0; /* sizeof(float) */
  int IO_word_size  = 4; /* (4 bytes) */

  /* ======================================== */
  /* Create an empty exodus file             */
  /* ====================================== */
  int exoid = ex_create("test.exo", EX_CLOBBER, &CPU_word_size, &IO_word_size);
  ex_close(exoid);

  /* ======================================== */
  /* Now try to open and read the empty file */
  /* ====================================== */

  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  printf("test.exo exoid = %d\n", exoid);

  {
    char title[MAX_LINE_LENGTH + 1];
    int  num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
    int  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk,
                            &num_node_sets, &num_side_sets);
    printf("after ex_get_init, error = %3d\n", error);
    if (error) {
      exit(-1);
    }
    else {
      exit(0);
    }
  }
}
