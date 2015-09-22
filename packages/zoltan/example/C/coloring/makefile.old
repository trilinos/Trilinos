# @HEADER
#
########################################################################
#
#  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
#                  Copyright 2012 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Karen Devine	kddevin@sandia.gov
#                    Erik Boman	        egboman@sandia.gov
#
########################################################################
#
# @HEADER
CC = mpicc
CFLAGS = -Wall -Wstrict-prototypes -O3 -m64
# CFLAGS = -Wall -Wstrict-prototypes -g -m64 -DZOLTANV31
# CFLAGS = -Wall -Wstrict-prototypes -g -m64 -D_DEBUG

BUILD = linux

ZOLTAN31 = ../../../Zoltan_v3.1_Release
ZOLTAN32 = ../../../Zoltan_v3.2_buggy

INC31 = -I$(ZOLTAN31)/src/include
LIBP31 = $(ZOLTAN31)/src/Obj_$(BUILD)
LIBS31 = -L$(LIBP31) -lzoltan -lm
INC32 = -I$(ZOLTAN32)/src/include
LIBP32 = $(ZOLTAN32)/src/Obj_$(BUILD)
LIBS32 = -L$(LIBP32) -lzoltan -lm

all: zcol zcoldriver zcoldriver-v3.1

zcol : zcol.c $(LIBP31)/libzoltan.a makefile
	$(CC) $(CFLAGS) -DZOLTANV31 -o $@ zcol.c $(INC31) $(LIBS31)


zcoldriver-v3.1 : zcoldriver.c $(LIBP31)/libzoltan.a makefile
	$(CC) $(CFLAGS) -DZOLTANV31 -o $@ zcoldriver.c $(INC31) $(LIBS31)

zcoldriver : zcoldriver.c $(LIBP32)/libzoltan.a makefile
	$(CC) $(CFLAGS) -o $@ zcoldriver.c $(INC32) $(LIBS32)

clean:
	rm -f *.o  zcoldriver zcoldriver-v3.1 zcol



