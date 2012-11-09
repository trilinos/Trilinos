#
#@HEADER
# ************************************************************************
#
#                   Trios: Trilinos I/O Support
#                 Copyright 2011 Sandia Corporation
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
#Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
#
# *************************************************************************
#@HEADER

#CC = cc
#FC = ftn
CC = mpicc
FC = mpif90
AR = ar
RM = rm -f
CFLAGS = -O2
OBJS = init.o logging.o wrappers.o
LIB = libcommsplitter.a
TEST = ctest.exe ftest.exe

CFLAGS = -O0 -g

all: $(LIB) $(TEST)

$(LIB): $(OBJS)
	$(AR) rcs $(LIB) $(OBJS)

test: $(TEST)
	@echo "******************************************************************"
	@echo "Each app should report that it belongs to a communicator of size 2"
	@echo "******************************************************************"
	mpirun -n 2 ./ctest.exe : -n 2 ./ftest.exe

debug: all

$(TEST): $(OBJS) ftest.f90 ctest.c
	$(FC) -g ftest.f90 -o ftest.exe -L. -lcommsplitter -lz
	$(CC) -g ctest.c -o ctest.exe -L. -lcommsplitter -lz

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(LIB) $(TEST) *~
