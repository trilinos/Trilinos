.PHONY: all

include ../make.include

#INCLUDE = $(shell perl $(TRILINOS_INCLUDE)/strip_dup_incl_paths.pl -I$(TRILINOS_INCLUDE))
INCLUDE = -I$(PB_INCLUDE) -I$(TRILINOS_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = PB_Utilities.o PB_BlockPreconditionerFactory.o \
       PB_JacobiPreconditionerFactory.o PB_LU2x2PreconditionerFactory.o \
       PB_BlockImplicitLinearOp.o PB_LU2x2InverseOp.o PB_BlockInvDiagonalStrategy.o \
       PB_BlockUpperTriInverseOp.o PB_GaussSeidelPreconditionerFactory.o \
       PB_BlockLowerTriInverseOp.o PB_BlockedReordering.o \
       PB_InverseFactory.o PB_AddPreconditionerFactory.o PB_MultPreconditionerFactory.o \
       PB_InverseLibrary.o

all: $(OBJS)
	make -C ./NS
	make -C ./Epetra

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@
	cp $@ $(PB_DIR)/objs/$@

clean:
	-rm -f $(OBJS) $(OBJS:%=$(PB_DIR)/objs/%)
	make -C ./NS clean
	make -C ./Epetra clean
