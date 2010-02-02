.PHONY: all

include ../make.include

#INCLUDE = $(shell perl $(TRILINOS_INCLUDE)/strip_dup_incl_paths.pl -I$(TRILINOS_INCLUDE))
INCLUDE = -I$(Teko_INCLUDE) -I$(TRILINOS_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = Teko_Utilities.o Teko_BlockPreconditionerFactory.o \
       Teko_JacobiPreconditionerFactory.o Teko_LU2x2PreconditionerFactory.o \
       Teko_BlockImplicitLinearOp.o Teko_LU2x2InverseOp.o Teko_BlockInvDiagonalStrategy.o \
       Teko_BlockUpperTriInverseOp.o Teko_GaussSeidelPreconditionerFactory.o \
       Teko_BlockLowerTriInverseOp.o Teko_BlockedReordering.o \
       Teko_InverseFactory.o Teko_AddPreconditionerFactory.o Teko_MultPreconditionerFactory.o \
       Teko_InverseLibrary.o Teko_LU2x2DiagonalStrategy.o Teko_ImplicitLinearOp.o \
       Teko_PreconditionerFactory.o Teko_PreconditionerState.o

all: $(OBJS)
	make -C ./NS
	make -C ./Epetra

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@
	cp $@ $(Teko_DIR)/objs/$@

clean:
	-rm -f $(OBJS) $(OBJS:%=$(Teko_DIR)/objs/%)
	make -C ./NS clean
	make -C ./Epetra clean
