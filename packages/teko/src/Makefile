.PHONY: all

include ../make.include

#INCLUDE = $(shell perl $(TRILINOS_INCLUDE)/strip_dup_incl_paths.pl -I$(TRILINOS_INCLUDE))
INCLUDE = -I$(TRILINOS_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = PB_Block2x2PreconditionerFactory.o PB_BlockJacobiPreconditionerFactory.o \
       PB_Helpers.o PB_SchurSolveLinearOp.o

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
