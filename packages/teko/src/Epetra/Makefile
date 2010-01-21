.PHONY: all

include ../../make.include

#INCLUDE = $(shell perl $(TRILINOS_INCLUDE)/strip_dup_incl_paths.pl -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE))
INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = Teko_EpetraHelpers.o \
       Teko_EpetraInverseOpWrapper.o \
       Teko_EpetraOperatorWrapper.o \
       Teko_EpetraThyraConverter.o \
       Teko_InterlacedEpetra.o \
       Teko_BlockingEpetra.o \
       Teko_StridedMappingStrategy.o \
       Teko_StridedEpetraOperator.o \
       Teko_ReorderedMappingStrategy.o \
       Teko_EpetraBlockPreconditioner.o \
       Teko_BlockedMappingStrategy.o \
       Teko_BlockedEpetraOperator.o \
       Teko_BasicMappingStrategy.o

all: $(OBJS)

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@
	cp $@ $(PB_DIR)/objs/$@

clean:
	-rm -f $(OBJS) $(OBJS:%=$(PB_DIR)/objs/%)
