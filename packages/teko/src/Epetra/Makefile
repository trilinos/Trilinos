.PHONY: all

include ../../make.include

#INCLUDE = $(shell perl $(TRILINOS_INCLUDE)/strip_dup_incl_paths.pl -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE))
INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = PB_EpetraHelpers.o \
       PB_EpetraInverseOpWrapper.o \
       PB_EpetraOperatorWrapper.o \
       PB_EpetraThyraConverter.o \
       PB_InterlacedEpetra.o \
       PB_BlockingEpetra.o \
       PB_StridedMappingStrategy.o \
       PB_StridedEpetraOperator.o \
       PB_ReorderedMappingStrategy.o \
       PB_EpetraBlockPreconditioner.o \
       PB_BlockedMappingStrategy.o \
       PB_BlockedEpetraOperator.o \
       Teko_BasicMappingStrategy.o

all: $(OBJS)

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@
	cp $@ $(PB_DIR)/objs/$@

clean:
	-rm -f $(OBJS) $(OBJS:%=$(PB_DIR)/objs/%)
