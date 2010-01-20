.PHONY: all

include ../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = Teko_LSCPreconditionerFactory.o Teko_SIMPLEPreconditionerFactory.o Teko_StaticLSCStrategy.o Teko_InvLSCStrategy.o \
       Teko_PCDStrategy.o

all: $(OBJS)

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@
	cp $@ $(PB_DIR)/objs/$@

clean:
	-rm -f $(OBJS) $(OBJS:%=$(PB_DIR)/objs/%)
