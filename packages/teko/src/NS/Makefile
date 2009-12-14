.PHONY: all

include ../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = PB_LSCPreconditionerFactory.o PB_SIMPLEPreconditionerFactory.o PB_StaticLSCStrategy.o PB_InvLSCStrategy.o \
       Teko_PCDStrategy.o

all: $(OBJS)

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@
	cp $@ $(PB_DIR)/objs/$@

clean:
	-rm -f $(OBJS) $(OBJS:%=$(PB_DIR)/objs/%)
