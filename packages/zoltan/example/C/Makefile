# $Id$
#

CEXAMPLES = simpleRCB zoltanSimple zoltanExample1 subDirs phgExample

DOTDOTINCLUDE= $(subst -I../, -I../../, $(INCLUDE))
DOTDOTLIBDIR= $(subst -L../, -L../../, $(LIBDIR))

all:    $(CEXAMPLES)

% : %.c ../lib/libexzoltan.a $(ZOLTAN_LIB_DIR)/libzoltan.a
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

../lib/libexzoltan.a:
	@cd ..; $(MAKE) exampleLibrary

subDirs:
	@cd sparse_matrix; \
        $(MAKE) all "CC=$(CC)" "CFLAGS=$(CFLAGS)" \
                    "INCLUDE=$(DOTDOTINCLUDE)" \
                    "LIBDIR=$(DOTDOTLIBDIR)" "LIBS=$(LIBS)"
clean:
	rm -f *.o  $(CEXAMPLES)
	@cd sparse_matrix; $(MAKE) clean



