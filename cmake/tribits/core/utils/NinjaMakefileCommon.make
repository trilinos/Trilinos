default_target: all
.NOTPARALLEL:

NINJA := ninja

NINJA_FLAGS :=
ifdef VERBOSE
NINJA_FLAGS += -v
endif
ifdef NP
NINJA_FLAGS += -j $(NP)
endif

STANDARD_TARGETS := install test package package_source edit_cache rebuild_cache

all $(STANDARD_TARGETS):
	$(NINJA) -C $(TOPDIR) $(NINJA_FLAGS) $(SUBDIR)/$@
$(TARGETS):
	$(NINJA) -C $(TOPDIR) $(NINJA_FLAGS) $@
clean:
	$(NINJA) -C $(TOPDIR) $(NINJA_FLAGS) -t clean $(SUBDIR)/all
help:
	@echo "This Makefile supports the following standard targets:"
	@echo ""
	@for t in "all (default)" clean help $(STANDARD_TARGETS); do echo "  $$t"; done
	@echo ""
	@echo "and the following project targets:"
	@echo ""
	@for t in $(sort $(TARGETS)); do echo "  $$t"; done

.PHONY: all clean help $(STANDARD_TARGETS) $(TARGETS)
