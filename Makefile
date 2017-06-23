# A dispatch makefile for Sorbentfit running make in each of it's component subdirs
SUBDIRS := mpi omp win_omp
CLEAN_TARGETS := $(SUBDIRS:%=clean-%)

.PHONY: all $(SUBDIRS) $(CLEAN_TARGETS)

all: $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) -sC $@


clean: $(CLEAN_TARGETS)

$(CLEAN_TARGETS):
	@$(MAKE) -sC $(@:clean-%=%) clean
