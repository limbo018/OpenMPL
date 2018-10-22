UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S), linux)
	-include src/Unix.mk
endif
ifeq ($(UNAME_S), Darwin)
	-include src/Darwin.mk
endif

all: 
	$(MAKE) -C sat_glucose/simp
	$(MAKE) -C sat_glucose/parallel
	$(MAKE) -C mace
	$(MAKE) -C src/mpl

clean:
	-$(MAKE) clean -C sat_glucose/simp
	-$(MAKE) clean -C sat_glucose/parallel
	-$(MAKE) clean -C src/mpl
	-$(MAKE) clean -C mace 
	
.PHONY:
	all clean
