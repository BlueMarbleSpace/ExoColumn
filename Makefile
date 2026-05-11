# Forwarding Makefile — delegates to build/Makefile.
# Usage from the project root:
#   make              build
#   make clean
#   make run
#   make USER_FC=gfortran

.PHONY: all clean run

all:
	$(MAKE) -C build $(if $(USER_FC),USER_FC=$(USER_FC))

clean:
	$(MAKE) -C build clean

run:
	run/exocol.exe
