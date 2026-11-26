# make          chooses the PETSc (if LocalPreferences.toml does not exist yet) and builds bin/mandyocjll: the C shim
#               and the compiled image (a few minutes). The PETSc is the first compatible one of: PETSC_DIR/PETSC_ARCH
#               of the environment, then a PETSc built in ~/opt/petsc and the conda environment tapioca (in the
#               opposite order on macOS); see tools/use_petsc.jl
# make PETSC_DIR=... [PETSC_ARCH=...]
#               the same with this PETSc, e.g. on Hypatia:
#               make PETSC_DIR=$HOME/opt/petsc PETSC_ARCH=optimized-v3.24.1-openmpi
# make petsc    chooses the PETSc again (or uses PETSC_DIR=... PETSC_ARCH=...), without building
# make shim     builds only the C shim
# make image    builds only the compiled image
# make clean    removes what was built
JULIA ?= julia
# PETSC_DIR/PETSC_ARCH of the environment (e.g. exported in ~/.bashrc) are only the first candidates of --auto

all: preferences packages shim image

preferences:
ifeq ($(origin PETSC_DIR),command line)
	$(JULIA) --startup-file=no tools/use_petsc.jl "$(PETSC_DIR)" "$(PETSC_ARCH)"
else
	@if [ -f LocalPreferences.toml ]; then echo "LocalPreferences.toml exists: keeping its PETSc (make petsc chooses again)"; \
	else $(JULIA) --startup-file=no tools/use_petsc.jl --auto; fi
endif

petsc:
ifeq ($(origin PETSC_DIR),command line)
	$(JULIA) --startup-file=no tools/use_petsc.jl "$(PETSC_DIR)" "$(PETSC_ARCH)"
else
	$(JULIA) --startup-file=no tools/use_petsc.jl --auto
endif

packages:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'

shim: packages
	$(JULIA) --project=. csrc/build.jl

image: packages
	$(JULIA) --project=build -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=build build/build_image.jl

clean:
	rm -f bin/mandyocjll.so bin/mandyocjll.dylib bin/mandyocjll.stamp csrc/libmandyoc_shim.so csrc/libmandyoc_shim.dylib

.PHONY: all preferences packages petsc shim image clean
