# make          builds bin/mandyocjll (shim and compiled image; a few minutes)
# make PETSC_DIR=... [PETSC_ARCH=...]
#               the same with a PETSc built from source and its MPI (writes LocalPreferences.toml), e.g. on Hypatia:
#               make PETSC_DIR=$HOME/opt/petsc PETSC_ARCH=optimized-v3.24.1-openmpi
# make petsc PETSC_DIR=... [PETSC_ARCH=...]   only writes LocalPreferences.toml
# make shim     builds only the C shim
# make clean    removes what was built
JULIA ?= julia

all: $(if $(PETSC_DIR),petsc) shim image

petsc:
	$(if $(PETSC_DIR),,$(error give PETSC_DIR (and PETSC_ARCH)))
	$(JULIA) tools/use_petsc.jl "$(PETSC_DIR)" "$(PETSC_ARCH)"

shim:
	$(JULIA) --project=. csrc/build.jl

image:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=build -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=build build/build_image.jl

clean:
	rm -f bin/mandyocjll.so bin/mandyocjll.dylib bin/mandyocjll.stamp csrc/libmandyoc_shim.so csrc/libmandyoc_shim.dylib

.PHONY: all petsc shim image clean
