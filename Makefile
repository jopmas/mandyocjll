# make          builds bin/mandyocjll (shim and compiled image; a few minutes)
# make shim     builds only the C shim
# make clean    removes what was built
JULIA ?= julia

all: shim image

shim:
	$(JULIA) --project=. csrc/build.jl

image:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=build -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=build build/build_image.jl

clean:
	rm -f bin/mandyocjll.so bin/mandyocjll.dylib bin/mandyocjll.stamp csrc/libmandyoc_shim.so csrc/libmandyoc_shim.dylib

.PHONY: all shim image clean
