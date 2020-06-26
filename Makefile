TENSOR_CONFIG_FLAGS ?= -DRANK=7 -DDEBUG=0
FCFLAGS ?= -g
FYPPFLAGS += $(TENSOR_CONFIG_FLAGS)

OBJS = tensor_lib.o
all: tensor_example

.SUFFIXES: .f90 .o

%.o: %.mod

%.f90: %.fpp
	fypp $(FYPPFLAGS) $< $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

tensor_example: tensor_example.f90 $(OBJS)
	$(FC) -g $(FCFLAGS) $< $(OBJS) -o $@

tensor_lib.o: tensor_lib.f90

.PHONY:
clean:
	rm -f *.o *.mod tensor_example tensor_lib.f90
