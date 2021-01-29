# fortran-einsum-example

This project goes back to a talk I gave at the [FortranCon conference 2020](https://tcevents.chem.uzh.ch/event/12/). It demonstrates how to develop generic APIs in Fortran by making use of polymorphism and generic programming techniques. This specific example implements a tensor contraction / Einstein summation API. To learn more have a look at the [presentation slides](presentation.pdf).


## Examples
sum(k) A(ijk) B(kj) = C(i)
```Fortran
class(tensor), allocatable :: a, b, c

real, dimension(:,:,:), allocatable :: data_a
real, dimension(:,:), allocatable :: data_b

! ... allocate and assign data_a, data_b ...

a = tensor(data_a)
b = tensor(data_b)

c = tensor_einsum( &
   a, [1,2,3], b, [3,2], [1])
```
sum(ij) A(ijkl) B(jim) = C(mkl)
```Fortran
class(tensor), allocatable :: a, b, c

integer, dimension(:,:,:,:), allocatable :: data_a
integer, dimension(:,:,:), allocatable :: data_b

! ... allocate and assign data_a, data_b ...

a = tensor(data_a)
b = tensor(data_b)

c = tensor_einsum( &
   a, [1,2,3,4], b, [2,1,5], [5,3,4])
```

A working example that can be compiled and run is [included](tensor_example.f90).

## Dependencies
* [Fypp](https://github.com/aradi/fypp) preprocessor
