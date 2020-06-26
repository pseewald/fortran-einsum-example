program tensor_example
   use tensor_lib

   implicit none

   class(tensor), allocatable :: a
   class(tensor), allocatable :: b
   class(tensor), allocatable :: c
   real, dimension(:, :, :), allocatable :: data_a_3d, data_b_3d
   real, dimension(:, :), allocatable :: data_a_2d, data_b_2d
   real, dimension(:), allocatable :: data_a_1d, data_b_1d
   integer :: i

   data_a_3d = reshape([(i, i=0, 59)], [3, 4, 5], order=[3, 2, 1])
   data_b_3d = reshape([(i, i=0, 23)], [4, 3, 2], order=[3, 2, 1])

   a = tensor(data_a_3d)
   b = tensor(data_b_3d)

   c = tensor_einsum(a, [1, 2, 3], b, [2, 1, 4], [3, 4])

   call tensor_print(c)
   deallocate (a, b, c)

   data_a_2d = reshape([(i, i=0, 9)], [5, 2], order=[2, 1])
   data_b_3d = reshape([(i, i=0, 29)], [3, 5, 2], order=[3, 2, 1])

   a = tensor(data_a_2d)
   b = tensor(data_b_3d)

   c = tensor_einsum(a, [1, 2], b, [3, 1, 4], [4, 2, 3])

   call tensor_print(c)
   deallocate (a, b, c)

   data_a_3d = reshape([(i, i=0, 23)], [3, 2, 4], order=[3, 2, 1])
   data_b_2d = reshape([(i, i=0, 7)], [4, 2], order=[2, 1])

   a = tensor(data_a_3d)
   b = tensor(data_b_2d)

   c = tensor_einsum(a, [1, 2, 3], b, [3, 2], [1])
   call tensor_print(c)
   deallocate (a, b, c)

   data_a_2d = reshape([(i, i=0, 5)], [3, 2], order=[2, 1])
   data_b_1d = reshape([(i, i=0, 4)], [5], order=[1])

   a = tensor(data_a_2d)
   b = tensor(data_b_1d)

   c = tensor_einsum(a, [1, 2], b, [3], [1, 3, 2])
   call tensor_print(c)
   deallocate (a, b, c)

   data_a_2d = reshape([(i, i=0, 5)], [3, 2], order=[2, 1])
   data_b_2d = reshape([(i, i=0, 5)], [2, 3], order=[2, 1])

   a = tensor(data_a_2d)
   b = tensor(data_b_2d)

   c = tensor_einsum(a, [1, 2], b, [2, 1])
   call tensor_print(c)
   deallocate (a, b, c)

end program
