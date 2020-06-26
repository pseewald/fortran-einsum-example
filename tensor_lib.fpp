module tensor_lib

#:include 'macros.fypp'
@:import_kinds()

   implicit none

   ! vector shape: unknown, row vector or column vector
   integer, parameter :: unknown_vec = 0, row_vec = 1, col_vec = 2

   ! abstract tensor type
   type, abstract :: tensor
      integer, dimension(:), allocatable :: shape
      class(tensor_data), allocatable :: data
   end type

   ! rank-specific tensor types
#:for rank in ranks
   type, extends(tensor) :: tensor_${rank}$d
#:if rank == 1
      integer :: vector_type = unknown_vec
#:endif
   end type
#:endfor

   ! abstract data type
   type, abstract :: tensor_data
   end type

   ! rank- and type-specific data types
#:for rank in ranks
#:for name, type, kind in data_params
   type, extends(tensor_data) :: data_${rank}$d_${name}$
      ${type}$, allocatable :: d${shape(rank)}$
   end type
#:endfor
#:endfor

   ! constructor
   interface tensor
#:for rank in ranks
#:for name in data_name
      module procedure tensor_${rank}$d_${name}$
#:endfor
#:endfor
   end interface

   ! outer product missing in the fortran standard
   interface outer_product
#:for name in data_name
      module procedure outer_product_${name}$
#:endfor
   end interface

contains

!--------------------------------------------------------------------------------------------------!
! CONSTRUCTORS                                                                                     !
!--------------------------------------------------------------------------------------------------!

   ! create tensor with shape but no data
   function tensor_nodata(shape) result(t)
      integer, dimension(:), intent(in) :: shape
      class(tensor), allocatable :: t
#:for rank in ranks
      type(tensor_${rank}$d), allocatable :: t_${rank}$d
#:endfor

      select case (size(shape))
#:for rank in ranks
      case (${rank}$)
         allocate (t_${rank}$d)
         allocate (t_${rank}$d%shape, source=shape)
         call move_alloc(t_${rank}$d, t)
#:endfor
@:assert_select_case_default()
      end select
   end function

   ! create tensor with data
#:for rank in ranks
#:for name, type, kind in data_params
   function tensor_${rank}$d_${name}$ (data) result(t)
      ${type}$, intent(in) :: data${shape(rank)}$
      integer, dimension(${rank}$) :: sh
      type(tensor_${rank}$d), allocatable :: t_${rank}$d
      class(tensor), allocatable :: t
      type(data_${rank}$d_${name}$), allocatable :: t_data

#:if rank > 0
      sh = shape(data)
#:endif
      allocate (t_${rank}$d)
      allocate (t_${rank}$d%shape(${rank}$), source=sh)

      allocate (t_data)
      allocate (t_data%d, source=data)
      call move_alloc(t_data, t_${rank}$d%data)

      call move_alloc(t_${rank}$d, t)
   end function
#:endfor
#:endfor

!--------------------------------------------------------------------------------------------------!
! TENSOR CONTRACTION (EINSTEIN SUMMATION)                                                          !
!--------------------------------------------------------------------------------------------------!

   ! generic tensor contraction routine using einstein summation convention
   ! examples:
   ! contraction [1,2,3][2,4]=[4,3,1] (sum over index '2')
   ! ind_1 = [1,2,3], ind_2 = [2,4], ind_3 = [4,3,1]
   function tensor_einsum(tensor_1, ind_1, tensor_2, ind_2, ind_3) result(tensor_3)
      class(tensor), intent(in) :: tensor_1
      integer, dimension(:), intent(in) :: ind_1
      class(tensor), intent(in) :: tensor_2
      integer, dimension(:), intent(in) :: ind_2
      class(tensor), allocatable :: tensor_3
      integer, dimension(:), intent(in), optional :: ind_3
      integer, dimension(:), allocatable :: &
         ind_1_l, ind_1_r, ind_2_l, ind_2_r, ind_3_l, ind_3_r, t3_shape
      class(tensor), allocatable :: matrix_1, matrix_2, matrix_3
      integer :: i

      call index_einstein_to_matrix_product(ind_1, ind_2, ind_3, ind_1_l, ind_1_r, ind_2_l, ind_2_r, ind_3_l, ind_3_r)

      if (.not. present(ind_3)) then
         ind_3_l = [(i, i=1, size(ind_1_l))]
         ind_3_r = [(i, i=1, size(ind_2_r))] + size(ind_1_l)
      endif

      matrix_1 = tensor_to_matrix(tensor_1, ind_1_l, ind_1_r)
      matrix_2 = tensor_to_matrix(tensor_2, ind_2_l, ind_2_r)

      matrix_3 = matrix_product(matrix_1, matrix_2)

      allocate (t3_shape(size(ind_3_l) + size(ind_3_r)))
      t3_shape([ind_3_l, ind_3_r]) = [tensor_1%shape(ind_1_l), tensor_2%shape(ind_2_r)]

      tensor_3 = tensor_from_matrix(matrix_3, t3_shape, ind_3_l, ind_3_r)

   end function

   ! translate einstein summation indices to matrix indices s.t. matrix multiplication is equivalent to tensor contraction
   !
   ! example:
   ! contraction [1,2,3][2,4]=[4,3,1] (sum over index '2')
   !   -> ind_1 = [1,2,3], ind_2 = [2,4], ind_3 = [4,3,1]
   ! this contraction translates to a matrix multiplication
   !      [1,3|2]x[1|2]=[3,2|1]
   !      notation [1,3|2]: 1st and 3rd tensor indices are mapped to 1st matrix index
   !      and 2nd tensor index is mapped to 2nd matrix index.
   !   -> ind_1_l = [1,3], ind_1_r = [2], ind_2_l = [1], ind_2_r = [2], ind_3_l = [3, 2], ind_3_r = [1]
   subroutine index_einstein_to_matrix_product(ind_1, ind_2, ind_3, ind_1_l, ind_1_r, ind_2_l, ind_2_r, ind_3_l, ind_3_r)
      integer, dimension(:), intent(in) :: ind_1, ind_2
      integer, dimension(:), intent(in), optional :: ind_3
      integer :: n_ind, ind, match1, match2, match3, i1l, i1r, i2l, i2r, i3l, i3r
      integer, dimension(:), allocatable :: ind_1_l, ind_1_r, ind_2_l, ind_2_r, ind_3_l, ind_3_r

      n_ind = maxval([ind_1, ind_2])

      allocate (ind_1_l(n_ind), ind_1_r(n_ind), ind_2_l(n_ind), ind_2_r(n_ind), ind_3_l(n_ind), ind_3_r(n_ind))
      ind_1_l = 0; ind_1_r = 0; ind_2_l = 0; ind_2_r = 0; ind_3_l = 0; ind_3_r = 0
      i1l = 0; i1r = 0; i2l = 0; i2r = 0; i3l = 0; i3r = 0

      do ind = 1, n_ind
         match1 = findloc(ind_1, ind, dim=1)
         match2 = findloc(ind_2, ind, dim=1)
         if (present(ind_3)) then
            match3 = findloc(ind_3, ind, dim=1)
         else
            match3 = 0
         endif

         if (match1 == 0) then
            i2r = i2r + 1
            ind_2_r(i2r) = match2
            i3r = i3r + 1
            ind_3_r(i3r) = match3
         elseif (match2 == 0) then
            i1l = i1l + 1
            ind_1_l(i1l) = match1
            i3l = i3l + 1
            ind_3_l(i3l) = match3
         elseif (match1 /= 0 .and. match2 /= 0) then
            i1r = i1r + 1
            ind_1_r(i1r) = match1
            i2l = i2l + 1
            ind_2_l(i2l) = match2
         endif
      enddo

      call remove_trailing_zeros(ind_1_l)
      call remove_trailing_zeros(ind_1_r)
      call remove_trailing_zeros(ind_2_l)
      call remove_trailing_zeros(ind_2_r)
      call remove_trailing_zeros(ind_3_l)
      call remove_trailing_zeros(ind_3_r)

   end subroutine

   ! remove trailing zero elements from an array
   subroutine remove_trailing_zeros(arr)
      integer, intent(inout), dimension(:), allocatable :: arr
      integer, dimension(:), allocatable :: tmp
      integer :: last

      call move_alloc(arr, tmp)
      last = findloc(tmp, 0, dim=1)
      if (last == 0) last = size(arr) + 1
      allocate (arr, source=tmp(:last - 1))

   end subroutine

!--------------------------------------------------------------------------------------------------!
! CONVERSION TENSOR <--> MATRIX, VECTOR, SCALAR                                                    !
!--------------------------------------------------------------------------------------------------!

   ! convert tensor to matrix or vector
   ! e.g. [1,2,3] -> [1,3|2] (matrix)
   !   -> ind_l = [1,3], ind_2_r = [2]
   ! or [1,2,3] -> [1,3,2|] (column vector)
   !   -> ind_l = [1,3,2], ind_2_r = []
   function tensor_to_matrix(t, ind_l, ind_r) result(matrix)
      class(tensor), intent(in) :: t
      integer, dimension(:), intent(in) :: ind_l, ind_r
      class(tensor), allocatable :: matrix

      if (size(ind_l) > 0 .and. size(ind_r) > 0) then
         matrix = tensor_to_2d(t, ind_l, ind_r)
      elseif (size(ind_l) == 0) then
         matrix = tensor_to_1d(t, ind_r, row_vec)
      elseif (size(ind_r) == 0) then
         matrix = tensor_to_1d(t, ind_l, col_vec)
      endif

   end function

   ! convert scalar or vector or matrix to tensor
   function tensor_from_matrix(matrix, t_shape, ind_l, ind_r) result(t)
      class(tensor), intent(in) :: matrix
      integer, dimension(:), intent(in) :: t_shape
      integer, dimension(:), intent(in) :: ind_l, ind_r
      class(tensor), allocatable :: t

      select type (matrix)
      type is (tensor_0d)
         t = tensor_from_0d(matrix)
      type is (tensor_1d)
         select case (matrix%vector_type)
         case (row_vec)
            t = tensor_from_1d(matrix, t_shape, ind_r)
         case (col_vec)
            t = tensor_from_1d(matrix, t_shape, ind_l)
@:assert_select_case_default()
         end select
      type is (tensor_2d)
         t = tensor_from_2d(matrix, t_shape, ind_l, ind_r)
@:assert_select_type_default()
      end select

   end function

   ! convert tensor to vector
   function tensor_to_1d(t, ind, vector_type) result(vector)
      class(tensor), intent(in) :: t
      integer, dimension(:), intent(in) :: ind
      integer, intent(in), optional :: vector_type
      class(tensor_1d), allocatable :: vector
      integer, dimension(size(t%shape)) :: order
      integer :: i
      integer, dimension(1) :: vector_shape

      vector_shape = [product(t%shape(ind))]
      vector = tensor_nodata(vector_shape)

      order(ind) = [(i, i=1, size(t%shape))]
      vector%data = reshape_data(t%data, order, vector%shape, 1)

      if (present(vector_type)) vector%vector_type = vector_type

   end function

   ! convert tensor to matrix
   function tensor_to_2d(t, ind_1, ind_2) result(matrix)
      class(tensor), intent(in) :: t
      integer, dimension(:), intent(in) :: ind_1, ind_2
      class(tensor_2d), allocatable :: matrix
      integer, dimension(size(t%shape)) :: order
      integer :: i
      integer, dimension(2) :: matrix_shape

      matrix_shape = [product(t%shape(ind_1)), product(t%shape(ind_2))]
      matrix = tensor_nodata(matrix_shape)

      order([ind_1, ind_2]) = [(i, i=1, size(t%shape))]
      matrix%data = reshape_data(t%data, order, matrix%shape, 1)

   end function

   ! convert matrix to tensor
   function tensor_from_2d(matrix, shape, ind_1, ind_2) result(t)
      class(tensor_2d), intent(in) :: matrix
      integer, dimension(:), intent(in) :: shape
      integer, dimension(:), intent(in) :: ind_1, ind_2
      class(tensor), allocatable :: t
      integer, dimension(size(shape)) :: order

      t = tensor_nodata(shape)
      order = [ind_1, ind_2]

      t%data = reshape_data(matrix%data, order, t%shape, 2)

   end function

   ! convert vector to tensor
   function tensor_from_1d(vector, shape, ind) result(t)
      class(tensor_1d), intent(in) :: vector
      integer, dimension(:), intent(in) :: shape
      integer, dimension(:), intent(in) :: ind
      class(tensor), allocatable :: t
      integer, dimension(size(shape))  :: order

      t = tensor_nodata(shape)

      order = ind
      t%data = reshape_data(vector%data, order, t%shape, 2)
   end function

   ! convert scalar to tensor
   function tensor_from_0d(scalar) result(t)
      class(tensor_0d), intent(in) :: scalar
      class(tensor), allocatable :: t

      select type (data => scalar%data)
#:for name in data_name
      type is (data_0d_${name}$)
         t = tensor(data%d)
#:endfor
@:assert_select_type_default()
      end select

   end function

   ! reshape data for conversion tensor <-> matrix / vector
   ! direction == 1: tensor -> matrix
   ! direction == 2: matrix -> tensor
   function reshape_data(data, order, shape_out, direction) result(data_out)
      class(tensor_data), intent(in) :: data
      integer, dimension(:), intent(in) :: order, shape_out
      integer, intent(in) :: direction
      class(tensor_data), allocatable :: data_out
      integer, dimension(:), allocatable :: shape_in

      select type (data)
#:for rank in ranks[1:]
#:for name in data_name
      type is (data_${rank}$d_${name}$)
         allocate (shape_in(${rank}$))
         select case (direction)
         case (1)
            shape_in(order) = shape(data%d)
         case (2)
            shape_in = shape(data%d)
@:assert_select_case_default()
         end select
#:for rank2 in ranks[1:]
         if (${rank2}$ == size(shape_out)) then
            allocate (data_${rank2}$d_${name}$ :: data_out)
            select type (data_out)
            type is (data_${rank2}$d_${name}$)
               select case (direction)
               case (1)
                  data_out%d = reshape( &
                               reshape(data%d, shape_in(1:${rank}$), order=order), &
                               shape_out(1:${rank2}$))
               case (2)
                  data_out%d = reshape(data%d, shape_out(1:${rank2}$), order=order)
@:assert_select_case_default()
               end select
@:assert_select_type_default()
            end select
         endif
#:endfor
#:endfor
#:endfor
@:assert_select_type_default()
      end select

   end function

!--------------------------------------------------------------------------------------------------!
! MATRIX / VECTOR PRODUCTS                                                                         !
!--------------------------------------------------------------------------------------------------!

   ! generic product involving matrices (tensor_2d) and/or vectors (tensor_1d):
   ! row vector x column vector (inner product) = scalar (tensor_0d)
   ! column vector x row vector (outer product) = matrix
   ! matrix x vector = vector
   ! vector x matrix = vector
   ! matrix x matrix = matrix
   function matrix_product(matrix_1, matrix_2) result(matrix_3)
      class(tensor), intent(in) :: matrix_1, matrix_2 ! dynamic type tensor_1d or tensor_2d
      class(tensor), allocatable :: matrix_3 ! dynamic type tensor_0d, tensor_1d or tensor_2d

      select type (matrix_1)
      type is (tensor_1d)
         select type (matrix_2)
         type is (tensor_1d)
            if (matrix_1%vector_type == row_vec .and. matrix_2%vector_type == col_vec) then
               select type (data_1 => matrix_1%data)
#:for name in data_name
               type is (data_1d_${name}$)
                  select type (data_2 => matrix_2%data)
                  type is (data_1d_${name}$)
                     matrix_3 = tensor(dot_product(data_1%d, data_2%d))
@:assert_select_type_default()
                  end select
#:endfor
@:assert_select_type_default()
               end select
            elseif (matrix_1%vector_type == col_vec .and. matrix_2%vector_type == row_vec) then
               select type (data_1 => matrix_1%data)
#:for name in data_name
               type is (data_1d_${name}$)
                  select type (data_2 => matrix_2%data)
                  type is (data_1d_${name}$)
                     matrix_3 = tensor(outer_product(data_1%d, data_2%d))
@:assert_select_type_default()
                  end select
#:endfor
@:assert_select_type_default()
               end select
            else
               print *, "unknown type", matrix_1%vector_type, matrix_2%vector_type
               call abort
            endif
         type is (tensor_2d)
            select type (data_1 => matrix_1%data)
#:for name in data_name
            type is (data_1d_${name}$)
               select type (data_2 => matrix_2%data)
               type is (data_2d_${name}$)
                  matrix_3 = tensor(matmul(data_1%d, data_2%d))
@:assert_select_type_default()
               end select
#:endfor
@:assert_select_type_default()
            end select
@:assert_select_type_default()
         end select
      type is (tensor_2d)
         select type (matrix_2)
         type is (tensor_1d)
            select type (data_1 => matrix_1%data)
#:for name in data_name
            type is (data_2d_${name}$)
               select type (data_2 => matrix_2%data)
               type is (data_1d_${name}$)
                  matrix_3 = tensor(matmul(data_1%d, data_2%d))
@:assert_select_type_default()
               end select
#:endfor
@:assert_select_type_default()
            end select
         type is (tensor_2d)
            select type (data_1 => matrix_1%data)
#:for name in data_name
            type is (data_2d_${name}$)
               select type (data_2 => matrix_2%data)
               type is (data_2d_${name}$)
                  matrix_3 = tensor(matmul(data_1%d, data_2%d))
@:assert_select_type_default()
               end select
#:endfor
@:assert_select_type_default()
            end select
@:assert_select_type_default()
         end select
@:assert_select_type_default()
      end select

      select type (matrix_3)
      type is (tensor_1d)
         select type (matrix_1)
         type is (tensor_1d)
            matrix_3%vector_type = matrix_1%vector_type
         end select
         select type (matrix_2)
         type is (tensor_1d)
            matrix_3%vector_type = matrix_2%vector_type
         end select
      end select

   end function

   ! outer vector product (missing in fortran standard)
#:for name, type, kind in data_params
   function outer_product_${name}$ (vector_1, vector_2) result(matrix)
      ${type}$, dimension(:), intent(in) :: vector_1, vector_2
      integer :: k, l
      ${type}$, dimension(:, :), allocatable :: matrix

      allocate (matrix(size(vector_1), size(vector_2)))
      do k = 1, size(vector_1)
         do l = 1, size(vector_2)
            matrix(k, l) = vector_1(k)*vector_2(l)
         enddo
      enddo
   end function
#:endfor

!--------------------------------------------------------------------------------------------------!
! OUTPUT                                                                                           !
!--------------------------------------------------------------------------------------------------!

   subroutine tensor_print(t)
      class(tensor), intent(in) :: t
      integer :: ${varlist('i', RANK)}$
      select type (t)
      type is (tensor_0d)
         select type (data => t%data)
#:for name in data_name
         type is (data_0d_${name}$)
            print *, data%d
#:endfor
         end select
#:for rank in ranks[1:]
      type is (tensor_${rank}$d)
#:for dim in range(1, rank+1)
         do i${dim}$ = 1, t%shape(${dim}$)
#:endfor
            select type (data => t%data)
#:for name in data_name
            type is (data_${rank}$d_${name}$)
               print *, "(", ${varlist('i', rank)}$, ")", data%d(${varlist('i', rank)}$)
#:endfor
@:assert_select_type_default()
            end select
#:for dim in range(1, rank+1)
         enddo
#:endfor
#:endfor
@:assert_select_type_default()
      end select
   end subroutine
end module
