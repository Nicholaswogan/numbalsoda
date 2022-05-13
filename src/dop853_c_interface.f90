module dop853_c_interface
  use iso_c_binding
  use dop853_module, only: dop853_class
  implicit none
  
  type, extends(dop853_class) :: dop853_custom
    procedure(rhs_function), pointer, nopass :: rhs_fcn
    type(c_ptr) :: data_ptr
    
    integer :: nt
    real(c_double), pointer :: usol(:,:)
    real(c_double), pointer :: teval(:)
  end type
  
  abstract interface
    subroutine rhs_function(t, u, du, data_ptr)
      use iso_c_binding
      real(c_double), value, intent(in) :: t
      real(c_double), intent(in) :: u(*)
      real(c_double), intent(out) :: du(*)
      type(c_ptr), value, intent(in) :: data_ptr
    end subroutine
  end interface

contains
  
  subroutine rhs(me, t, u, du)
    class(dop853_class), intent(inout) :: me
    real(c_double), intent(in) :: t
    real(c_double), intent(in), dimension(:) :: u
    real(c_double), intent(out), dimension(:) :: du
    
    select type (me)
    class is (dop853_custom)
      call me%rhs_fcn(t, u, du, me%data_ptr)
    end select
    
  end subroutine
  
  subroutine solout(me,nr,xold,x,y,irtrn,xout)
    use dop853_module, only: dop853_class
    
    class(dop853_class),intent(inout) :: me
    integer,intent(in)                :: nr
    real(c_double),intent(in)               :: xold
    real(c_double),intent(in)               :: x
    real(c_double),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(c_double),intent(out)              :: xout
    
    integer :: i
    
    select type (me)
    class is (dop853_custom)
      
      if (nr == 1) then
        xout = me%teval(1)
        me%nt = 1
      else
        do while (me%teval(me%nt) <= x .and. me%nt <= size(me%teval))
          do i = 1,size(y)
            me%usol(i,me%nt) = me%contd8(i, me%teval(me%nt))
          enddo
          me%nt = me%nt + 1
        enddo
        if (me%nt <= size(me%teval)) then
          xout = me%teval(me%nt)
        endif
      endif
    
    end select
    
  end subroutine
  
  subroutine dop853_wrapper(rhs_fcn_c, neq, u0, data_ptr, nt, &
                            teval, usol, rtol, atol, mxstep, success) bind(c,name="dop853_wrapper")
    type(c_funptr), value, intent(in) :: rhs_fcn_c
    integer(c_int), value, intent(in) :: neq
    real(c_double), intent(in) :: u0(neq)
    type(c_ptr), value, intent(in) :: data_ptr
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in), target :: teval(nt)
    real(c_double), intent(out), target :: usol(neq,nt)
    real(c_double), value, intent(in) :: rtol
    real(c_double), value, intent(in) :: atol
    integer(c_int), value, intent(in) :: mxstep
    integer(c_int), intent(out) :: success

    integer :: i, j
    integer :: idid
    integer :: icomp(neq)
    logical :: status_ok
    real(c_double) :: u(neq)
    real(c_double) :: t, tout
    
    type(dop853_custom) :: dop
    
    success = 1
    
    do i = 1,neq
      icomp(i) = i
    enddo
    
    ! initialize
    call dop%initialize(fcn=rhs, solout=solout, n=neq, nmax=mxstep, icomp=icomp, status_ok=status_ok)
    if (.not. status_ok) then
      success = 0
      return
    endif
    
    call c_f_procpointer(rhs_fcn_c, dop%rhs_fcn)
    dop%data_ptr = data_ptr
    dop%teval => teval
    dop%usol => usol
  
    ! load initial conditions
    do i = 1,neq
      u(i) = u0(i)
    enddo
  
    ! perform integration
    t = teval(1)
    tout = teval(nt)
    call dop%integrate(t, u, tout, [rtol], [atol], iout=3, idid=idid)
    if (idid < 0) then
      success = 0
      return
    endif
  
  end subroutine
  
end module