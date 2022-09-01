module dop853_solve_ivp
  use iso_c_binding
  use iso_fortran_env, only: dp => real64
  use dop853_module, only: dop853_class
  implicit none
  private

  public :: dop853_solve_ivp_wrapper
  public :: dop853_solve_ivp_result

  type :: dop853_data
    ! input
    procedure(rhs_function), pointer, nopass :: rhs_fcn => NULL()
    real(dp), pointer :: t_span(:)
    logical :: t_eval_given
    logical :: events_given
    integer :: n_events
    procedure(event_function), pointer, nopass :: events_fcn => NULL()
    type(c_ptr) :: data_ptr
    integer :: neq

    ! work
    integer :: nt
    real(dp), allocatable :: events_old(:)
    real(dp), allocatable :: events_new(:)

    ! output
    character(:), allocatable :: message
    integer :: nfev = 0
    logical :: success

    real(dp), allocatable :: t(:)
    real(dp), allocatable :: y(:,:)

    real(dp), allocatable :: t_event
    real(dp), allocatable :: y_event(:)
    integer, allocatable :: ind_event
  end type

  type, extends(dop853_class) :: dop853_custom
    type(dop853_data), pointer :: d => NULL()
  end type
  
  abstract interface
    subroutine rhs_function(t, u, du, data_ptr) bind(c)
      use iso_c_binding
      real(c_double), value, intent(in) :: t
      real(c_double), intent(in) :: u(*)
      real(c_double), intent(out) :: du(*)
      type(c_ptr), value, intent(in) :: data_ptr
    end subroutine

    subroutine event_function(t, u, event, data_ptr) bind(c)
      use iso_c_binding
      real(c_double), value, intent(in) :: t
      real(c_double), intent(in) :: u(*)
      real(c_double), intent(out) :: event(*)
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
      call me%d%rhs_fcn(t, u, du, me%d%data_ptr)
    end select
    
  end subroutine

  subroutine find_root(dop, d, xold, x, ind)
    use brent_mod, only: brent_class
    type(dop853_class), intent(inout) :: dop
    type(dop853_data), intent(inout) :: d
    real(dp), intent(in) :: xold
    real(dp), intent(in) :: x
    integer, intent(in) :: ind

    real(dp), parameter :: tol = 1.0e-8_dp
    real(dp) :: xzero, fzero
    integer :: info
    integer :: i
    real(dp), allocatable :: y_tmp(:)
    real(dp), allocatable :: events_tmp(:)

    type(brent_class) :: b

    allocate(y_tmp(size(d%y,1)))
    allocate(events_tmp(size(d%events_new)))

    call b%set_function(fcn)
    call b%find_zero(xold, x, tol, xzero, fzero, info)
    if (info /= 0) then
      print*,'brent failed'
      stop 1
    endif

    allocate(d%t_event)
    d%t_event = xzero
    allocate(d%y_event(size(d%y,1)))
    do i = 1,size(d%y_event)
      d%y_event(i) = dop%contd8(i, xzero)
    enddo
    allocate(d%ind_event)
    d%ind_event = ind
    
  contains
    function fcn(me, x_) result(f)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x_
      real(dp) :: f

      do i = 1,size(y_tmp)
        y_tmp(i) = dop%contd8(i, x_)
      enddo

      call d%events_fcn(x_, y_tmp, events_tmp, d%data_ptr)

      f = events_tmp(ind)

    end function
  end subroutine
 
  subroutine solout(me,nr,xold,x,y,irtrn,xout)
    use dop853_module, only: dop853_class
    
    class(dop853_class),intent(inout) :: me
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout

    type(dop853_data), pointer :: d
    integer :: i
    real(dp) :: xx

    select type (me)
    class is (dop853_custom)
      d => me%d
    end select

    ! search for events
    if (d%events_given) then
      call d%events_fcn(x, y, d%events_new, d%data_ptr)
      if (nr == 1) then
        d%events_old(:) = d%events_new(:)
      endif

      ! check for a sign change since last step
      do i = 1,size(d%events_old)
        if ((d%events_old(i) < 0 .and. d%events_new(i) > 0) .or. &
            (d%events_old(i) > 0 .and. d%events_new(i) < 0)) then
          call find_root(me, d, xold, x, i)
          irtrn = -1 ! exit integration
          exit
        endif
      enddo

      d%events_old(:) = d%events_new(:)
    endif

    ! save the solution
    if (d%t_eval_given) then
      if (nr == 1) then
        xout = d%t(1)
        d%nt = 1
      else

        if (allocated(d%t_event)) then
          xx = d%t_event
        else
          xx = x
        endif

        do while (d%t(d%nt) <= xx .and. d%nt <= size(d%t))
          do i = 1,size(y)
            d%y(i,d%nt) = me%contd8(i, d%t(d%nt))
          enddo
          d%nt = d%nt + 1
        enddo
        if (d%nt <= size(d%t) .and. .not.d%events_given) then
          xout = d%t(d%nt)
        endif
      endif
    else
      ! save every step
      if (.not. allocated(d%t_event) .and. (x >= d%t_span(1) .and. x <= d%t_span(2))) then
        d%t = [d%t, x]
        d%y = reshape(d%y,shape=[size(y),size(d%y,2)+1], pad=y)
      endif
    endif
    
  end subroutine

  subroutine dop853_solve_ivp_wrapper(rhs_fcn_c, t_span, neq, u0, &
                                      nt, teval, n_events, events_fcn, data_ptr, &
                                      first_step, max_step, rtol, atol, &
                                      len_message, nt_out, n_events_found, d_p) bind(c,name="dop853_solve_ivp_wrapper")
    type(c_funptr), value, intent(in) :: rhs_fcn_c
    real(c_double), target, intent(in) :: t_span(2)
    integer(c_int), value, intent(in) :: neq
    real(c_double), intent(in) :: u0(neq)

    ! optional times to evaluate the solution
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in), target :: teval(nt)

    ! events stuff
    integer(c_int), value, intent(in) :: n_events
    type(c_funptr), value, intent(in) :: events_fcn

    ! data to passed to rhs
    type(c_ptr), value, intent(in) :: data_ptr

    ! options
    real(c_double), value, intent(in) :: first_step
    real(c_double), value, intent(in) :: max_step
    real(c_double), value, intent(in) :: rtol
    real(c_double), value, intent(in) :: atol
    
    ! output
    integer(c_int), intent(out) :: len_message
    integer(c_int), intent(out) :: nt_out
    integer(c_int), intent(out) :: n_events_found
    type(c_ptr), intent(out) :: d_p

    ! LOCAL VARIABLES
    integer :: i, j
    integer :: idid, iout
    integer :: icomp(neq)
    logical :: status_ok
    real(dp) :: u(neq)
    real(dp) :: tn, tout
    
    type(dop853_data), pointer :: d
    type(dop853_custom) :: dop

    allocate(d)

    ! outputs
    len_message = 0
    nt_out = 0
    n_events_found = 0
    d_p = c_loc(d)
    
    d%success = .true.

    do i = 1,neq
      icomp(i) = i
    enddo
    
    ! initialize
    call dop%initialize(fcn=rhs, solout=solout, n=neq, icomp=icomp, hinitial=first_step, &
                        hmax=max_step, status_ok=status_ok, iprint=0)
    if (.not. status_ok) then
      d%success = .false.
      d%message = 'Integrator failed to initialize.'
      len_message = len(d%message)
      nt_out = 0
      n_events_found = 0
      return
    endif

    dop%d => d

    d%t_span => t_span
    call c_f_procpointer(rhs_fcn_c, d%rhs_fcn)
    d%neq = neq

    if (nt > 0) then
      d%t_eval_given = .true.
      allocate(d%t(nt))
      allocate(d%y(neq,nt))
      d%t(:) = teval(:)
    else
      d%t_eval_given = .false.
      allocate(d%t(0))
      allocate(d%y(neq,0))
    endif

    d%n_events = n_events
    if (n_events > 0 .and. c_associated(events_fcn)) then
      d%events_given = .true.
      call c_f_procpointer(events_fcn, d%events_fcn)
    else
      d%events_given = .false.
    endif

    if (d%t_eval_given .and. d%events_given) then
      iout = 2 ! call solout every step
    elseif (d%t_eval_given .and. .not.d%events_given) then
      iout = 3 ! call solout at a specified time
    elseif (.not.d%t_eval_given .and. d%events_given) then
      iout = 2 ! call solout every step
    elseif (.not.d%t_eval_given .and. .not.d%events_given) then
      iout = 2 ! call solout every step
    endif

    d%data_ptr = data_ptr

    ! work
    d%nt = 1
    allocate(d%events_old(n_events))
    allocate(d%events_new(n_events))

    ! initial conditions
    u(:) = u0(:)

    ! perform integration
    tn = t_span(1)
    tout = t_span(2)
    call dop%integrate(tn, u, tout, [rtol], [atol], iout=iout, idid=idid)
    if (idid < 0) then
      d%success = .false.
      d%message = 'Failed to complete integration.'
      len_message = len(d%message)
      nt_out = size(d%t)
      n_events_found = 0
      return
    endif

    ! get nfev
    call dop%info(nfcn=d%nfev)

    if (d%events_given .and. allocated(d%t_event)) then
      ! an event was found. 
      n_events_found = 1
      d%message = 'A termination event occurred.'

      if (d%t_eval_given) then; block
        real(dp), allocatable :: t_tmp(:), y_tmp(:,:)
        ! trim the results
        nt_out = d%nt - 1
        allocate(t_tmp(nt_out))
        allocate(y_tmp(neq,nt_out))

        do i = 1,nt_out
          t_tmp(i) = d%t(i)
          y_tmp(:,i) = d%y(:,i)
        enddo

        deallocate(d%t,d%y)
        allocate(d%t(nt_out))
        allocate(d%y(neq,nt_out))
        d%t(:) = t_tmp(:)
        d%y(:,:) = y_tmp(:,:)

      endblock; endif
      
    else
      n_events_found = 0
      d%message = 'The solver successfully reached the end of the integration interval.'
    endif

    len_message = len(d%message)
    nt_out = size(d%t)

  end subroutine

  subroutine dop853_solve_ivp_result(d_p, neq, len_message, nt_out, message, success, nfev, t, y, &
    event_found, ind_event, t_event, y_event) bind(c,name="dop853_solve_ivp_result")
    type(c_ptr), value, intent(in) :: d_p
    integer(c_int), value, intent(in) :: neq
    integer(c_int), value, intent(in) :: len_message
    integer(c_int), value, intent(in) :: nt_out
    character(c_char), intent(out) :: message(len_message+1)
    logical(c_bool), intent(out) :: success
    integer(c_int), intent(out) :: nfev
    real(c_double), intent(out) :: t(nt_out)
    real(c_double), intent(out) :: y(neq,nt_out)
    logical(c_bool), intent(out) :: event_found
    integer(c_int), intent(out) :: ind_event
    real(c_double), intent(out) :: t_event
    real(c_double), intent(out) :: y_event(neq)

    type(dop853_data), pointer :: d

    call c_f_pointer(d_p, d)

    call copy_string_ftoc(d%message, message)

    success = d%success
    nfev = d%nfev
    t(:) = d%t(:)
    y(:,:) = d%y(:,:)
    event_found = .false.
    if (allocated(d%t_event)) then
      event_found = .true.
      ind_event = d%ind_event - 1 ! convert to python indexing
      t_event = d%t_event
      y_event(:) = d%y_event(:)
    endif

    deallocate(d)

  end subroutine

  subroutine copy_string_ftoc(stringf,stringc)
    ! utility function to convert c string to fortran string
    character(len=*), intent(in) :: stringf
    character(c_char), intent(out) :: stringc(:)
    integer j, n, n1, n2
    n1 = len_trim(stringf)  
    n2 = size(stringc) - 1
    n = min(n1, n2)
    do j=1,n    
      stringc(j) = stringf(j:j)   
    end do
    stringc(n+1) = c_null_char
  end subroutine copy_string_ftoc
  
end module