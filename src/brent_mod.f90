module brent_mod
  use iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: brent_class

  type :: brent_class
    !! the main class
    procedure(func), pointer :: f => null() !! function to be minimized
  contains
    procedure :: set_function
    procedure :: minimize => fmin
    procedure :: find_zero => zeroin
  end type brent_class

  abstract interface
    function func(me, x) result(f)
      !! Interface to the function to be minimized.
      !! It should evaluate f(x) for any x in the interval (ax,bx)
      import :: dp, brent_class
      implicit none
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x
      real(dp) :: f
    end function func
  end interface

contains

  subroutine set_function(me, f)
    class(brent_class), intent(inout) :: me
    procedure(func) :: f

    me%f => f

  end subroutine set_function


!>
!  An approximation x to the point where f attains a minimum on
!  the interval (ax,bx) is determined.
!
!  the method used is a combination of golden section search and
!  successive parabolic interpolation. convergence is never much slower
!  than that for a fibonacci search. if f has a continuous second
!  derivative which is positive at the minimum (which is not at ax or
!  bx), then convergence is superlinear, and usually of the order of
!  about 1.324.
!
!  the function f is never evaluated at two points closer together
!  than eps*abs(fmin) + (tol/3), where eps is approximately the square
!  root of the relative machine precision. if f is a unimodal
!  function and the computed values of f are always unimodal when
!  separated by at least eps*abs(x) + (tol/3), then fmin approximates
!  the abcissa of the global minimum of f on the interval ax,bx with
!  an error less than 3*eps*abs(fmin) + tol. if f is not unimodal,
!  then fmin may approximate a local, but perhaps non-global, minimum to
!  the same accuracy.
!
!  this function subprogram is a slightly modified version of the
!  algol 60 procedure localmin given in richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!### See also
!  [1] http://www.netlib.org/fmm/fmin.f


  function fmin(me, ax, bx, tol) result(xmin)
    class(brent_class), intent(inout) :: me
    real(dp), intent(in) :: ax !! left endpoint of initial interval
    real(dp), intent(in) :: bx !! right endpoint of initial interval
    real(dp), intent(in) :: tol !! desired length of the interval of
    !! uncertainty of the final result (>=0)
    real(dp) :: xmin !! abcissa approximating the point where
    !! f attains a minimum

    real(dp) :: a, b, d, e, xm, p, q, r, tol1, tol2, u, v, w
    real(dp) :: fu, fv, fw, fx, x
    real(dp) :: abs, sqrt, sign

    real(dp), parameter :: c = (3.0_dp - sqrt(5.0_dp))/2.0_dp !! squared inverse of golden ratio
    real(dp), parameter :: half = 0.5_dp
    real(dp), parameter :: eps = sqrt(epsilon(1.0_dp))

    ! initialization

    a = ax
    b = bx
    v = a + c*(b - a)
    w = v
    x = v
    e = 0.0_dp
    fx = me%f(x)
    fv = fx
    fw = fx

    do !  main loop starts here

      xm = half*(a + b)
      tol1 = eps*abs(x) + tol/3.0_dp
      tol2 = 2.0_dp*tol1

    !  check stopping criterion

      if (abs(x - xm) <= tol2 - half*(b - a)) then
        exit
      end if

      ! is golden-section necessary

      if (abs(e) <= tol1) then

        !  a golden-section step

        if (x >= xm) then
          e = a - x
        else
          e = b - x
        end if
        d = c*e

      else

        !  fit parabola

        r = (x - w)*(fx - fv)
        q = (x - v)*(fx - fw)
        p = (x - v)*q - (x - w)*r
        q = 2.0_dp*(q - r)
        if (q > 0.0_dp) then
          p = -p
        end if
        q = abs(q)
        r = e
        e = d

        !  is parabola acceptable

        if (abs(p) >= abs(half*q*r) .or. p <= q*(a - x) .or. p >= q*(b - x)) then

          !  a golden-section step

          if (x >= xm) then
            e = a - x
          else
            e = b - x
          end if
          d = c*e

        else

          !  a parabolic interpolation step

          d = p/q
          u = x + d

          !  f must not be evaluated too close to ax or bx

          if (u - a < tol2 .or. b - u < tol2) then
            d = sign(tol1, xm - x)
          end if

        end if

      end if

      !  f must not be evaluated too close to x

      if (abs(d) >= tol1) then
        u = x + d
      else
        u = x + sign(tol1, d)
      end if
      fu = me%f(u)

      !  update a, b, v, w, and x

      if (fu <= fx) then
        if (u >= x) then
          a = x
        end if
        if (u < x) then
          b = x
        end if
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
      else
        if (u < x) then
          a = u
        end if
        if (u >= x) then
          b = u
        end if
        if (fu <= fw .or. w == x) then
          v = w
          fv = fw
          w = u
          fw = fu
        else
          if (fu <= fv .or. v == x .or. v == w) then
            v = u
            fv = fu
          end if
        end if
      end if

    end do !  end of main loop

    xmin = x

  end function fmin

!>
!  Find a zero of the function \( f(x) \) in the given interval
!  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
!  where \( \epsilon \) is the relative machine precision defined as
!  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
!
!  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
!
!### References
!  * R. P. Brent, "[An algorithm with guaranteed convergence for
!    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
!    The Computer Journal, Vol 14, No. 4., 1971.
!  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
!    Prentice-Hall, Inc., 1973.
!
!### See also
!  1. [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib


  subroutine zeroin(me, ax, bx, tol, xzero, fzero, iflag, fax, fbx)

    use iso_fortran_env, only: error_unit

    class(brent_class), intent(inout) :: me
    real(dp), intent(in) :: ax !! left endpoint of initial interval
    real(dp), intent(in) :: bx !! right endpoint of initial interval
    real(dp), intent(in) :: tol !! desired length of the interval of uncertainty of the final result (>=0)
    real(dp), intent(out) :: xzero !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(dp), intent(out) :: fzero !! value of `f` at the root (`f(xzero)`)
    integer, intent(out) :: iflag !! status flag (`-1`=error, `0`=root found)
    real(dp), intent(in), optional :: fax !! if `f(ax)` is already known, it can be input here
    real(dp), intent(in), optional :: fbx !! if `f(ax)` is already known, it can be input here

    real(dp), parameter :: eps = epsilon(1.0_dp) !! original code had d1mach(4)
    real(dp) :: a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s

    tol1 = eps + 1.0_dp

    a = ax
    b = bx

    if (present(fax)) then
      fa = fax
    else
      fa = me%f(a)
    end if
    if (present(fbx)) then
      fb = fbx
    else
      fb = me%f(b)
    end if

    !check trivial cases first:
    if (fa == 0.0_dp) then

      iflag = 0
      xzero = a
      fzero = fa

    else
      if (fb == 0.0_dp) then

        iflag = 0
        xzero = b
        fzero = fb

      else
        if (fa*fb/abs(fb) < 0.0_dp) then ! check that f(ax) and f(bx) have different signs

          c = a
          fc = fa
          d = b - a
          e = d

          do

            if (abs(fc) < abs(fb)) then
              a = b
              b = c
              c = a
              fa = fb
              fb = fc
              fc = fa
            end if

            tol1 = 2.0_dp*eps*abs(b) + 0.5_dp*tol
            xm = 0.5_dp*(c - b)
            if (abs(xm) <= tol1 .or. fb == 0.0_dp) then
              exit
            end if

          ! see if a bisection is forced
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
              s = fb/fa
              if (a /= c) then
            ! inverse quadratic interpolation
                q = fa/fc
                r = fb/fc
                p = s*(2.0_dp*xm*q*(q - r) - (b - a)*(r - 1.0_dp))
                q = (q - 1.0_dp)*(r - 1.0_dp)*(s - 1.0_dp)
              else
                ! linear interpolation
                p = 2.0_dp*xm*s
                q = 1.0_dp - s
              end if
              if (p <= 0.0_dp) then
                p = -p
              else
                q = -q
              end if
              s = e
              e = d
              if (2.0_dp*p >= 3.0_dp*xm*q - abs(tol1*q) .or. p >= abs(0.5_dp*s*q)) then
                d = xm
                e = d
              else
                d = p/q
              end if
            else
              d = xm
              e = d
            end if

            a = b
            fa = fb
            if (abs(d) <= tol1) then
              if (xm <= 0.0_dp) then
                b = b - tol1
              else
                b = b + tol1
              end if
            else
              b = b + d
            end if
            fb = me%f(b)
            if (fb*fc/abs(fc) > 0.0_dp) then
              c = a
              fc = fa
              d = b - a
              e = d
            end if

          end do

          iflag = 0
          xzero = b
          fzero = fb

        else

          iflag = -1

        end if
      end if
    end if

  end subroutine zeroin

end module

