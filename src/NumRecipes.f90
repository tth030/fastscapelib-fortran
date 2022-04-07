#include "Error.fpp"
module NumRecipes

use FastScapeErrorCodes

contains

  subroutine spline(x,y,yp1,ypn,y2,ierr)
    !use nrtype
    !use nrutil, only : assert_eq
    !use nr, only : tridag

    implicit none
    real ( kind = 8 ), dimension(:), intent(in) :: x,y
    real ( kind = 8 ), intent(in) :: yp1,ypn
    real ( kind = 8 ), dimension(:), intent(OUT) :: y2
    integer,intent(inout):: ierr
    !Given arrays x and y of length N containing a tabulated function, i.e., yi = f (xi), with x1 <
    !x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
    !function at points 1 and N , respectively, this routine returns an array y2 of length N
    !that contains the second derivatives of the interpolating function at the tabulated points
    !xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
    !corresponding boundary condition for a natural spline, with zero second derivative on that
    !boundary.
    integer :: n
    real ( kind = 8 ), dimension(size(x)) :: a,b,c,r

    n=assert_eq(size(x),size(y),size(y2),'spline',ierr)
    c(1:n-1)=x(2:n)-x(1:n-1)                       !Set up the tridiagonal equations.
    r(1:n-1)=6.0d0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0d0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99d30) then                     !The lower boundary condition is set either to be natural
      r(1)=0.0
      c(1)=0.0
    else                                           !or else to have a specified first derivative.
      r(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      c(1)=0.5
    end if                                         !The upper boundary condition is set either to be natural
    if (ypn > 0.99d30) then
      r(n)=0.0
      a(n)=0.0
    else                                           !or else to have a specified first derivative.
      r(n)=(-3.0d0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
      a(n)=0.5
    end if
    call tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n),ierr)

  end subroutine spline

  function splint(xa,ya,y2a,x,ierr)

    !use nrtype
    !use nrutil, only : assert_eq,nrerror
    !use nr, only: locate
    implicit none
    real ( kind = 8 ), dimension(:), intent(in) :: xa,ya,y2a
    real ( kind = 8 ), intent(in) :: x
    integer,intent(inout):: ierr
    real ( kind = 8 ) :: splint
    !Given the arrays xa and ya, which tabulate a function (with the xai in increasing or
    !decreasing order), and given the array y2a, which is the output from spline above, and
    !given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
    !and y2a are all of the same size.
    integer :: khi,klo,n
    real ( kind = 8 ) :: a,b,h

    n=assert_eq(size(xa),size(ya),size(y2a),'splint',ierr)
    klo=max(min(locate(xa,x),n-1),1)
    !We will find the right place in the table by means of locate bisection algorithm. This is
    !optimal if sequential calls to this routine are at random values of x. If sequential calls are in
    !order, and closely spaced, one would do better to store previous values of klo and khi and
    !test if they remain appropriate on the next call.
    khi=klo+1                                               !klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint',ierr)    !The xa must be distinct.
    a=(xa(khi)-x)/h                                         !Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

  end function splint

  subroutine splie2(x1a,x2a,ya,y2a,ierr)
    !USE nrtype
    !USE nrutil, ONLY : assert_eq
    !USE nr, ONLY : spline
    IMPLICIT NONE
    real ( kind = 8 ), dimension(:), intent(IN) :: x1a,x2a
    real ( kind = 8 ), dimension(:,:), intent(IN) :: ya
    real ( kind = 8 ), dimension(:,:), intent(OUT) :: y2a
    integer,intent(inout):: ierr
    !Given an M x N tabulated function ya, and N tabulated independent variables x2a, this
    !routine constructs one-dimensional natural cubic splines of the rows of ya and returns the
    !second derivatives in the M x N array y2a. (The array x1a is included in the argument
    !list merely for consistency with routine splin2.)
    integer :: j,m,ndum
    m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m',ierr)
    ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum',ierr)
    do j=1,m
      call spline(x2a,ya(j,:),1.0d30,1.0d30,y2a(j,:),ierr)
      !Values 1x10^30 signal a natural spline.
    end do
  end subroutine splie2

  function splin2(x1a,x2a,ya,y2a,x1,x2,ierr)
    !USE nrtype
    !USE nrutil, ONLY : assert_eq
    !USE nr, ONLY : spline,splint

    implicit none
    real ( kind = 8 ), dimension(:), intent(IN) :: x1a,x2a
    real ( kind = 8 ), dimension(:,:), intent(IN) :: ya,y2a
    real ( kind = 8 ), intent(IN) :: x1,x2
    integer,intent(inout):: ierr
    real ( kind = 8 ) :: splin2
    !Given x1a, x2a, ya as described in splie2 and y2a as produced by that routine; and given
    !a desired interpolating point x1,x2; this routine returns an interpolated function value by
    !bicubic spline interpolation.
    integer :: j,m,ndum
    real ( kind = 8 ), dimension(size(x1a)) :: yytmp,y2tmp2

    m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m',ierr)
    ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum',ierr)
    do j=1,m
      yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2,ierr)
      !Performm evaluations of the row splines constructed by splie2, using the one-dimensional
      !spline evaluator splint.
    end do
    call spline(x1a,yytmp,1.0d30,1.0d30,y2tmp2,ierr)
    !Construct the one-dimensional column spline and evaluate it.
    splin2=splint(x1a,yytmp,y2tmp2,x1,ierr)
  end function splin2

  function locate(xx,x)
    !use nrtype
    implicit none
    real ( kind = 8 ), dimension(:), intent(in) :: xx
    real ( kind = 8 ), intent(in) :: x
    integer :: locate
    !Given an array xx(1:N ), and given a value x, returns a value j such that x is between
    !xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
    !j = N is returned to indicate that x is out of range.
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))           !True if ascending order of table, false otherwise.
    jl=0                               !Initialize lower
    ju=n+1                             !and upper limits.
    do
      if (ju-jl <= 1) exit               !Repeat until this condition is satisfied.
      jm=(ju+jl)/2                       !Compute a midpoint,
      if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm                            !and replace either the lower limit
      else
        ju=jm                            !or the upper limit, as appropriate.
      end if
    end do
    if (x == xx(1)) then                 !Then set the output, being careful with the endpoints.
      locate=1
    else if (x == xx(n)) then
      locate=n-1
    else
      locate=jl
    end if

  end function locate

  function assert_eq(n1,n2,n3,string,ierr)
    character(len=*), intent(in) :: string
    character(len=len(string)+35) :: message
    integer, intent(in) :: n1,n2,n3
    integer,intent(inout):: ierr
    integer :: assert_eq
    if (n1 == n2 .and. n2 == n3) then
      assert_eq=n1
    else
      message='an assert_eq failed with this tag: '//string
      call nrerror(message,ierr)
      !write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      !string
      !stop 'program terminated by assert_eq'
    end if
  end function assert_eq

  function assert_eqn(nn,string,ierr)
    character(len=*), intent(in) :: string
    character(len=len(string)+35) :: message
    integer, dimension(:), intent(in) :: nn
    integer,intent(inout):: ierr
    integer :: assert_eqn
    if (all(nn(2:) == nn(1))) then
      assert_eqn=nn(1)
    else
      message='an assert_eq failed with this tag: '//string
      call nrerror(message,ierr)
      !write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      !  string
      !STOP 'program terminated by assert_eqn'
    end if
  end function assert_eqn

  subroutine nrerror(string,ierr)
    !Report a message, then die.
    character(len=*), intent(in) :: string
    integer,intent(inout):: ierr
    character(len=len(string)+8) :: message
    !write (*,*) 'nrerror:',string
    !stop 'program terminated by nrerror'
    message='nrerror:'//string
    FSCAPE_RAISE_MESSAGE(message,ERR_NrerrorNumRecipes,ierr);FSCAPE_CHKERR(ierr)
  end subroutine nrerror

  subroutine tridag_ser(a,b,c,r,u,ierr)
    !use nrtype
    !use nrutil, ONLY : assert_eq,nrerror

    implicit none
    real ( kind = 8 ), dimension(:), intent(in) :: a,b,c,r
    real ( kind = 8 ), dimension(:), intent(out) :: u
    integer,intent(inout):: ierr
    !Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
    !serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N ,
    !while a and c (off-diagonal elements) are size N − 1.
    real ( kind = 8 ), dimension(size(b)) :: gam    !One vector of workspace, gam is needed.
    integer :: n,j
    real ( kind = 8 ) :: bet
    n=assert_eqn((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser',ierr)
    bet=b(1)
    if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1',ierr)
    !If this happens then you should rewrite your equations as a set of order N − 1, with u2
    !trivially eliminated.
    u(1)=r(1)/bet
    do j=2,n                !Decomposition and forward substitution.
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      if (bet == 0.0) &     !Algorithm fails; see below routine in Vol. 1.
        call nrerror('tridag_ser: Error at code stage',ierr)
      u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1              !Backsubstitution.
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do

  end subroutine tridag_ser

end module NumRecipes
