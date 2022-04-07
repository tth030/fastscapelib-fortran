#include "Error.fpp"
subroutine Advect3d_2 (ierr)

  use FastScapeContext
  use Qshep2d_mod

  ! routine to advect the topography (h), basement (b) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy
  integer i,j
  integer ierr

  ! 3d advection scheme
  double precision, dimension(:), allocatable :: xtemp,ytemp,ztemp
  double precision, dimension(:), allocatable :: zres,bres,etotres
  double precision :: xres,yres
  integer :: nq, nr, nw, counter
  integer, dimension(:), allocatable :: lnext
  integer, dimension(:,:), allocatable :: lcell
  real ( kind = 8 ), dimension(:), allocatable :: rsq
  real ( kind = 8 ), dimension(:,:), allocatable :: aaa
  real ( kind = 8 ) xmin,ymin,rmax

  !print*,'Advect3d'

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  !===============================
  !=====[advect free surface]=====
  !===============================

  allocate(xtemp(nn),ytemp(nn),ztemp(nn),zres(nn),bres(nn),etotres(nn))
  counter = 0
  do j=1,ny
    do i=1,nx
      counter=counter+1
      xtemp(counter)=(i-1)*dx
      ytemp(counter)=(j-1)*dy
      ztemp(counter)=h2(i,j)
    enddo
  enddo
 
  xtemp=xtemp + vx * dt
  ytemp=ytemp + vy * dt
  ztemp=ztemp +  u * dt

  !============================
  !=====[resample surface]=====
  !============================

  nq = 13
  nw = 19
  nr = int(sqrt(dble(nn)/3.d0))

  allocate(lcell(nr,nr))
  allocate(lnext(nn))
  allocate(rsq(nn))
  allocate(aaa(5,nn))

  call qshep2 (nn,xtemp,ytemp,ztemp,nq,nw,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,aaa,ierr)

  counter=0
  do j=1,ny
  do i=1,nx
     counter=counter+1
     xres=(i-1)*dx
     yres=(j-1)*dy
     zres(counter) = qs2val ( xres, yres, nn, xtemp, ytemp, ztemp, nr, lcell, lnext, xmin, &
                   ymin, dx, dy, rmax, rsq, aaa )
     bres(counter) = qs2val ( xres, yres, nn, xtemp, ytemp, b, nr, lcell, lnext, xmin, &
                   ymin, dx, dy, rmax, rsq, aaa )
     etotres(counter) = qs2val ( xres, yres, nn, xtemp, ytemp, etot, nr, lcell, lnext, xmin, &
                   ymin, dx, dy, rmax, rsq, aaa )
  end do
  end do

  deallocate(lcell)
  deallocate(lnext)
  deallocate(rsq)
  deallocate(aaa)

  h    = zres
  b    = bres
  etot = etotres
  b    = min(b,h)

  deallocate(xtemp,ytemp,ztemp,zres,bres,etotres)
  return
  end subroutine Advect3d_2 
