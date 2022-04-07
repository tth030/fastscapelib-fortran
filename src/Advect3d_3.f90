#include "Error.fpp"
subroutine Advect3d_3 (ierr)

  use FastScapeContext
  use NumRecipes

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

  double precision, dimension(:,:), allocatable :: h2_2d,b2_2d,etot2_2d
  double precision, dimension(:), allocatable :: xcoord, ycoord
  integer :: counter
  double precision :: xres,yres

  !print*,'Advect3d'

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  !===============================
  !=====[advect free surface]=====
  !===============================

  allocate(xtemp(nn),ytemp(nn),ztemp(nn),zres(nn),bres(nn),etotres(nn))
  allocate(xcoord(nx),ycoord(ny))
  counter = 0
  do j=1,ny
    do i=1,nx
      counter=counter+1
      xtemp(counter)=(i-1)*dx
      ytemp(counter)=(j-1)*dy
      ztemp(counter)=h2(i,j)
      if (j==1) then
        xcoord(i)=(i-1)*dx
      endif
    enddo
    ycoord(j)=(j-1)*dy
  enddo
 
  xtemp=xtemp + vx * dt
  ytemp=ytemp + vy * dt
  ztemp=ztemp +  u * dt

  !============================
  !=====[resample surface]=====
  !============================

  allocate(h2_2d(nx,ny),b2_2d(nx,ny),etot2_2d(nx,ny))

  call splie2(xcoord,ycoord,h2,h2_2d,ierr)
  call splie2(xcoord,ycoord,b2,b2_2d,ierr)
  call splie2(xcoord,ycoord,etot2,etot2_2d,ierr)

  counter=0
  do j=1,ny
  do i=1,nx
     counter=counter+1
     xres=(i-1)*dx
     yres=(j-1)*dy
     zres(counter)    = splin2(xcoord,ycoord,h2,h2_2d,xres,yres,ierr)
     bres(counter)    = splin2(xcoord,ycoord,b2,b2_2d,xres,yres,ierr)
     etotres(counter) = splin2(xcoord,ycoord,etot2,etot2_2d,xres,yres,ierr)
  end do
  end do

  deallocate(h2_2d,b2_2d,etot2_2d,xcoord,ycoord)

  h    = zres
  b    = bres
  etot = etotres
  b    = min(b,h)

  deallocate(xtemp,ytemp,ztemp,zres,bres,etotres)
  return
  end subroutine Advect3d_3 
