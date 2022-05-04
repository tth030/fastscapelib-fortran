#include "Error.fpp"
subroutine Advect_leapfrog (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy,alpha,cA,cB,cC
  double precision,dimension(:),allocatable :: hres, bres, etotres
  double precision,dimension(:),allocatable :: hprev_saved, bprev_saved, etotprev_saved
  integer i,j
  integer,intent(out) :: ierr

  !print*,'Advect'

  ierr = 0
  dx=xl/(nx-1)
  dy=yl/(ny-1)
  allocate(hprev_saved(nn), bprev_saved(nn), etotprev_saved(nn))
  hprev_saved=h
  bprev_saved=b
  etotprev_saved=etot

  ! x-advection using an explicit, second-order scheme to solve advection eq.
  allocate(hres(nx), bres(nx), etotres(nx))
  do j=1,ny
    do i=2,nx-1
      alpha   = (vx2(i,j)**2)*(dt**2)/(dx**2)

!debug
!      write(*,*) 'alpha ',alpha
!      write(*,*) 'h2prev(i,j) ', h2prev(i,j)
!      write(*,*) 'h2(i+1,j) h2(i-1,j) ', h2(i+1,j), h2(i-1,j)

      cA       = h2(i+1,j) + h2(i-1,j) 
      cB       = 2.d0*(1-alpha)*h2(i,j)
      cC       = h2prev(i,j)
      hres(i)  = alpha*cA + cB - cC

      cA       = b2(i+1,j) + b2(i-1,j) 
      cB       = 2.d0*(1-alpha)*b2(i,j)
      cC       = b2prev(i,j)
      bres(i)  = alpha*cA + cB - cC

      cA         = etot2(i+1,j) + etot2(i-1,j) 
      cB         = 2.d0*(1-alpha)*etot2(i,j)
      cC         = etot2prev(i,j)
      etotres(i) = alpha*cA + cB - cC
    enddo
    hres(1)     = hres(2)
    hres(nx)    = hres(nx-1)
    bres(1)     = bres(2)
    bres(nx)    = bres(nx-1)
    etotres(1)  = etotres(2)
    etotres(nx) = etotres(nx-1)
    h2(1:nx,j)    = hres
    b2(1:nx,j)    = bres
    etot2(1:nx,j) = etotres

!debug
!    do i=1,nx
!      alpha   = (vx2(i,j)**2)*(dt**2)/(dx**2)
!      write(*,'(2I7,F13.3)') i, j, h2(i,j)
!      if (h2(i,j) > 1.e8 ) then
!        write(*,*) 'voila pb '
!        write(*,*) 'alpha ',alpha
!        write(*,*) 'h2prev(i,j) ', h2prev(i,j)
!        write(*,*) 'h2(i+1,j) h2(i-1,j) ', h2(i+1,j), h2(i-1,j)
!        stop 'h2(i,j) > 1.e8'
!      endif
!    enddo

  enddo
  deallocate(hres, bres, etotres)

  ! y-advection using an explicit, second-order scheme to solve advection eq.
  allocate(hres(ny), bres(ny), etotres(ny))
  do i=1,nx
    do j=2,ny-1
      alpha   = (vy2(i,j)**2)*(dt**2)/(dx**2)

      cA       = h2(i,j+1) + h2(i,j-1) 
      cB       = 2.d0*(1-alpha)*h2(i,j)
      cC       = h2prev(i,j)
      hres(j)  = alpha*cA + cB - cC

      cA       = b2(i,j+1) + b2(i,j-1) 
      cB       = 2.d0*(1-alpha)*b2(i,j)
      cC       = b2prev(i,j)
      bres(j)  = alpha*cA + cB - cC

      cA         = etot2(i,j+1) + etot2(i,j-1) 
      cB         = 2.d0*(1-alpha)*etot2(i,j)
      cC         = etot2prev(i,j)
      etotres(j) = alpha*cA + cB - cC
    enddo
    hres(1)     = hres(2)
    hres(ny)    = hres(ny-1)
    bres(1)     = bres(2)
    bres(ny)    = bres(ny-1)
    etotres(1)  = etotres(2)
    etotres(ny) = etotres(ny-1)
    h2(i,1:ny)    = hres
    b2(i,1:ny)    = bres
    etot2(i,1:ny) = etotres
  enddo
  deallocate(hres, bres, etotres)

  b=min(b,h)
  hprev=hprev_saved
  bprev=bprev_saved
  etotprev=etotprev_saved
  deallocate(hprev_saved, bprev_saved, etotprev_saved)

  return
  end subroutine Advect_leapfrog
