#include "Error.fpp"
subroutine Advect_leapfrog (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy,alpha,cA,cB
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

  if (step<=2) then
    call Advect_laxwendroff(ierr);FSCAPE_CHKERR(ierr)
  else
    ! x-advection using an explicit, second-order scheme to solve advection eq.
    allocate(hres(nx), bres(nx), etotres(nx))
    do j=1,ny
      do i=2,nx-1
        alpha   = vx2(i,j)*dt/dx
  
        cA       = h2(i+1,j) - h2(i-1,j) 
        cB       = h2prev(i,j)
        hres(i)  = cB - alpha*cA

        cA       = b2(i+1,j) - b2(i-1,j) 
        cB       = b2prev(i,j)
        bres(i)  = cB - alpha*cA
  
        cA         = etot2(i+1,j) - etot2(i-1,j) 
        cB         = etot2prev(i,j)
        etotres(i) = cB - alpha*cA
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
    enddo
    deallocate(hres, bres, etotres)
  
    ! y-advection using an explicit, second-order scheme to solve advection eq.
    allocate(hres(ny), bres(ny), etotres(ny))
    do i=1,nx
      do j=2,ny-1
        alpha   = vy2(i,j)*dt/dy
  
        cA       = h2(i,j+1) - h2(i,j-1) 
        cB       = h2prev(i,j)
        hres(j)  = cB - alpha*cA
  
        cA       = b2(i,j+1) - b2(i,j-1) 
        cB       = b2prev(i,j)
        bres(j)  = cB - alpha*cA
  
        cA         = etot2(i,j+1) - etot2(i,j-1) 
        cB         = etot2prev(i,j)
        etotres(j) = cB - alpha*cA
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
  endif

  hprev=hprev_saved
  bprev=bprev_saved
  etotprev=etotprev_saved
  deallocate(hprev_saved, bprev_saved, etotprev_saved)

  return
  end subroutine Advect_leapfrog
