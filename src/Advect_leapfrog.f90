#include "Error.fpp"
subroutine Advect_leapfrog (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy,alpha,cA,cB
  double precision,dimension(:),allocatable :: hres, bres, etotres
  double precision,dimension(:),allocatable :: hprev_x_saved, bprev_x_saved, etotprev_x_saved
  double precision,dimension(:),allocatable :: hprev_y_saved, bprev_y_saved, etotprev_y_saved
  integer i,j
  integer,intent(out) :: ierr

  !print*,'Advect'

  ierr = 0
  dx=xl/(nx-1)
  dy=yl/(ny-1)

  if (step<=2) then

    hprev_x = h
    bprev_x = b 
    etotprev_x = etot

    hprev_y = h
    bprev_y = b
    etotprev_y = etot

    call Advect_laxwendroff(ierr);FSCAPE_CHKERR(ierr)

  else

    ! y-advection using an explicit, second-order scheme to solve advection eq.
    allocate(hprev_y_saved(nn), bprev_y_saved(nn), etotprev_y_saved(nn))
    hprev_y_saved=h
    bprev_y_saved=b
    etotprev_y_saved=etot

    allocate(hres(ny), bres(ny), etotres(ny))
    do i=1,nx
      do j=2,ny-1
        alpha   = vy2(i,j)*dt/dy
  
        cA       = h2(i,j+1) - h2(i,j-1) 
        cB       = h2prev_y(i,j)
        hres(j)  = cB - alpha*cA
  
        cA       = b2(i,j+1) - b2(i,j-1) 
        cB       = b2prev_y(i,j)
        bres(j)  = cB - alpha*cA
  
        cA         = etot2(i,j+1) - etot2(i,j-1) 
        cB         = etot2prev_y(i,j)
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

 
    ! x-advection using an explicit, second-order scheme to solve advection eq.
    allocate(hprev_x_saved(nn), bprev_x_saved(nn), etotprev_x_saved(nn))
    hprev_x_saved=h
    bprev_x_saved=b
    etotprev_x_saved=etot

    allocate(hres(nx), bres(nx), etotres(nx))
    do j=1,ny
      do i=2,nx-1
        alpha   = vx2(i,j)*dt/dx
  
        cA       = h2(i+1,j) - h2(i-1,j) 
        cB       = h2prev_x(i,j)
        hres(i)  = cB - alpha*cA

        cA       = b2(i+1,j) - b2(i-1,j) 
        cB       = b2prev_x(i,j)
        bres(i)  = cB - alpha*cA
  
        cA         = etot2(i+1,j) - etot2(i-1,j) 
        cB         = etot2prev_x(i,j)
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
 
    b=min(b,h)

    hprev_x=hprev_x_saved
    bprev_x=bprev_x_saved
    etotprev_x=etotprev_x_saved

    hprev_y=hprev_y_saved
    bprev_y=bprev_y_saved
    etotprev_y=etotprev_y_saved
    deallocate(hprev_x_saved, bprev_x_saved, etotprev_x_saved)
    deallocate(hprev_y_saved, bprev_y_saved, etotprev_y_saved)
  endif

  return
  end subroutine Advect_leapfrog
