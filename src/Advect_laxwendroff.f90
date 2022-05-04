#include "Error.fpp"
subroutine Advect_laxwendroff (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy,alpha,cA,cB,cC,cD
  double precision,dimension(:),allocatable :: hres, bres, etotres
  integer i,j
  integer,intent(out) :: ierr

  !print*,'Advect'

  ierr = 0
  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! x-advection using an explicit, second-order scheme to solve advection eq.
  allocate(hres(nx), bres(nx), etotres(nx))
  do j=1,ny
    do i=2,nx-1
      alpha   = vx2(i,j)*dt/dx

      cA       = h2(i+1,j) + h2(i,j) 
      cB       = h2(i+1,j) - h2(i,j)
      cC       = h2(i,j) + h2(i-1,j)
      cD       = h2(i,j) - h2(i-1,j)
      hres(i) = h2(i,j) - alpha*( (cA/2) - (cB/2)*alpha - (cC/2) + (cD/2)*alpha )
      
      cA       = b2(i+1,j) + b2(i,j)
      cB       = b2(i+1,j) - b2(i,j)
      cC       = b2(i,j) + b2(i-1,j)
      cD       = b2(i,j) - b2(i-1,j)
      bres(i) = b2(i,j) - alpha*( (cA/2) - (cB/2)*alpha - (cC/2) + (cD/2)*alpha )

      cA          = etot2(i+1,j) + etot2(i,j)
      cB          = etot2(i+1,j) - etot2(i,j)
      cC          = etot2(i,j) + etot2(i-1,j)
      cD          = etot2(i,j) - etot2(i-1,j)
      etotres(i) = etot2(i,j) - alpha*( (cA/2) - (cB/2)*alpha - (cC/2) + (cD/2)*alpha )
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

      cA       = h2(i,j+1) + h2(i,j) 
      cB       = h2(i,j+1) - h2(i,j)
      cC       = h2(i,j) + h2(i,j-1)
      cD       = h2(i,j) - h2(i,j-1)
      hres(j) = h2(i,j) - alpha*( (cA/2) - (cB/2)*alpha - (cC/2) + (cD/2)*alpha )
      
      cA       = b2(i,j+1) + b2(i,j)
      cB       = b2(i,j+1) - b2(i,j)
      cC       = b2(i,j) + b2(i,j-1)
      cD       = b2(i,j) - b2(i,j-1)
      bres(j) = b2(i,j) - alpha*( (cA/2) - (cB/2)*alpha - (cC/2) + (cD/2)*alpha )

      cA          = etot2(i,j+1) + etot2(i,j)
      cB          = etot2(i,j+1) - etot2(i,j)
      cC          = etot2(i,j) + etot2(i,j-1)
      cD          = etot2(i,j) - etot2(i,j-1)
      etotres(j) = etot2(i,j) - alpha*( (cA/2) - (cB/2)*alpha - (cC/2) + (cD/2)*alpha )
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

  return
  end subroutine Advect_laxwendroff
