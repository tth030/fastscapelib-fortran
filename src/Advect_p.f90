#include "Error.fpp"
subroutine Advect_p (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res
  double precision, dimension(:,:), allocatable :: h_used, b_used, etot_used, vx_used, vy_used
  double precision dx,dy
  integer i,j, nx_used, ny_used, multiplicator
  integer ierr

  !bilinear interpolation
  integer counter, ii, jj
  double precision dxx, dyy, dx_used, dy_used
  double precision deltafx, deltafy, deltafxy

  !print*,'Advect'

  multiplicator = 1

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  if ( multiplicator <= 1 ) then
    nx_used = nx
    ny_used = ny
    allocate(h_used(nx,ny),b_used(nx,ny),etot_used(nx,ny),vx_used(nx,ny),vy_used(nx,ny))
    dx_used   = xl/(nx_used-1)
    dy_used   = yl/(ny_used-1)
    h_used    = h2
    b_used    = b2
    etot_used = etot2
    vx_used   = vx2
    vy_used   = vy2
  else
    nx_used = (nx-1)*multiplicator
    ny_used = (ny-1)*multiplicator
    allocate(h_used(nx_used,ny_used),b_used(nx_used,ny_used),etot_used(nx_used,ny_used),vx_used(nx_used,ny_used),vy_used(nx_used,ny_used))
    dx_used = xl/(nx_used-1)
    dy_used = yl/(ny_used-1)
    ! bilinear interpolation
    !$omp parallel do shared(nx_used,ny_used,multiplicator,dx_used,dx,dy_used,dy,h2,b2,etot2,vx2,vy2,h_used,b_used,etot_used,vx_used,vy_used) private(counter,i,j,ii,jj,dxx,dyy,deltafx,deltafy,deltafxy)
    do counter=1,nx_used*ny_used
      j  = ceiling(real(counter)/real(nx_used))
      i  = counter - (j-1)*nx_used
      ii = ceiling(real(i)/real(multiplicator))
      jj = ceiling(real(j)/real(multiplicator))
      dxx = (i-1)*dx_used - (ii-1)*dx
      dyy = (j-1)*dy_used - (jj-1)*dy

      deltafx        = h2(ii+1,jj) - h2(ii,jj) ; deltafy = h2(ii,jj+1) - h2(ii,jj) ; deltafxy = h2(ii,jj) + h2(ii+1,jj+1) - h2(ii+1,jj) - h2(ii,jj+1)
      h_used(i,j)    = deltafx*(dxx/dx_used) + deltafy*(dyy/dy_used) + deltafxy*(dxx*dyy/(dx_used*dy_used)) + h2(ii,jj)

      deltafx        = b2(ii+1,jj) - b2(ii,jj) ; deltafy = b2(ii,jj+1) - b2(ii,jj) ; deltafxy = b2(ii,jj) + b2(ii+1,jj+1) - b2(ii+1,jj) - b2(ii,jj+1)
      b_used(i,j)    = deltafx*(dxx/dx_used) + deltafy*(dyy/dy_used) + deltafxy*(dxx*dyy/(dx_used*dy_used)) + b2(ii,jj)

      deltafx        = etot2(ii+1,jj) - etot2(ii,jj) ; deltafy = etot2(ii,jj+1) - etot2(ii,jj) ; deltafxy = etot2(ii,jj) + etot2(ii+1,jj+1) - etot2(ii+1,jj) - etot2(ii,jj+1)
      etot_used(i,j) = deltafx*(dxx/dx_used) + deltafy*(dyy/dy_used) + deltafxy*(dxx*dyy/(dx_used*dy_used)) + etot2(ii,jj)

      deltafx        = vx2(ii+1,jj) - vx2(ii,jj) ; deltafy = vx2(ii,jj+1) - vx2(ii,jj) ; deltafxy = vx2(ii,jj) + vx2(ii+1,jj+1) - vx2(ii+1,jj) - vx2(ii,jj+1)
      vx_used(i,j)   = deltafx*(dxx/dx_used) + deltafy*(dyy/dy_used) + deltafxy*(dxx*dyy/(dx_used*dy_used)) + vx2(ii,jj)

      deltafx        = vy2(ii+1,jj) - vy2(ii,jj) ; deltafy = vy2(ii,jj+1) - vy2(ii,jj) ; deltafxy = vy2(ii,jj) + vy2(ii+1,jj+1) - vy2(ii+1,jj) - vy2(ii,jj+1)
      vy_used(i,j)   = deltafx*(dxx/dx_used) + deltafy*(dyy/dy_used) + deltafxy*(dxx*dyy/(dx_used*dy_used)) + vy2(ii,jj)
    enddo
    !$omp end parallel do
  endif
 
  ! x-advection using an implicit, second-order scheme to solve advection eq.

  allocate (diag(nx_used),sup(nx_used),inf(nx_used),rhs(nx_used),res(nx_used))

  do j=1,ny_used

    diag=1.d0
    sup=0.d0
    inf=0.d0

    do i=1,nx_used
      if (vx_used(i,j).gt.0.d0) then
        diag(i)=1.d0+vx_used(i,j)*dt/dx_used
        inf(i)=-vx_used(i,j)*dt/dx_used
      elseif (vx_used(i,j).lt.0.d0) then
        diag(i)=1.d0-vx_used(i,j)*dt/dx_used
        sup(i)=vx_used(i,j)*dt/dx_used
      endif
    enddo
    sup(1)=0.d0
    diag(1)=1.d0
    diag(nx_used)=1.d0
    inf(nx_used)=0.d0

    rhs=h_used(:,j)
    call tridag (inf,diag,sup,rhs,res,nx_used,ierr);FSCAPE_CHKERR(ierr)
    h_used(:,j)=res

    rhs=b_used(:,j)
    call tridag (inf,diag,sup,rhs,res,nx_used,ierr);FSCAPE_CHKERR(ierr)
    b_used(:,j)=res

    rhs=etot_used(:,j)
    call tridag (inf,diag,sup,rhs,res,nx_used,ierr);FSCAPE_CHKERR(ierr)
    etot_used(:,j)=res

  enddo

  deallocate (diag,sup,inf,rhs,res)

  ! y-advection using an implicit, second-order scheme to solve advection eq.

  allocate (diag(ny_used),sup(ny_used),inf(ny_used),rhs(ny_used),res(ny_used))

  do i=1,nx_used

    diag=1.d0
    sup=0.d0
    inf=0.d0

    do j=1,ny_used
      if (vy_used(i,j).gt.0.d0) then
        diag(j)=1.d0+vy_used(i,j)*dt/dy_used
        inf(j)=-vy_used(i,j)*dt/dy_used
      elseif (vy_used(i,j).lt.0.d0) then
        diag(j)=1.d0-vy_used(i,j)*dt/dy_used
        sup(j)=vy_used(i,j)*dt/dy_used
      endif
    enddo
    sup(1)=0.d0
    diag(1)=1.d0
    diag(ny_used)=1.d0
    inf(ny_used)=0.d0

    rhs=h_used(i,:)
    call tridag (inf,diag,sup,rhs,res,ny_used,ierr);FSCAPE_CHKERR(ierr)
    h_used(i,:)=res

    rhs=b_used(i,:)
    call tridag (inf,diag,sup,rhs,res,ny_used,ierr);FSCAPE_CHKERR(ierr)
    b_used(i,:)=res

    rhs=etot_used(i,:)
    call tridag (inf,diag,sup,rhs,res,ny_used,ierr);FSCAPE_CHKERR(ierr)
    etot_used(i,:)=res

  enddo

  deallocate (diag,sup,inf,rhs,res)

  if ( multiplicator <= 1 ) then
    h2    = h_used
    b2    = b_used
    etot2 = etot_used
  else
    !$omp parallel do shared(nn,nx,h2,b2,etot2,h_used,b_used,etot_used,multiplicator) private(counter,i,j)
    do counter=1,nn
      j  = ceiling(real(counter)/real(nx))
      i  = counter - (j-1)*nx
      h2(i,j)    = h_used((i-1)*multiplicator+1,(j-1)*multiplicator+1)
      b2(i,j)    = b_used((i-1)*multiplicator+1,(j-1)*multiplicator+1)
      etot2(i,j) = etot_used((i-1)*multiplicator+1,(j-1)*multiplicator+1)
    enddo
    !$omp end parallel do
  endif

  b=min(b,h)

  deallocate(h_used,b_used,etot_used,vx_used,vy_used)

  return
  end subroutine Advect_p
