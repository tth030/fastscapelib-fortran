#include "Error.fpp"
subroutine Advect_leapfrog (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy,Cx,Cy,alpha,PI
!  double precision,dimension(:),allocatable :: hres, bres, etotres
!  double precision,dimension(:),allocatable :: hprev_x_saved, bprev_x_saved, etotprev_x_saved
!  double precision,dimension(:),allocatable :: hprev_y_saved, bprev_y_saved, etotprev_y_saved
  
  double precision, target, dimension(:), allocatable   :: hres, bres, etotres
  double precision, dimension(:,:), pointer, contiguous :: hres2, bres2, etotres2
  double precision, dimension(:), allocatable           :: hprev_saved, bprev_saved,etotprev_saved
  integer i,j
  integer,intent(out) :: ierr

  !print*,'Advect'

  ierr = 0
  dx=xl/(nx-1)
  dy=yl/(ny-1)

  if (step<=2) then

    hprev = h
    bprev = b 
    etotprev = etot

    call Advect_laxwendroff(ierr);FSCAPE_CHKERR(ierr)

  else

    allocate(hres(nn),bres(nn),etotres(nn))
    allocate(hprev_saved(nn),bprev_saved(nn),etotprev_saved(nn))
    hres2(1:nx,1:ny) => hres
    bres2(1:nx,1:ny) => bres
    etotres2(1:nx,1:ny) => etotres
    hprev_saved = h
    bprev_saved = b
    etotprev_saved  = etot

    PI = 4.d0*datan(1.d0)

    do i=2,nx-1
      do j=2,ny-1
        !direction
        if (vx2(i,j)==0.d0) then
          alpha = PI/2.d0
          if (vy2(i,j)<0.d0) then
          alpha = alpha*-1.d0
          endif
        else
          alpha = atan(vy2(i,j)/vx2(i,j))
        endif
        !if (vx2(i,j)>0 .and. vy2(i,j)>=0) then
        !  pass
        if (vx2(i,j)<0 .and. vy2(i,j)>=0) then
          alpha = alpha + PI
        else if (vx2(i,j)<0 .and. vy2(i,j)<0) then
          alpha = alpha + PI
        else if (vx2(i,j)>0 .and. vy2(i,j)<0) then
          if (alpha<-1.d0*PI/4.d0) alpha = 2*PI + alpha
        endif

        !Courant numbers
        Cx = vx2(i,j)*dt/dx
        Cy = vy2(i,j)*dt/dy

        !Advection
        call compute_advection_point(hres2(i,j),h2,hprev2,nx,ny,i,j,alpha,Cx,Cy)
        call compute_advection_point(bres2(i,j),b2,bprev2,nx,ny,i,j,alpha,Cx,Cy)
        call compute_advection_point(etotres2(i,j),etot2,etotprev2,nx,ny,i,j,alpha,Cx,Cy)    
    
      enddo
    enddo

    !Advection on boundaries
    call advect_boundaries(hres2,h2,hprev2,vx2,vy2,nx,ny,dt,dx,dy)
    call advect_boundaries(bres2,b2,bprev2,vx2,vy2,nx,ny,dt,dx,dy)
    call advect_boundaries(etotres2,etot2,etotprev2,vx2,vy2,nx,ny,dt,dx,dy)

    !Transfer values and deallocation
    h        = hres
    b        = bres
    b        = min(b,h)
    etot     = etotres
    hprev    = hprev_saved
    bprev    = bprev_saved
    etotprev = etotprev_saved
    deallocate(hres,bres,etotres,hprev_saved,bprev_saved,etotprev_saved)

  endif

  return

end subroutine Advect_leapfrog

!--------------------------------------------------------------------------

subroutine compute_advection_point (advected_value, array,array_prev,nx,ny,i,j,alpha, Cx, Cy)

  implicit none

  double precision, intent(out) :: advected_value
  double precision, dimension(nx,ny), intent(in) :: array, array_prev
  double precision, intent(in) :: alpha, Cx, Cy
  integer, intent(in) :: nx, ny, i, j
  double precision PI, cA, cB

  PI = 4.d0*datan(1.d0)
  
  if      (alpha>=-1.d0*PI/4.d0 .and. alpha<=     PI/4.d0) then
    cA         = array(i,j) - array(i-1,j) 
    cB         = array(i,j+1) - array(i,j-1) + array(i-1,j+1) - array(i-1,j-1)
    advected_value = array_prev(i-1,j) + array(i,j) - array(i-1,j)  - 2.d0*Cx*cA - (Cy/2.d0)*cB
  else if (alpha>=      PI/4.d0 .and. alpha<=3.d0*PI/4.d0) then
    cA         = array(i+1,j) - array(i-1,j) + array(i+1,j-1) - array(i-1,j-1)
    cB         = array(i,j) - array(i,j-1)
    advected_value = array_prev(i,j-1) + array(i,j) - array(i,j-1) - (Cx/2.d0)*cA - 2.d0*Cy*cB
  else if (alpha>= 3.d0*PI/4.d0 .and. alpha<=5.d0*PI/4.d0) then
    cA         = array(i+1,j) - array(i,j)
    cB         = array(i,j+1) - array(i,j-1) + array(i+1,j+1) - array(i+1,j-1)
    advected_value = array_prev(i+1,j) + array(i,j) - array(i+1,j) - 2.d0*Cx*cA - (Cy/2.d0)*cB
  else if (alpha>= 5.d0*PI/4.d0 .and. alpha<=7.d0*PI/4.d0) then
    cA         = array(i+1,j) - array(i-1,j) + array(i+1,j+1) - array(i-1,j+1)
    cB         = array(i,j+1) - array(i,j)
    advected_value = array_prev(i,j+1) + array(i,j) - array(i,j+1) - (Cx/2.d0)*cA - 2.d0*Cy*cB
  endif

end subroutine compute_advection_point

!--------------------------------------------------------------------------

subroutine advect_boundaries(final_array,array,array_prev,vx2,vy2,nx,ny,dt,dx,dy)

  implicit none
  
  integer, intent(in) :: nx, ny
  double precision, dimension(nx,ny), intent(inout) :: final_array
  double precision, dimension(nx,ny), intent(in) :: array, array_prev, vx2, vy2
  double precision, intent(in) :: dt,dx,dy

  double precision Cx, Cy, cA, cB
  integer i, j 
 
  !Right and top sides
  i = nx
  do j=2,ny
    Cx         = vx2(i,j)*dt/dx
    Cy         = vy2(i,j)*dt/dy
    cA         = array(i,j) - array(i-1,j) + array(i,j-1) - array(i-1,j-1)
    cB         = array(i,j) - array(i,j-1) + array(i-1,j) - array(i-1,j-1)
    final_array(i,j) = array_prev(i-1,j-1) + array(i,j) - array(i-1,j-1) - Cx*cA - Cy*cB
  enddo
  j = ny
  do i=2,nx
    Cx         = vx2(i,j)*dt/dx
    Cy         = vy2(i,j)*dt/dy
    cA         = array(i,j) - array(i-1,j) + array(i,j-1) - array(i-1,j-1)
    cB         = array(i,j) - array(i,j-1) + array(i-1,j) - array(i-1,j-1)
    final_array(i,j) = array_prev(i-1,j-1) + array(i,j) - array(i-1,j-1) - Cx*cA - Cy*cB
  enddo
  !Left and botttom sides
  i = 1
  do j=1,ny-1
    Cx         = vx2(i,j)*dt/dx
    Cy         = vy2(i,j)*dt/dy
    cA         = array(i,j) - array(i+1,j) + array(i,j+1) - array(i+1,j+1)
    cB         = array(i,j) - array(i,j+1) + array(i+1,j) - array(i+1,j+1)
    final_array(i,j) = array_prev(i+1,j+1) + array(i,j) - array(i+1,j+1) - Cx*cA - Cy*cB
  enddo
  j = 1
  do i=1,nx-1
    Cx         = vx2(i,j)*dt/dx
    Cy         = vy2(i,j)*dt/dy
    cA         = array(i,j) - array(i+1,j) + array(i,j+1) - array(i+1,j+1)
    cB         = array(i,j) - array(i,j+1) + array(i+1,j) - array(i+1,j+1)
    final_array(i,j) = array_prev(i+1,j+1) + array(i,j) - array(i+1,j+1) - Cx*cA - Cy*cB
  enddo
  !Bottom right
  j = 1 ; i = nx
  Cx         = vx2(i,j)*dt/dx
  Cy         = vy2(i,j)*dt/dy
  cA         = array(i,j) - array(i-1,j) + array(i,j+1) - array(i-1,j+1)
  cB         = array(i,j) - array(i,j+1) + array(i-1,j) - array(i-1,j+1)
  final_array(i,j) = array_prev(i-1,j+1) + array(i,j) - array(i-1,j+1) - Cx*cA - Cy*cB
  !Top left
  j = ny ; i = 1
  Cx         = vx2(i,j)*dt/dx
  Cy         = vy2(i,j)*dt/dy
  cA         = array(i,j) - array(i+1,j) + array(i,j-1) - array(i+1,j-1)
  cB         = array(i,j) - array(i,j-1) + array(i+1,j) - array(i+1,j-1)
  final_array(i,j) = array_prev(i+1,j-1) + array(i,j) - array(i+1,j-1) - Cx*cA - Cy*cB


end subroutine advect_boundaries
