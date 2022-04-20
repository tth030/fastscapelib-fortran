#include "Error.fpp"
subroutine Advect3d (ierr)

  use FastScapeContext

  ! routine to advect the topography (h), basement (b) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy
  integer i,j
  integer, intent(out) ::ierr

  ! 3d advection scheme
  double precision, dimension(:), allocatable :: xtemp,ytemp,ztemp
  double precision, dimension(:), allocatable :: xres,yres,zres,bres,etotres
  double precision rcut, dist
  integer next, ip, nneighbours, k, counter, il, jk, l
  integer, dimension(:), allocatable :: pair
  double precision, dimension(:), allocatable :: Nnneighbours
  double precision, dimension(:), allocatable :: z_nneighbours,etot_nneighbours,b_nneighbours

  !print*,'Advect3d'

  ierr = 0

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  !===============================
  !=====[advect free surface]=====
  !===============================

  allocate(xtemp(nn),ytemp(nn),ztemp(nn),xres(nn),yres(nn),zres(nn),bres(nn),etotres(nn))
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

  rcut=2.3d0*dx

  next = 2

  allocate(pair(21))
  allocate(Nnneighbours(21))
  allocate(z_nneighbours(21))
  allocate(etot_nneighbours(21))
  allocate(b_nneighbours(21))

  counter=0
  do j=1,ny
  do i=1,nx
     counter=counter+1
     xres(counter)=(i-1)*dx
     yres(counter)=(j-1)*dy

     ! find neighbouring nodes
     nneighbours=0
     do k=-next,next
        do l=-next,next
           jk=j+k
           il=i+l
           if (il >= 1 .and. il <= nx .and. jk >= 1 .and. jk <= ny) then
              ip=nx*(jk-1)+il

              dist = sqrt((xres(counter)-xtemp(ip))**2 + (yres(counter)-ytemp(ip))**2)

              if (dist <= rcut) then
                 nneighbours=nneighbours+1
                 pair(nneighbours)=ip
              end if

           end if
        end do
     end do

     ! store their height in array
     do k=1,nneighbours
        z_nneighbours(k)   =ztemp(pair(k))
        etot_nneighbours(k)=etot(pair(k))
        b_nneighbours(k)   =b(pair(k))
     end do

     call compute_SF3 (nneighbours,pair(1:nneighbours),nn,xtemp,ytemp,xres(counter),yres(counter),rcut,Nnneighbours(1:nneighbours))

     zres(counter)    = dot_product(Nnneighbours(1:nneighbours),z_nneighbours(1:nneighbours))
     etotres(counter) = dot_product(Nnneighbours(1:nneighbours),etot_nneighbours(1:nneighbours))
     bres(counter)    = dot_product(Nnneighbours(1:nneighbours),b_nneighbours(1:nneighbours))

  end do
  end do

  deallocate(pair)
  deallocate(Nnneighbours)
  deallocate(z_nneighbours)
  deallocate(etot_nneighbours)
  deallocate(b_nneighbours)

  h    = zres
  b    = bres
  etot = etotres
  b    = min(b,h)

  deallocate(xtemp,ytemp,ztemp,xres,yres,zres,bres,etotres)
  return
  end subroutine Advect3d 
