#include "Error.fpp"
subroutine Advect3d_p (ierr)

  use FastScapeContext
  use omp_lib

  ! routine to advect the topography (h), basement (b) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision dx,dy
  integer i,j
  integer ierr

  ! 3d advection scheme
  double precision, dimension(:), allocatable :: xtemp,ytemp,ztemp
  double precision, dimension(:), allocatable :: xres,yres,zres,bres,etotres
  double precision rcut, dist
  integer next, ip, nneighbours, k, counter, il, jk, l
  integer, dimension(:), allocatable :: pair
  double precision, dimension(:), allocatable :: Nnneighbours
  double precision, dimension(:), allocatable :: z_nneighbours,etot_nneighbours,b_nneighbours

  !integer :: thread_id

  !print*,'Advect3d'
  !!$omp parallel private(thread_id)
  !thread_id = omp_get_thread_num()
  !print *, 'Hello from process: ', thread_id
  !!$omp end parallel

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  !===============================
  !=====[advect free surface]=====
  !===============================

  allocate(xtemp(nn),ytemp(nn),ztemp(nn),xres(nn),yres(nn),zres(nn),bres(nn),etotres(nn))
 
  ! Initialisation
  ! not needed for omp parallel do shared(ny,nx,nn,xtemp,ytemp,ztemp,dx,dy,h2) private(counter,i,j)
  do counter=1,nn
      j = ceiling(real(counter)/real(nx))
      i = counter - (j-1)*nx
      xtemp(counter)=(i-1)*dx
      ytemp(counter)=(j-1)*dy
      ztemp(counter)=h2(i,j)
      xres(counter)=xtemp(counter)
      yres(counter)=ytemp(counter)
  enddo
  ! not needed end parallel do
 
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

  !$omp parallel do shared(ny,nx,nn,xtemp,ytemp,ztemp,etot,b,xres,yres,dx,dy,next,rcut) private(counter,i,j,Nnneighbours,nneighbours,k,l,jk,il,ip,dist,pair,z_nneighbours,etot_nneighbours,b_nneighbours)
  do counter=1,nn
      j = ceiling(real(counter)/real(nx))
      i = counter - (j-1)*nx

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
  !$omp end parallel do

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
  end subroutine Advect3d_p 
