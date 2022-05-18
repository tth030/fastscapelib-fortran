#include "Error.fpp"
program FastScapeRUN

implicit none

integer :: nx,ny,istep,nstep,nn,nfreq
double precision, dimension(:), allocatable :: u,ux,uy,h,h_init,b,etot,erate,a,chi,catchment,sedflux,sedflux_shore,kf1,kd1
!double precision, dimension(:,:,:), allocatable :: field
real :: time_in,time_out
double precision :: m,n,g1,g2,preci_rate,p_flow_dir_exp,kf2,kd2
double precision xl,yl,dx,dy,dt
double precision sealevel, poro1, poro2, z1, z2, ratio, L, kds1, kds2
double precision vexp, xpos, ypos, vel

integer i,j,ij
integer ierr

integer :: nx_used,ny_used,use_advection,advect_every_step
integer :: useNoise,runSPL,runMarine,useAdvect3d_lag
double precision :: kf1_used,kf2_used,kd1_used,kd2_used,kds1_used,kds2_used
character(len=50) :: inputtopofile

open (1, file = 'inputfile', status = 'old')
read(1,*) nx_used
read(1,*) ny_used
read(1,*) kf1_used
read(1,*) kf2_used
read(1,*) kd1_used
read(1,*) kd2_used
read(1,*) kds1_used
read(1,*) kds2_used
read(1,*) use_advection
read(1,*) advect_every_step
read(1,'(a)') inputtopofile
close(1)

useAdvect3d_lag = 1
runSPL    = 1
runMarine = 1
useNoise  = 1
write(*,*) 'useAdvect3d_lag : ',useAdvect3d_lag
write(*,*) 'useNoise  : ',useNoise
write(*,*) 'runMarine : ',runMarine
write(*,*) 'runSPL    : ',runSPL

write(*,*) 'Input topofile is ',inputtopofile

nx=nx_used  !1025  !4097
ny=ny_used  !513   !2049
nn=nx*ny
allocate (h(nn),h_init(nn),b(nn),u(nn),ux(nn),uy(nn),etot(nn),erate(nn),a(nn),chi(nn),catchment(nn),sedflux(nn),sedflux_shore(nn))

call FastScape_Init(ierr);FSCAPE_CHKERR_ABORT(ierr)
call FastScape_Set_NX_NY (nx,ny,ierr);FSCAPE_CHKERR_ABORT(ierr)
call FastScape_Set_ADVECT_EVERY_STEP (advect_every_step,ierr);FSCAPE_CHKERR_ABORT(ierr)
call FastScape_Setup(ierr);FSCAPE_CHKERR_ABORT(ierr)

!call FastScape_Use_Marine_Aggradation(.true.,ierr);FSCAPE_CHKERR_ABORT(ierr)
!call FastScape_Set_Marine_Aggradation_rate(-0.0001d0,ierr);FSCAPE_CHKERR_ABORT(ierr)

call FastScape_Set_Enforce_Marine_Mass_cons(.true.,ierr)
call FastScape_Set_Enforce_Marine_No_Erosion(.true.,ierr)
call FastScape_Set_Enforce_Marine_Sed_Below_Sealevel(.true.,ierr)

xl=2400.d3
yl=1200.d3
dx=xl/(nx-1)
dy=yl/(ny-1)
call FastScape_Set_XL_YL (xl,yl,ierr);FSCAPE_CHKERR_ABORT(ierr)
dt=1.d3
call FastScape_Set_DT (dt,ierr);FSCAPE_CHKERR_ABORT(ierr)
allocate (kf1(nn),kd1(nn))
kf1=kf1_used
kf2=kf2_used
!kf2=kf1
m=0.4d0
n=1.d0
p_flow_dir_exp = -2.d0
kd1 = kd1_used
kd2 = kd2_used
!kd2=kd1
g1=1.d0
g2=1.d0
preci_rate = 1.d0 ! precipitation rate
if (runSPL) then
  call FastScape_Set_Erosional_Parameters (kf1,kf2,m,n,kd1,kd2,g1,g2,p_flow_dir_exp,ierr);FSCAPE_CHKERR_ABORT(ierr)
  !call FastScape_Set_Precipitation_Rate (preci_rate)
endif
sealevel = 0.d0
!poro1 = 0.63d0
!poro2 = 0.49d0
poro1 = 0.0d0
poro2 = 0.0d0
ratio = 0.5d0
L = 0.5d2
kds1 = kds1_used
kds2 = kds2_used
z1 = 2.d3
z2 = 4.d3
if (runMarine) then
  call FastScape_Set_Marine_Parameters (sealevel, poro1, poro2, z1, z2, ratio, L, kds1, kds2,ierr);FSCAPE_CHKERR_ABORT(ierr)
endif

call FastScape_Set_BC (1111,ierr);FSCAPE_CHKERR_ABORT(ierr)

! Read initial topography
call random_number (h)
open (1, file = inputtopofile, status = 'old')
do i = 1,nn  
  read(1,*) h_init(i)
  h(i) = (h(i)-0.5d0)*0.25d0 + h_init(i)
end do
close(1)


!allocate (field(nx,ny,2))
!call random_number (h)
!  do j=1,ny
!    do i=1,nx
!    ij=(j-1)*nx+i
!      if (j.lt.ny/2) then
!      h(ij)=h(ij)-200.d0
!      !u(ij)=-5.d-4
!      elseif (j.gt.ny/2) then
!      h(ij)=h(ij)+1000.d0
!      !u(ij)=5.d-4
!      endif
!    enddo
!  enddo

call FastScape_Init_H (h,ierr);FSCAPE_CHKERR_ABORT(ierr)

!do j=1,ny
!    do i=1,nx
!        ij=(j-1)*nx+i
!        u(ij)=5.d-4
!        if (j.lt.ny/4) then
!            u(ij)=-5.d-4*float(j-1)/(ny/4-1)
!        elseif (j.gt.3*ny/4) then
!            u(ij)=0.d0
!        endif
!    enddo
!enddo


if (use_advection==1) then
  u = 0.d0
  call FastScape_Set_U (u,ierr);FSCAPE_CHKERR_ABORT(ierr)
 
  ! E-W extension 
  vel = 0.5e-2 
  uy  = 0.d0
  ux  = -vel
  do j=1,ny
    do i=1,nx
      ij=(j-1)*nx+i
      xpos = (i-1)*dx 
      ypos = (j-1)*dy
      if (xpos > 2100e3 .and. ypos <= 295e3) then
          ux(ij) = vel
      endif
      if (xpos > 1530e3 .and. ypos > 295e3 .and. ypos <= 600e3) then
          ux(ij) = vel
      endif
      if (xpos >  870e3 .and. ypos > 600e3 .and. ypos <= 900e3) then
          ux(ij) = vel
      endif
      if (xpos >  300e3 .and. ypos > 900e3 .and. ypos <= 1200e3) then
          ux(ij) = vel
      endif
    enddo
  enddo

!  ! oblique extension
!  vel = 0.35356e-2
!  uy  = -vel
!  ux  = -vel
!  do j=1,ny
!    do i=1,nx
!      ij=(j-1)*nx+i
!      xpos = (i-1)*dx
!      ypos = (j-1)*dy
!      if (xpos > 2100e3 .and. ypos <= 295e3) then
!          ux(ij) = vel
!          uy(ij) = vel
!      endif
!      if (xpos > 1530e3 .and. ypos > 295e3 .and. ypos <= 600e3) then
!          ux(ij) = vel
!          uy(ij) = vel
!      endif
!      if (xpos >  870e3 .and. ypos > 600e3 .and. ypos <= 900e3) then
!          ux(ij) = vel
!          uy(ij) = vel
!      endif
!      if (xpos >  300e3 .and. ypos > 900e3 .and. ypos <= 1200e3) then
!          ux(ij) = vel
!          uy(ij) = vel
!      endif
!    enddo
!  enddo

  call FastScape_Set_V (ux,uy,ierr);FSCAPE_CHKERR_ABORT(ierr)
 
  if (useAdvect3d_lag) call FastScape_Set_RunLagToEul(.true.,ierr)
endif

nstep=30e3
nfreq=200 ! frequency of output
call FastScape_View(ierr);FSCAPE_CHKERR_ABORT(ierr)
istep=0

! initial state
call FastScape_Copy_Total_Erosion (etot,ierr);FSCAPE_CHKERR_ABORT(ierr)
vexp = 50.d0
call FastScape_VTK (etot, vexp,ierr);FSCAPE_CHKERR_ABORT(ierr)

call cpu_time (time_in)
  do while (istep.le.nstep)
  call FastScape_Execute_Step (ierr);FSCAPE_CHKERR_ABORT(ierr)
  call FastScape_Get_Step (istep,ierr);FSCAPE_CHKERR_ABORT(ierr)
    if (mod(istep,nfreq)==0) then
    call FastScape_Copy_H (h,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Basement (b,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Total_Erosion (etot,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Erosion_Rate (erate,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Drainage_Area (a,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Chi (chi,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Catchment (catchment,ierr);FSCAPE_CHKERR_ABORT(ierr)
    !call FastScape_Copy_Sediment_Flux(sedflux)
    !call FastScape_Copy_Sediment_Flux_Shore(sedflux_shore)
    print*,istep
    print*,'topo',minval(h),sum(h)/nn,maxval(h)
    print*,'basement',minval(b),sum(b)/nn,maxval(b)
    print*,'etot',minval(etot),sum(etot)/nn,maxval(etot)
    print*,'erate',minval(erate),sum(erate)/nn,maxval(erate)
    print*,'a',minval(a),sum(a)/nn,maxval(a)
    print*,'chi',minval(chi),sum(chi)/nn,maxval(chi)
    print*,'catchment',minval(catchment),sum(catchment)/nn,maxval(catchment)
    call FastScape_Debug(ierr);FSCAPE_CHKERR_ABORT(ierr)
    !field(:,1)=etot
    !field(:,2)=erate
    !vexp = 2.d0
    !call FastScape_VTK (chi, vexp,ierr);FSCAPE_CHKERR_ABORT(ierr) ! if value > 0, elevation plus value is written to file. If value <0, basement and sealevel + value is written to file.
    vexp = 50.d0
    !call FastScape_VTK (etot, vexp,ierr);FSCAPE_CHKERR_ABORT(ierr) ! if value > 0, elevation plus value is written to file. If value <0, basement and sealevel + value is written to file.
    call FastScape_VTK (b, vexp,ierr);FSCAPE_CHKERR_ABORT(ierr) ! if value > 0, elevation plus value is written to file. If value <0, basement and sealevel + value is written to file.
    !call FastScape_VTK (catchment, vexp,ierr);FSCAPE_CHKERR_ABORT(ierr)
    !call FastScape_VTK (etot, vexp,ierr);FSCAPE_CHKERR_ABORT(ierr)
    !field(:,3)=sedflux
    !field(:,4)=sedflux_shore
    !call VTK (h,b,2,field,nx,ny,dx,dy,istep)
    !VTK (h,name,nf,f,fname,nx,ny,dx,dy,istep,vex)
    endif

    if (useNoise) then
      ! apply additional noise
      call FastScape_Copy_H (h,ierr);FSCAPE_CHKERR_ABORT(ierr)
      call random_number (h_init)
      where (h>sealevel+0.125d0)  h = h + (h_init-0.5d0)*0.25d0
      where (h<=sealevel+0.125d0 .and. h>=sealevel) h = h + h_init*0.125d0
      call FastScape_Set_H(h,ierr);FSCAPE_CHKERR_ABORT(ierr)
    endif
    
  enddo
call cpu_time (time_out)
print*,'Total run time',time_out-time_in

call FastScape_Destroy (ierr);FSCAPE_CHKERR_ABORT(ierr)

deallocate (u,ux,uy,h,h_init,b,a,etot,erate,chi,catchment,sedflux,sedflux_shore)
deallocate (kf1,kd1)

end program FastScapeRun
