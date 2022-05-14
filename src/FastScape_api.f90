! FastScape API

! -----------------------------------------------------------------------------------------

! contains a series of subroutines that can be accessed from outside FastScape
! to setup, run and close FastScape when used as a subroutine from another program

! subroutines and their functions

! FastScape_Init ()
! Must be called before any other routine to initialize nx, ny and step

! FastScape_SetUp ()
! Must be called to allocate memory for all internal arrays
! can only be called once FastScapeSetNXNY has been used to set up nx and ny

! FastScape_Execute_Step ()
! Executes a single step solving the SPL and diffusion equations

! FastScape_Destroy ()
! Must be called to deallocate memory

! FastScape_View ()
! prints to standard output the value of nx,ny,nn,step,xl,yl,dt,K,m,n,kd,ibc
! as well as min, mean and max values of h and u

! FastScape_Set_ADVECT_EVERY_STEP
! sets the value of advect_every_step

! FastScape_Set_NX_NY (nx,ny)
! sets the value of nx and ny, the rectangular grid dimensions
! nx and ny are integer

! FastScape_Set_XL_YL (xl,yl)
! sets the value of xl,yl, the rectangular grid extent (in m)
! xl and yl are double precision

! FastScape_Set_Erosional_Parameters (k1,k2,m,n,kd1,kd2,g1,g2)
! sets the value of the erosional parameters
! k1,k2 are rate coefficient in the stream power law (in m^(1-2*m)/yr) for bedrock and sediment respectively
! m is the area exponent in the stream power law
! kd1, kd2 are the hillslope transport coefficient or diffusivity (in m^2/yr) for bedrock and sediment respectively
! g1, g2 are the sediment fluvial transport/deposition coefficients (dimensionless) for bedrock and sediment respectively
! all parameters are double precision

! FastScape_Set_Marine_Parameters (sealevel, poro1, poro2, zporo1, zporo2, ratio, length, kds1, kds2)
! sets the value of the marine transport parameters
! sl is sea level (in m)
! poro1 is surface porosity for shale (dimensionless)
! poro2 is surface porosity for sand (dimensionless)
! zporo1 is e-folding porosity depth for shale (in m)
! zporo2 is e-folding porosity depth for sand (in m)
! ratio is the ratio of sand in the incoming flux from the continent (dimensionless)
! length is the thickness of the "mixed" surface layer (in m) at the bottom of the ocean
! kds1 and kds2 are the marine transport coefficients (diffusivities) for shale and sand respectively (in m^2/yr)

! FastScape_Set_DT (dt)
! sets the time step length (in yr)
! dt is double precision

! FastScape_Set_BC (ibc)
! sets the boundary conditions
! two types are allowed (0 is reflective bc and 1 is fixed base level)
! ibc should be an integer made of 0 and 1 corresponding to the four boundaries in the
! following order: bottom, right, top and left
! ibc is integer

! FastScape_Set_U (u)
! sets the uplift velocity/rate (in m/yr)
! an array of dimension nn(=nx*ny) should be passed
! u is double precision of size nn

! FastScape_Set_V (vx,vy)
! sets the x- and y-direction advection velocities/rates (in m/yr)
! two array of dimension nn(=nx*ny) should be passed
! vx and vy are double precision of size nn

! FastScape_Init_H (h)
! sets the initial topography (in m) as well as the basement heigh to h
! an array of dimension nn(=nx*ny) should be passed
! h is double precision of size nn

! FastScape_Init_F (F)
! sets the initial sand to shale ratio to F
! an array of dimension nn(=nx*ny) should be passed
! F is double precision of size nn

! FastScape_Copy_H (h)
! returns the current topographic height (in m)
! as an array of dimension nn(=nx*ny)
! h is double precision of size nn

! FastScape_Copy_Basement (b)
! returns the current basement height (in m)
! as an array of dimension nn(=nx*ny)
! b is double precision of size nn

! FastScape_Copy_F (F)
! returns the current surface sand-to-shale ratio (in m)
! as an array of dimension nn(=nx*ny)
! F is double precision of size nn

! FastScape_Copy_dh (dh)
! returns pure silt and sand during deposition.
! as an array of dimension nn(=nx*ny)
! dh is double precision of size nn

! FastScape_Copy_Etot (etot)
! returns the current cumulative erosion (in m)
! as an array of dimension nn(=nx*ny)
! etot is double precision of size nn

! FastScape_Reset_Cumulative_Erosion ()
! resets current cumulative erosion

! FastScape_Copy_Area (area)
! returns the drainage area at each point (in m^2)
! as an array of dimension nn(=nx*ny)
! area is double precision of size nn

! FastScape_Copy_Sediment_Flux_Shoreline (sedflux)
! returns the sediment flux at the shoreline as given to the marine diffusion routine
! as an array of dimension nn(=nx*ny)
! sedflux is double precision of size nn

! FastScape_Copy_Erate (erate)
! returns the current erosion rate (in m/yr)
! as an array of dimension nn(=nx*ny)
! erate is double precision of size nn

! FastScape_Get_Sizes (nx,ny)
! returns the value of the grid size
! nx and ny are integer

! FastScape_Get_Step (step)
! returns the value of the current time step
! step is integer

! FastScape_Set_H (h)
! resets the surface topography (in m)
! an array of dimension nn(=nx*ny) should be passed
! h is double precision of size nn

! FastScape_Set_Basement (b)
! resets the basement topography (in m)
! an array of dimension nn(=nx*ny) should be passed
! b is double precision of size nn

! FastScape_Set_Precip (p)
! resets the precipitation rate (in m/yr)
! an array of dimension nn(=nx*ny) should be passed
! p is double precision of size nn

! FastScape_Debug()
! writes debugging information to the default output

! -----------------------------------------------------------------------------------------
#include "Error.fpp"

subroutine FastScape_Init(ierr)

  use FastScapeContext
  implicit none

  integer, intent(out):: ierr

  ierr=0

  call Init()

  return

end subroutine FastScape_Init

!--------------------------------------------------------------------------

subroutine FastScape_Setup(ierr)

  use FastScapeContext
  implicit none

  integer, intent(out):: ierr

  ierr=0

  if (nx.eq.0) then
    FSCAPE_RAISE_MESSAGE('FastScape_Setup(): nx cannot be zero',ERR_ParameterInvalid,ierr)
  end if
  if (nx.le.0) then
     FSCAPE_RAISE_MESSAGE('FastScape_Setup(): nx cannot be negative',ERR_ParameterOutOfRange,ierr)
  end if
  if (ny.eq.0) then
     FSCAPE_RAISE_MESSAGE('FastScape_Setup(): ny cannot be zero',ERR_ParameterInvalid,ierr)
  end if
  if (ny.le.0) then
    FSCAPE_RAISE_MESSAGE('FastScape_Setup(): ny cannot be negative',ERR_ParameterOutOfRange,ierr)
  end if
  FSCAPE_CHKERR(ierr) ! Call FSCAPE_CHKERR() so that all possible exceptions above will be displayed

  call SetUp()

  return

end subroutine FastScape_Setup

!--------------------------------------------------------------------------

subroutine FastScape_Destroy(ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr

  ierr=0

  call Destroy()

  return

end subroutine FastScape_Destroy

!--------------------------------------------------------------------------

subroutine FastScape_View(ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr

  ierr=0

  call View()

  return

end subroutine FastScape_View

!--------------------------------------------------------------------------
subroutine FastScape_Execute_Step(ierr)

  use FastScapeContext
  use omp_lib

  implicit none

  integer, intent(out):: ierr
  !real :: time_in, time_out
  double precision :: dtime_in, dtime_out
  double precision, dimension(:), allocatable :: h_before_sp, b_before_sp, etot_before_sp, erate_before_sp

  ierr=0

  totaltime = totaltime + dt

  if (runAdvect3d) then
    !call cpu_time (time_in)
    !call Advect3d (ierr)
    !!call Advect3d_2 (ierr)
    !!call Advect3d_3 (ierr)
    !call cpu_time (time_out)
    dtime_in = omp_get_wtime()
    if (runLagToEul) then
      call Advect3d_lag (ierr)
    else
      call Advect3d_p (ierr)
    endif
    dtime_out = omp_get_wtime()
    timeAdvect3d = timeAdvect3d + dtime_out-dtime_in
  else
    if (runAdvect) then
      !call cpu_time (time_in)
      !call Advect (ierr)
      !call cpu_time (time_out)
      dtime_in = omp_get_wtime()
      !call Advect_p (ierr)
      call Advect (ierr)
      !call Advect_laxwendroff (ierr)
      !call Advect_leapfrog (ierr)
      dtime_out = omp_get_wtime()
      timeAdvect = timeAdvect + dtime_out-dtime_in
    endif
  
    if (runUplift) then
      !call cpu_time (time_in)
      dtime_in = omp_get_wtime()
      call Uplift()
      dtime_out = omp_get_wtime()
      !call cpu_time (time_out)
      timeUplift = timeUplift + dtime_out-dtime_in
    endif
  endif
  
  if (runAdvect3d .and. runLagToEul) then
    allocate(h_before_sp(nn),b_before_sp(nn),etot_before_sp(nn),erate_before_sp(nn))
    h_before_sp     = h
    b_before_sp     = b
    etot_before_sp  = etot
    erate_before_sp = erate
  endif

  if (runSPL) then
    !call cpu_time (time_in)
    dtime_in = omp_get_wtime()
    if (SingleFlowDirection) then
      call FlowRoutingSingleFlowDirection (ierr);FSCAPE_CHKERR(ierr)
      call FlowAccumulationSingleFlowDirection ()
      call StreamPowerLawSingleFlowDirection ()
    else
      call FlowRouting (ierr);FSCAPE_CHKERR(ierr)
      call FlowAccumulation ()
      call StreamPowerLaw ()
    endif
    dtime_out = omp_get_wtime()
    !call cpu_time (time_out)
    timeSPL = timeSPL + dtime_out-dtime_in
  endif

  if (runDiffusion) then
    !call cpu_time (time_in)
    dtime_in = omp_get_wtime()
    call Diffusion (ierr);FSCAPE_CHKERR(ierr)
    dtime_out = omp_get_wtime()
    !call cpu_time (time_out)
    timeDiffusion = timeDiffusion + dtime_out-dtime_in
  endif

  if (runMarine) then
     !call cpu_time (time_in)
     dtime_in = omp_get_wtime()
     if (.not. use_marine_aggradation) then
       call Marine (ierr);FSCAPE_CHKERR(ierr)
     else
       call MarineAggradation(ierr);FSCAPE_CHKERR(ierr)
     end if
     dtime_out = omp_get_wtime()
     !call cpu_time (time_out)
     timeMarine = timeMarine + dtime_out-dtime_in
  endif

  if (runStrati) then
     !call cpu_time (time_in)
     dtime_in = omp_get_wtime()
     call Run_Strati ()
     dtime_out = omp_get_wtime()
     !call cpu_time (time_out)
     timeStrati = timeStrati + dtime_out-dtime_in
  endif

  if (runAdvect3d .and. runLagToEul)  then
     dtime_in = omp_get_wtime()
     call EulToLag (h_before_sp,b_before_sp,etot_before_sp,erate_before_sp,ierr)
     deallocate(h_before_sp,b_before_sp,etot_before_sp,erate_before_sp)
     dtime_out = omp_get_wtime()
     timeEulToLag = timeEulToLag + dtime_out-dtime_in
  endif

  step=step+1

  return

end subroutine FastScape_Execute_Step

!--------------------------------------------------------------------------

subroutine FastScape_Init_H(hp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: hp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call InitH(hp)

  return

end subroutine FastScape_Init_H

!--------------------------------------------------------------------------

subroutine FastScape_Init_F(Fmixp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: Fmixp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call InitF (Fmixp)

  return

end subroutine FastScape_Init_F

!--------------------------------------------------------------------------

subroutine FastScape_Copy_H(hp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: hp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyH(hp)

  return

end subroutine FastScape_Copy_H

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Basement(bp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: bp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyBasement(bp)

  return

end subroutine FastScape_Copy_Basement

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Total_Erosion (etotp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: etotp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyEtot(etotp)

  return

end subroutine FastScape_Copy_Total_Erosion

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Drainage_Area (ap,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: ap

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyArea(ap)

  return

end subroutine FastScape_Copy_Drainage_Area

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Erosion_Rate (eratep,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: eratep

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyERate(eratep)

  return

end subroutine FastScape_Copy_Erosion_Rate

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Chi (chip,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: chip

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyChi(chip)

  return

end subroutine FastScape_Copy_Chi

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Slope (slopep,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: slopep

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopySlope(slopep)

  return

end subroutine FastScape_Copy_Slope

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Curvature (curvaturep,ierr)

  use FastScapeContext

  implicit none


  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: curvaturep

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyCurvature(curvaturep)

  return

end subroutine FastScape_Copy_Curvature

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Catchment (catchp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: catchp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyCatchment (catchp)

  return

end subroutine FastScape_Copy_Catchment

!--------------------------------------------------------------------------

subroutine FastScape_Copy_F(Fmixp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: Fmixp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyF(Fmixp)

  return

end subroutine FastScape_Copy_F

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Lake_Depth(Lp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: Lp

  ierr=0
  if (.not.setup_has_been_run) then
    FSCAPE_RAISE(ERR_SetupNotRun,ierr);FSCAPE_CHKERR(ierr)
  end if

  call CopyLakeDepth(Lp)

  return

end subroutine FastScape_Copy_Lake_Depth

!--------------------------------------------------------------------------

subroutine FastScape_Set_ADVECT_EVERY_STEP (every_step,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  integer, intent(in) :: every_step

  ierr=0

  call SetADVECTEVERYSTEP (every_step)

  return

end subroutine FastScape_Set_ADVECT_EVERY_STEP

!--------------------------------------------------------------------------

subroutine FastScape_Set_NX_NY (nnx,nny,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  integer, intent(in) :: nnx,nny

  ierr=0

  call SetNXNY (nnx,nny)

  return

end subroutine FastScape_Set_NX_NY

!--------------------------------------------------------------------------

subroutine FastScape_Set_XL_YL (xxl,yyl,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(in) :: xxl,yyl

  ierr=0

  call SetXLYL (xxl,yyl)

  return

end subroutine FastScape_Set_XL_YL

!--------------------------------------------------------------------------

subroutine FastScape_Set_DT (dtt,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(in) :: dtt

  ierr=0

  call SetDT (dtt)

  return

end subroutine FastScape_Set_DT

!--------------------------------------------------------------------------

subroutine FastScape_Set_Erosional_Parameters (kkf,kkfsed,mm,nnn,kkd,kkdsed,gg1,gg2,pp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(in), dimension(*) :: kkf,kkd
  double precision, intent(in) :: kkfsed,mm,nnn,kkdsed,gg1,gg2,pp

  ierr=0

  call SetErosionalParam (kkf,kkfsed,mm,nnn,kkd,kkdsed,gg1,gg2,pp)

  return

end subroutine FastScape_Set_Erosional_Parameters

!--------------------------------------------------------------------------

subroutine FastScape_Set_Marine_Parameters (sl, p1, p2, z1, z2, r, l, kds1, kds2,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(in) :: sl, p1, p2, z1, z2, r, l, kds1, kds2

ierr=0

call SetMarineParam (sl, p1, p2, z1, z2, r, l, kds1, kds2)

return

end subroutine FastScape_Set_Marine_Parameters

!--------------------------------------------------------------------------

subroutine FastScape_Get_Sizes (nnx,nny,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  integer, intent(out) :: nnx,nny

  ierr=0

  call GetSizes (nnx,nny)

  return

end subroutine FastScape_Get_Sizes

!--------------------------------------------------------------------------

subroutine FastScape_Get_Step (sstep,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  integer, intent(out) :: sstep

  ierr=0

  call GetStep (sstep)

  return

end subroutine FastScape_Get_Step

!--------------------------------------------------------------------------

subroutine FastScape_Debug(ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr

  ierr=0

  call Debug()

  return

end subroutine FastScape_Debug

!--------------------------------------------------------------------------

subroutine FastScape_Set_BC(jbc,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  integer, intent(in) :: jbc

  ierr=0

  call SetBC (jbc)

  return

end subroutine FastScape_Set_BC

!--------------------------------------------------------------------------

subroutine FastScape_Set_U (up,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(in), dimension(*) :: up

  ierr=0

  call SetU(up)

  return

end subroutine FastScape_Set_U

!--------------------------------------------------------------------------

subroutine FastScape_Set_V (ux,uy,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(in), dimension(*) :: ux,uy

  ierr=0

  call SetV(ux,uy)

  return

end subroutine FastScape_Set_V

!--------------------------------------------------------------------------

subroutine FastScape_Reset_Cumulative_Erosion (ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr

  ierr=0

  call ResetCumulativeErosion ()

  return

end subroutine FastScape_Reset_Cumulative_Erosion

!--------------------------------------------------------------------------

subroutine FastScape_Set_H(hp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: hp

  ierr=0

  call SetH(hp)

  return

end subroutine FastScape_Set_H

!--------------------------------------------------------------------------

subroutine FastScape_Set_All_Layers (dhp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: dhp

  ierr=0

  call SetAllLayers(dhp)

  return

end subroutine FastScape_Set_All_Layers

!--------------------------------------------------------------------------

subroutine FastScape_Set_Basement(bp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: bp

  ierr=0

  call SetBasement(bp)

  return

end subroutine FastScape_Set_Basement

!--------------------------------------------------------------------------

subroutine FastScape_Set_Precip (precipp,ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: precipp

  ierr=0

  call SetPrecip (precipp)

  return

end subroutine FastScape_Set_Precip

!--------------------------------------------------------------------------

subroutine FastScape_VTK (fp, vexp, ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(inout), dimension(*) :: fp
  double precision, intent(inout) :: vexp

  ierr=0

  call Make_VTK (fp, vexp)

  return

end subroutine FastScape_VTK

!--------------------------------------------------------------------------

subroutine FastScape_Strati (nstepp, nreflectorp, nfreqp, vexp, ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  integer, intent(inout) :: nstepp, nreflectorp, nfreqp
  double precision, intent(inout) :: vexp

  ierr=0

  call Activate_Strati (nstepp, nreflectorp, nfreqp, vexp)

  return

end subroutine FastScape_Strati

!--------------------------------------------------------------------------

subroutine FastScape_Get_Fluxes (ttectonic_flux, eerosion_flux, bboundary_flux, ierr)

  use FastScapeContext

  implicit none

  integer, intent(out):: ierr
  double precision, intent(out) :: ttectonic_flux, eerosion_flux, bboundary_flux

  ierr=0

  call compute_fluxes (ttectonic_flux, eerosion_flux, bboundary_flux)

  return

end subroutine FastScape_Get_Fluxes

!--------------------------------------------------------------------------

subroutine FastScape_Copy_dh(dhp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(inout), dimension(*) :: dhp

ierr=0

call Copydh(dhp)

return

end subroutine FastScape_Copy_dh

!--------------------------------------------------------------------------

subroutine FastScape_Copy_RockType(rocktype,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
integer, intent(inout), dimension(*) :: rocktype

ierr=0

call CopyRockType(rocktype)

return

end subroutine FastScape_Copy_RockType

!--------------------------------------------------------------------------

subroutine FastScape_Set_RockType(rocktype,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
integer, intent(in), dimension(*) :: rocktype

ierr=0

call SetRockType(rocktype)

return

end subroutine FastScape_Set_RockType

!--------------------------------------------------------------------------

subroutine FastScape_Copy_Sediment_Flux_Shoreline(sedfluxshorep,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(inout), dimension(*) :: sedfluxshorep

ierr=0

call CopySedimentFluxShore(sedfluxshorep)

return

end subroutine FastScape_Copy_Sediment_Flux_Shoreline

!--------------------------------------------------------------------------

subroutine FastScape_Set_Cumulative_Erosion (etotp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(in) :: etotp(*)

ierr=0

call SetCumulativeErosion (etotp)

return

end subroutine FastScape_Set_Cumulative_Erosion

!--------------------------------------------------------------------------

subroutine FastScape_Set_Enforce_Marine_Mass_cons (enforce_marine_mass_consp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
logical, intent(in) :: enforce_marine_mass_consp

ierr=0

call SetEnforceMarineMassCons (enforce_marine_mass_consp)

return

end subroutine FastScape_Set_Enforce_Marine_Mass_cons

!--------------------------------------------------------------------------

subroutine FastScape_Set_Enforce_Marine_Sed_Below_Sealevel (enforce_marine_sed_below_sealevelp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
logical, intent(in) :: enforce_marine_sed_below_sealevelp

ierr=0

call SetEnforceMarineSedBelowSealevel (enforce_marine_sed_below_sealevelp)

return

end subroutine FastScape_Set_Enforce_Marine_Sed_Below_Sealevel

!--------------------------------------------------------------------------

subroutine FastScape_Set_Enforce_Marine_No_Erosion (enforce_marine_no_erosionp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
logical, intent(in) :: enforce_marine_no_erosionp

ierr=0

call SetEnforceMarineNoErosion (enforce_marine_no_erosionp)

return

end subroutine FastScape_Set_Enforce_Marine_No_Erosion

!--------------------------------------------------------------------------

subroutine FastScape_Set_RunLagToEul (runLagToEulp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
logical, intent(in) :: runLagToEulp

ierr=0

call SetRunLagToEul (runLagToEulp)

return

end subroutine FastScape_Set_RunLagToEul

!--------------------------------------------------------------------------

subroutine FastScape_Set_Correct_Shallow_Sealevel (low_sealevel_at_shallow_seap,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
logical, intent(in) :: low_sealevel_at_shallow_seap

ierr=0

call SetCorrectShallowSealevel (low_sealevel_at_shallow_seap)

return

end subroutine FastScape_Set_Correct_Shallow_Sealevel

!--------------------------------------------------------------------------

subroutine FastScape_Set_Sealevel (sealevelp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(in) :: sealevelp

ierr=0

call SetSealevel (sealevelp)

return

end subroutine FastScape_Set_Sealevel

!--------------------------------------------------------------------------
subroutine FastScape_Set_Tolerance_SPL (atol_SPLp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(in) :: atol_SPLp

ierr=0

call SetAtolSPL (atol_SPLp)

return

end subroutine FastScape_Set_Tolerance_SPL

!--------------------------------------------------------------------------

subroutine FastScape_Use_Marine_Aggradation (use_marine_aggp,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
logical, intent(in) :: use_marine_aggp

ierr=0

call UseMarineAggradation (use_marine_aggp)

return

end subroutine FastScape_Use_Marine_Aggradation

!--------------------------------------------------------------------------

subroutine FastScape_Set_Marine_Aggradation_rate (marine_agg_ratep,ierr)

use FastScapeContext

implicit none

integer, intent(out):: ierr
double precision, intent(in) :: marine_agg_ratep

ierr=0

call SetMarineAggradationRate (marine_agg_ratep)

return

end subroutine FastScape_Set_Marine_Aggradation_rate
