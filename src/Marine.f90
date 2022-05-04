#include "Error.fpp"
subroutine Marine(ierr)

  ! Marine transport component
  ! using silt and sand coupling diffusion solver
  ! developed by Xiaoping Yuan (2017-2018)

  use FastScapeContext

  implicit none

  double precision, dimension(:), allocatable :: flux,shelfdepth,ht,Fs,dh,dh1,dh2,Fmixt,mwater
  double precision, dimension(:), allocatable :: dhs, dhs1, F1, F2, zi, zo
  integer, dimension(:), allocatable :: flag,mmnrec,mmstack
  integer, dimension(:,:), allocatable :: mmrec
  double precision, dimension(:,:), allocatable :: mmwrec,mmlrec
  double precision shelfslope,ratio1,ratio2,dx,dy,f_mass_cons,f_mass_cons_deep,sed_above_sealevel
  integer ij,ijr,ijk,k
  integer, intent(inout):: ierr
  
  double precision, dimension(:), allocatable :: dh_above, dh_below
  !TT make the marine diffusion sequential ------
  character cbc*4
  ! define other parameters
  integer i,j,ipj,imj,ijp,ijm
  double precision, dimension(:), allocatable :: Q1,Q2,hp,fp,ft,hhalf,fhalf,fhalfp,tint
  double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res
  double precision K1,K2,tol,err1,err2
  double precision Ap,Bp,Cp,Dp,Ep,Mp,Np
  !TT -------------------------------------------

  allocate (flux(nn),shelfdepth(nn),ht(nn),Fs(nn),dh(nn),dh1(nn),dh2(nn),Fmixt(nn),flag(nn))
  allocate (dhs(nn),dhs1(nn),F1(nn),F2(nn),zi(nn),zo(nn))

  ! set nodes at transition between ocean and continent
  flag=0

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! computing flux from continental erosion
  flux=0.d0
  where (h.gt.sealevel) flux=Sedflux

  if ( .not. runSPL ) then
    do ij=nn,1,-1
      ijk=stack(ij)
      ijr=rec(ijk)
      if (ijr.ne.ijk.and.h(ijk).gt.sealevel) then
        flux(ijr)=flux(ijr)+flux(ijk)
      endif
    enddo
    ! here the integral of erosion/deposition has been done
    ! and distributed as flux to ocean
    where (h.gt.sealevel) flux=0.d0
  
    ! set nodes at transition between ocean and continent
    !where (flux.gt.tiny(flux)) flag=1
  
    ! decompact volume of pure solid phase (silt and sand) from onshore
    ratio1=ratio/(1.d0-poro1)
    ratio2=(1.d0-ratio)/(1.d0-poro2)
    ! total volume of silt and sand after decompaction
    flux=flux*(ratio1+ratio2)
  
    ! modifications made by Jean for multiple flow to distribute continental flux to ocean on the shelf
    ! Dec 2018
    !write(*,*)'before find_mult_rec'
  
    allocate (mmrec(8,nn),mmnrec(nn),mmwrec(8,nn),mmlrec(8,nn),mmstack(nn),mwater(nn))
  
    call find_mult_rec (h,rec,stack,mwater,mmrec,mmnrec,mmwrec,mmlrec,mmstack,nx,ny,dx,dy,0.d0,p_mfd_exp, &
      bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic, ierr)
  
    !write(*,*)'after find_mult_rec'
    !print*,count(flux>0.and.mmnrec==0),count(flux>0),count(mmstack==0)
  
    ! modifications made by Jean
    ! to compute shelf depth
    shelfdepth=sealevel
    shelfslope=-1.d-4
    do ij=1,nn
      ijk=mmstack(ij)
      do k=1,mmnrec(ijk)
        ijr=mmrec(k,ijk)
        if (h(ijk).lt.sealevel) then
          shelfdepth(ijr)=min(shelfdepth(ijr),shelfdepth(ijk)+mmlrec(k,ijk)*shelfslope)
          shelfdepth(ijr)=max(shelfdepth(ijr),h(ijr))
        endif
      enddo
    enddo
    ! end modifications
  
    ! passes the flux across the shelf
    ! modifications made by Jean
  
    where (h.lt.sealevel) flux=flux+(h-shelfdepth)
    do ij=1,nn
      ijk=mmstack(ij)
      do k=1,mmnrec(ijk)
        ijr=mmrec(k,ijk)
        flux(ijr)=flux(ijr)+max(0.d0,flux(ijk)*mmwrec(k,ijk))
      enddo
    enddo
    ! modifications made by Jean
  
    deallocate (mmrec,mmnrec,mmwrec,mmlrec,mmstack,mwater)
  
    where (flux.gt.0.d0.and.h.lt.sealevel) flux=-(h-shelfdepth)
    where (flux.le.0.d0.and.h.lt.sealevel) flux=flux-(h-shelfdepth)
    where (h.ge.sealevel) flux=0.d0
    flux=max(flux,0.d0)
  
  endif

  ! silt fraction (after decompaction) in shelf
  Fs=0.d0
  where (flux.gt.0.d0) Fs=ratio1/(ratio1+ratio2)

  ! scales flux by time step
  flux=flux/dt

  ! store flux at shoreline in additional array for visualisation purposes,
  ! and to compute timestepping criterion based on the maximum sediment flux at the shoreline.
  sedflux_shore = flux

  ! stores initial height and fraction
  ht=h
  Fmixt=Fmix

  !print*,'flux',minval(flux),sum(flux)/nx/ny,maxval(flux)
  !print*,'Fmix',minval(Fmix),sum(Fmix)/nx/ny,maxval(Fmix)

  !write(*,*)'before silt and sand coupling diffusion in ocean'

!  ! silt and sand coupling diffusion in ocean
!  call SiltSandCouplingDiffusion (h,Fmix,flux*Fs,flux*(1.d0-Fs), &
!  nx,ny,dx,dy,dt,sealevel,layer,kdsea1,kdsea2,nGSMarine,flag,bounds_ibc,ierr);FSCAPE_CHKERR(ierr)

!TT ========================================================================================================
! sequential SiltSandCouplingDiffusion 

!  if (any(h<sealevel)) then

  K1 = kdsea1
  K2 = kdsea2
  write (cbc,'(i4)') bounds_ibc
  allocate (hp(nn),fp(nn),ft(nn),hhalf(nn),fhalf(nn),fhalfp(nn),tint(nn))
  allocate (Q1(nn), Q2(nn))
  Q1 = flux*Fs
  Q2 = flux*(1.d0-Fs)

  ! initilize the elevation and silt fraction at time t and t+dt/2
  ft=Fmix
  hhalf=h
  fhalf=Fmix

  ! tolerance is in m
  tol=1.d0
  err1=2*tol
  err2=2*tol
  nGSMarine=0

  ! iteration until convergence is reached
  do while (err1.gt.tol)
    ! update the elevation and silt fraction during each iteration
    hp=h
    fp=Fmix
    fhalfp=fhalf
    nGSMarine=nGSMarine+1
    ! calculate the elevation h in x-direction
    allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))
    do j=2,ny-1
      do i=1,nx
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          if (i.eq.1) then
            if (cbc(4:4).eq.'1') then
              diag(i)=1.d0
              sup(i)=0.d0
              rhs(i)=ht(ij)
            else
              Ap=dt/2.d0*(K2+(K1-K2)*(fhalfp(ipj)+fhalfp(ij))/2.d0)/dx**2
              diag(i)=1.d0+Ap
              sup(i)=-Ap
              Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
              -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(i)=Cp+ht(ij)
            endif
          elseif (i.eq.nx) then
            if (cbc(2:2).eq.'1') then
              diag(i)=1.d0
              inf(i)=0.d0
              rhs(i)=ht(ij)
            else
              Bp=-dt/2.d0*(K2+(K1-K2)*(fhalfp(ij)+fhalfp(imj))/2.d0)/dx**2
              diag(i)=1.d0-Bp
              inf(i)=Bp
              Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
              -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(i)=Cp+ht(ij)
            endif
          else
            Ap=dt/2.d0*(K2+(K1-K2)*(fhalfp(ipj)+fhalfp(ij))/2.d0)/dx**2
            Bp=-dt/2.d0*(K2+(K1-K2)*(fhalfp(ij)+fhalfp(imj))/2.d0)/dx**2
            diag(i)=1.d0+Ap-Bp
            sup(i)=-Ap
            inf(i)=Bp
            Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
            -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
            +(Q1(ij)+Q2(ij))*dt/2.d0
            rhs(i)=Cp+ht(ij)
          endif
          ! in continent
        else
          diag(i)=1.d0
          sup(i)=0.d0
          inf(i)=0.d0
          rhs(i)=ht(ij)
        endif
      enddo
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,nx,ierr);FSCAPE_CHKERR(ierr)
      do i=1,nx
        ij=(j-1)*nx+i
        hhalf(ij)=res(i)
      enddo
    enddo
    tint=hhalf
    ! the corner nodes (1,1) and (1,ny)
    hhalf(1)=hhalf(2)
    hhalf((ny-1)*nx+1)=hhalf((ny-1)*nx+2)
    ! the corner nodes (nx,1) and (nx,ny)
    hhalf(nx)=hhalf(nx-1)
    hhalf(nx*ny)=hhalf(nx*ny-1)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the silt fraction F in x-direction
    allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))
    do j=2,ny-1
      do i=2,nx-1
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          ! deposition
          if (hhalf(ij).ge.(1.d0+1.d-6)*ht(ij)) then
            Dp=(hhalf(ij)-ht(ij))/dt
            Ep=K1/2.d0*(hhalf(ipj)-hhalf(ij))/dx**2
            Mp=-K1/2.d0*(hhalf(ij)-hhalf(imj))/dx**2
            Np=K1/2.d0*(ft(ijp)+ft(ij))*(ht(ijp)-ht(ij))/dy**2 &
            -K1/2.d0*(ft(ij)+ft(ijm))*(ht(ij)-ht(ijm))/dy**2 &
            +Q1(ij)
            diag(i)=2.d0*layer/dt+Dp-Mp-Ep
            sup(i)=-Ep
            inf(i)=-Mp
            rhs(i)=Np-Dp*ft(ij)+2.d0*layer*ft(ij)/dt
            ! erosion
          else
            diag(i)=1.d0
            sup(i)=0.d0
            inf(i)=0.d0
            rhs(i)=ft(ij)
          endif
          ! in continent
        else
          diag(i)=1.d0
          sup(i)=0.d0
          inf(i)=0.d0
          rhs(i)=ft(ij)
        endif
      enddo
      ! bc on i=1
      diag(1)=1.d0
      sup(1)=-1.d0
      rhs(1)=0.d0
      ! bc on i=nx
      diag(nx)=1.d0
      inf(nx)=-1.d0
      rhs(nx)=0.d0
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,nx,ierr);FSCAPE_CHKERR(ierr)
      do i=1,nx
        ij=(j-1)*nx+i
        fhalf(ij)=res(i)
      enddo
    enddo
    fhalf=max(0.d0,fhalf)
    fhalf=min(1.d0,fhalf)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the elevation h in y-direction
    allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))
    do i=2,nx-1
      do j=1,ny
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          if (j.eq.1) then
            if (cbc(1:1).eq.'1') then
              diag(j)=1.d0
              sup(j)=0.d0
              rhs(j)=hhalf(ij)
            else
              Ap=dt/2.d0*(K2+(K1-K2)*(fp(ijp)+fp(ij))/2.d0)/dy**2
              diag(j)=1.d0+Ap
              sup(j)=-Ap
              Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
              -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(j)=Cp+hhalf(ij)
            endif
          elseif (j.eq.ny) then
            if (cbc(3:3).eq.'1') then
              diag(j)=1.d0
              inf(j)=0.d0
              rhs(j)=hhalf(ij)
            else
              Bp=-dt/2.d0*(K2+(K1-K2)*(fp(ij)+fp(ijm))/2.d0)/dy**2
              diag(j)=1.d0-Bp
              inf(j)=Bp
              Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
              -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(j)=Cp+hhalf(ij)
            endif
          else
            Ap=dt/2.d0*(K2+(K1-K2)*(fp(ijp)+fp(ij))/2.d0)/dy**2
            Bp=-dt/2.d0*(K2+(K1-K2)*(fp(ij)+fp(ijm))/2.d0)/dy**2
            diag(j)=1.d0+Ap-Bp
            sup(j)=-Ap
            inf(j)=Bp
            Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
            -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
            +(Q1(ij)+Q2(ij))*dt/2.d0
            rhs(j)=Cp+hhalf(ij)
          endif
          ! in continent
        else
          diag(j)=1.d0
          sup(j)=0.d0
          inf(j)=0.d0
          rhs(j)=hhalf(ij)
        endif
      enddo
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,ny,ierr);FSCAPE_CHKERR(ierr)
      do j=1,ny
        ij=(j-1)*nx+i
        tint(ij)=res(j)
      enddo
    enddo
    h=tint
    ! the corner nodes (1,1) and (1,ny)
    h(1)=h(2)
    h((ny-1)*nx+1)=h((ny-1)*nx+2)
    ! the corner nodes (nx,1) and (nx,ny)
    h(nx)=h(nx-1)
    h(nx*ny)=h(nx*ny-1)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the silt fraction F in y-direction
    allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))
    do i=2,nx-1
      do j=2,ny-1
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          ! deposition
          if (h(ij).ge.(1.d0+1.d-6)*hhalf(ij)) then
            Dp=(h(ij)-hhalf(ij))/dt
            Ep=K1/2.d0*(h(ijp)-h(ij))/dy**2
            Mp=-K1/2.d0*(h(ij)-h(ijm))/dy**2
            Np=K1/2.d0*(fhalf(ipj)+fhalf(ij))*(hhalf(ipj)-hhalf(ij))/dx**2 &
            -K1/2.d0*(fhalf(ij)+fhalf(imj))*(hhalf(ij)-hhalf(imj))/dx**2 &
            +Q1(ij)
            diag(j)=2.d0*layer/dt+Dp-Mp-Ep
            sup(j)=-Ep
            inf(j)=-Mp
            rhs(j)=Np-Dp*fhalf(ij)+2.d0*layer*fhalf(ij)/dt
            ! erosion
          else
            diag(j)=1.d0
            sup(j)=0.d0
            inf(j)=0.d0
            rhs(j)=fhalf(ij)
          endif
          ! in continent
        else
          diag(j)=1.d0
          sup(j)=0.d0
          inf(j)=0.d0
          rhs(j)=fhalf(ij)
        endif
      enddo
      ! bc on j=1
      diag(1)=1.d0
      sup(1)=-1.d0
      rhs(1)=0.d0
      ! bc on j=ny
      diag(ny)=1.d0
      inf(ny)=-1.d0
      rhs(ny)=0.d0
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,ny,ierr);FSCAPE_CHKERR(ierr)
      do j=1,ny
        ij=(j-1)*nx+i
        Fmix(ij)=res(j)
      enddo
    enddo
    Fmix=max(0.d0,Fmix)
    Fmix=min(1.d0,Fmix)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the errors in each iteration
    err1=maxval(abs(h-hp))
    err2=maxval(abs(h-hp)/(1.d0+abs(h)))

!    print*,'nGSMarine',nGSMarine,minval(h-hp),sum(h-hp)/nn,maxval(h-hp),err1

    if (nGSMarine.gt.1000) then
      print*,'Marine error: Multi-lithology diffusion not converging; decrease time step'
      !FSCAPE_RAISE_MESSAGE('Marine error: Multi-lithology diffusion not converging; decrease time step',ERR_NotConverged,ierr)
      !FSCAPE_CHKERR(ierr)
      write(*,*) 'WARNING: dt_spm too high - cannot converge marine diffusion err > tol (m)', err1, tol
      tol = err1 + 1.d0
    endif

    ! end of iteration
  enddo

  deallocate (hp,fp,ft,hhalf,fhalf,fhalfp,tint)
  deallocate (Q1,Q2)

!  endif ! seas exist (h<sealevel)

  ! set the silt fraction for continent
  where (h.ge.sealevel+1.d-3) Fmix=0.d0


!TT ========================================================================================================


  ! pure silt and sand during deposition/erosion
  dh1=((h-ht)*Fmix+layer*(Fmix-Fmixt))*(1.d0-poro1)
  dh2=((h-ht)*(1.d0-Fmix)+layer*(Fmixt-Fmix))*(1.d0-poro2)
  dh=dh1+dh2

  if (enforce_marine_no_erosion .and. any(dh<0.d0)) then
    !debug
    !allocate(dh_above(nn))
    !dh_above = dh
    where (dh.lt.0.d0 .and. (ht+dh).lt.sealevel-50) dh = 0.d0
    where (dh1.lt.0.d0 .and. (ht+dh).lt.sealevel-50) dh1 = 0.d0
    where (dh2.lt.0.d0 .and. (ht+dh).lt.sealevel-50) dh2 = 0.d0
    !totaleroded_check = sum(dh-dh_above)
    !write(*,'(a,F13.4)') 'totaleroded_check    = ',totaleroded_check
    !deallocate(dh_above)
    h = ht + dh
    dh1=((h-ht)*Fmix+layer*(Fmix-Fmixt))*(1.d0-poro1)
    dh2=((h-ht)*(1.d0-Fmix)+layer*(Fmixt-Fmix))*(1.d0-poro2)
  endif

  if (enforce_marine_sed_below_sealevel) then
    !TT Should not have sediments above sealevel after diffusion
    allocate(dh_above(nn))
    dh_above = dh
    where ( ht.le.sealevel .and. (ht+dh) .gt. sealevel ) dh = sealevel - ht
    sed_above_sealevel = sum(dh_above-dh)
    write(*,'(a,F13.4,a)') '  --> correction of sediments above sealevel, volume ~ ',sed_above_sealevel,' m'
    deallocate(dh_above)
    h = ht + dh
    dh1=((h-ht)*Fmix+layer*(Fmix-Fmixt))*(1.d0-poro1)
    dh2=((h-ht)*(1.d0-Fmix)+layer*(Fmixt-Fmix))*(1.d0-poro2)
  endif

  ! compute whether Marine diffusion was mass conserving. This might not be the case in basins
  ! with a lot of "underwater" topgraphy and little sediment input
  ! printing is only done if option enforce_marine_mass_cons is used.
  if (enforce_marine_mass_cons) then
     f_mass_cons = sum(dh)/(sum(flux/(ratio1+ratio2))*dt)
     if ((sum(flux)*dt) > 0.d0 .and. f_mass_cons > 1.0d0) then
!         write(*,*)'Too much marine deposits with a factor f_mass_cons compared to sed flux = ',f_mass_cons
         allocate(dh_above(nn),dh_below(nn))
         dh_above = dh
         dh_below = dh
         !we mostly remove sediments deeper then 50 m
         where (h.le.(sealevel-50)) dh_above=0.d0
         where (h.gt.(sealevel-50)) dh_below=0.d0
         f_mass_cons_deep = sum(dh_below) / ((sum(dh)/f_mass_cons) - sum(dh_above))
         if (f_mass_cons_deep>0) then
           where (h.le.(sealevel-50)) dh=dh/f_mass_cons_deep
           write(*,'(a,F14.3,a)')'Scaled marine deposits down according to f_mass_cons = ',f_mass_cons_deep, ' for sediments deeper than 50 m'
         else
           write(*,'(a,F14.3)') '50 m f_mass_cons_deep is < 0 ',f_mass_cons_deep
           dh_above = dh
           dh_below = dh
           !we mostly remove sediments deeper then 25 m
           where (h.le.(sealevel-25)) dh_above=0.d0
           where (h.gt.(sealevel-25)) dh_below=0.d0
           f_mass_cons_deep = sum(dh_below) / ((sum(dh)/f_mass_cons) - sum(dh_above))
           if (f_mass_cons_deep>0) then
             where (h.le.(sealevel-25)) dh=dh/f_mass_cons_deep
             write(*,'(a,F13.4,a)')'Scaled marine deposits down according to f_mass_cons = ',f_mass_cons_deep, ' for sediments deeper than 25 m'
           else
             write(*,'(a,F14.3)') '25 m f_mass_cons_deep is < 0 ',f_mass_cons_deep
             dh = dh/f_mass_cons
             write(*,'(a,F14.3)')'Scaled marine deposits down according to f_mass_cons = ',f_mass_cons
           endif
         endif
         deallocate(dh_above,dh_below)
         h = ht + dh
         dh1=((h-ht)*Fmix+layer*(Fmix-Fmixt))*(1.d0-poro1)
         dh2=((h-ht)*(1.d0-Fmix)+layer*(Fmixt-Fmix))*(1.d0-poro2)
         dh=dh1+dh2
         f_mass_cons = sum(dh)/(sum(flux/(ratio1+ratio2))*dt)
         write(*,'(a,F7.3)')'After correction f_mass_cons =',f_mass_cons
     else if ((sum(flux)*dt) > 0.d0 .and. f_mass_cons < 1.d0) then
         write(*,'(a,F14.3)')'Mass lost during marine diffusion: f_mass_cons = ',f_mass_cons
     endif
  else
    f_mass_cons = sum(dh)/(sum(flux/(ratio1+ratio2))*dt)
    if ((sum(flux)*dt) > 0.d0 .and. f_mass_cons > 1.0d0) then
      write(*,'(a,F14.3)') 'Mass gain during marine diffusion: f_mass_cons = ',f_mass_cons
    else if ((sum(flux)*dt) > 0.d0 .and. f_mass_cons < 1.d0) then
      write(*,'(a,F14.3)') 'Mass lost during marine diffusion: f_mass_cons = ',f_mass_cons
    endif
  endif

  ! store deposited material in dh_dep to be able to transfer this value to coupled thermo-mechanical code
  ! where compaction can be applied if necessary.
  dh_dep = dh

  ! >>>>>>>> compaction starts added by Jean (Dec 2018)

  ! sum of pure silt and solid phase
  if (step.eq.0) then
    dhs1=dh1
    dhs=dh
  else
    dhs1=dhs1+dh1
    dhs=dhs+dh
  endif
  where (dhs1.lt.0.d0) dhs1=0.d0
  where (dhs.lt.0.d0) dhs=0.d0

  ! calculate the average silt (and sand) fraction in ocean part
  F1=0.d0;F2=0.d0
  where (h.le.sealevel.and.dhs.gt.0.d0) F1=dhs1/dhs
  F1=max(0.d0,F1);F1=min(1.d0,F1)
  where (h.le.sealevel.and.dhs.gt.0.d0) F2=1.d0-F1

  ! calculate the thickness after compaction, initial thickness of sediments
  !zi=ht-b
  !call compaction (F1,F2,poro1,poro2,zporo1,zporo2,nn,dh,zi,zo)
  ! update the elevation
  !h=b+zo

  ! >>>>>>>> compaction ends

  ! update the elevation
  etot=etot+ht-h
  erate=erate+(ht-h)/dt
  where (h.lt.sealevel) Sedflux=0.d0
  where (h.lt.sealevel) etot=0.d0
  where (h.lt.sealevel) erate=0.d0

  ! set the silt fraction in continent
  where (h.ge.sealevel+1.d-3) Fmix=0.d-1
  Fmix=max(0.d0,Fmix)
  Fmix=min(1.d0,Fmix)

  ! updates basement
  b=min(h,b)

  ! update rock_type
  where (dh > 0.0d0 .and. h <= sealevel ) rock_type = 3 ! marine sed.

  deallocate (flux,shelfdepth,ht,Fs,dh,dh1,dh2,Fmixt)

  return

end subroutine Marine

!----------------------------------------------------------------------------------

subroutine SiltSandCouplingDiffusion (h,f,Q1,Q2,nx,ny,dx,dy,dt, &
  sealevel,L,kdsea1,kdsea2,niter,flag,ibc,ierr)

  use FastScapeErrorCodes

  implicit none

  ! arguments
  double precision,dimension(:),intent(inout) :: h, f
  double precision,dimension(:),intent(in)    :: Q1, Q2
  integer nx,ny
  double precision dx,dy,dt,sealevel,L,kdsea1,kdsea2
  integer niter
  integer,dimension(:),intent(in) :: flag
  integer ibc
  integer, intent(inout) :: ierr

  ! define other parameters
  integer i,j,ij,ipj,imj,ijp,ijm,nn
  double precision, dimension(:), allocatable :: hp,fp,ht,ft,hhalf,fhalf,fhalfp,tint
  double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res
  double precision K1,K2,tol,err1,err2
  double precision Ap,Bp,Cp,Dp,Ep,Mp,Np
  character cbc*4

  write(*,*) 'SiltSandCouplingDiffusion before checking seas exist'

  if (any(h<sealevel)) then  

  write(*,*) 'SiltSandCouplingDiffusion there are seas, sea level = ',sealevel

  write (cbc,'(i4)') ibc

  K1=kdsea1
  K2=kdsea2

  nn=nx*ny

  write(*,*) 'before allocation'

  allocate (hp(nn),fp(nn),ht(nn),ft(nn),hhalf(nn),fhalf(nn),fhalfp(nn),tint(nn))

  write(*,*) 'after alllocation'

  ! initilize the elevation and silt fraction at time t and t+dt/2
  ht=h
  ft=f
  hhalf=h
  fhalf=f

  ! tolerance is in m
  tol=1.d0
  err1=2*tol
  err2=2*tol
  niter=0

  ! iteration until convergence is reached
  do while (err1.gt.tol)
    ! update the elevation and silt fraction during each iteration
    hp=h
    fp=f
    fhalfp=fhalf
    niter=niter+1
    ! calculate the elevation h in x-direction
    allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))
    do j=2,ny-1
      do i=1,nx
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          if (i.eq.1) then
            if (cbc(4:4).eq.'1') then
              diag(i)=1.d0
              sup(i)=0.d0
              rhs(i)=ht(ij)
            else
              Ap=dt/2.d0*(K2+(K1-K2)*(fhalfp(ipj)+fhalfp(ij))/2.d0)/dx**2
              diag(i)=1.d0+Ap
              sup(i)=-Ap
              Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
              -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(i)=Cp+ht(ij)
            endif
          elseif (i.eq.nx) then
            if (cbc(2:2).eq.'1') then
              diag(i)=1.d0
              inf(i)=0.d0
              rhs(i)=ht(ij)
            else
              Bp=-dt/2.d0*(K2+(K1-K2)*(fhalfp(ij)+fhalfp(imj))/2.d0)/dx**2
              diag(i)=1.d0-Bp
              inf(i)=Bp
              Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
              -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(i)=Cp+ht(ij)
            endif
          else
            Ap=dt/2.d0*(K2+(K1-K2)*(fhalfp(ipj)+fhalfp(ij))/2.d0)/dx**2
            Bp=-dt/2.d0*(K2+(K1-K2)*(fhalfp(ij)+fhalfp(imj))/2.d0)/dx**2
            diag(i)=1.d0+Ap-Bp
            sup(i)=-Ap
            inf(i)=Bp
            Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
            -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
            +(Q1(ij)+Q2(ij))*dt/2.d0
            rhs(i)=Cp+ht(ij)
          endif
          ! in continent
        else
          diag(i)=1.d0
          sup(i)=0.d0
          inf(i)=0.d0
          rhs(i)=ht(ij)
        endif
      enddo
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,nx,ierr);FSCAPE_CHKERR(ierr)
      do i=1,nx
        ij=(j-1)*nx+i
        hhalf(ij)=res(i)
      enddo
    enddo
    tint=hhalf
    ! the corner nodes (1,1) and (1,ny)
    hhalf(1)=hhalf(2)
    hhalf((ny-1)*nx+1)=hhalf((ny-1)*nx+2)
    ! the corner nodes (nx,1) and (nx,ny)
    hhalf(nx)=hhalf(nx-1)
    hhalf(nx*ny)=hhalf(nx*ny-1)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the silt fraction F in x-direction
    allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))
    do j=2,ny-1
      do i=2,nx-1
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          ! deposition
          if (hhalf(ij).ge.(1.d0+1.d-6)*ht(ij)) then
            Dp=(hhalf(ij)-ht(ij))/dt
            Ep=K1/2.d0*(hhalf(ipj)-hhalf(ij))/dx**2
            Mp=-K1/2.d0*(hhalf(ij)-hhalf(imj))/dx**2
            Np=K1/2.d0*(ft(ijp)+ft(ij))*(ht(ijp)-ht(ij))/dy**2 &
            -K1/2.d0*(ft(ij)+ft(ijm))*(ht(ij)-ht(ijm))/dy**2 &
            +Q1(ij)
            diag(i)=2.d0*L/dt+Dp-Mp-Ep
            sup(i)=-Ep
            inf(i)=-Mp
            rhs(i)=Np-Dp*ft(ij)+2.d0*L*ft(ij)/dt
            ! erosion
          else
            diag(i)=1.d0
            sup(i)=0.d0
            inf(i)=0.d0
            rhs(i)=ft(ij)
          endif
          ! in continent
        else
          diag(i)=1.d0
          sup(i)=0.d0
          inf(i)=0.d0
          rhs(i)=ft(ij)
        endif
      enddo
      ! bc on i=1
      diag(1)=1.d0
      sup(1)=-1.d0
      rhs(1)=0.d0
      ! bc on i=nx
      diag(nx)=1.d0
      inf(nx)=-1.d0
      rhs(nx)=0.d0
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,nx,ierr);FSCAPE_CHKERR(ierr)
      do i=1,nx
        ij=(j-1)*nx+i
        fhalf(ij)=res(i)
      enddo
    enddo
    fhalf=max(0.d0,fhalf)
    fhalf=min(1.d0,fhalf)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the elevation h in y-direction
    allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))
    do i=2,nx-1
      do j=1,ny
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          if (j.eq.1) then
            if (cbc(1:1).eq.'1') then
              diag(j)=1.d0
              sup(j)=0.d0
              rhs(j)=hhalf(ij)
            else
              Ap=dt/2.d0*(K2+(K1-K2)*(fp(ijp)+fp(ij))/2.d0)/dy**2
              diag(j)=1.d0+Ap
              sup(j)=-Ap
              Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
              -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(j)=Cp+hhalf(ij)
            endif
          elseif (j.eq.ny) then
            if (cbc(3:3).eq.'1') then
              diag(j)=1.d0
              inf(j)=0.d0
              rhs(j)=hhalf(ij)
            else
              Bp=-dt/2.d0*(K2+(K1-K2)*(fp(ij)+fp(ijm))/2.d0)/dy**2
              diag(j)=1.d0-Bp
              inf(j)=Bp
              Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
              -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(j)=Cp+hhalf(ij)
            endif
          else
            Ap=dt/2.d0*(K2+(K1-K2)*(fp(ijp)+fp(ij))/2.d0)/dy**2
            Bp=-dt/2.d0*(K2+(K1-K2)*(fp(ij)+fp(ijm))/2.d0)/dy**2
            diag(j)=1.d0+Ap-Bp
            sup(j)=-Ap
            inf(j)=Bp
            Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
            -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
            +(Q1(ij)+Q2(ij))*dt/2.d0
            rhs(j)=Cp+hhalf(ij)
          endif
          ! in continent
        else
          diag(j)=1.d0
          sup(j)=0.d0
          inf(j)=0.d0
          rhs(j)=hhalf(ij)
        endif
      enddo
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,ny,ierr);FSCAPE_CHKERR(ierr)
      do j=1,ny
        ij=(j-1)*nx+i
        tint(ij)=res(j)
      enddo
    enddo
    h=tint
    ! the corner nodes (1,1) and (1,ny)
    h(1)=h(2)
    h((ny-1)*nx+1)=h((ny-1)*nx+2)
    ! the corner nodes (nx,1) and (nx,ny)
    h(nx)=h(nx-1)
    h(nx*ny)=h(nx*ny-1)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the silt fraction F in y-direction
    allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))
    do i=2,nx-1
      do j=2,ny-1
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          ! deposition
          if (h(ij).ge.(1.d0+1.d-6)*hhalf(ij)) then
            Dp=(h(ij)-hhalf(ij))/dt
            Ep=K1/2.d0*(h(ijp)-h(ij))/dy**2
            Mp=-K1/2.d0*(h(ij)-h(ijm))/dy**2
            Np=K1/2.d0*(fhalf(ipj)+fhalf(ij))*(hhalf(ipj)-hhalf(ij))/dx**2 &
            -K1/2.d0*(fhalf(ij)+fhalf(imj))*(hhalf(ij)-hhalf(imj))/dx**2 &
            +Q1(ij)
            diag(j)=2.d0*L/dt+Dp-Mp-Ep
            sup(j)=-Ep
            inf(j)=-Mp
            rhs(j)=Np-Dp*fhalf(ij)+2.d0*L*fhalf(ij)/dt
            ! erosion
          else
            diag(j)=1.d0
            sup(j)=0.d0
            inf(j)=0.d0
            rhs(j)=fhalf(ij)
          endif
          ! in continent
        else
          diag(j)=1.d0
          sup(j)=0.d0
          inf(j)=0.d0
          rhs(j)=fhalf(ij)
        endif
      enddo
      ! bc on j=1
      diag(1)=1.d0
      sup(1)=-1.d0
      rhs(1)=0.d0
      ! bc on j=ny
      diag(ny)=1.d0
      inf(ny)=-1.d0
      rhs(ny)=0.d0
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,ny,ierr);FSCAPE_CHKERR(ierr)
      do j=1,ny
        ij=(j-1)*nx+i
        f(ij)=res(j)
      enddo
    enddo
    f=max(0.d0,f)
    f=min(1.d0,f)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the errors in each iteration
    err1=maxval(abs(h-hp))
    err2=maxval(abs(h-hp)/(1.d0+abs(h)))

    !print*,'niter',niter,minval(h-hp),sum(h-hp)/nn,maxval(h-hp),err1

    if (niter.gt.1000) then
      FSCAPE_RAISE_MESSAGE('Marine error: Multi-lithology diffusion not converging; decrease time step',ERR_NotConverged,ierr)
      FSCAPE_CHKERR(ierr)
    endif

    ! end of iteration
  enddo

  deallocate (hp,fp,ht,ft,hhalf,fhalf,fhalfp,tint)

  endif ! seas exist (h<sealevel)

  ! set the silt fraction for continent
  where (h.ge.sealevel+1.d-3) f=0.d0

  ! end of the subroutine
end subroutine SiltSandCouplingDiffusion

!-----------------------------------------------------

subroutine compaction (F1,F2,poro1,poro2,z1,z2,nn,dh,zi,zo)

  ! Newton iteration to calculate the thickness after compaction

  implicit none

  integer k,nn
  double precision poro1,poro2,z1,z2,fx,dfx
  double precision F1(nn),F2(nn),dh(nn),zi(nn),zo(nn)

  ! initial guess on zo
  zo=zi
  ! iteration process
  do k=1,nn
    1000 continue
    fx=zo(k)-zi(k)+F1(k)*poro1*z1*(exp(-zo(k)/z1)-exp(-zi(k)/z1)) &
    +F2(k)*poro2*z2*(exp(-zo(k)/z2)-exp(-zi(k)/z2))-dh(k)
    dfx=1.d0-F1(k)*poro1*exp(-zo(k)/z1)-F2(k)*poro2*exp(-zo(k)/z2)
    zo(k)=zo(k)-fx/dfx
    if (abs(fx/dfx).gt.1.d-6) goto 1000
  enddo

  return

end subroutine compaction

!-----------------------------------------------------
! added computation of timestep for marine deposition
! this functionality allows to run a loop on SiltSandCouplingDiffusion,
! so that the timestep can be maximized. This function is not used atm.
subroutine compute_dt_marine(flux,dt,dt_crit_marine,dt_marine,nstep_marine,layer,nn)

implicit none
!==============================================================================!

! arguments
double precision, intent(in) :: dt_crit_marine, dt, layer
integer, intent(in) :: nn
double precision,dimension(nn),intent(in) :: flux
integer, intent(out) :: nstep_marine
double precision, intent(out) :: dt_marine
double precision :: MaxSedflux, dt_max

MaxSedflux = maxval(flux)
if (MaxSedflux <= 0.d0) then
    dt_max = dt
else
    dt_max = dt_crit_marine * layer / MaxSedflux
endif

if (dt_max >= dt) then
    dt_marine = dt
    nstep_marine = 1
else
    nstep_marine = ceiling(dt/dt_max)
    dt_marine    = dt/dble(nstep_marine)
endif
write(*,*) 'Fastscape marine timestepping criterion used.'
write(*,*) 'timestepping factor = ',dt_crit_marine,'computed dt used = ', dt_marine,'times ', nstep_marine

end subroutine compute_dt_marine


!-----------------------------------------------------
! Mass conserving aggradation of marine domain with sediments arriving at the shoreline.
subroutine MarineAggradation(ierr)

use FastScapeContext

implicit none
!==============================================================================!
double precision, dimension(:), allocatable :: flux,shelfdepth,ht,Fs,dh,Fmixt
double precision ratio1,ratio2,dx,dy,fill_level
integer ij,ijr,ijk,i,iter_counter
double precision bottom,top,error,tol,vol,vol_tester
integer, intent(inout):: ierr

allocate (flux(nn),shelfdepth(nn),ht(nn),Fs(nn),dh(nn),Fmixt(nn))

! arguments
dx=xl/(nx-1)
dy=yl/(ny-1)

! computing flux from continental erosion
flux=0.d0
where (h.gt.sealevel) flux=Sedflux
do ij=nn,1,-1
  ijk=stack(ij)
  ijr=rec(ijk)
  if (ijr.ne.ijk.and.h(ijk).gt.sealevel) then
    flux(ijr)=flux(ijr)+flux(ijk)
  endif
enddo
! here the integral of erosion/deposition has been done
! and distributed as flux to ocean
where (h.gt.sealevel) flux=0.d0

! set nodes at transition between ocean and continent
!where (flux.gt.tiny(flux)) flag=1

! decompact volume of pure solid phase (silt and sand) from onshore
ratio1=ratio/(1.d0-poro1)
ratio2=(1.d0-ratio)/(1.d0-poro2)
! total volume of silt and sand after decompaction
flux=flux*(ratio1+ratio2)

! silt fraction (after decompaction) in shelf
Fs=0.d0
where (flux.gt.0.d0) Fs=ratio1/(ratio1+ratio2)

! scales flux by time step
flux=flux/dt

! store flux at shoreline in additional array for visualisation purposes,
! and to compute timestepping criterion based on the maximum sediment flux at the shoreline.
sedflux_shore = flux
ht  = h
if (marine_aggradation_rate < 0.d0) then
  if (sum(sedflux_shore) > 0.d0) then

    vol           = sum(sedflux_shore*dt)*dx*dy
    tol           = 1e-4*dx*dy      ! 0.1 mm on average
    fill_level    = sealevel        ! initialise fill_level to be sealevel
    top           = maxval(h)       ! initially top and bottom are the model max and min
    bottom        = minval(h)
    iter_counter  = 0
    error         = 2*tol

    while_loop: do while (abs(error) >= abs(tol))
      iter_counter = iter_counter + 1
      ht  = h

      ! Fill the tester to fill_level
      do i=1,nn
        if (ht(i) <= fill_level) then
          ht(i)  = fill_level
        end if
      end do

      vol_tester  = sum(ht-h)*dx*dy
      error       = vol_tester - vol

      ! testing error and adjusting fill level
      if (error >= 0.d0) then ! fill level was too high
        top = fill_level
        bottom = bottom
        fill_level = (top+bottom)/2
      else ! fill level was too low
        top = top
        bottom = fill_level
        fill_level = (top+bottom)/2
      end if

      ! save guard to prevent filling above sealevel
      if (fill_level > sealevel) then
        fill_level = sealevel
        write(*,*)'Basin filled up to sealevel. Stopping aggradation at sealevel.'
        exit while_loop
      endif

      if (iter_counter >=10000) then
        write(*,'(a)')'------ FastScape mass conserving marine aggradation failed ------'
        write(*,'(a,es15.4,a,es15.4)')'average error [mm]:',abs(error)/(dx*dy),' average tol [mm]:',abs(tol)/(dx*dy)
        write(*,'(a,i6,a)')'sed routine needed',iter_counter,' iterations'
        write(*,'(a,f8.4)')'percent of deposited sed: ',vol_tester/vol * 100
        write(*,'(a,f13.3,a,f13.3)')'fill_level = ',fill_level,', sealevel = ',sealevel
        write(*,'(a,f8.4)')'sum sedflux: ',sum(sedflux_shore)
        FSCAPE_RAISE_MESSAGE('Marine error: Mass conserving filling with available sediment not converged',ERR_NotConverged,ierr)
        FSCAPE_CHKERR(ierr)
      end if

    end do while_loop

    write(*,'(a)')'------ FastScape mass conserving marine aggradation used ------'
    write(*,'(a,es15.4,a,es15.4)')'average error [mm]:',abs(error)/(dx*dy),' average tol [mm]:',abs(tol)/(dx*dy)
    write(*,'(a,i6,a)')'sed routine needed',iter_counter,' iterations'
    write(*,'(a,f8.4)')'percent of deposited sed: ',vol_tester/vol * 100
    write(*,'(a,f13.3,a,f13.3)')'fill_level = ',fill_level,', sealevel = ',sealevel
  end if
else if (marine_aggradation_rate > 0.d0) then
  ! Filling with constant aggradation rate up to sealevel
  do i=1,nn
    if (ht(i) <= sealevel) then
      ht(i)  = ht(i) + marine_aggradation_rate*dt
      if (ht(i) > sealevel) ht(i) = sealevel
    end if
  end do
  write(*,'(a)')'------ FastScape marine filling with constant aggradation rate used ------'
  write(*,'(a,f13.6)')'Aggradation rate = ',marine_aggradation_rate
end if

h = ht

dh = ht-h

dh_dep = dh

etot=etot+ht-h
erate=erate+(ht-h)/dt

! updates basement
b=min(h,b)

! update rock_type
where (dh > 0.0d0 .and. h <= sealevel ) rock_type = 3 ! marine sed.

deallocate (flux,shelfdepth,ht,Fs,dh,Fmixt)

end subroutine MarineAggradation
