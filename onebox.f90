!#########################################
!           Matthew Copper
! One box model of io plasma torus
! Based on IDL program by Andrew Steffl 
! Started on 2/25/2013
!#########################################

PROGRAM Onebox

  USE DEBUG
  USE FUNCTIONS
  USE TIMESTEP 
  USE ReadEmis
  USE INPUTS
  USE FTMIX
  USE PARALLELVARIABLES
  USE OUTPUTS
  USE DIFFUSION
  USE MPI
  USE DIMENSIONS

  IMPLICIT NONE

  character(len=8)    ::x1
  real                ::thing1
  double precision    ::thing2  

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
  
  lnggrid=mod(mype,LNG_GRID)+1
  radgrid=(mype/LNG_GRID)+1

if(.true.) then
  write (x1, '(I3.3)') mype !write the integer 'mygrid' to non-existent file
 
  num_char=trim(x1)  !trim non-existent file and store as num_char
  !num_char can be used to name output files by grid space using "output"//num_char//".dat" 

  if ( npes .ne. LNG_GRID*RAD_GRID ) then
    print *, "The current version only supports ", LNG_GRID*RAD_GRID, " processors."   
  else 
   call model()
  endif
endif
call MPI_FINALIZE(ierr)

CONTAINS 

subroutine model()
  integer             ::nit
  real                ::lontemp, day, source_tot, trans_tot
  real                ::tm, tm0, source_mult
  type(density)       ::n, ni, np, nl2, nl2e
  real                ::const
  type(temp)          ::T, Ti, Tp
  real                ::Te0, Ti0, Teh0
  type(height)        ::h, hi
  type(r_ind)         ::ind
  type(nT)            ::nrg, nTi, nTp
  integer             ::i, j, k
  real                ::var, n_height
  type(nu)            ::v, vi
  type(r_dep)         ::dep, depi
  type(lat_dist)      ::lat, lati
  type(ft_int)        ::ft
  type(energy)        ::nrgy
  type(ft_mix)        ::plot
  real                ::output_it
  character(len=8)    ::x1
  character(len=4)    ::day_char
  integer             ::file_num
  real                ::elecHot_multiplier, intensity, n_ave, T_ave, test_multiplier, volume, smooth
  logical             ::isNaN

!  call initNu(v)
  longitude = (mype * 360.0 / LNG_GRID)

  do i=1, LAT_SIZE
    lat%z(i)= (i-1) * Rj / 10.0  !Initializing lat%z
  end do

  call readInputs()  !call to input.f90 module to read initial variables from 'input.dat'
!print *, source
  call read_rec_tables()

!set trans_ variables (user prompt or formatted file migh be used in the future)
  trans_exp=1.0
  trans_type=.false.

!set dt (2000)
  dt=dt_input
!  source = source *2000.0/dt

!set run time
  runt=run_days*8.64e4 !one day = 86400 seconds

  nit=(runt/dt)+1 ! number of iterations to reach run_days

!set radial distance
!  rdist= 6   !in Rj
  dr=((OUT_LIM-IN_LIM)/(RAD_GRID-1))*Rj
  rdist= IN_LIM+dr*(radgrid-1)/Rj   !in Rj
  torus_circumference = Rj * rdist * 2.0 * PI
  dx = torus_circumference / LNG_GRID
  volume=PI*((((rdist*Rj+dr/2.0)*1.0e5))**2 - ((rdist*Rj-dr/2.0)*1.0e5)**2)*0.5*ROOTPI*Rj*1.0e5/LNG_GRID
!  volume=dr*dx*1.0e10*ROOTPI*Rj*.5*1.0e5
!  source=source/volume
  net_source = source * (rdist/6.0)**source_exp
!  source = source*((6.0**source_exp)*(source_exp+1.0))/((9.0**(source_exp+1.0))-(6.0**(source_exp+1.0)))
!  net_source= (source/((6.0**source_exp)*(source_exp+1.0)*(dr/Rj)))*((rdist+(dr/Rj))**(source_exp+1.0) - (rdist**(source_exp+1.0)))
  call MPI_ALLREDUCE(net_source, source_tot, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  var=net_source
  if (mype .eq. 0) print *, "Source correction coefficient : ", source/source_tot
  net_source=net_source*(source/source_tot)
!  if (mype .eq. 0) then
!    print*, "TOTAL SOURCE", var, net_source, source_tot, volume, net_source/volume
!  endif
  call MPI_ALLREDUCE(net_source, source_tot, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
!  if (mype .eq. 0) then
!    print*, "TOTAL SOURCE",  source_tot, source_tot*((o_to_s*8.0+16.0)/(1.0+o_to_s))*mp
!  endif
  net_source=net_source/volume
!  print *, mype, "is at the radial distance", rdist, "Rj, with volumetric source ",  net_source, " /cm3s"
!  print *, source/(ROOTPI*Rj*1e5*.5)

!  print *, net_source, volume, source

  numerical_c_neutral = v_neutral*dt/dx
  numerical_c_ion = v_ion*dt/dx

  if( vrad ) v_ion=1.0-abs(rdist-6.8)
  if( .not. vrad .and. .not. vmass) v_ion=1.0
  if( vrad .and. v_ion .lt. 0.0 ) v_ion=0.0

!set sys3 longitude of box
  lon3=200

!set zoff
  zoff= abs(6.0* cos((lon3-longitude) * dTOr) * dTOr * rdist * Rj) !in km
  n_height = Rj/2.0

  tm0=0.00

!set density values
  const=1800.0

  test_multiplier=1.0
  if( test_pattern ) then
!!    if( rdist .gt. 7.0 .and. rdist .lt. 9.0) test_multiplier=4.0
    test_multiplier=1.0+0.2*cos((longitude-110.0)*PI/180.0)
  endif

!  n%sp = 0.060 * const * test_multiplier* (rdist/6.0)**(-8.0)
!  n%s2p= 0.212 * const * test_multiplier* (rdist/6.0)**(-8.0)
!  n%s3p= 0.034 * const * test_multiplier* (rdist/6.0)**(-3.0)
!  n%op = 0.242 * const * test_multiplier* (rdist/6.0)**(-8.0)
!  n%o2p= 0.242 * const * test_multiplier* (rdist/6.0)**(-3.0)

  n%sp = 150.0 * test_multiplier* (rdist/6.0)**(-10.0)
  n%s2p= 600.0 * test_multiplier* (rdist/6.0)**(-8.0)
  n%s3p=  100.0 * test_multiplier* (rdist/6.0)**(-5.0)
  n%op = 400.0 * test_multiplier* (rdist/6.0)**(-8.0)
  n%o2p=  40.0 * test_multiplier* (rdist/6.0)**(-2.0)

  n%s=0.0!50.0 * test_multiplier* (rdist/6.0)**source_exp
  n%o=0.0!100.0 * test_multiplier* (rdist/6.0)**source_exp


  Te0 = 5.0
  Ti0 = 70.0
  Teh0= tehot*(rdist/6.0)**tehot_alpha
  if(Teh0 .gt. 400.0) Teh0=400.0
  fehot_const= fehot_const*(rdist/6.0)**fehot_exp
  n%fh=fehot_const
!  trans = 4.62963e-7
!  trans = 1.0/((v_r0/dr)*86400.0) 
!  net_source = source*(6.0/rdist)**20 ! ~6.3e6 fix FIX
!  if (radgrid .eq. 1 ) net_source = source

!  do i=1, RAD_GRID
!    source_mult=source_mult+1.0/(10.0**i)
!  end do 

!  net_source=(source/((10.0**radgrid)*source_mult))

  n%elec = (n%sp + n%op + 2 * (n%s2p + n%o2p) + 3 * n%s3p) !* (1.0 - n%protons)
  n%elecHot = n%fh * n%elec! / (1.0-n%fh)
!  n%elecHot = 0.01 * n%elec

  n%fc = 1.0 - n%fh

!set temp values
  T%sp      = Ti0
  T%s2p     = Ti0
  T%s3p     = Ti0
  T%op      = Ti0
  T%o2p     = Ti0
  T%elec    = Te0
  T%elecHot = Teh0
  

!get scale heights 
  call get_scale_heights(h, T, n)

  if (protons > 0.0) then
    n%protons = 0.0!protons
  endif

  ind%o_to_s= o_to_s
  ind%o2s_spike=2.0

!  v_r0=v_r0*(rdist/6.0)**expv_r0
!  v_r0= (v_r0/((6.0**expv_r0)*(expv_r0+1.0)*(dr/Rj)))*((rdist+(dr/Rj))**(expv_r0+1.0) - (rdist**(expv_r0+1.0)))
!  call MPI_ALLREDUCE(dr/v_r0, trans_tot, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  if(mype .eq. 0) then
!    print *, "TOTAL TRANSPORT TIME:::", trans_tot
  endif 
!  net_source= (source/((6.0**source_exp)*(source_exp+1.0)*dr))*((rdist+dr)**(source_exp+1.0) - (rdist**(source_exp+1.0)))
  numerical_c_r=v_r0*dt/dr
 
  if(mype-LNG_GRID.ge.0) call MPI_SEND(numerical_c_r, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype+LNG_GRID<npes) call MPI_RECV(numerical_c_rout, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)

  if(mype+LNG_GRID<npes) call MPI_SEND(numerical_c_r, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype-LNG_GRID.ge.0) call MPI_RECV(numerical_c_rin, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)

!do i=1, npes-1
!  if(mype .eq. i) then
!    open(unit=320, file='AllSource.dat', status='unknown', position='append')
!      write(320,*) rdist, net_source/(ROOTPI*Rj*1e5*.5)  
!    close(320)
!    open(unit=330, file='AllTrans.dat', status='unknown', position='append')
!      write(330,*) rdist, v_r0*dt/dr
!    close(330)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!enddo

!  transport = (v_r0/dr)*86400.0 
!  transport = transport * (6.0/rdist)**5.6

!  tau0=transport !1.0/(trans*8.64e4)
  net_source0=net_source 
  !fh0 = fehot_const

  h%s=n_height
  h%o=n_height

  call InitIndependentRates(ind)

  T%pu_s = Tpu(32.0, rdist*1.0)
  T%pu_o = Tpu(16.0, rdist*1.0)

  T%elecHot=Teh0

  call independent_rates(ind, T, h)

  n%fc= 1.0 - n%fh   

  n%elec = ((n%sp + n%op) + 2.0*(n%s2p + n%o2p) + 3.0 * n%s3p)!/(1.0-n%protons)
  n%elecHot = n%elec * n%fh!/n%fc
  nrg%elec = n%elec * T%elec
  nrg%elecHot = n%elecHot * T%elecHot
  nrg%sp = n%sp * T%sp
  nrg%s2p = n%s2p * T%s2p
  nrg%s3p = n%s3p * T%s3p
  nrg%op = n%op * T%op
  nrg%o2p = n%o2p * T%o2p

  ni=n
  np=n

  Ti=T
  Tp=T

  hi=h
 
  nTi=nrg
  nTp=nrg

  vi=v

  lati=lat

  call get_scale_heights(h, T, n)

  output_it = 0 !This variable determine when data is output. 
  Io_loc=0      !Io's location in the torus
  sys4_loc=(110.0/360.0)*torus_circumference    !The location of the sys4 hot electron population
  file_num=0    !Output files are numbered so they can be assembled as a animated visualization (refer to scripts)

  nl2=NLsquared(n, T, nl2e, h)

  mass_loading(mype+1)=1.0
  ave_loading=1.0

  do i=1, nit
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    if( mype .eq. 0 ) then
!      print *, "((((((((((((((((((( i = ", i, " )))))))))))))))))))"
!    endif
    tm = tm0 + (i-1) * dt / 86400.0
!    if( mype .eq. 0) print*, tm
!if(mype .eq. 6) print *, n%s2p
    var =exp(-((tm-neutral_t0)/neutral_width)**2)

  !  net_source = net_source0*(1.0 + neutral_amp*var) !Ubiquitous source
    if( moving_Io ) then
      if( mype .eq. int(Io_loc*LNG_GRID/torus_circumference) )then
        net_source = LNG_GRID*net_source0*(1.0+neutral_amp*var)
      else
        if( i .eq. 1 ) then
          net_source = net_source0*(1.0+neutral_amp*var)
        else
          net_source=0
        endif
      endif
    endif

    if( .not. moving_Io ) then
      net_source = (net_source0*(1.0 + neutral_amp*var))!/LNG_GRID !ubiquitous
    endif

    ind%o_to_s = o_to_s
!    ind%o_to_s = (otos + o2s_spike * neutral_amp * var) & !o2s_spike
!               / (1.0 + neutral_amp * var)
    n%fh  = fehot_const * (1.0 + hote_amp * var)

    elecHot_multiplier=1.0

    if( sys3hot ) then
      elecHot_multiplier=elecHot_multiplier+sys3_amp*(cos((290.0-longitude)*dTOr))
    endif

    if( sys4hot ) then
      elecHot_multiplier=elecHot_multiplier&
             +sys4_amp*cos(((mype/(LNG_GRID-1.0))-(sys4_loc/torus_circumference))*2.0*PI)
    endif

 !   elecHot_multiplier=elecHot_multiplier*(1.0+0.5*((mass_loading(mype+1)/ave_loading)-1.0))

    n%fh  = fehot_const * (1.0 + hote_amp * var)*elecHot_multiplier

    ni%fh = n%fh
    np%fh = n%fh

    n%fc  = 1.0 - n%fh
    ni%fc = n%fc
    np%fc = n%fc

    n%elecHot = n%elec * n%fh!/n%fc
    nrg%elecHot = n%elecHot * T%elecHot

    do j=1, LAT_SIZE
      lat%z(j)= (j-1) * h%elec / 10.0  !Initializing lat%z
      lat%elec(j) = n%elec*exp(-(lat%z(j)/h%elec)**2)
      lati%elecHot(j) = n%elecHot!*exp(-(lat%z(j)/h%elec)**2)
    end do

    if ( DEBUG_GEN ) then !this variable set in debug.f90
      call DebugOutput(i, n, h, T, v, nrg)
    endif
    
    if( rdist < reac_off_dist ) then
      call cm3_latavg_model(n, T, nrg, h, v, ni, Ti, nTi, hi &
                         ,vi, np, Tp, nTp, ind, dep, depi, lat, lati, ft, zoff) 

      call update_temp(n, nrg, T)
    
    else

      mass_loading(mype+1)=0.0

    end if

!    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    do j=0, npes
!      if(mype .eq. j) print *, rdist, h%sp, T%sp
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    enddo
!    if(mype .eq. 0) print *, ""
 
    call MPI_ALLREDUCE(mass_loading, all_loading, NUMPES, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
!    if( mype .eq. 8 ) print *, all_loading
!    call MPI_ALLREDUCE(mass_loading(mype+1), ave_loading, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    ave_loading=0.0
    do j=(mype/LNG_GRID)*LNG_GRID, (((mype/LNG_GRID)+1)*LNG_GRID)-1
      ave_loading=ave_loading+all_loading(j+1)
    end do 
    ave_loading=ave_loading/(LNG_GRID*1.0)
    call MPI_ALLREDUCE(mass_loading(mype+1), min_loading, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(mass_loading(mype+1), max_loading, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
    if( rdist < reac_off_dist ) then
       smooth=all_loading(int((mype+1)/LNG_GRID)*LNG_GRID+MOD(mype,LNG_GRID)+1)
!       if( mype .eq. 8) print *, smooth, int((mype+1)/LNG_GRID)*LNG_GRID+MOD(mype,LNG_GRID)+1
       do k=1, 1
         smooth=smooth+all_loading(int((mype+1)/LNG_GRID)*LNG_GRID+MOD(mype+k,LNG_GRID)+1)
         smooth=smooth+all_loading(int((mype+1)/LNG_GRID)*LNG_GRID+MODULO(mype-k,LNG_GRID)+1)
!         if( mype .eq. 8) print *, smooth, all_loading(abs(int((mype+1)/LNG_GRID))*LNG_GRID+MODULO(MOD(mype,LNG_GRID)-k,LNG_GRID)+1)&
!                                 ,all_loading(int(mype+1/LNG_GRID)*LNG_GRID+MOD(mype+k,LNG_GRID)+1)
       end do
       smooth=smooth/(1.0+2.0*(k-1))
!       if( mype .eq. 8) print *, smooth, mass_loading(mype+1)
!       print *, mype, smooth, mass_loading(mype+1)
!      v_ion=1.0!+((mass_loading(mype+1)/ave_loading)**1.5)*2.05 !ion lag in km/s where v_ion is scaled on mass_loading
!      n%fh=n%fh*((mass_loading(mype+1)/ave_loading)**(2.0)) !Should replace sys4 population.
!      v_ion=0.0+(((mass_loading(mype+1)*volume*LNG_GRID)*57.0*(rdist)**5)/(0.8*sqrt(1.0-(1.0/(rdist)))*4.0*PI*1.5*((Rj*1.0e3)**2)*(4.2e-4)**2))
      if( vmass ) then
        v_ion=0.1+(((mass_loading(mype+1)*volume*LNG_GRID)*57.0*(rdist)**5)/(0.012*sqrt(1.0-(1.0/(rdist)))*4.0*PI*1.5*((Rj*1.0e3)**2)*(4.2e-4)**2))
        n%fh=n%fh*((mass_loading(mype+1)/ave_loading)**(15.0)) !Should replace sys4 population.
      endif
    end if
    numerical_c_ion = v_ion*dt/dx
!    if(mype .eq. 0) print *, v_ion, ave_loading*volume*LNG_GRID, mass_loading(mype+1)/ave_loading

    if( i .eq. nit .and. .true. ) then
       open(unit=15, file='sub.dat', status='unknown',position='append')
       do j=1, npes
         if(mype .eq. j) write(15,*) MOD(longitude, 360.0), v_ion,  rdist
         if( mod(j-1,LNG_GRID) .eq. 0) write (15,*) ""
       enddo
       close(15)
    end if

    call get_scale_heights(h, T, n)

    call energyBudget(n, h, T, dep, ind, ft, lat, v, nrgy)

    if (nint(output_it)+1 .eq. i .and. (OUTPUT_MIXR .or. OUTPUT_DENS .or. OUTPUT_TEMP .or. OUTPUT_INTS)) then !Output at set intervals when OUTPUT_MIX is true (from debug.f90)
        miscOutput=mass_loading(mype+1)
        day = (i-1.0)*dt/86400
        write (x1, '(I4.4)') file_num
        day_char=trim(x1)  !trim non-existent file and store as day_char
        !if( mype .eq. 0 ) then
        !endif
        call dens_ave(n_ave, n) 
        call temp_ave(T_ave, T) 
        do k=0, RAD_GRID-1
          do j=0, LNG_GRID
            if( mype .eq. mod(j,LNG_GRID)+(k*LNG_GRID)) then
              if(OUTPUT_DENS) then
                 call IonElecOutput(n%sp, n%s2p, n%s3p, n%op, n%o2p, n%elec,&
                  longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'DENS')
                 call IonElecOutput3D(n%sp, n%s2p, n%s3p, n%op, n%o2p, n%elec,&
                  longitude, rdist, day_char, 'DENS')
!                 if( rdist < reac_off_dist ) call OtherOutput(200.0*mass_loading(mype+1)/ave_loading, longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'LOAD')
!                 if( rdist < reac_off_dist ) call OtherOutput3D(200.0*mass_loading(mype+1)/ave_loading, longitude+((j+1)/(LNG_GRID+1))*360.0, rdist, day_char, 'LOAD')
                 call OtherOutput(mass_loading(mype+1)/ave_loading, longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'LOAD')
                 call OtherOutput3D(mass_loading(mype+1), longitude+((j+1)/(LNG_GRID+1))*360.0, rdist, day_char, 'LOAD')
                 call OtherOutput3D(v_ion, longitude+((j+1)/(LNG_GRID+1))*360.0, rdist, day_char, 'VSUB')
                 call OtherOutput3D(nrgy%puv, longitude+((j+1)/(LNG_GRID+1))*360.0, rdist, day_char, 'PUV_')
                 call OtherOutput(real(miscOutput), longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'MOUT')
                 if( j .eq. LNG_GRID ) call spacer3D(day_char, 'DENS')
                 if( j .eq. LNG_GRID ) call otherspacer3D(day_char, 'LOAD')
                 if( j .eq. LNG_GRID ) call otherspacer3D(day_char, 'VSUB')
                 if( j .eq. LNG_GRID ) call otherspacer3D(day_char, 'PUV_')
              endif 
              if(OUTPUT_NL2) then
                call IonElecOutput(nl2%sp, nl2%s2p, nl2%s3p, nl2%op, nl2%o2p,nl2%sp+nl2%s2p+nl2%s3p+nl2%op+nl2%o2p,&
                  longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'NL2_')
                call IonElecOutput3D(nl2%sp, nl2%s2p, nl2%s3p, nl2%op, nl2%o2p,nl2%sp+nl2%s2p+nl2%s3p+nl2%op+nl2%o2p,&
                  longitude, rdist, day_char, 'NL2_')
                if( j .eq. LNG_GRID ) call spacer3D(day_char, 'NL2_')
                !call IonElecOutput(nl2e%sp, nl2e%s2p, nl2e%s3p, nl2e%op, nl2e%o2p,nl2e%sp+nl2e%s2p+nl2e%s3p+nl2e%op+nl2e%o2p,&
                !#  longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'NL2e')
                !call IonElecOutput3D(nl2e%sp, nl2e%s2p, nl2e%s3p, nl2e%op, nl2e%o2p,nl2e%sp+nl2e%s2p+nl2e%s3p+nl2e%op+nl2e%o2p,&
                !  longitude, rdist, day_char, 'NL2e')
                !if( j .eq. LNG_GRID ) call spacer3D(day_char, 'NL2e')
              endif
              if(OUTPUT_MIXR) then  
                plot = ftint_mix(n, h) !calculate values to be plotted
                call IonOutput(plot%sp, plot%s2p, plot%s3p, plot%op, plot%o2p, &
                  longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'MIXR')
                 call IonOutput3D(plot%sp, plot%s2p, plot%s3p, plot%op, plot%o2p, &
                  longitude, rdist, day_char, 'MIXR')
                if( j .eq. LNG_GRID ) call spacer3D(day_char, 'MIXR')
              endif
              if(OUTPUT_TEMP) then
                call IonElecOutput(T%sp, T%s2p, T%s3p, T%op, T%o2p, T%elec, & !longitude, day_char, 'TEMP')
                  longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'TEMP')
                call IonElecOutput3D(T%sp, T%s2p, T%s3p, T%op, T%o2p, T%elec, & !longitude, day_char, 'TEMP')
                  longitude, rdist, day_char, 'TEMP')
                if( j .eq. LNG_GRID ) call spacer3D(day_char, 'TEMP')
              end if
              if(OUTPUT_INTS) then !Intensity
                call IonOutput(n%sp*T%sp, n%s2p*T%s2p, n%s3p*T%s3p, n%op*T%op, n%o2p*T%o2p,&
                  longitude+((j+1)/(LNG_GRID+1))*360.0, day_char, 'INTS')
                 intensity= n%sp*T_ave/(n_ave*T%sp)
                 open(unit=120, file='intensity'//day_char//'.dat', status='unknown', position='append')
                 write(120,*) longitude, intensity 
                 close(120)
              end if
!              open(unit=200, file='feh'//day_char//'.dat', status='unknown', position='append')
!              open(unit=210, file='vr'//day_char//'.dat', status='unknown', position='append')
!              open(unit=220, file='source'//day_char//'.dat', status='unknown', position='append')
!                 write(200,*) longitude, n%fh
!                 write(210,*) longitude, v_r0, numerical_c_r
!                 write(220,*) longitude, net_source
!              close(200)
!              close(210)
!              close(220)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          end do
        end do
        output_it=output_it + (86400.0/(dt*per_day)) !Determines when data is output. Set for once each run day (86400/dt).
        file_num = file_num + 1
    endif        
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    isNaN=NaNcatch(n%sp, 1000, mype) !fix
    call Grid_transport(n, T, nrg, dep, h, nl2, nl2e)
!    isNaN=NaNcatch(n%sp, 100, mype) !fix
    Io_loc = mod(Io_loc+(dt*v_Io), torus_circumference)  
    sys4_loc = mod(sys4_loc+(dt*v_sys4), torus_circumference)  

  end do

call FinalOutput(nrgy)

end subroutine model

subroutine dens_ave(n_ave, n)!, i)
  real                ::n_tot, n_ave, space_ave
  type(density)       ::n
!  integer             ::i

  n_tot=n%sp !+n%s2p+n%s3p+n%op+n%o2p
  call MPI_REDUCE(n_tot, n_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  n_ave=n_ave/LNG_GRID
  call MPI_BCAST(n_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

!  n_tot=n%sp+n%s2p+n%s3p+n%op+n%o2p
!  call MPI_REDUCE(n_tot, space_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!  space_ave=space_ave/LNG_GRID
!  n_ave=(n_ave*(i-1.0)+space_ave)/(i*1.0)
!  call MPI_BCAST(n_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

end subroutine dens_ave

subroutine temp_ave(T_ave, T)!, i)
  real                ::T_tot, T_ave, space_ave
  type(temp)          ::T
!  integer             ::i

  T_tot=T%sp !+T%s2p+T%s3p+T%op+T%o2p
  call MPI_REDUCE(T_tot, T_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  T_ave=T_ave/LNG_GRID
  call MPI_BCAST(T_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

!  T_tot=T%sp+T%s2p+T%s3p+T%op+T%o2p
!  call MPI_REDUCE(T_tot, space_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!  space_ave=space_ave/LNG_GRID
!  T_ave=(T_ave*(i-1.0)+space_ave)/(i*1.0)
!  call MPI_BCAST(T_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

end subroutine temp_ave

subroutine FinalOutput(nrgy)
  type(energy)        ::nrgy, avg
  integer             ::j

  call MPI_REDUCE(nrgy%s_ion, avg%s_ion, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%s_ion=avg%s_ion/LNG_GRID

  call MPI_REDUCE(nrgy%s_cx, avg%s_cx, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%s_cx=avg%s_cx/LNG_GRID

  call MPI_REDUCE(nrgy%o_ion, avg%o_ion, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%o_ion=avg%o_ion/LNG_GRID

  call MPI_REDUCE(nrgy%o_cx, avg%o_cx, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%o_cx=avg%o_cx/LNG_GRID

  call MPI_REDUCE(nrgy%elecHot_eq, avg%elecHot_eq, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%elecHot_eq=avg%elecHot_eq/LNG_GRID

  call MPI_REDUCE(nrgy%tot_eq, avg%tot_eq, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%tot_eq=avg%tot_eq/LNG_GRID

  call MPI_REDUCE(nrgy%P_in, avg%P_in, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%P_in=avg%P_in/LNG_GRID

  if( nrgy%Puv .ne. nrgy%Puv ) nrgy%Puv =0.0
!  print *, mype, nrgy%Puv
  call MPI_REDUCE(nrgy%Puv, avg%Puv, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Puv=avg%Puv

  call MPI_REDUCE(nrgy%Pfast, avg%Pfast, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Pfast=avg%Pfast/LNG_GRID

  call MPI_REDUCE(nrgy%Ptrans, avg%Ptrans, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Ptrans=avg%Ptrans/LNG_GRID

  call MPI_REDUCE(nrgy%Ptrans_elecHot, avg%Ptrans_elecHot, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Ptrans_elecHot=avg%Ptrans_elecHot/LNG_GRID

  call MPI_REDUCE(nrgy%P_out, avg%P_out, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%P_out=avg%P_out/LNG_GRID

  if( mype .eq. 0 ) then
    print *, "AVERAGE VALUES"
    call FinalTable(avg)
    print *, ""
  endif

!  do j=1, LNG_GRID
!    if( mygrid .eq. j ) then
!      print *, "mygrid = ", mygrid
!      call FinalTable(nrgy)
!      print *, ""
!    endif
!    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!  enddo

end subroutine FinalOutput

subroutine FinalTable(nrgy)
  type(energy)        ::nrgy

  print *, '$$--------------------------------'
  print *, '$$ INPUT PARAMETERS'
  print *, '$$--------------------------------'
  print *, '$$ Source Rate..........', source 
  print *, '$$ Hot Elec Fraction....', fehot_const 
  print *, '$$ Transport............', transport 
  print *, '$$ Hot Elec Temp........', tehot  
  print *, '$$ Subcorotation........', v_ion
  print *, '$$ Sys 3 Hot Elec Frac..', sys3_amp 
  print *, '$$ Sys 4 Hot Elec Amp...', sys4_amp
  print *, '$$ Sys 4 Speed..........', v_sys4

  print *, ''
  print *, '$$ GAUSSIAN SOURCE CHANGE VARIABLES'  
  print *, '$$ Neutral Variation....', neutral_amp 
  print *, '$$ Neutral Temp.........', neutral_t0 
  print *, '$$ Variation Duration...', neutral_width
  print *, '$$ Hot Elec Variation...', hote_amp  
  print *, '$$ Hot Elec Temp........', hote_t0
  print *, '$$ Variation Duration...', hote_width
 
  print *, ''
  print *, '$$ Run Length(days)....',  run_days
  print *, '$$ Outputs per day.....',  per_day

  print *, ''
  print *, '$$--------------------------------'
  print *, '$$ IN-CODE ENERGY BUDGET'
  print *, '$$--------------------------------'
  print *, '$$ ionized S............', nrgy%s_ion
  print *, '$$ ionized O............', nrgy%o_ion
  print *, '$$ charge exchange S....', nrgy%s_cx
  print *, '$$ charge exchange O....', nrgy%o_cx
  print *, '$$ equil with ehot......', nrgy%elecHot_eq + nrgy%tot_eq
  print *, '$$ total in.............', nrgy%P_in + nrgy%tot_eq
  print *, '$$ puv..................', nrgy%Puv
  print *, '$$ fast/ena.............', nrgy%pfast - nrgy%tot_eq
  print *, '$$ transport............', nrgy%ptrans + nrgy%ptrans_elecHot
  print *, '$$ total out............', nrgy%P_out - nrgy%tot_eq
  print *, '$$ in/out...............', (nrgy%P_in + nrgy%tot_eq )/(nrgy%P_out - nrgy%tot_eq )
!  print *, ""
!  print *, '++++++++++++++++++++++++++++++++++++'
!  print *, 'Final Variable Values'
!  print *, '++++++++++++++++++++++++++++++++++++'
!  print *, 'O/S.........................', o_to_s
!  print *, 'Fraction of Hot Electrons...', fehot_const
!  print *, 'Transport...................', transport
!  print *, 'Hot Electron Temp...........', tehot
!  print *, 'Lag Constant................', lag_const
!  print *, 'Neutral Amplitude...........', neutral_amp
!  print *, 'Inital Neutral Temperature..', neutral_t0
!  print *, 'Neutral Width...............', neutral_width
!  print *, 'Hot Electron Amplitude......', hote_amp
!  print *, 'Hot Electron Initial Temp...', hote_t0
!  print *, 'Hot Electron Width..........', hote_width

end subroutine FinalTable

subroutine DebugOutput(i, n, h, T, v, nrg)
  integer             ::i
  type(density)       ::n
  type(height)        ::h
  type(temp)          ::T
  type(nu)            ::v
  type(nT)            ::nrg

  print *,  "||||||||||||||||||||||||||||||||||||||||||||||"
  print *,  "lnggrid = ", lnggrid
  print *,  "radgrid = ", radgrid
  print *,  "i = ", i-1
  print *,  "||||||||||||||||||||||||||||||||||||||||||||||"
  print *, "~~~~~~~~~~~~~DENSITY~~~~~~~~~~~~~"
  call output(n)
  print *, "~~~~~~~~~~~~~HEIGHT~~~~~~~~~~~~~~"
  call output(h)
  print *, "~~~~~~~~~~~TEMPERATURE~~~~~~~~~~~"
  call output(T)
  print *, "~~~~~~~~~~~~~~~NU~~~~~~~~~~~~~~~~"
  call output(v)
  print *, "~~~~~~~~~~~~~ENERGY~~~~~~~~~~~~~~"
  call output(nrg)
 

end subroutine DebugOutput

subroutine Grid_transport(n, T, nrg, dep, h, nl2, nl2e)
  type(height)        ::h
  type(temp)          ::T
  type(density)       ::n, dens_source, nl2e, nl2
  type(nT)            ::nrg, nrg_source
  type(r_dep)         ::dep
!  real                ::dll0, dlla
  logical             ::isNaN
  integer             ::i

!  dll0=4.2E-8 !4.2E-7
!  dlla=4.6

  do i=1, aztrans_it
    call az_transport(n, nrg)
  enddo
  call update_temp(n, nrg, T)
  isNaN=NaNcatch(n%sp, -1, mype)
  nl2=NLsquared(n, T, nl2e, h)
!  if(mype .eq. 0) print*, n%sp, n%s2p, n%s3p, n%op, n%o2p
!  if(mype .eq. 0) print*, nl2%sp, nl2%s2p, nl2%s3p, nl2%op, nl2%o2p
  isNaN=NaNcatch(nl2%sp, 0, mype)
  do i=1, radtrans_it
    call transport_nl2(nl2, nl2e, dll0, dlla)
!    if(mype .eq. 0) print *, i 
  enddo

  call iterate_NL2(nl2, nl2e, n, T, h)

  T%sp=(nl2e%sp/(nl2%sp*rdist**2))**(3.0/4.0)
  nrg%sp=n%sp*T%sp
  T%s2p=(nl2e%s2p/(nl2%s2p*rdist**2))**(3.0/4.0)
  nrg%s2p=n%s2p*T%s2p
  T%s3p=(nl2e%s3p/(nl2%s3p*rdist**2))**(3.0/4.0)
  nrg%s3p=n%s3p*T%s3p
  T%op=(nl2e%op/(nl2%op*rdist**2))**(3.0/4.0)
  nrg%op=n%op*T%op
  T%o2p=(nl2e%o2p/(nl2%o2p*rdist**2))**(3.0/4.0)
  nrg%o2p=n%o2p*T%o2p

!  T%sp=(nl2e%sp/(nl2%sp*rdist**2))**(3/4)
!  T%s2p=(nl2e%s2p/(nl2%s2p*rdist**2))**(3/4)
!  T%s3p=(nl2e%s3p/(nl2%s3p*rdist**2))**(3/4)
!  T%op=(nl2e%op/(nl2%op*rdist**2))**(3/4)
!  T%o2p=(nl2e%o2p/(nl2%o2p*rdist**2))**(3/4)

  isNaN=NaNcatch(nl2%sp, 10, mype)
!  if(mype .eq. 0) print*, n%sp, n%s2p, n%s3p, n%op, n%o2p
!  if(mype .eq. 0) print*, ""
  n%elec=(n%sp + 2.0*n%s2p + 3.0*n%s3p + n%op + 2.0*n%o2p)/(1.0-n%protons) 
!  if(rdist .lt. reac_off_dist) T%elec=nrg%elec/n%elec
  n%elecHot=n%fh*n%elec/(1.0-n%fh)
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine Grid_transport

subroutine az_transport(n, nrg)
  type(density)     ::n
  type(nT)          ::nrg
  real              ::cleft, cright, cn, c
  cn=numerical_c_neutral/aztrans_it
  c=numerical_c_ion/aztrans_it

  if(UseLaxWendroff) then
    call GetAzNeighbors(n, nrg, cleft, cright)
    cleft=cleft/aztrans_it
    cright=cright/aztrans_it
    n%s=LaxWendroff(nleft%s, n%s, nright%s, cn, cn, cn)
    n%sp=LaxWendroff(nleft%sp, n%sp, nright%sp, cleft, c, cright)
    n%s2p=LaxWendroff(nleft%s2p, n%s2p, nright%s2p, cleft, c, cright)
    n%s3p=LaxWendroff(nleft%s3p, n%s3p, nright%s3p, cleft, c, cright)
    n%o=LaxWendroff(nleft%o, n%o, nright%o, cn, cn, cn)
    n%op=LaxWendroff(nleft%op, n%op, nright%op, cleft, c, cright)
    n%o2p=LaxWendroff(nleft%o2p, n%o2p, nright%o2p, cleft, c, cright)
!double precision function LaxWendroff(left, center, right, uleft, u, uright)

    nrg%sp=LaxWendroff(nTleft%sp, nrg%sp, nTright%sp, cleft, c, cright)
    nrg%s2p=LaxWendroff(nTleft%s2p, nrg%s2p, nTright%s2p, cleft, c, cright)
    nrg%s3p=LaxWendroff(nTleft%s3p, nrg%s3p, nTright%s3p, cleft, c, cright)
    nrg%op=LaxWendroff(nTleft%op, nrg%op, nTright%op, cleft, c, cright)
    nrg%o2p=LaxWendroff(nTleft%o2p, nrg%o2p, nTright%o2p, cleft, c, cright)
  endif
    
  if( Upwind ) then
    call GetAzNeighbors(n, nrg, cleft, cright)
!    print *, mype, numerical_c_ion, cleft, "mype, c, cleft"
!    if(mype .eq. 0) then
!      print *, "AZ Transport:::", n%sp-UpwindTransport(nleft%sp , n%sp ,numerical_c_ion)
!    endif 
    n%s   =UpwindTransport(nleft%s  , n%s  ,numerical_c_neutral, numerical_c_neutral)
    n%sp  =UpwindTransport(nleft%sp , n%sp ,cleft, numerical_c_ion)
    n%s2p =UpwindTransport(nleft%s2p, n%s2p,cleft, numerical_c_ion)
    n%s3p =UpwindTransport(nleft%s3p, n%s3p,cleft, numerical_c_ion)
    n%o   =UpwindTransport(nleft%o  , n%o  ,numerical_c_neutral, numerical_c_neutral)
    n%op  =UpwindTransport(nleft%op , n%op ,cleft, numerical_c_ion)
    n%o2p =UpwindTransport(nleft%o2p, n%o2p,cleft, numerical_c_ion)
  
    nrg%sp  =UpwindTransport(nTleft%sp , nrg%sp ,cleft, numerical_c_ion)
    nrg%s2p =UpwindTransport(nTleft%s2p, nrg%s2p,cleft, numerical_c_ion)
    nrg%s3p =UpwindTransport(nTleft%s3p, nrg%s3p,cleft, numerical_c_ion)
    nrg%op  =UpwindTransport(nTleft%op , nrg%op ,cleft, numerical_c_ion)
    nrg%o2p =UpwindTransport(nTleft%o2p, nrg%o2p,cleft, numerical_c_ion)
  endif

  if( Euler ) then
    n%s   = EulerTransport(n%s  , v_neutral)
    n%sp  = EulerTransport(n%sp , v_ion)
    n%s2p = EulerTransport(n%s2p, v_ion)
    n%s3p = EulerTransport(n%s3p, v_ion)
    n%o   = EulerTransport(n%o  , v_neutral)
    n%op  = EulerTransport(n%op , v_ion)
    n%o2p = EulerTransport(n%o2p, v_ion)

    nrg%sp  = EulerTransport(nrg%sp , v_ion)
    nrg%s2p = EulerTransport(nrg%s2p, v_ion)
    nrg%s3p = EulerTransport(nrg%s3p, v_ion)
    nrg%op  = EulerTransport(nrg%op , v_ion)
    nrg%o2p = EulerTransport(nrg%o2p, v_ion)
  endif

end subroutine az_transport

subroutine rad_transport(n, nrg, dep, h)
  type(density)     ::n, nempty, ni,ninold
  type(nT)          ::nrg, nTempty, nrgi,nTinold
  type(r_dep)       ::dep
  type(height)      ::h
  double precision  ::loss
  
  call GetRadLeftNeighbors(n, nrg, h)
  call GetRadRightNeighbors(n, nrg, h)

!  nlast=n
!  nTlast=nrg

  if (mype .lt. LNG_GRID) then
    nin%sp=n%sp
    nin%s2p=n%s2p
    nin%s3p=n%s3p
    nin%op=n%op
    nin%o2p=n%o2p
    nin%elec=n%elec
!    call output(nin) 
!    call output(n) 
    nTin%sp=70.0*n%sp
    nTin%s2p=70.0*n%s2p
    nTin%s3p=70.0*n%s3p
    nTin%op=70.0*n%op
    nTin%o2p=70.0*n%o2p
    nTin%elec=5.0*n%elec
    numerical_c_rin=numerical_c_r*((6.0-dr)/6.0)**expv_r0
  endif
!    if (mype .eq. 0) then
      !print*, n%sp-UpwindTransport(nin%sp ,  n%sp , numerical_c_r)
!    endif
!  if (mype .eq. 0) then
!    print *, "NCR DEBUG", numerical_c_r, v_r0, dr, dt
!  endif 
!    if( mype .eq. 0) then
!      print *, "Transport loss::" , n%sp-UpwindRadTransport(nin%sp ,  n%sp , numerical_c_r)
!
!    endif
   if( .false. ) then
    n%sp  =UpwindRadTransport(nin%sp ,  n%sp , numerical_c_r, numerical_c_rin)
    n%s2p =UpwindRadTransport(nin%s2p,  n%s2p, numerical_c_r, numerical_c_rin)
    n%s3p =UpwindRadTransport(nin%s3p,  n%s3p, numerical_c_r, numerical_c_rin)
    n%op  =UpwindRadTransport(nin%op ,  n%op , numerical_c_r, numerical_c_rin)
    n%o2p =UpwindRadTransport(nin%o2p,  n%o2p, numerical_c_r, numerical_c_rin)
    n%elec=UpwindRadTransport(nin%elec, n%elec,numerical_c_r, numerical_c_rin)
 
    nrg%sp  =UpwindRadTransport(nTin%sp ,  nrg%sp , numerical_c_r, numerical_c_rin)
    nrg%s2p =UpwindRadTransport(nTin%s2p,  nrg%s2p, numerical_c_r, numerical_c_rin)
    nrg%s3p =UpwindRadTransport(nTin%s3p,  nrg%s3p, numerical_c_r, numerical_c_rin)
    nrg%op  =UpwindRadTransport(nTin%op ,  nrg%op , numerical_c_r, numerical_c_rin)
    nrg%o2p =UpwindRadTransport(nTin%o2p,  nrg%o2p, numerical_c_r, numerical_c_rin)
    nrg%elec=UpwindRadTransport(nTin%elec, nrg%elec,numerical_c_r, numerical_c_rin)
 !   nrg%elecHot =UpwindTransport(nTin%elecHot, nrg%elecHot,numerical_c_r)
  else
!double precision function LWRadTransport(inside, center, outside, cin, c, cout)
    if (mype .ge. npes-LNG_GRID ) then 
      nTout%sp=100.0*n%sp
      nTout%s2p=100.0*n%s2p
      nTout%s3p=100.0*n%s3p
      nTout%op=100.0*n%op
      nTout%o2p=100.0*n%o2p
!      nTout%elec=100.0*n%elec
      nout%sp=n%sp
      nout%s2p=n%s2p
      nout%s3p=n%s3p
      nout%op=n%op
      nout%o2p=n%o2p
      nout%elec=n%elec
      numerical_c_rout=numerical_c_r*((9.0+dr)/6.0)**expv_r0
    endif 
    n%sp  =LWRadTransport(nin%sp  ,  n%sp  , nout%sp  , numerical_c_rin, numerical_c_r, numerical_c_rout)
    n%s2p =LWRadTransport(nin%s2p ,  n%s2p , nout%s2p , numerical_c_rin, numerical_c_r, numerical_c_rout)
    n%s3p =LWRadTransport(nin%s3p ,  n%s3p , nout%s3p , numerical_c_rin, numerical_c_r, numerical_c_rout)
    n%op  =LWRadTransport(nin%op  ,  n%op  , nout%op  , numerical_c_rin, numerical_c_r, numerical_c_rout)
    n%o2p =LWRadTransport(nin%o2p ,  n%o2p , nout%o2p , numerical_c_rin, numerical_c_r, numerical_c_rout)
    n%elec=LWRadTransport(nin%elec,  n%elec, nout%elec, numerical_c_rin, numerical_c_r, numerical_c_rout)
  
!    nrg%sp  =LWRadTransport(nTin%sp  ,  nrg%sp  , nTout%sp  , numerical_c_rin, numerical_c_r, numerical_c_rout)
!    nrg%s2p =LWRadTransport(nTin%s2p ,  nrg%s2p , nTout%s2p , numerical_c_rin, numerical_c_r, numerical_c_rout)
!    nrg%s3p =LWRadTransport(nTin%s3p ,  nrg%s3p , nTout%s3p , numerical_c_rin, numerical_c_r, numerical_c_rout)
!    nrg%op  =LWRadTransport(nTin%op  ,  nrg%op  , nTout%op  , numerical_c_rin, numerical_c_r, numerical_c_rout)
!    nrg%o2p =LWRadTransport(nTin%o2p ,  nrg%o2p , nTout%o2p , numerical_c_rin, numerical_c_r, numerical_c_rout)
!    if (mype .eq. 5) then
!      print *, nrg%elec, LWRadTransport(nTin%elec,  nrg%elec, nTout%elec, numerical_c_rin, numerical_c_r, numerical_c_rout)
!    endif
!    nrg%elec=LWRadTransport(nTin%elec,  nrg%elec, nTout%elec, numerical_c_rin, numerical_c_r, numerical_c_rout)
  endif
!Improved Euler method
!  ni%sp   = n%sp  + (nin%sp  - n%sp )*numerical_c_r
!  ni%s2p  = n%s2p + (nin%s2p - n%s2p)*numerical_c_r 
!  ni%s3p  = n%s3p + (nin%s3p - n%s3p)*numerical_c_r 
!  ni%op   = n%op  + (nin%op  - n%op )*numerical_c_r 
!  ni%o2p  = n%o2p + (nin%o2p - n%o2p)*numerical_c_r 
!
!  nrgi%sp   = nrg%sp   + (nTin%sp   - nrg%sp )*numerical_c_r 
!  nrgi%s2p  = nrg%s2p  + (nTin%s2p  - nrg%s2p)*numerical_c_r 
!  nrgi%s3p  = nrg%s3p  + (nTin%s3p  - nrg%s3p)*numerical_c_r 
!  nrgi%op   = nrg%op   + (nTin%op   - nrg%op )*numerical_c_r 
!!  nrgi%o2p  = nrg%o2p  + (nTin%o2p  - nrg%o2p)*numerical_c_r 

!  ninold=nin
!  nTinold=nTin
!
!  call GetRadNeighbors(ni, nrgi)
!
!  n%sp   = n%sp  + .5 * (nin%sp  + ninold%sp - ni%sp - n%sp )*numerical_c_r
!  n%s2p  = n%s2p + .5 * (nin%s2p + ninold%s2p- ni%s2p- n%s2p)*numerical_c_r
!  n%s3p  = n%s3p + .5 * (nin%s3p + ninold%s3p- ni%s3p- n%s3p)*numerical_c_r
!  n%op   = n%op  + .5 * (nin%op  + ninold%op - ni%op - n%op )*numerical_c_r
!  n%o2p  = n%o2p + .5 * (nin%o2p + ninold%o2p- ni%o2p- n%o2p)*numerical_c_r
!
!  nrg%sp   = nrg%sp  + .5 * (nTin%sp  + nTinold%sp  - nrgi%sp - nrg%sp )*numerical_c_r 
!  nrg%s2p  = nrg%s2p + .5 * (nTin%s2p + nTinold%s2p - nrgi%s2p- nrg%s2p)*numerical_c_r 
!!  nrg%s3p  = nrg%s3p + .5 * (nTin%s3p + nTinold%s3p - nrgi%s3p- nrg%s3p)*numerical_c_r 
!  nrg%op   = nrg%op  + .5 * (nTin%op  + nTinold%op  - nrgi%op - nrg%op )*numerical_c_r 
!  nrg%o2p  = nrg%o2p + .5 * (nTin%o2p + nTinold%o2p - nrgi%o2p- nrg%o2p)*numerical_c_r 
    !handles all radial transport 
    !must remove radial loss from F_* and EF_* functions in functions.f90

end subroutine rad_transport


double precision function UpwindTransport(left, center, cleft, c)
  double precision    ::left, center
  real                ::c, cleft

!  UpwindTransport = (numerical_s + c)*left + (1 - 2*numerical_s - c)*center + numerical_s*right
  UpwindTransport =  cleft*left + (1.0 - c)*center

end function UpwindTransport

double precision function LaxWendroff(left, center, right, cleft, c, cright)
  double precision    :: left, center, right
  real                :: cleft, c, cright
  
  LaxWendroff=center+0.5*(cleft*left - cright*right)+0.5*((cleft**2.0)*left-2.0*c*c*center+(cright**2.0)*right)

end function LaxWendroff

double precision function UpwindRadTransport(left, center, c, cin)
  double precision    ::left, center
  real                ::c, cin

!  UpwindTransport = (numerical_s + c)*left + (1 - 2*numerical_s - c)*center + numerical_s*right
  UpwindRadTransport = cin*left + (1.0 - c)*center

end function UpwindRadTransport

double precision function LWRadTransport(inside, center, outside, cin, c, cout)
  double precision    ::inside, center, outside
  real                ::cin, c, cout
  double precision    ::qplushalf, qminushalf, qND, qNDout, qNDin
  double precision    ::fDphalf, fADphalf, fDmhalf, fADmhalf
  double precision    ::qD, qDin, qDout, qD2in, qD2out
  double precision    ::delqDout, delqDin, delqD2out, delqD2in
  double precision    ::FCTout, FCTin, sigout, sigin
  double precision    ::minout, minin, maxout, maxin
!  double precision    ::

!calculate half time steps
  qplushalf=.5*(center+outside) - .5*(cout*outside - c*center)
  qminushalf=.5*(center+inside) - .5*(cin*inside - c*center)

!calculate non-diffusive updated value
  qND=center-(cout*qplushalf - cin*qminushalf)

  qNDout=qND
  qNDin=qND

  if(mype-LNG_GRID.ge.0) call MPI_SEND(qND, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype+LNG_GRID<npes) call MPI_RECV(qNDout, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)
  
  if(mype+LNG_GRID<npes) call MPI_SEND(qND, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype-LNG_GRID.ge.0) call MPI_RECV(qNDin, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)
  
!calculate diffusive and antidiffusive flux
  fDphalf=((1.0/6.0)+cout*cout/3.0)*(outside-center)      
  fADphalf=(1.0-cout*cout)*(outside-center)/6.0      
  fDmhalf=-((1.0/6.0)+cin*cin/3.0)*(inside-center)      
  fADmhalf=-(1.0-cin*cin)*(inside-center)/6.0      

!value with second order diffusion
  qD=qND+fDphalf-fDmhalf

  qDout=qD
  qDin=qD

  if(mype-LNG_GRID.ge.0) call MPI_SEND(qD, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype+LNG_GRID<npes) call MPI_RECV(qDout, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)
  
  if(mype+LNG_GRID<npes) call MPI_SEND(qD, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype-LNG_GRID.ge.0) call MPI_RECV(qDin, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)

  qD2out=qDout
  qD2in=qDin

  if(mype-LNG_GRID.ge.0) call MPI_SEND(qDout, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype+LNG_GRID<npes) call MPI_RECV(qD2out, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)
  
  if(mype+LNG_GRID<npes) call MPI_SEND(qDin, 1, MPI_DOUBLE_PRECISION, mype+LNG_GRID, 22, MPI_COMM_WORLD, ierr)
  if(mype-LNG_GRID.ge.0) call MPI_RECV(qD2in, 1, MPI_DOUBLE_PRECISION, mype-LNG_GRID, 22, MPI_COMM_WORLD, stat, ierr)


!spatial variation of diffusive value
  delqDout=qDout-qD
  delqDin=qDin-qD

  delqD2out=qD2out-qDout
  delqD2in=qD2in-qDin

  sigout=sign(1.0,real(fADphalf))
  sigin=sign(1.0,real(fADmhalf))

  if(sigout*delqDin < fADphalf) then
    minout=sigout*delqDin
  else 
    minout=fADphalf
  endif
  if(sigout*delqD2out<minout) minout=sigout*delqD2out
 
  if(sigin*delqD2in < fADmhalf) then
    minin=sigin*delqD2in
  else 
    minin=fADmhalf
  endif
  if(sigin*delqDout<minin) minin=sigin*delqDout

  maxout=minout
  if(maxout<0.0) maxout = 0.0

  maxin=minin
  if(maxin<0.0) maxin = 0.0

  FCTout=sigout*maxout
  FCTin=sigin*maxin

  LWRadTransport=center+FCTout-FCTin

end function LWRadTransport

double precision function EulerTransport(old, v) !improved euler method applied to azimuthal transport
  double precision    ::old, intermediate, loss
  real                ::v

  loss = getLoss(v, old)

  intermediate = old - loss

  EulerTransport = old - .5 * (loss + getLoss(v, intermediate))

end function EulerTransport

double precision function getLoss(v, val)
  real                ::v
  double precision    ::val, source

  getLoss = val * dt * v * LNG_GRID / torus_circumference
  
  call MPI_SEND(getLoss, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(source, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  getLoss = getLoss - source

end function getLoss

subroutine GetAzNeighbors(n, nrg, cleft, cright)
  type(density)       ::n
  type(nT)            ::nrg
  integer             ::left, right, i
  real                ::cleft, cright

! ALGORITHM
! Send to right
! Receive from left
! Send to left X           upwind only needs left
! Receive from right X
!!!!!!!!!!!!!!!!!!!!DENSITY!!!!!!!!!!!!!!!!!!!!!
  do i=0, RAD_GRID-1
    if(mype .ge. i*LNG_GRID .and. mype < (i+1)*LNG_GRID) then
      left = mype - 1 
      right= mype + 1

      if( left < (i*LNG_GRID) )           left = left+LNG_GRID
      if( right .ge. ((i+1)*LNG_GRID) ) right=right-LNG_GRID

      call MPI_SENDRECV(numerical_c_ion, 1, MPI_REAL, right, 22, cleft, 1, MPI_REAL, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%s, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%s, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)
!      call MPI_RECV(nleft%s, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%sp, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%sp, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%s2p, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%s2p, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%s3p, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%s3p, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%o, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%o, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%op, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%op, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%o2p, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%o2p, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%elec, 1, MPI_DOUBLE_PRECISION, right, 22, nleft%elec, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(numerical_c_ion, 1, MPI_REAL, left, 22, cright, 1, MPI_REAL, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%s, 1, MPI_DOUBLE_PRECISION, left, 22, nright%s, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)
!      call MPI_RECV(nleft%s, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%sp, 1, MPI_DOUBLE_PRECISION, left, 22, nright%sp, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%s2p, 1, MPI_DOUBLE_PRECISION, left, 22, nright%s2p, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%s3p, 1, MPI_DOUBLE_PRECISION, left, 22, nright%s3p, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%o, 1, MPI_DOUBLE_PRECISION, left, 22, nright%o, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%op, 1, MPI_DOUBLE_PRECISION, left, 22, nright%op, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%o2p, 1, MPI_DOUBLE_PRECISION, left, 22, nright%o2p, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(n%elec, 1, MPI_DOUBLE_PRECISION, left, 22, nright%elec, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)


!      call MPI_SEND(n%elecHot, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
!      call MPI_RECV(nleft%elecHot, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

!!!!!!!!!!!!!!!!!!!!ENERGY!!!!!!!!!!!!!!!!!!!!!!
      call MPI_SENDRECV(nrg%sp, 1, MPI_DOUBLE_PRECISION, right, 22, nTleft%sp, 1, MPI_DOUBLE_PRECISION, left, 22,& 
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%s2p, 1, MPI_DOUBLE_PRECISION, right, 22, nTleft%s2p, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%s3p, 1, MPI_DOUBLE_PRECISION, right, 22, nTleft%s3p, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%op, 1, MPI_DOUBLE_PRECISION, right, 22, nTleft%op, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%o2p, 1, MPI_DOUBLE_PRECISION, right, 22, nTleft%o2p, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%elec, 1, MPI_DOUBLE_PRECISION, right, 22, nTleft%elec, 1, MPI_DOUBLE_PRECISION, left, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%sp, 1, MPI_DOUBLE_PRECISION, left, 22, nTright%sp, 1, MPI_DOUBLE_PRECISION, right, 22,& 
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%s2p, 1, MPI_DOUBLE_PRECISION, left, 22, nTright%s2p, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%s3p, 1, MPI_DOUBLE_PRECISION, left, 22, nTright%s3p, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%op, 1, MPI_DOUBLE_PRECISION, left, 22, nTright%op, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%o2p, 1, MPI_DOUBLE_PRECISION, left, 22, nTright%o2p, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)

      call MPI_SENDRECV(nrg%elec, 1, MPI_DOUBLE_PRECISION, left, 22, nTright%elec, 1, MPI_DOUBLE_PRECISION, right, 22,&
           MPI_COMM_WORLD, stat, ierr)


!      call MPI_SEND(nrg%elecHot, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
!      call MPI_RECV(nTleft%elecHot, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

    endif
  end do
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine GetAzNeighbors

subroutine GetRadLeftNeighbors(n, nrg, h)
  type(height)        ::h
  type(density)       ::n
  type(nT)            ::nrg
  real                ::rado, nT_exp
  integer             ::inside, outside, i

  nT_exp=5.0/3.0

! ALGORITHM
! Send to outside
! Receive from inside
! Send to inside X           upwind only needs inside
! Receive from outside X
!!!!!!!!!!!!!!!!!!!!DENSITY!!!!!!!!!!!!!!!!!!!!!
 
  outside = mype+LNG_GRID
  inside  = mype-LNG_GRID

  if(inside.ge.0) call MPI_SEND(rdist, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(rado, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(n%sp*rdist/rado, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nin%sp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(n%s2p*rdist/rado, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nin%s2p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(n%s3p*rdist/rado, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nin%s3p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(n%op*rdist/rado, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nin%op, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(n%o2p*rdist/rado, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nin%o2p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(n%elec*rdist/rado, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nin%elec, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

!!!!!!!!!!!!!!!!!!!!HEIGHT!!!!!!!!!!!!!!!!!!!!!!
  if(inside.ge.0) call MPI_SEND(h%sp, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(hout%sp, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(h%s2p, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(hout%s2p, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(h%s3p, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(hout%s3p, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(h%op, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(hout%op, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(h%o2p, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(hout%o2p, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(h%elec, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(hout%elec, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, stat, ierr)

!*(rdist/rado)**(2/3)
!!!!!!!!!!!!!!!!!!!!ENERGY!!!!!!!!!!!!!!!!!!!!!!
  if(outside<npes) call MPI_SEND(nrg%sp*&
    ((rdist*h%sp)/(rado*hout%sp))**nT_exp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nTin%sp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(nrg%s2p*&
    ((rdist*h%s2p)/(rado*hout%s2p))**nT_exp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nTin%s2p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(nrg%s3p*&
    ((rdist*h%s3p)/(rado*hout%s3p))**nT_exp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nTin%s3p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(nrg%op*&
    ((rdist*h%op)/(rado*hout%op))**nT_exp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nTin%op, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(nrg%o2p*&
    ((rdist*h%o2p)/(rado*hout%o2p))**nT_exp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nTin%o2p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(nrg%elec*&
    ((rdist*h%elec)/(rado*hout%elec))**nT_exp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(nTin%elec, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine GetRadLeftNeighbors

subroutine GetRadRightNeighbors(n, nrg, h)
  type(density)       ::n
  type(height)        ::h
  type(nT)            ::nrg
  real                ::radi, nT_exp
  integer             ::inside, outside, i

  nt_exp=5.0/3.0

! ALGORITHM
! Send to inside
! Receive from outside
!!!!!!!!!!!!!!!!!!!!DENSITY!!!!!!!!!!!!!!!!!!!!!
 
  outside = mype+LNG_GRID
  inside  = mype-LNG_GRID

  if(outside<npes) call MPI_SEND(rdist, 1, MPI_REAL, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(radi, 1, MPI_REAL, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(n%sp*rdist/radi, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nout%sp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(n%s2p*rdist/radi, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nout%s2p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(n%s3p*rdist/radi, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nout%s3p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(n%op*rdist/radi, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nout%op, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(n%o2p*rdist/radi, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nout%o2p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(n%elec*rdist/radi, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nout%elec, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

!!!!!!!!!!!!!!!!!!!!HEIGHT!!!!!!!!!!!!!!!!!!!!!!
  if(outside<npes) call MPI_SEND(h%sp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(hin%sp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(h%s2p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(hin%s2p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(h%s3p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(hin%s3p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(h%op, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(hin%op, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(h%o2p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(hin%o2p, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)

  if(outside<npes) call MPI_SEND(h%elec, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, ierr)
  if(inside.ge.0) call MPI_RECV(hin%elec, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, stat, ierr)
!*(rdist/radi)**(2/3)
!!!!!!!!!!!!!!!!!!!!ENERGY!!!!!!!!!!!!!!!!!!!!!!
  if(inside.ge.0) call MPI_SEND(nrg%sp*&
    ((rdist*h%sp)/(radi*hin%sp))**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nTout%sp, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(nrg%s2p*&
    ((rdist*h%s2p)/(radi*hin%s2p))**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
!  if(inside.ge.0) call MPI_SEND(nrg%s2p*(rdist/radi)**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nTout%s2p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(nrg%s3p*&
    ((rdist*h%s3p)/(radi*hin%s3p))**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
!  if(inside.ge.0) call MPI_SEND(nrg%s3p*(rdist/radi)**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nTout%s3p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(nrg%op*&
    ((rdist*h%op)/(radi*hin%op))**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
!  if(inside.ge.0) call MPI_SEND(nrg%op*(rdist/radi)**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nTout%op, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(nrg%o2p*&
    ((rdist*h%o2p)/(radi*hin%o2p))**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
!  if(inside.ge.0) call MPI_SEND(nrg%o2p*(rdist/radi)**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nTout%o2p, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)

  if(inside.ge.0) call MPI_SEND(nrg%elec*&
    ((rdist*h%elec)/(radi*hin%elec))**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
!  if(inside.ge.0) call MPI_SEND(nrg%elec*(rdist/radi)**nT_exp, 1, MPI_DOUBLE_PRECISION, inside, 22, MPI_COMM_WORLD, ierr)
  if(outside<npes) call MPI_RECV(nTout%elec, 1, MPI_DOUBLE_PRECISION, outside, 22, MPI_COMM_WORLD, stat, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine GetRadRightNeighbors

END PROGRAM Onebox

