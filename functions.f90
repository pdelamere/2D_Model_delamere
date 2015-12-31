MODULE FUNCTIONS
 
  USE REACTIONS
  USE varTypes
  USE GLOBAL
  USE FITC
  USE FITR
  USE dielectronic
  USE DEBUG
  USE PARALLELVARIABLES

  IMPLICIT NONE


  CONTAINS

  SUBROUTINE cm3_reactions(r_ind2, r_dep2, h2, n2, ft_int2, zoff)

    type(height)   ::h2
    type(density)  ::n2
    type(r_dep)    ::r_dep2
    type(r_ind)    ::r_ind2
    type(ft_int)   ::ft_int2
    real           ::zoff

    mass_loading(mype+1)=0.0
    
    call ion_neutral_reactions(h2, zoff)
    call ion_ion_reactions(h2)
    call electron_neutral_reactions(n2)
    call electron_ion_reactions(n2)
    call exchange_reactions(ft_int2, r_ind2, n2, h2)
    call impact_ionization(ft_int2, r_dep2, r_ind2, n2, h2)
    call recombination(ft_int2, r_dep2, n2)

  END SUBROUTINE cm3_reactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ion_neutral_reactions(h1, z0)
    type(height)   :: h1
    real           :: z0

    s_sp  = ion_neutral(h1%sp,  h1%s, z0)
    s_s2p = ion_neutral(h1%s2p, h1%s, z0)
    s_s3p = ion_neutral(h1%s3p, h1%s, z0)
    s_op  = ion_neutral(h1%op,  h1%s, z0)
    s_o2p = ion_neutral(h1%o2p, h1%s, z0)

    o_sp  = ion_neutral(h1%sp,  h1%o, z0)
    o_s2p = ion_neutral(h1%s2p, h1%o, z0)
    o_s3p = ion_neutral(h1%s3p, h1%o, z0)
    o_op  = ion_neutral(h1%op,  h1%o, z0)
    o_o2p = ion_neutral(h1%o2p, h1%o, z0)

  END SUBROUTINE ion_neutral_reactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ion_ion_reactions(h1)
    type(height)   :: h1

    sp_sp   = ion_ion(h1%sp , h1%sp )
    sp_s2p  = ion_ion(h1%sp , h1%s2p)
    sp_s3p  = ion_ion(h1%sp , h1%s3p)
    sp_op   = ion_ion(h1%sp , h1%op )
    sp_o2p  = ion_ion(h1%sp , h1%o2p)

    s2p_s2p = ion_ion(h1%s2p , h1%s2p)
    s2p_s3p = ion_ion(h1%s2p , h1%s3p)
    s2p_op  = ion_ion(h1%s2p , h1%op )
    s2p_o2p = ion_ion(h1%s2p , h1%o2p)

    s3p_s3p = ion_ion(h1%s3p , h1%s3p)
    s3p_op  = ion_ion(h1%s3p , h1%op )
    s3p_o2p = ion_ion(h1%s3p , h1%o2p)

    op_op   = ion_ion(h1%op , h1%op )
    op_o2p  = ion_ion(h1%op , h1%o2p)

    o2p_o2p = ion_ion(h1%o2p , h1%o2p)

  END SUBROUTINE ion_ion_reactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL FUNCTION ion_neutral(hi, hn, z0)
    REAL              ::z0
    DOUBLE PRECISION  ::hi, hn
    REAL              ::a, b, c
    a = (hi*hi+hn*hn)/(hi*hi*hn*hn)
    b = 2.0 * z0 / (hn * hn)
    c = (z0*z0)/(hn*hn)
    ion_neutral = sqrt(1.0/a) * exp((b*b - 4.0*a*c)/(4.0*a))
  END FUNCTION ion_neutral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL FUNCTION ion_ion(ha, hb)
    DOUBLE PRECISION  ::ha, hb
    ion_ion = ha * hb /(sqrt(ha*ha + hb*hb))

  END FUNCTION ion_ion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE electron_neutral_reactions(n1)
    type(density)  :: n1

    e_s = (1.0 * n1%sp  * s_sp  + &
           2.0 * n1%s2p * s_s2p + & 
           3.0 * n1%s3p * s_s3p + & 
           1.0 * n1%op  * s_op  + & 
           2.0 * n1%o2p * s_o2p   ) 

    e_o = (1.0 * n1%sp  * o_sp  + &
           2.0 * n1%s2p * o_s2p + & 
           3.0 * n1%s3p * o_s3p + & 
           1.0 * n1%op  * o_op  + & 
           2.0 * n1%o2p * o_o2p   ) 

  END SUBROUTINE electron_neutral_reactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE electron_ion_reactions(n1)
    type(density)  :: n1
    
    e_sp  = (1.0 * n1%sp  * sp_sp   + &
             2.0 * n1%s2p * sp_s2p  + & 
             3.0 * n1%s3p * sp_s3p  + & 
             1.0 * n1%op  * sp_op   + & 
             2.0 * n1%o2p * sp_o2p    ) 
  
    e_s2p = (1.0 * n1%sp  * sp_s2p  + &
             2.0 * n1%s2p * s2p_s2p + & 
             3.0 * n1%s3p * s2p_s3p + & 
             1.0 * n1%op  * s2p_op  + & 
             2.0 * n1%o2p * s2p_o2p   ) 

    e_s3p = (1.0 * n1%sp  * sp_s3p  + &
             2.0 * n1%s2p * s2p_s3p + & 
             3.0 * n1%s3p * s3p_s3p + & 
             1.0 * n1%op  * s3p_op  + & 
             2.0 * n1%o2p * s3p_o2p   ) 

    e_op  = (1.0 * n1%sp  * sp_op   + &
             2.0 * n1%s2p * s2p_op  + & 
             3.0 * n1%s3p * s3p_op  + & 
             1.0 * n1%op  * op_op   + & 
             2.0 * n1%o2p * op_o2p    ) 

    e_o2p = (1.0 * n1%sp  * sp_o2p  + &
             2.0 * n1%s2p * s2p_o2p + & 
             3.0 * n1%s3p * s3p_o2p + & 
             1.0 * n1%op  * op_o2p  + & 
             2.0 * n1%o2p * o2p_o2p   ) 

  END SUBROUTINE electron_ion_reactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE exchange_reactions(ft_int1, r_ind1, n1, h)
    type(ft_int)   :: ft_int1
    type(r_ind)   :: r_ind1
    type(density)  :: n1
    type(height)   ::h
    integer        :: i

    ft_int1%cx(1)  = rootpi * r_ind1%cx(1)  * n1%sp  * n1%s2p  * sp_s2p
    ft_int1%cx(2)  = rootpi * r_ind1%cx(2)  * n1%s   * n1%sp   * s_sp
    ft_int1%cx(3)  = rootpi * r_ind1%cx(3)  * n1%s   * n1%s2p  * s_s2p
    ft_int1%cx(4)  = rootpi * r_ind1%cx(4)  * n1%s   * n1%s2p  * s_s2p
    ft_int1%cx(5)  = rootpi * r_ind1%cx(5)  * n1%s   * n1%s3p  * s_s3p
    ft_int1%cx(6)  = rootpi * r_ind1%cx(6)  * n1%o   * n1%op   * o_op
    ft_int1%cx(7)  = rootpi * r_ind1%cx(7)  * n1%o   * n1%o2p  * o_o2p
    ft_int1%cx(8)  = rootpi * r_ind1%cx(8)  * n1%o   * n1%o2p  * o_o2p
    ft_int1%cx(9)  = rootpi * r_ind1%cx(9)  * n1%o   * n1%sp   * o_sp
    ft_int1%cx(10) = rootpi * r_ind1%cx(10) * n1%s   * n1%op   * s_op
    ft_int1%cx(11) = rootpi * r_ind1%cx(11) * n1%s   * n1%o2p  * s_o2p
    ft_int1%cx(12) = rootpi * r_ind1%cx(12) * n1%s   * n1%o2p  * s_o2p
    ft_int1%cx(13) = rootpi * r_ind1%cx(13) * n1%o   * n1%s2p  * o_s2p
    ft_int1%cx(14) = rootpi * r_ind1%cx(14) * n1%o2p * n1%sp   * sp_o2p
    ft_int1%cx(15) = rootpi * r_ind1%cx(15) * n1%o   * n1%s3p  * o_s3p
    ft_int1%cx(16) = rootpi * r_ind1%cx(16) * n1%o2p * n1%s2p  * s2p_o2p
    ft_int1%cx(17) = rootpi * r_ind1%cx(17) * n1%s3p * n1%sp   * sp_s3p
!    if(mype .eq. 0) print *, r_ind1%cx(1), sp_s2p
!    if(mype .eq. 0) print *, r_ind1%cx(10), s_op
  END SUBROUTINE exchange_reactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE impact_ionization(ft_int1, r_dep1, r_ind1, n1, h1)
    type(ft_int)   :: ft_int1
    type(r_dep)    :: r_dep1
    type(r_ind)    :: r_ind1
    type(density)  :: n1
    type(height)   :: h1

    ft_int1%is   =  rootpi * r_dep1%is   * n1%s   * e_s  
    ft_int1%isp  =  rootpi * r_dep1%isp  * n1%sp  * e_sp  
    ft_int1%is2p =  rootpi * r_dep1%is2p * n1%s2p * e_s2p 
    ft_int1%is3p =  rootpi * r_dep1%is3p * n1%s3p * e_s3p  

    ft_int1%io   =  rootpi * r_dep1%io   * n1%o   * e_o  
    ft_int1%iop  =  rootpi * r_dep1%iop  * n1%op  * e_op  
    ft_int1%io2p =  rootpi * r_dep1%io2p * n1%o2p * e_o2p  

    ft_int1%ish   =  rootpi * r_ind1%ish   * n1%s   * n1%elecHot * h1%s
    ft_int1%isph  =  rootpi * r_ind1%isph  * n1%sp  * n1%elecHot * h1%sp
    ft_int1%is2ph =  rootpi * r_ind1%is2ph * n1%s2p * n1%elecHot * h1%s2p
    ft_int1%is3ph =  rootpi * r_ind1%is3ph * n1%s3p * n1%elecHot * h1%s3p

    ft_int1%ioh   =  rootpi * r_ind1%ioh   * n1%o   * n1%elecHot * h1%o
    ft_int1%ioph  =  rootpi * r_ind1%ioph  * n1%op  * n1%elecHot * h1%op
    ft_int1%io2ph =  rootpi * r_ind1%io2ph * n1%o2p * n1%elecHot * h1%o2p
!    if( rdist .lt. 7.0) print *, "R Dep", r_dep1%is, r_dep1%io, r_ind1%ish
  END SUBROUTINE impact_ionization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE recombination(ft_int1, r_dep1, n1)
    type(ft_int)   :: ft_int1
    type(r_dep)    :: r_dep1
    type(density)  :: n1

    ft_int1%rsp  = rootpi * r_dep1%rsp  * n1%sp  * e_sp
    ft_int1%rs2p = rootpi * r_dep1%rs2p * n1%s2p * e_s2p
    ft_int1%rs3p = rootpi * r_dep1%rs3p * n1%s3p * e_s3p

    ft_int1%rop  = rootpi * r_dep1%rop  * n1%op  * e_op
    ft_int1%ro2p = rootpi * r_dep1%ro2p * n1%o2p * e_o2p

  END SUBROUTINE recombination
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_scale_heights(h, T, n)
    type(height)          ::h
    type(temp)            ::T
    type(density)         ::n

    double precision      ::v, v_corot, omega, ms, mo, sc_const

    v_corot = 2.0 * PI * rdist * Rj / (9.925 *3600)
    v = v_corot - v_ion + lag_amp * cos((lon3-lag_phase)*dTOr) !v_ion has been used instead of lag_const
    omega = v/(rdist*Rj)

    ms = 32.0 * mp
    mo = 16.0 * mp
 
    sc_const = 1.6e-19

    h%sp  = sqrt(2.0*T%sp *sc_const*(1.0+1.0 *T%elec/T%sp )/(3.0*ms))/(1000.0*omega)
    h%s2p = sqrt(2.0*T%s2p *sc_const*(1.0+2.0 *T%elec/T%s2p )/(3.0*ms))/(1000.0*omega)
    h%s3p = sqrt(2.0*T%s3p *sc_const*(1.0+3.0 *T%elec/T%s3p )/(3.0*ms))/(1000.0*omega)
!    h%s4p = sqrt(2.0*T%s4p *sc_const*(1+4 *T%elec/T%s4p )/(3.0*ms))/(1000.0*omega)
    h%op  = sqrt(2.0*T%op *sc_const*(1.0+1.0 *T%elec/T%op )/(3.0*mo))/(1000.0*omega)
    h%o2p = sqrt(2.0*T%o2p *sc_const*(1.0+2.0 *T%elec/T%o2p )/(3.0*mo))/(1000.0*omega)

    h%elec=(n%sp*h%sp + 2.0*n%s2p*h%s2p + 3.0*n%s3p*h%s3p + n%op*h%op + 2.0*n%o2p*h%o2p) &
           /(n%sp + 2.0*n%s2p + 3.0*n%s3p + n%op + 2.0*n%o2p)
   ! print *, h%op, rdist

  END SUBROUTINE get_scale_heights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lat_distribution(n, h, lat)
    type(density)   ::n
    type(height)    ::h
    type(lat_dist)  ::lat
    integer         ::i

    do i=1, LAT_SIZE
      lat%sp(i)  = n%sp *exp(-(lat%z(i)/h%sp)**2)
      lat%s2p(i) = n%s2p *exp(-(lat%z(i)/h%s2p)**2)
      lat%s3p(i) = n%s3p *exp(-(lat%z(i)/h%s3p)**2)
!      lat%s4p(i) = n%s4p *exp(-(lat%z(i)/h%s4p)**2)
      lat%op(i)  = n%op *exp(-(lat%z(i)/h%op)**2)
      lat%o2p(i) = n%o2p *exp(-(lat%z(i)/h%o2p)**2)

      lat%elec(i) = lat%sp(i) + lat%op(i) + 2.0*(lat%s2p(i)+lat%o2p(i)) + 3.0*lat%s3p(i)! + 4*lat%s4p(i) !fix s4pa
!      if(mype .eq. 0) print *, lat%elec(i),  lat%sp(i),lat%s2p(i),lat%s3p(i),lat%op(i),lat%o2p(i)
    end do


  END SUBROUTINE lat_distribution 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lat_print(lat) !for debugging purposes
    double precision  ::lat(LAT_SIZE)
    integer         ::i

    do i = 1, LAT_SIZE
      print *, lat(i)
    end do

  END SUBROUTINE lat_print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dependent_rates(dep, ind, n, T, h)
!    USE INPUTS
    type(r_dep)       ::dep
    type(r_ind)       ::ind
    type(density)     ::n
    type(temp)        ::T
    type(height)      ::h

    double precision  ::tkelv, kTOev
    real              ::logtemp, frac_s, frac_o, iondens, tau

! NEVER LOSE THIS
! http://www.pa.uky.edu/~verner/atom.html
! cfit.f77 and rrfit.f77 origin
! FROM DIMA VERNER

    kTOev= 8.617385e-05
    tkelv = T%elec/kTOev
    logtemp = log10(tkelv)

    call cfit(16,16, tkelv, dep%is)   
    call cfit(16,15, tkelv, dep%isp)   
    call cfit(16,14, tkelv, dep%is2p)   
    call cfit(16,13, tkelv, dep%is3p)   
    call cfit(8,8, tkelv, dep%io)   
    call cfit(8,7, tkelv, dep%iop)   
    call cfit(8,6, tkelv, dep%io2p)

    call rrfit(16,16, tkelv, dep%rsp)  
!    call rrfit(16,14, tkelv, dep%rs2p)  
!    call rrfit(16,13, tkelv, dep%rs3p)  
!    call rrfit(8,7, tkelv, dep%rop)  
!    call rrfit(8,6, tkelv, dep%ro2p) 

    frac_s = 1.0/(1.0 + ind%o_to_s)
    frac_o = ind%o_to_s/(1.0+ind%o_to_s)

    dep%rsp =dep%rsp+dielectronic_rate(16,16,T%elec) 
    dep%rs2p=S_chart_interpolate(logtemp, 1)
    dep%rs3p=S_chart_interpolate(logtemp, 2)

    dep%rop =O_chart_interpolate(logtemp, 1)
    dep%ro2p=O_chart_interpolate(logtemp, 2)

    if (trans_type) then
      iondens = n%sp+n%s2p+n%s3p+n%op+n%o2p
    endif

!    tau = tau0 * (net_source0/net_source)**trans_exp
    
!    dep%transport = 1.0/(tau*86400.0)
!    dep%transport = 1.0/(transport*86400.0)
    h%s=(Rj/2.0)!*(rdist/6.0)
    h%o=(Rj/2.0)!*(rdist/6.0) 
    ind%S_production = frac_s*net_source!/(ROOTPI*h%s*(1e5))
!    if( mype .eq. 0 ) print *, frac_s, net_source 
    ind%O_production = frac_o*net_source!/(ROOTPI*h%o*(1e5))
  END SUBROUTINE dependent_rates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE independent_rates(ind, T, h)
    type(r_ind)       ::ind
    type(temp)        ::T
    type(height)      ::h
    double precision  ::tkelv, kTOev

    kTOev= 8.617385e-05
    tkelv=T%elecHot/kTOev

    call cfit(16,16, tkelv, ind%ish)
    call cfit(16,15, tkelv, ind%isph)
    call cfit(16,14, tkelv, ind%is2ph)
    call cfit(16,13, tkelv, ind%is3ph)
    call cfit(8,8, tkelv, ind%ioh)
    call cfit(8,7, tkelv, ind%ioph)
    call cfit(8,6, tkelv, ind%io2ph)

    !ind%cx=[1.2e-8,2.4e-8,3.0e-10,7.8e-9,1.32e-8,1.32e-8,5.2e-10   &
    ind%cx=[8.1e-9,2.4e-8,3.0e-10,7.8e-9,1.32e-8,1.32e-8,5.2e-10   &
           ,5.4e-9,6.0e-11,3.12e-9,2.34e-8,1.62e-8,2.28e-9,1.38e-9 &
           ,1.92e-8,9.0e-10,3.6e-10]

  END SUBROUTINE independent_rates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function az_loss(v, val)
  real                ::v
  double precision    ::val, source

  az_loss = val * dt * v * LNG_GRID / torus_circumference

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_SEND(az_loss, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(source, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  az_loss = az_loss - source

end function az_loss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_s(n,h,ind,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_ind)       ::ind
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    integer           ::i
    real              ::s,l !source and loss
    
    s = ind%S_production
!    print *, "should be ~10^-4", ind%S_production
    l = (ft%is + ft%ish) 
    do i=2, 5
      l=l+ft%cx(i)
    end do

    do i=10, 12
      l=l+ft%cx(i)
    end do
!    if( mype .eq. 0 ) then
!      print *, "SOURCES S"
!      print *, "Sulfur Production      : ", ind%S_production
!      print *, "LOSSES S"
!      print *, "Sulfur Impact ionization      : ", ft%is/l
!      print *, "Sulfur Impact ionization(hot) : ", ft%ish/l
!      print *, "Charge Exchange #2            : ", ft%cx(2)/l
!      print *, "Charge Exchange #3            : ", ft%cx(3)/l
!      print *, "Charge Exchange #4            : ", ft%cx(4)/l
!      print *, "Charge Exchange #5            : ", ft%cx(5)/l
!      print *, "Charge Exchange #10           : ", ft%cx(10)/l
!      print *, "Charge Exchange #11           : ", ft%cx(11)/l
!      print *, "Charge Exchange #12           : ", ft%cx(12)/l
!    end if

    if( Euler ) l=l+az_loss(v_neutral,n%s)
    l=l/(ROOTPI*h%s)
    mass_loading(mype+1)=l*32.0*mp+mass_loading(mype+1)
    F_s = s-l

  end function F_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_sp(n,h,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    real              ::s,l !source and loss

    s= ft%is + ft%ish + ft%rs2p + (ft%cx(3) + ft%cx(3)) + ft%cx(5) &
     + ft%cx(10) + ft%cx(11) + ft%cx(13)
    l= ft%isp + ft%isph + ft%rsp + ft%cx(9) + ft%cx(14) + ft%cx(17) 
    
!    if( mype .eq. 0 ) then
!      print *, "SOURCES S+"
!      print *, "Sulfur Impact ionization      : ", ft%is/s
!      print *, "Sulfur Impact ionization(hot) : ", ft%ish/s
!      print *, "Sulfur ++ Recombination       : ", ft%rs2p/s
!      print *, "Charge Exchange #3            : ", 2.0*ft%cx(3)/s
!      print *, "Charge Exchange #5            : ", ft%cx(5)/s
!      print *, "Charge Exchange #10           : ", ft%cx(10)/s
!      print *, "Charge Exchange #11           : ", ft%cx(11)/s
!      print *, "Charge Exchange #13           : ", ft%cx(13)/s
!      print *, "LOSSES S+"
!      print *, "Sulfur Plus Impact ionization      : ", ft%isp/l
!      print *, "Sulfur Plus Impact ionization(hot) : ", ft%isph/l
!      print *, "Sulfur Plus Recombination          : ", ft%rsp/l
!      print *, "Charge Exchange #9                 : ", ft%cx(9)/l
!      print *, "Charge Exchange #14                : ", ft%cx(14)/l
!      print *, "Charge Exchange #17                : ", ft%cx(17)/l
!    end if

    F_sp=(s-l)/(ROOTPI*h%sp) 
!    if ( mype .eq. 1 .or. mype .eq. 6 .or. mype .eq. 7 ) then
!      print *, "SOURCE ::: ", s, ":::LOSS:::", l, ft%cx(3) , mype
!    endif

  end function F_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_s2p(n,h,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    real              ::s,l !source and loss

    s= ft%isp + ft%isph + ft%rs3p +ft%cx(5) + ft%cx(12) + ft%cx(14) &
     + ft%cx(15) + ft%cx(17) + ft%cx(17)
    l= ft%is2p + ft%is2ph + ft%rs2p + ft%cx(3) + ft%cx(13) + ft%cx(16)
    F_s2p=(s-l)/(ROOTPI*h%s2p) 

  end function F_s2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_s3p(n,h,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    real              ::s,l !source and loss

    s= ft%is2p + ft%is2ph + ft%cx(16)
    l= ft%is3p + ft%is3ph + ft%rs3p + ft%cx(5) + ft%cx(15) + ft%cx(17) 
!    if( Euler ) l=l+az_loss(v_ion,n%s3p)

    F_s3p = (s-l)/(ROOTPI*h%s3p) 

  end function F_s3p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_o(n,h,ind,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_ind)       ::ind
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    real              ::s,l !source and loss

    s = ind%O_production
    l= ft%io + ft%ioh + ft%cx(6) + ft%cx(7) + ft%cx(8) + ft%cx(9) &
     + ft%cx(13) + ft%cx(15) 
    l= l/(ROOTPI*h%o)
    mass_loading(mype+1)=mass_loading(mype+1)+l*16.0*mp
    F_o= (s-l)

  end function F_o
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_op(n,h,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    real              ::s,l !source and loss

    s= ft%io + ft%ioh + ft%ro2p + ft%cx(7) + ft%cx(7)              &
     + ft%cx(9) + ft%cx(11) + ft%cx(12) + ft%cx(13) + ft%cx(14)    &
     + ft%cx(15) + ft%cx(16)
    l= ft%iop +ft%ioph + ft%rop + ft%cx(10) 

    F_op= (s-l)/(ROOTPI*h%op)

  end function F_op
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function F_o2p(n,h,dep,ft)
    type(density)     ::n
    type(height)      ::h
    type(r_dep)       ::dep
    type(ft_int)      ::ft
    real              ::s,l !source and loss

    s= ft%iop + ft%ioph
    l= ft%io2p + ft%io2ph + ft%ro2p + ft%cx(7) + ft%cx(11)  &
     + ft%cx(12) + ft%cx(14) + ft%cx(16)

    F_o2p= (s-l)/(ROOTPI*h%o2p)

  end function F_o2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function Tpu(m, L)
    real    ::m, L, mp, vrel !mass as an atomic number and (L radial distance in Rj?)

    mp= 1.67262158e-27   !in kg
    
    vrel= (1e3)*(12.57*L - (42.0/sqrt(L)))  !km/s
      
    Tpu=(m*mp*vrel*vrel)/(3.0*(1.60217646e-19))  !eV

  end function Tpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!From NRL plasma formulary rev 2007 pg. 34
  real function lambda_ee(ne, Te)
    double precision  ::ne, Te
!    real              ::ne, Te
    lambda_ee = 23.5-.5*log(ne*sqrt(Te**(-5.0)))-sqrt((1.0e-5) + (((log(Te)-2.0)**2.0)/16.0))
  end function lambda_ee
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!From NRL plasma formulary rev 2007 pg. 34
  real function lambda_ei(ne,Te,ni,Ti,zi,mui)
    double precision  ::ne, ni
    double precision  ::Te, Ti
    real              ::zi, mui, cond1, cond2

     cond1=Ti*me/(mp*mui)
     cond2=10.0*zi*zi

     if(cond1<Te .and. Te<cond2) then
       lambda_ei= 23 - log(ne*zi*zi/(Te**3))/2   
!       print *, "CASE: 1"  , "    Lambda:", lambda_ei
     elseif(cond1<cond2 .and. cond2<Te) then
       lambda_ei= 24 - log(sqrt(ne)/Te)     
!       print *, "CASE: 2"  , "    Lambda:", lambda_ei
     elseif(Te<(cond1*zi)) then
       lambda_ei= 30 - log((ni*(zi**4))/((Ti**3)*mui*mui))/2     
!       print *, "CASE: 3"  , "    Lambda:", lambda_ei
     else
!       print *, "Coulomb case not found: lambda_ei in functions.f90" !fix
       lambda_ei= 23 - log(ne*zi*zi/(Te**3))/2  !default to first case
     endif

  end function lambda_ei
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!From NRL plasma formulary rev 2007 pg. 34
  real function lambda_ii(n1,T1,z1,mu1,n2,T2,z2,mu2)
    double precision  ::n1, n2
    double precision  ::T1, T2
    real              ::z1, mu1, z2, mu2, arg ! <- like a pirate 
!    real              ::n1, n2, T1, T2
    
    arg= z1*z2*(mu1+mu2)/(mu1*T2+mu2*T1)
    arg= arg*arg*((n1*z1*z1/T1)+(n2*z2*z2/T2)) 
    lambda_ii= 23 - log(arg)/2

  end function lambda_ii
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function nu_ei(ne,Te,ni,Ti,z,mu)
    double precision  ::ni(LAT_SIZE), Ti, ne(LAT_SIZE), Te
!    real              ::ni, Ti, ne, Te
    real              ::mu, z !mu is atomic mass and z is charge number
    double precision  ::lambda, kgTOg, mi, top, bottom
    double precision  ::nu_arr(LAT_SIZE), top_tot, bot_tot
    integer           ::i

    mi=mu*mpg
    bottom=(sqrt((mi*Te+meg*Ti)**3))   
    top_tot=0.0
    bot_tot=0.0

    do i=1, LAT_SIZE
      lambda= lambda_ei(ne(i), Te, ni(i), Ti, z, mu)
      top= (1.8e-19)*sqrt(meg*mi)*z*z*ne(i)*ni(i)*lambda
      nu_arr(i)=top
    enddo

    do i=1, LAT_SIZE
      top_tot=top_tot+(nu_arr(i)*ne(i)*ni(i))
      bot_tot=bot_tot+(ni(i)*ne(i))
    end do

    nu_ei=top_tot/(bot_tot*bottom)
  end function nu_ei
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function nu_ehi(neh,Te,ni,Ti,z,mu)
    double precision  ::ni(LAT_SIZE), Ti, neh, Te
!    real              ::ni, Ti, ne, Te
    real              ::mu, z !mu is atomic mass and z is charge number
    double precision  ::lambda, kgTOg, mi, top, bottom
    double precision  ::nu_arr(LAT_SIZE), top_tot, bot_tot
    integer           ::i

    mi=mu*mpg
    bottom=(sqrt((mi*Te+meg*Ti)**3))   
    top_tot=0.0
    bot_tot=0.0

    do i=1, LAT_SIZE
      lambda= lambda_ei(neh, Te, ni(i), Ti, z, mu)
      top= (1.8e-19)*sqrt(meg*mi)*z*z*neh*ni(i)*lambda
      nu_arr(i)=top
    enddo

    do i=1, LAT_SIZE
      top_tot=top_tot+(nu_arr(i)*neh*ni(i))
      bot_tot=bot_tot+(ni(i)*neh)
    end do

    nu_ehi=top_tot/(bot_tot*bottom)
  end function nu_ehi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real function nu_ee(ne, Te, neh, Teh)
    double precision  ::ne(LAT_SIZE), Te, neh, Teh
!    real              ::n1, T1, n2, T2
    double precision  ::lambda, nu_arr(LAT_SIZE),top,bot, const
    integer           ::i
    
    const=1.8e-19

    do i=1, LAT_SIZE
      lambda= lambda_ee(ne(i), Te)    
      nu_arr(i)= const*meg*ne(i)*neh*lambda/sqrt((meg*(Te+Teh))**3)
    end do

    top=0.0
    bot=0.0
    do i=1, LAT_SIZE
      top=top+nu_arr(i)*ne(i)*neh
      bot=bot+ne(i)*neh
    end do

    nu_ee=top/bot
!    print *, rdist, nu_ee
!    if( mype .eq. 0) print *, "ne: ", ne
!    if( mype .eq. 0) print *, "Te: ", Te
!    if( mype .eq. 0) print *, "neh: ", neh
!    if( mype .eq. 0) print *, "Teh: ", Teh
!    if( mype .eq. 0) print *, "Lambda: ", Lambda
!    if( mype .eq. 0) print *, "nu_arr: ", nu_arr
!    if( mype .eq. 0) print *, "nu_ee: ", nu_ee
!    if( mype .eq. 0) print *, ""
  end function nu_ee  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function nu_ii(n1, T1, z1, mu1, n2, T2, z2, mu2)
    double precision  ::n1(LAT_SIZE), n2(LAT_SIZE)
    double precision  ::T1, T2
!    real              ::n1, n2, T1, T2
    real              ::z1, mu1, z2, mu2
    double precision  ::m1, m2, lambda
    double precision  ::nu_arr(LAT_SIZE), top, bot
    integer           ::i

    m1=mu1*mpg
    m2=mu2*mpg
    do i=1, LAT_SIZE
      lambda=lambda_ii(n1(i), T1, z1, mu1, n2(i), T2, z2, mu2)
      nu_arr(i)= (1.8e-19)*sqrt(m1*m2)*z1*z1*z2*z2*n1(i)*n2(i)*lambda/sqrt((m1*T1+m2*T2)**3)
    end do

    top=0.0
    bot=0.0
    do i=1, LAT_SIZE
      top=top+nu_arr(i)*n1(i)*n2(i)
      bot=bot+n1(i)*n2(i)
    end do

    nu_ii=top/bot

  end function nu_ii
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function nu_ii_dens(n1, T1, z1, mu1, n2, T2, z2, mu2)
    double precision  ::n1, n2
    double precision  ::T1, T2
!    real              ::n1, n2, T1, T2
    real              ::z1, mu1, z2, mu2
    double precision  ::m1, m2, lambda
    double precision  ::nu_arr, top, bot
    integer           ::i

    m1=mu1*mp*1000.0
    m2=mu2*mp*1000.0
    lambda=lambda_ii(n1, T1, z1, mu1, n2, T2, z2, mu2)
    nu_arr= (1.8e-19)*sqrt(m1*m2)*z1*z1*z2*z2*n1*n2*lambda/sqrt((m1*T1+m2*T2)**3)

    top=nu_arr*n1*n2
    bot=n1*n2

    nu_ii_dens=top/bot

  end function nu_ii_dens
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_temp(n,nrg, T)
    type(density)     ::n
    type(nT)          ::nrg
    type(temp)        ::T

    T%sp = nrg%sp / n%sp
    T%s2p = nrg%s2p / n%s2p
    T%s3p = nrg%s3p / n%s3p
!    T%s4p = nrg%s4p / n%s4p
    T%op = nrg%op / n%op
    T%o2p = nrg%o2p / n%o2p
    T%elec = nrg%elec / n%elec

  end subroutine update_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function ft_rad(lat, T, ind, h)
    type(lat_dist)    ::lat
    type(height)      ::h
    type(temp)        ::T
    type(r_ind)       ::ind
    real              ::rad_sp(LAT_SIZE), rad_s2p(LAT_SIZE), rad_s3p(LAT_SIZE)
    real              ::rad_op(LAT_SIZE), rad_o2p(LAT_SIZE), rad_tot(LAT_SIZE)
    real              ::x, yarr(LAT_SIZE), elec_tot
    integer           ::i 
    real              ::const,stuff

    ft_rad=0.0
    elec_tot=0.0
    const=3.69897 !log10(5000)
    x =(100.0*(1.0+ LOG10(T%elec))/const)
!    x = T%elec
    do i=1, LAT_SIZE
      yarr(i) = 100.0*log10(lat%elec(i))/const
!      yarr(i) = lat%elec(i)
    end do

    do i=1, LAT_SIZE
!      rad_sp(i)=interpolate(ind, ind%emis_sp, x, yarr(i), T%elec, lat%elec(i), i)*lat%sp(i)
!      rad_s2p(i)=interpolate(ind, ind%emis_s2p, x, yarr(i), T%elec, lat%elec(i), i)*lat%s2p(i)
!      rad_s3p(i)=interpolate(ind, ind%emis_s3p, x, yarr(i), T%elec, lat%elec(i), i)*lat%s3p(i)
!      rad_op(i)=interpolate(ind, ind%emis_op, x, yarr(i), T%elec, lat%elec(i), i)*lat%op(i)
!      rad_o2p(i)=interpolate(ind, ind%emis_o2p, x, yarr(i), T%elec, lat%elec(i), i)*lat%o2p(i)
      rad_sp(i)=bilinearInterpolate(ind, ind%emis_sp, x, yarr(i), T%elec, lat%elec(i), i)*lat%sp(i)
      rad_s2p(i)=bilinearInterpolate(ind, ind%emis_s2p, x, yarr(i), T%elec, lat%elec(i), i)*lat%s2p(i)
      rad_s3p(i)=bilinearInterpolate(ind, ind%emis_s3p, x, yarr(i), T%elec, lat%elec(i), i)*lat%s3p(i)
      rad_op(i)=bilinearInterpolate(ind, ind%emis_op, x, yarr(i), T%elec, lat%elec(i), i)*lat%op(i)
      rad_o2p(i)=bilinearInterpolate(ind, ind%emis_o2p, x, yarr(i), T%elec, lat%elec(i), i)*lat%o2p(i)
!      rad_sp(i)=interpolate_II(ind, ind%emis_sp, x, yarr(i))*lat%sp(i)
!      if(mype .eq. 0) print *, i, rad_sp(i), lat%z(i), T%sp, T%elec, lat%elec(i), lat%sp(i) 
!      rad_s2p(i)=interpolate_II(ind, ind%emis_s2p, x, yarr(i))*lat%s2p(i)
!      rad_s3p(i)=interpolate_II(ind, ind%emis_s3p, x, yarr(i))*lat%s3p(i)
!      rad_op(i)=interpolate_II(ind, ind%emis_op, x, yarr(i))*lat%op(i)
!      rad_o2p(i)=interpolate_II(ind, ind%emis_o2p, x, yarr(i))*lat%o2p(i)
      rad_tot(i)=rad_sp(i) + rad_s2p(i) + rad_s3p(i) + rad_op(i) +rad_o2p(i)
!    print *, i, "<====>", rad_tot(i)
      ft_rad=ft_rad+rad_tot(i)*lat%elec(i)
      elec_tot=elec_tot + lat%elec(i)
!      if(mype .eq. 0) print *, i, rad_tot(i)
    end do
 
    ft_rad=ft_rad/elec_tot
!    if(mype .eq. 0) print *, ft_rad, h%elec/Rj, h%sp/Rj, h%s2p/Rj, h%s3p/Rj, h%op/Rj, h%o2p/Rj
!    if(mype .eq. 0) print *, ""
  end function ft_rad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine emisPrint(xmin, ymin, xmax, ymax, emis)
    integer           ::i, j, xmin, ymin, xmax, ymax
    real              ::emis(EMIS_SIZE, EMIS_SIZE)
    
    do i=ymin, ymax
      do j=xmin, xmax
        write(*,"(e11.4, 3x)", advance='no') emis(j,i)
      end do
      print *, ""
    end do

  end subroutine emisPrint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function bilinearInterpolate(ind, table, tl, nl, t, n, i)
    type(r_ind)       ::ind
    real              ::tgt, tlt, ngt, nlt, table(EMIS_SIZE, EMIS_SIZE) !temperature greater than, less than, same for density....etc.
    double precision  ::t, n
    real              ::tl, nl
    integer           ::ti, ni, i !x and y indices for table search
    real              ::dt(3), dn(3), A(5)

   ti=INT(tl)
   ni=INT(nl)

   if(ti > EMIS_SIZE) ti=EMIS_SIZE-1
   if(ni > EMIS_SIZE) ni=EMIS_SIZE-1
   if(ti < 1) ti=1
   if(ni < 1) ni=1

   tgt=ind%emis_temp(ti+1)
   tlt=ind%emis_temp(ti)
   dt(1)=ind%emis_temp(ti+1)-ind%emis_temp(ti)
   ngt=ind%emis_dens(ni+1)
   nlt=ind%emis_dens(ni)
   dn(1)=ind%emis_dens(ni+1)-ind%emis_dens(ni)

   A(1)=dt(1)*dn(1) 

   dt(2)=t-tlt
   dt(3)=tgt-t
   dn(2)=n-nlt
   dn(3)=ngt-n

   A(2)=dt(2)*dn(2)/A(1)
   A(3)=dt(3)*dn(2)/A(1)
   A(4)=dt(2)*dn(3)/A(1)
   A(5)=dt(3)*dn(3)/A(1)

   if ( t > ind%emis_temp(EMIS_SIZE) .or. n > ind%emis_dens(EMIS_SIZE) .or. t < ind%emis_temp(1) .or. n < ind%emis_dens(1) ) then
     if ( .not. HUSH .and. HUSH) then !THIS VARIABLE IS SET IN THE DEBUG MODULE
       print *, "interpolate function out of bounds(called by ft_rad)"
       print *, "t : ",t
       print *, "n : ",n
     endif
   endif

   bilinearInterpolate=table(ti,ni)*A(5)+table(ti+1,ni)*A(4)+table(ti,ni+1)*A(3)+table(ti+1,ni+1)*A(2)
!   if(mype .eq. 0 .and. i .eq. 1) print *, A(2) + A(3) + A(4) + A(5)
  end function bilinearInterpolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function interpolate_II(ind, table, t, n)
    type(r_ind)       ::ind
    real              ::t, n, gt, lt
    real              ::table(EMIS_SIZE, EMIS_SIZE)
    integer           ::ti, ni

    ti=NINT(t)
    ni=NINT(n)

    gt=ind%emis_temp(1)
    do while ( ti<=EMIS_SIZE .and. gt<t )
      ti=ti+1
      lt=gt
      gt=ind%emis_temp(ti)
    end do

    gt=ind%emis_dens(1)
    do while ( ni<=EMIS_SIZE .and. gt<n )
      ni=ni+1
      lt=gt
      gt=ind%emis_dens(ni)
    end do

    if( ni .le. 0 ) ni = 1
    if( ti .le. 0 ) ti = 1
    if( ni .ge. EMIS_SIZE ) ni = EMIS_SIZE
    if( ti .ge. EMIS_SIZE ) ti = EMIS_SIZE

    interpolate_II=table(ti,ni)

  end function interpolate_II
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function interpolate(ind, table, tl, nl, t, n, i)
    type(r_ind)       ::ind
    double precision  ::t, n
    integer           ::ti, ni !x and y indices for table search
    real              ::tgt, tlt, ngt, nlt, table(EMIS_SIZE, EMIS_SIZE) !temperature greater than, less than, same for density....etc.
    double precision  ::tslope, nslope
    real              ::tl, nl
    integer           ::i

   ti=INT(tl)
   ni=INT(nl)

   if(ti > EMIS_SIZE) ti=EMIS_SIZE-1
   if(ni > EMIS_SIZE) ni=EMIS_SIZE-1
   if(ti < 1) ti=1
   if(ni < 1) ni=1

   tgt=ind%emis_temp(ti+1)
   tlt=ind%emis_temp(ti)
   ngt=ind%emis_dens(ni+1)
   nlt=ind%emis_dens(ni)

   if ( t > ind%emis_temp(EMIS_SIZE) .or. n > ind%emis_dens(EMIS_SIZE) .or. t < ind%emis_temp(1) .or. n < ind%emis_dens(1) ) then
     if ( .not. HUSH .and. HUSH) then !THIS VARIABLE IS SET IN THE DEBUG MODULE
       print *, "interpolate function out of bounds(called by ft_rad)"
       print *, "t : ",t
       print *, "n : ",n
     endif
   endif

    tslope=((table(ti+1,ni+1)-table(ti, ni+1))*(1.0-((ngt-n)/(ngt-nlt))) &
          +(table(ti+1,ni)-table(ti,ni))*(1.0-((n-nlt)/(ngt-nlt))))/(tgt-tlt)
    nslope=((table(ti+1,ni+1)-table(ti+1 ,(ni)))*(1.0-((tgt-t)/(tgt-tlt))) &
          +(table(ti,ni+1)-table(ti,ni))*(1.0-((t-tlt)/(tgt-tlt))))/(ngt-nlt)

    if ( (tgt-t >= t-tlt) .and. (ngt-n) >= (n-nlt) ) then
      interpolate= table(ti,ni) + (n-nlt)*nslope + (t-tlt)*tslope
    elseif ( (t-tlt) > (tgt-t) .and. (ngt-n) >= (n-nlt) ) then
      interpolate= table(ti+1,ni) + (n-nlt)*nslope - (tgt-t)*tslope
    elseif ( (tgt-t) >= (t-tlt) .and. (n-nlt) > (ngt-n) ) then
      interpolate= table(ti,ni+1) - (ngt-n)*nslope + (t-tlt)*tslope
    elseif ( (t-tlt) > (tgt-t) .and. (t-tlt) > (tgt-t) ) then
      interpolate= table(ti+1,ni+1) - (ngt-n)*nslope - (tlt-t)*tslope
    endif

  end function interpolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function EF_elec(n, T, h, ind, dep, lat, v, ft)
    type(density)      ::n
    type(temp)         ::T
    type(height)       ::h
    type(r_ind)        ::ind
    type(r_dep)        ::dep
    type(lat_dist)     ::lat
    type(nu)           ::v
    type(ft_int)       ::ft

    real               ::Teq!uila (Equilibrium temp)
    real               ::ipo, ipop, ipo2p, ips, ipsp, ips2p, ips3p!ionization potentials
    real               ::NeTot, ipPer, ipTot, rad, en_ip

    Teq= v%sp_elec *(T%sp  - T%elec)&
       + v%s2p_elec*(T%s2p - T%elec)&
       + v%s3p_elec*(T%s3p - T%elec)&
!       + v%s4p_elec*(T%s4p - T%elec)&
       + v%op_elec *(T%op  - T%elec)&
       + v%o2p_elec*(T%o2p - T%elec)&
       + v%elec_elecHot*(T%elecHot - T%elec)

    call lat_distribution(n, h, lat)

    rad = ft_rad(lat,T, ind, h) 

    ipo   = 13.61806
    ipop  = 35.1173
    ipo2p = 54.9355
    ips   = 10.36001
    ipsp  = 23.3379
    ips2p = 34.79
    ips3p = 47.222

    en_ip= ft%io * ipo + ft%iop * ipop + ft%io2p * ipo2p &
         + ft%is * ips + ft%isp * ipsp + ft%is2p * ips2p + ft%is3p * ips3p

    NeTot= (ROOTPI / (1.0 - n%protons)) &
          * (1.0 * n%sp  * h%sp          &
          +  2.0 * n%s2p * h%s2p         &
          +  3.0 * n%s3p * h%s3p         &
          +  1.0 * n%op  * h%op          &
          +  2.0 * n%o2p * h%o2p         )

    ipPer= en_ip / NeTot

    ipTot= 2.0 * ipPer * n%elec /3.0

    EF_elec= Teq - (2.0 * rad / 3.0)! - ipTot !- (v_r0/dr * n%elec * T%elec) 

!    if(mype .eq. 0) print *, rdist, EF_elec, Teq, (2.0*rad/3.0)!, "ipTot: ", ipTot
 
!    print *, rdist, EF_elec

  end function EF_elec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function EF_sp(n, T, h, ind, dep, v, ft)
    type(density)      ::n
    type(temp)         ::T
    type(r_ind)        ::ind
    type(r_dep)        ::dep
    type(height)       ::h
    type(nu)           ::v
    type(ft_int)       ::ft
    integer            ::i
    real               ::stuff

    double precision   ::Teq !(Equilibrium temp)
    double precision   ::S, L
   
    stuff=ROOTPI*h%sp

    S= T%pu_s * (ft%is + ft%ish + ft%cx(2) + ft%cx(3) + ft%cx(5) &
              +  ft%cx(10) + ft%cx(11))                          &
              + T%s2p  * (ft%rs2p + ft%cx(1) + ft%cx(3) + ft%cx(13))

    L= T%sp * (ft%isp + ft%isph + ft%rsp + ft%cx(1) + ft%cx(2) &
            +  ft%cx(9) + ft%cx(14) + ft%cx(17))                
!            +  dep%transport * n%sp * ROOTPI * h%sp)
!    if( Euler ) L=L+az_loss(v_ion,nrg%op)

    Teq= v%sp_s2p *(T%s2p  - T%sp)&
       + v%sp_s3p*(T%s3p   - T%sp)&
!       + v%sp_s4p*(T%s4p  - T%sp)&
       + v%sp_op *(T%op    - T%sp)&
       + v%sp_o2p*(T%o2p   - T%sp)&
       + v%sp_elec*(T%elec - T%sp)&
       + v%sp_elecHot*(T%elecHot - T%sp)
!  if (rdist .lt. 7.0) print *, rdist, S/(ROOTPI*h%sp), L/(ROOTPI*h%sp)
!    if (mype .eq. 0) print *, rdist  
!    if (mype .eq. 0) print *, "source", S/stuff
!    if (mype .eq. 0) print *, "is  ", ft%is*T%pu_s/stuff  
!    if (mype .eq. 0) print *, "ish ", ft%ish *T%pu_s /stuff 
!    if (mype .eq. 0) print *, "cx5 ", ft%cx(5)*T%pu_s  /stuff 
!    if (mype .eq. 0) print *, "cx11", ft%cx(11)*T%pu_s  /stuff 
!    if (mype .eq. 0) print *, "cx2 ", ft%cx(2)*T%pu_s  /stuff 
!    if (mype .eq. 0) print *, "cx3 ", ft%cx(3)*T%pu_s  /stuff 
!    if (mype .eq. 0) print *, "cx10", ft%cx(10)*T%s2p  /stuff 
!    if (mype .eq. 0) print *, "cx3 ", ft%cx(3)*T%s2p /stuff 
!    if (mype .eq. 0) print *, "cx13", ft%cx(13)*T%s2p  /stuff 
!    if (mype .eq. 0) print *, "cx1 ", ft%cx(1)*T%s2p  /stuff 
!    if (mype .eq. 0) print *, "rs2p", ft%rs2p*T%s2p  /stuff
!    if (mype .eq. 0) print *, "+++++++++++++++++++++++++++++++++++"
!    if (mype .eq. 0) print *, "loss", L/stuff
!    if (mype .eq. 0) print *, "isp ", ft%isp *T%sp  /stuff 
!    if (mype .eq. 0) print *, "isph", ft%isph*T%sp  /stuff 
!    if (mype .eq. 0) print *, "rsp ", ft%rsp*T%sp  /stuff 
!    if (mype .eq. 0) print *, "cx1 ", ft%cx(1 )*T%sp  /stuff 
!    if (mype .eq. 0) print *, "cx2 ", ft%cx(2 )*T%sp  /stuff 
!    if (mype .eq. 0) print *, "cx9 ", ft%cx(9 )*T%sp  /stuff 
!    if (mype .eq. 0) print *, "cx14", ft%cx(14)*T%sp  /stuff 
!    if (mype .eq. 0) print *, "cx17", ft%cx(17)*T%sp  /stuff 
!    if (mype .eq. 0) print *, "Teq", Teq
!    if (mype .eq. 0) print *, "||||||||||||||||||||||||||||||||||||"

  EF_sp = Teq + (S - L)/(ROOTPI * h%sp)
  miscOutput=EF_sp

  end function EF_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function EF_s2p(n, T, h, ind, dep, v, ft)
    type(density)      ::n
    type(temp)         ::T
    type(r_ind)        ::ind
    type(r_dep)        ::dep
    type(height)       ::h
    type(nu)           ::v
    type(ft_int)       ::ft

    double precision   ::Teq !(Equilibrium temp)
    double precision   ::S, L
   
    S= T%pu_s * (ft%cx(4) + ft%cx(12))                                &
     + T%sp   * (ft%isp + ft%isph + ft%cx(1) + ft%cx(14) + ft%cx(17)) &
     + T%s3p  * (ft%rs3p + ft%cx(5) + ft%cx(15) + ft%cx(17))

    L= T%s2p * (ft%is2p + ft%is2ph + ft%rs2p + ft%cx(1) + ft%cx(3) &
            +  ft%cx(4) + ft%cx(13) + ft%cx(16))                
!            +  dep%transport * n%s2p * ROOTPI * h%s2p)

    Teq= v%sp_s2p *(T%sp    - T%s2p)&
       + v%s2p_s3p*(T%s3p   - T%s2p)&
!       + v%s2p_s4p*(T%s4p  - T%s2p)&
       + v%s2p_op *(T%op    - T%s2p)&
       + v%s2p_o2p*(T%o2p   - T%s2p)&
       + v%s2p_elec*(T%elec - T%s2p)&
       + v%s2p_elecHot*(T%elecHot - T%s2p)

  EF_s2p = Teq + (S - L)/(ROOTPI * h%s2p)

  end function EF_s2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function EF_s3p(n, T, h, ind, dep, v, ft)
    type(density)      ::n
    type(temp)         ::T
    type(r_ind)        ::ind
    type(r_dep)        ::dep
    type(height)       ::h
    type(nu)           ::v
    type(ft_int)       ::ft

    double precision   ::Teq !(Equilibrium temp)
    double precision   ::S, L
   
    S= T%s2p * (ft%is2p + ft%is2ph + ft%cx(16))

    L= T%s3p * (ft%is3p + ft%is3ph + ft%rs3p + ft%cx(15) + ft%cx(17) )
!            +  dep%transport * n%s3p * ROOTPI * h%s3p)

    Teq= v%sp_s3p *(T%sp    - T%s3p)&
       + v%s2p_s3p*(T%s2p   - T%s3p)&
!       + v%s3p_s4p*(T%s4p  - T%s3p)&
       + v%s3p_op *(T%op    - T%s3p)&
       + v%s3p_o2p*(T%o2p   - T%s3p)&
       + v%s3p_elec*(T%elec - T%s3p)&
       + v%s3p_elecHot*(T%elecHot - T%s3p)

  EF_s3p = Teq + (S - L)/(ROOTPI * h%s3p)

  end function EF_s3p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function EF_op(n, T, h, ind, dep, v, ft)
    type(density)      ::n
    type(temp)         ::T
    type(r_ind)        ::ind
    type(r_dep)        ::dep
    type(height)       ::h
    type(nu)           ::v
    type(ft_int)       ::ft

    double precision   ::Teq !(Equilibrium temp)
    double precision   ::S, L
   
    S= T%pu_o * (ft%io + ft%ioh + ft%cx(6) + ft%cx(7) + ft%cx(9) &
              +  ft%cx(13) + ft%cx(15))                          &
     + T%o2p  * (ft%ro2p + ft%cx(7) + ft%cx(11) + ft%cx(12) + ft%cx(14) + ft%cx(16))

    L= T%op * (ft%iop + ft%ioph + ft%rop + ft%cx(6) + ft%cx(10) )
!            +  dep%transport * n%op * ROOTPI * h%op)

    Teq= v%sp_op *(T%sp    - T%op)&
       + v%s2p_op*(T%s2p   - T%op)&
       + v%s3p_op*(T%s3p   - T%op)&
!       + v%s4p_op *(T%4p  - T%op)&
       + v%op_o2p*(T%o2p   - T%op)&
       + v%op_elec*(T%elec - T%op)&
       + v%op_elecHot*(T%elecHot - T%op)

  EF_op = Teq + (S - L)/(ROOTPI * h%op)

  end function EF_op
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function EF_o2p(n, T, h, ind, dep, v, ft)
    type(density)      ::n
    type(temp)         ::T
    type(r_ind)        ::ind
    type(r_dep)        ::dep
    type(height)       ::h
    type(nu)           ::v
    type(ft_int)       ::ft

    double precision   ::Teq !(Equilibrium temp)
    double precision   ::S, L
   
    S= T%pu_o * ft%cx(8) &
     + T%op   * (ft%iop + ft%ioph)

    L= T%o2p * (ft%io2p + ft%io2ph + ft%ro2p + ft%cx(7) + ft%cx(8) &
             +  ft%cx(11) + ft%cx(12) + ft%cx(14) + ft%cx(16))      
 !            +  dep%transport * n%o2p * ROOTPI * h%o2p)


    Teq= v%sp_o2p *(T%sp    - T%o2p)&
       + v%s2p_o2p*(T%s2p   - T%o2p)&
       + v%s3p_o2p*(T%s3p   - T%o2p)&
!       + v%s4p_o2p*(T%4p   - T%o2p)&
       + v%op_o2p*(T%op     - T%o2p)&
       + v%o2p_elec*(T%elec - T%o2p)&
       + v%o2p_elecHot*(T%elecHot - T%o2p)

  EF_o2p = Teq + (S - L)/(ROOTPI * h%o2p)

  end function EF_o2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energyBudget(n, h, T, dep, ind, ft, lat, v, nrgy)
  type(density)       ::n
  type(height)        ::h
  type(temp)          ::T
  type(r_dep)         ::dep
  type(r_ind)         ::ind
  type(ft_int)        ::ft
  type(lat_dist)      ::lat
  type(nu)            ::v
  type(energy)        ::nrgy 

  real                ::Teq_sp, Teq_s2p, Teq_s3p
  real                ::Teq_op, Teq_o2p, Teq_elh

  real                ::ion_sp, ionh_sp, cx_sp, fast_sp, trans_sp    
  real                ::ion_s2p, ionh_s2p, cx_s2p, fast_s2p, trans_s2p    
  real                ::ion_s3p, ionh_s3p, cx_s3p, fast_s3p, trans_s3p    
  real                ::ion_op, ionh_op, cx_op, fast_op, trans_op    
  real                ::ion_o2p, ionh_o2p, cx_o2p, fast_o2p, trans_o2p    

  Teq_sp = v%sp_s2p* (T%s2p- T%sp) + &
           v%sp_s3p* (T%s3p- T%sp) + &
           v%sp_op*  (T%op - T%sp) + &
           v%sp_o2p* (T%o2p- T%sp) + &
           v%sp_elec*(T%elec-T%sp)

  ion_sp   = T%pu_s*ft%is /(ROOTPI*h%sp)
  ionh_sp  = T%pu_s*ft%ish/(ROOTPI*h%sp)
  cx_sp    = T%pu_s*(ft%cx(2) + ft%cx(3) + ft%cx(5) + ft%cx(10) + ft%cx(11))/(ROOTPI*h%sp)
  fast_sp  = T%sp*(ft%cx(2) + ft%cx(9))/(ROOTPI*h%sp)
  trans_sp = T%sp*(v_r0/dr*n%sp)

  Teq_s2p = v%sp_s2p*  (T%sp-  T%s2p) + &
            v%s2p_s3p* (T%s3p- T%s2p) + &
            v%s2p_op*  (T%op-  T%s2p) + &
            v%s2p_o2p* (T%o2p- T%s2p) + &
            v%s2p_elec*(T%elec-T%s2p)

  ion_s2p   = T%sp*ft%isp /(ROOTPI*h%s2p)
  ionh_s2p  = T%sp*ft%isph/(ROOTPI*h%s2p)
  cx_s2p    = T%pu_s*(ft%cx(4) + ft%cx(12))/(ROOTPI*h%s2p)
  fast_s2p  = T%s2p*(ft%cx(4))/(ROOTPI*h%s2p)
  trans_s2p = T%s2p*(v_r0/dr*n%s2p)

    Teq_s3p = v%sp_s3p*  (T%sp - T%s3p) + &
              v%s2p_s3p* (T%s2p -T%s3p) + &
              v%s3p_op*  (T%op - T%s3p) + &
              v%s3p_o2p* (T%o2p -T%s3p) + &
              v%s3p_elec*(T%elec-T%s3p)

  ion_s3p   = T%s2p*ft%is2p /(ROOTPI*h%s3p)
  ionh_s3p  = T%s2p*ft%is2ph/(ROOTPI*h%s3p)
  !cx_s3p    = T%s2p*ft%cx_k15/(ROOTPI*h%s3p)
  trans_s3p = T%s3p*(v_r0/dr*n%s3p)

  Teq_op = v%sp_op*  (T%sp-  T%op) + &
           v%s2p_op* (T%s2p- T%op) + &
           v%s3p_op* (T%s3p- T%op) + &
           v%op_o2p* (T%o2p- T%op) + &
           v%op_elec*(T%elec-T%op)

  ion_op   = T%pu_o*ft%io /(ROOTPI*h%op)
  ionh_op  = T%pu_o*ft%ioh/(ROOTPI*h%op)
  cx_op    = T%pu_o*(ft%cx(6) + ft%cx(7) + ft%cx(9) + ft%cx(13) + ft%cx(15))/(ROOTPI*h%op)
  fast_op  = T%op*(ft%cx(7) + ft%cx(10))/(ROOTPI*h%op)
  trans_op = T%op*(v_r0/dr*n%op)

  Teq_o2p = v%sp_o2p*  (T%sp-  T%o2p) + &
            v%s2p_o2p* (T%s2p- T%o2p) + &
            v%s3p_o2p* (T%s3p- T%o2p) + &
            v%op_o2p*  (T%op-  T%o2p) + &
            v%o2p_elec*(T%elec-T%o2p)

  ion_o2p   = T%op*ft%iop /(ROOTPI*h%o2p)
  ionh_o2p  = T%op*ft%ioph/(ROOTPI*h%o2p)
  cx_o2p    = T%pu_o*(ft%cx(8))/(ROOTPI*h%o2p)
  fast_o2p  = T%o2p*(ft%cx(8))/(ROOTPI*h%o2p)
  trans_o2p = T%o2p*(v_r0/dr*n%o2p)

  Teq_elh = v%sp_elecHot *(T%elecHot - T%sp ) + &
            v%s2p_elecHot*(T%elecHot - T%s2p) + &
            v%s3p_elecHot*(T%elecHot - T%s3p) + &
            v%op_elecHot *(T%elecHot - T%op ) + &
            v%o2p_elecHot*(T%elecHot - T%o2p) + &
            v%elec_elecHot *(T%elecHot - T%elec )

  nrgy%elecHot_eq  = 1.5 * Teq_elh
  nrgy%tot_eq      = 1.5 * (Teq_sp + Teq_s2p + Teq_s3p + Teq_op + Teq_o2p)

  nrgy%s_ion = 1.5 * (ion_sp + ionh_sp + ion_s2p + ionh_s2p + ion_s3p + ionh_s3p)
  nrgy%o_ion = 1.5 * (ion_op + ionh_op + ion_o2p + ionh_o2p)
  nrgy%s_cx  = 1.5 * (cx_sp + cx_s2p) ! + cx_s3p)
  nrgy%o_cx  = 1.5 * (cx_op + cx_o2p)

  nrgy%P_in  = nrgy%s_ion + nrgy%o_ion + nrgy%s_cx + nrgy%o_cx + nrgy%elecHot_eq

  nrgy%Pfast = 1.5 * (fast_sp + fast_s2p + fast_op + fast_o2p)
  nrgy%Puv   = ft_rad(lat, T, ind, h)

  nrgy%Ptrans         = 1.5 * v_r0/dr * (n%sp*T%sp + n%s2p*T%s2p + n%s3p*T%s3p + n%op*T%op + n%o2p*T%o2p + n%elec*T%elec)
  nrgy%Ptrans_elecHot = 1.5 * v_r0/dr * n%elecHot * T%elecHot

  nrgy%P_out = nrgy%puv + nrgy%pfast + nrgy%ptrans + nrgy%ptrans_elecHot


end subroutine energyBudget
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPACER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




END MODULE FUNCTIONS


