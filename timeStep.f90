MODULE TIMESTEP

  USE FUNCTIONS
  USE GLOBAL

  IMPLICIT NONE

  CONTAINS

  subroutine updateNu(v, lat, T, elecHot)
!This subroutine calculates the new nu values for each time step. 
!this subroutine uses the nu_* functions from functions.f90
!StepNu is a near copy. These two should be consolidated if possible. 
    type(nu)          ::v
    type(lat_dist)    ::lat
    type(temp)        ::T
    double precision  ::elecHot
    integer           ::i

    v%sp_s2p= nu_ii(lat%sp, T%sp, 1.0, 32.0, lat%s2p, T%s2p, 2.0, 32.0)
    v%sp_s3p= nu_ii(lat%sp, T%sp, 1.0, 32.0, lat%s3p, T%s3p, 3.0, 32.0)
!    v%sp_s4p= nu_ii(lat%sp, T%sp, 1.0, 32.0, lat%s4p, T%s4p, 4.0, 32.0)   
    v%sp_op = nu_ii(lat%sp, T%sp, 1.0, 32.0, lat%op , T%op, 1.0, 16.0 )
    v%sp_o2p= nu_ii(lat%sp, T%sp, 1.0, 32.0, lat%o2p, T%o2p, 2.0, 16.0)

    v%s2p_s3p= nu_ii(lat%s2p, T%s2p, 2.0, 32.0, lat%s3p, T%s3p, 3.0, 32.0)
!    v%s2p_s4p= nu_ii(lat%s2p, T%s2p, 2.0, 32.0, lat%s4p, T%s4p, 4.0, 32.0)   
    v%s2p_op = nu_ii(lat%s2p, T%s2p, 2.0, 32.0, lat%op , T%op , 1.0, 16.0)
    v%s2p_o2p= nu_ii(lat%s2p, T%s2p, 2.0, 32.0, lat%o2p, T%o2p, 2.0, 16.0)

!    v%s3p_s4p= nu_ii(lat%s3p, T%s3p, 3.0, 32.0, lat%s4p, T%s4p, 4.0, 32.0)   
    v%s3p_op = nu_ii(lat%s3p, T%s3p, 3.0, 32.0, lat%op , T%op , 1.0, 16.0)
    v%s3p_o2p= nu_ii(lat%s3p, T%s3p, 3.0, 32.0, lat%o2p, T%o2p, 2.0, 16.0)

!    v%s4p_op = nu_ii(lat%s4p, T%s4p, 4.0, 32.0, lat%op , T%op , 1.0, 16.0)   
!    v%s4p_o2p= nu_ii(lat%s4p, T%s4p, 4.0, 32.0, lat%o2p, T%o2p, 2.0, 16.0)   

    v%op_o2p= nu_ii(lat%op, T%op, 1.0, 16.0, lat%o2p, T%o2p, 2.0, 16.0)

    v%sp_elec = nu_ei(lat%elec, T%elec, lat%sp , T%sp , 1.0, 32.0)
    v%s2p_elec= nu_ei(lat%elec, T%elec, lat%s2p, T%s2p, 2.0, 32.0)
    v%s3p_elec= nu_ei(lat%elec, T%elec, lat%s3p, T%s3p, 3.0, 32.0)
!    v%s4p_elec= nu_ei(lat%elec, T%elec, lat%s4p, T%s4p, 4.0, 32.0)   
    v%op_elec = nu_ei(lat%elec, T%elec, lat%op , T%op , 1.0, 16.0)
    v%o2p_elec= nu_ei(lat%elec, T%elec, lat%o2p, T%o2p, 2.0, 16.0)

!    v%sp_elecHot = nu_ei(lat%elecHot, T%elecHot, lat%sp , T%sp , 1.0, 32.0)
!    v%s2p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s2p, T%s2p, 2.0, 32.0)
!    v%s3p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s3p, T%s3p, 3.0, 32.0)
!    v%s4p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s4p, T%s4p, 4.0, 32.0)   
!    v%op_elecHot = nu_ei(lat%elecHot, T%elecHot, lat%op , T%op , 1.0, 16.0)
!    v%o2p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%o2p, T%o2p, 2.0, 16.0)
!
!    v%elec_elecHot= nu_ee(lat%elec, T%elec, lat%elecHot, T%elecHot)

    v%sp_elecHot = nu_ehi(elecHot, T%elecHot, lat%sp , T%sp , 1.0, 32.0)
    v%s2p_elecHot= nu_ehi(elecHot, T%elecHot, lat%s2p, T%s2p, 2.0, 32.0)
    v%s3p_elecHot= nu_ehi(elecHot, T%elecHot, lat%s3p, T%s3p, 3.0, 32.0)
!    v%s4p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s4p, T%s4p, 4.0, 32.0)   
    v%op_elecHot = nu_ehi(elecHot, T%elecHot, lat%op , T%op , 1.0, 16.0)
    v%o2p_elecHot= nu_ehi(elecHot, T%elecHot, lat%o2p, T%o2p, 2.0, 16.0)

    v%elec_elecHot= nu_ee(lat%elec, T%elec, elecHot, T%elecHot) 
    
!    print *, rdist, v%elec_elecHot, fehot_const

  end subroutine updateNu

  subroutine stepNu(v, n, lat, T)
!Probably redundant
    type(nu)          ::v
    type(density)     ::n
    type(lat_dist)    ::lat
    type(temp)        ::T

    integer           ::i

    v%sp_s2p= nu_ii_dens(n%sp, T%sp, 1.0, 32.0, n%s2p, T%s2p, 2.0, 32.0)
    v%sp_s3p= nu_ii_dens(n%sp, T%sp, 1.0, 32.0, n%s3p, T%s3p, 3.0, 32.0)
!    v%sp_s4p= nu_ii_dens(n%sp, T%sp, 1.0, 32.0, n%s4p, T%s4p, 4.0, 32.0)   
    v%sp_op = nu_ii_dens(n%sp, T%sp, 1.0, 32.0, n%op , T%op, 1.0, 16.0 )
    v%sp_o2p= nu_ii_dens(n%sp, T%sp, 1.0, 32.0, n%o2p, T%o2p, 2.0, 16.0)

    v%s2p_s3p= nu_ii_dens(n%s2p, T%s2p, 2.0, 32.0, n%s3p, T%s3p, 3.0, 32.0)
!    v%s2p_s4p= nu_ii_dens(n%s2p, T%s2p, 2.0, 32.0, n%s4p, T%s4p, 4.0, 32.0)   
    v%s2p_op = nu_ii_dens(n%s2p, T%s2p, 2.0, 32.0, n%op , T%op , 1.0, 16.0)
    v%s2p_o2p= nu_ii_dens(n%s2p, T%s2p, 2.0, 32.0, n%o2p, T%o2p, 2.0, 16.0)

!    v%s3p_s4p= nu_ii_dens(n%s3p, T%s3p, 3.0, 32.0, n%s4p, T%s4p, 4.0, 32.0)   
    v%s3p_op = nu_ii_dens(n%s3p, T%s3p, 3.0, 32.0, n%op , T%op , 1.0, 16.0)
    v%s3p_o2p= nu_ii_dens(n%s3p, T%s3p, 3.0, 32.0, n%o2p, T%o2p, 2.0, 16.0)

!    v%s4p_op = nu_ii_dens(n%s4p, T%s4p, 4.0, 32.0, n%op , T%op , 1.0, 16.0)   
!    v%s4p_o2p= nu_ii_dens(n%s4p, T%s4p, 4.0, 32.0, n%o2p, T%o2p, 2.0, 16.0)   

    v%op_o2p= nu_ii_dens(n%op, T%op, 1.0, 16.0, n%o2p, T%o2p, 2.0, 16.0)

    v%sp_elec = nu_ei(lat%elec, T%elec, lat%sp , T%sp , 1.0, 32.0)
    v%s2p_elec= nu_ei(lat%elec, T%elec, lat%s2p, T%s2p, 2.0, 32.0)
    v%s3p_elec= nu_ei(lat%elec, T%elec, lat%s3p, T%s3p, 3.0, 32.0)
!    v%s4p_elec= nu_ei(lat%elec, T%elec, lat%s4p, T%s4p, 4.0, 32.0)   
    v%op_elec = nu_ei(lat%elec, T%elec, lat%op , T%op , 1.0, 16.0)
    v%o2p_elec= nu_ei(lat%elec, T%elec, lat%o2p, T%o2p, 2.0, 16.0)

    v%sp_elecHot = nu_ei(lat%elecHot, T%elecHot, lat%sp , T%sp , 1.0, 32.0)
    v%s2p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s2p, T%s2p, 2.0, 32.0)
    v%s3p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s3p, T%s3p, 3.0, 32.0)
!    v%s4p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%s4p, T%s4p, 4.0, 32.0)   
    v%op_elecHot = nu_ei(lat%elecHot, T%elecHot, lat%op , T%op , 1.0, 16.0)
    v%o2p_elecHot= nu_ei(lat%elecHot, T%elecHot, lat%o2p, T%o2p, 2.0, 16.0)

    v%elec_elecHot= nu_ee(lat%elec, T%elec, n%elecHot, T%elecHot)

  end subroutine stepNu


  subroutine fluxCorrect(nold, n, h, ind, dep, ft)
!This subroutine handles transport in the torus. I relies on the Flux functions (F_*) located in functions.f90
    type(density)     ::nold, n   
    type(height)      ::h  
    type(r_ind)       ::ind   
    type(r_dep)       ::dep   
    type(ft_int)      ::ft

!EF__ and F__ are kept in the global varible module "global.f90"
    Fs   = F_s(nold,h,ind,dep,ft)
    n%s   = nold%s   + dt * Fs
    call gtzero(n%s)
    Fsp  = F_sp(nold,h,dep,ft)
!    if( mype .eq. 0 ) print *, "r_ind%ish: ", ind%ish
!    if( mype .eq. 0 ) print *, "r_ind%isph: ", ind%isph
    n%sp  = nold%sp  + dt * Fsp
    call gtzero(n%sp)
    Fs2p = F_s2p(nold,h,dep,ft)
    n%s2p = nold%s2p + dt * Fs2p
    call gtzero(n%s2p)
    Fs3p = F_s3p(nold,h,dep,ft)
    n%s3p = nold%s3p + dt * Fs3p
    call gtzero(n%s3p)
!    n%s4p = nold%s4p + dt * Fs4p(nold,h,dep,ft)
!    call gtzero(n%s4p)
    Fo   = F_o(nold,h,ind,dep,ft)
    n%o   = nold%o   + dt * Fo
    call gtzero(n%o)
    Fop  = F_op(nold,h,dep,ft)
    n%op  = nold%op  + dt * Fop
    call gtzero(n%op)
    Fo2p = F_o2p(nold,h,dep,ft)
    n%o2p = nold%o2p + dt * Fo2p
    call gtzero(n%o2p)
    n%elec=(n%sp + 2.0*n%s2p + 3.0*n%s3p + n%op + 2.0*n%o2p)/(1.0-n%protons) 
    n%elecHot=n%fh*n%elec/(1.0-n%fh)

  end subroutine fluxCorrect

  subroutine gtzero(num)
    double precision  ::num
    if (num<0.0) then 
      num=0.0 
    end if
  end subroutine gtzero

  function updateNT(n, T, h, ind, dep, v, ft, lat, nrg)
!Calculates energy values for the next time step
    type(density)    ::n
    type(temp)       ::T
    type(height)     ::h
    type(r_ind)      ::ind
    type(r_dep)      ::dep
    type(nu)         ::v
    type(ft_int)     ::ft
    type(lat_dist)   ::lat
    type(nT)         ::nrg
    type(nT)         ::updateNT
!EF__ and F__ are kept in the global variable module "global.f90"
    EFsp = EF_sp(n, T, h, ind, dep, v, ft)
    updateNT%sp  = nrg%sp  + dt *  EFsp
    call gtzero(updateNT%sp )
    EFs2p = EF_s2p(n, T, h, ind, dep, v, ft)
    updateNT%s2p = nrg%s2p + dt * EFs2p
    call gtzero(updateNT%s2p)
    EFs3p = EF_s3p(n, T, h, ind, dep, v, ft)
    updateNT%s3p = nrg%s3p + dt * EFs3p
    call gtzero(updateNT%s3p)
!    updateNT%s4p = nrg%s4p + dt * EF_s4p(n, T, h, ind, dep, v, ft)
!    call gtzero(updateNT%s4p)
    EFop = EF_op(n, T, h, ind, dep, v, ft)
    updateNT%op  = nrg%op  + dt *  EFop
    call gtzero(updateNT%op )
    EFo2p = EF_o2p(n, T, h, ind, dep, v, ft)
    updateNT%o2p = nrg%o2p + dt * EFo2p
    call gtzero(updateNT%o2p)
    EFelec = EF_elec(n, T, h, ind, dep, lat, v, ft)
    updateNT%elec=nrg%elec + dt * EFelec
    call gtzero(updateNT%elec)
!    print *, rdist, EFelec

  end function updateNT

  subroutine improvedEuler(np, n, n1, T1, h1, ind1, dep1, nu1, ft1, lat1, nTp, nrg)
!Uses an improved Euler approximation to calculate the density and energy for the next time step
    type(density)    ::np, n, n1
    type(temp)       ::T1
    type(height)     ::h1
    type(r_ind)      ::ind1
    type(r_dep)      ::dep1
    type(ft_int)     ::ft1
    type(nu)         ::nu1
    type(lat_dist)   ::lat1
    type(nT)         ::nTp, nrg
    real             ::frac_s, frac_o

!    frac_s = 1.0/(1.0 + ind1%o_to_s)
!    frac_o = ind1%o_to_s/(1.0+ind1%o_to_s)

!    ind1%S_production = frac_s*net_source/(ROOTPI*h1%s*(1e5))
!    ind1%O_production = frac_o*net_source/(ROOTPI*h1%o*(1e5))
!    if (mype .eq. 0) then
!      print *, ind1%S_production, ind1%o_to_s
!    endif

!EF__ and F__ are kept in the global varible module "global.f90"
    np%s   = (n%s   + dt * 0.5 * (Fs   + F_s(n1,h1,ind1,dep1, ft1)))  
    call gtzero(np%s )
    np%sp  = (n%sp  + dt * 0.5 * (Fsp  + F_sp(n1,h1,dep1, ft1)))
    call gtzero(np%sp )
    np%s2p = (n%s2p + dt * 0.5 * (Fs2p + F_s2p(n1,h1,dep1, ft1)))
    call gtzero(np%s2p )
    np%s3p = (n%s3p + dt * 0.5 * (Fs3p + F_s3p(n1,h1,dep1, ft1)))
    call gtzero(np%s3p )
!    np%s4p = (n%s4p + dt * 0.5 * (Fs4p + F_s4p(n1,h1,dep1, ft1)))
!    call gtzero(np% )
    np%o   = (n%o   + dt * 0.5 * (Fo   + F_o(n1,h1,ind1,dep1, ft1)))  
    call gtzero(np%o )
    np%op  = (n%op  + dt * 0.5 * (Fop  + F_op(n1,h1,dep1, ft1))) 
    call gtzero(np%op )
    np%o2p = (n%o2p + dt * 0.5 * (Fo2p + F_o2p(n1,h1,dep1, ft1)))
    call gtzero(np%o2p )

    np%elec = (np%sp + 2.0*np%s2p + 3.0*np%s3p + np%op + 2.0*np%o2p)!/(1.0-np%protons)
    np%elecHot = np%fh/(1.0-np%fh) * np%elec

    nTp%sp  = (nrg%sp  + dt * 0.5 * (EFsp  + EF_sp(n1,T1,h1,ind1,dep1,nu1, ft1))) 
    call gtzero(nTp%sp )
    nTp%s2p = (nrg%s2p + dt * 0.5 * (EFs2p + EF_s2p(n1,T1,h1,ind1,dep1,nu1, ft1)))
    call gtzero(nTp%s2p )
    nTp%s3p = (nrg%s3p + dt * 0.5 * (EFs3p + EF_s3p(n1,T1,h1,ind1,dep1,nu1, ft1)))
    call gtzero(nTp%s3p )
!    nTp%s4p = (nrg%s4p + dt * 0.5 * (EFs4p + EF_s4p(n1,T1,ind1,dep1,nu1, ft1)))
!    call gtzero(nTp%op )
    nTp%op  = (nrg%op  + dt * 0.5 * (EFop  + EF_op(n1,T1,h1,ind1,dep1,nu1, ft1))) 
    call gtzero(nTp%op )
    nTp%o2p = (nrg%o2p + dt * 0.5 * (EFo2p + EF_o2p(n1,T1,h1,ind1,dep1,nu1, ft1)))
    call gtzero(nTp%o2p )
    nTp%elec  = (nrg%elec  + dt * 0.5 * (EFelec  + EF_elec(n1,T1,h1, ind1,dep1,lat1,nu1, ft1)))
    call gtzero(nTp%elec )

!    if(mype .le. 2) print *, rdist, EFelec, T1%elec
!    if(mype .le. 2) print *, rdist, EF_elec(n1,T1,h1, ind1,dep1,lat1,nu1, ft1), nTp%elec/np%elec
!    if(mype .le. 2) print *, n%elec, n1%elec, np%elec
!    if(mype .le. 2) print *, nrg%elec, nTp%elec
!    if(mype .le. 2) print *, ""
 
  end subroutine improvedEuler

  subroutine cm3_latavg_model(n, T, nrg, h, v, n1, T1, nrg1, h1, v1, np, Tp, nTp, ind, dep, dep1, lat, lat1, ft, z)
    type(density)   ::n, n1, np
    type(temp)      ::T, T1, Tp
    type(nT)        ::nrg, nrg1, nTp
    type(height)    ::h, h1
    type(nu)        ::v, v1
    type(r_ind)     ::ind
    type(r_dep)     ::dep, dep1
    type(lat_dist)  ::lat, lat1
    type(ft_int)    ::ft
    real            ::z
    integer         ::i

    call dependent_rates(dep, ind, n, T, h)

    call cm3_reactions(ind, dep, h, n, ft, z)
 
    call lat_distribution(n, h, lat)
    
    call updateNu(v, lat, T, n%elecHot)
    mass_loading(mype+1)=0.0
    call fluxCorrect(n, n1, h, ind, dep, ft)

    do i=1, LAT_SIZE 
      lat%elecHot(i)=n%elecHot
    enddo

    nrg1 = updateNT(n, T, h, ind, dep, v, ft, lat, nrg)

    call update_temp(n1, nrg1, T1)

    call get_scale_heights(h1, T1, n1)

    call dependent_rates(dep1, ind, n1, T1, h1)

    call lat_distribution(n1, h1, lat1)

    call stepNu(v1, n1, lat1, T1)
    mass_loading(mype+1)=0.0
    call improvedEuler(np, n, n1, T1, h1, ind, dep1, v1, ft, lat1, nTp, nrg)
    call update_temp(np, nTp, Tp)
    !if(mype .le. 1) print *, rdist, nTp%elec, Tp%elec
!    print *, rdist, Tp%elecHot, n%fh
!    if(mype .eq. 0) print *, "+++++++++++++++++++++++++++"

!    print *, rdist, EFelec
    call get_scale_heights(h, Tp, np)

!    if(mype .eq. 0) print *, T%sp, T1%sp, Tp%sp

    n = np
    nrg = nTp
    T = Tp

!    if(mype .le. 2) print *, rdist, T%elec
!    call az_transport(n, nrg)

  end subroutine cm3_latavg_model

!  subroutine az_transport(n, nrg)
!    type(density)     ::n
!    type(nT)          ::nrg
!
!    dens_loss%s = n%s * v_neutral * dt /( torus_circumference / LNG_GRID )
!    dens_loss%o = n%o * v_neutral * dt /( torus_circumference / LNG_GRID )
!
!    dens_loss%sp  = n%sp  * v_ion * dt /( torus_circumference / LNG_GRID )
!    dens_loss%s2p = n%s2p * v_ion * dt /( torus_circumference / LNG_GRID )
!    dens_loss%s3p = n%s3p * v_ion * dt /( torus_circumference / LNG_GRID )
!    dens_loss%op  = n%op  * v_ion * dt /( torus_circumference / LNG_GRID )
!    dens_loss%o2p = n%o2p * v_ion * dt /( torus_circumference / LNG_GRID )
!
!    nrg_loss%sp = nrg%sp * v_ion * dt /( torus_circumference / LNG_GRID )
!    nrg_loss%s2p = nrg%s2p * v_ion * dt /( torus_circumference / LNG_GRID )
!    nrg_loss%s3p = nrg%s3p * v_ion * dt /( torus_circumference / LNG_GRID )
!    nrg_loss%op = nrg%op * v_ion * dt /( torus_circumference / LNG_GRID )
!    nrg_loss%o2p = nrg%o2p * v_ion * dt /( torus_circumference / LNG_GRID )
!
!  end subroutine az_transport

END MODULE TIMESTEP

