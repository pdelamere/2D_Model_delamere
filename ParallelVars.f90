MODULE PARALLELVARIABLES
  USE VARTYPES  
  USE MPI
  USE DIMENSIONS

  integer             ::ierr, mype, npes, lnggrid, radgrid 
  integer             ::stat(MPI_STATUS_SIZE)
  character(len=3)    ::num_char

  type(density)        ::nleft, nright, nin, nout
  type(nT)             ::nTleft, nTright, nTin, nTout
  type(height)         ::hin, hout

  type(density)       ::dens_loss
  type(nT)            ::nrg_loss
  real                ::Io_loc, v_Io, torus_circumference, v_ion, v_neutral, dx, dr
  real                ::sys4_loc, v_sys4, sys4_amp, sys3_amp
  real                ::numerical_s, numerical_c_ion, numerical_c_neutral
  real                ::numerical_c_r, numerical_c_rin, numerical_c_rout
  real                ::mass_loading(NUMPES), all_loading(NUMPES), ave_loading, min_loading, max_loading
  real                ::miscOutput 

  parameter(v_Io=57.0)                !in km/s
  parameter(v_neutral=57.0)           !in km/s 57 is default
! parameter(v_ion=1.0)                !in km/s
  parameter(v_sys4=1.05)              !in km/s (1.05 corresponds to 12.5 degrees/day at 6 Rj)
  parameter(sys4_amp=0.40)            !amplification factor for moving hot electron populatioN
  parameter(sys3_amp=0.30)            !amplification for stationary (sys 3) hot electrons

END MODULE


