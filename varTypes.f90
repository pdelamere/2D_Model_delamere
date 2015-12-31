!supporting file for onebox.f90
!sets constants and declares necessary types for onebox.f90

MODULE varTypes

  USE DIMENSIONS

  integer :: EMIS_SIZE, CX_SIZE, LAT_SIZE, REC_ROWS, REC_O, REC_S
  real    ::PI, rootpi, Rj, dTOr, mp, me, mpg, meg, charge
!Simulation grid dimensions (Must multiply to equal npes)
!  parameter(LNG_GRID=1)       !number of longitudinal slices of the torus
!  parameter(RAD_GRID=8)        !number of longitudinal slices of the torus

  parameter(EMIS_SIZE=101)    !size of chianti emmision tables
  parameter(CX_SIZE=17)       !charge exchange reactions
  parameter(LAT_SIZE=31)      !latitudinal grid size
  parameter(REC_ROWS=81)      !Recombination table length
  parameter(REC_O=9)          !Number of O reaction in recombination_O.dat
  parameter(REC_S=3)          !Number of columns in recombination_S.dat
  parameter(PI=3.1415927)     ! PI
  parameter(rootpi=sqrt(PI))  !sqrt(pi) (reduces calculations)
  parameter(Rj=71492.0)       !radius of jupiter in km
  parameter(Rjm=71492000.0)   !radius of jupiter in m
  parameter(dTOr=(PI/180.0))  !conversion from degree to radian
  parameter(mp=1.672621e-27)  !kg
  parameter(mpg=1.672621e-24) !grams 
  parameter(me=9.109381e-31)  !kg
  parameter(meg=9.109381e-28) !grams
  parameter(charge=1.6e-19)   !elementary charge
  parameter(omega=1.76e-4)    !angular velocity of Jupiter (rad/s)
!  parameter()      


!! Make necessary structures to mimic Steffl Onebox model
!!//////////////////////////////////////////////////
  TYPE ::  density  !n in Steffl Onebox model
    double precision  :: sp=0.0, s2p=0.0, s3p=0.0, s4p=0.0, op=0.0, o2p=0.0, elec=0.0, elecHot=0.0, s=0.0, o=0.0 
    double precision  :: fc=0.0, fh=0.0, protons=0.0
  END TYPE density

  TYPE :: nT 
    double precision  :: sp=0.0, s2p=0.0, s3p=0.0, s4p=0.0, op=0.0, o2p=0.0, elec=0.0, elecHot=0.0
  END TYPE nT

  TYPE :: temp !T in Steffl Onebox model
    double precision  :: sp=0.0, s2p=0.0, s3p=0.0, s4p=0.0, op=0.0, o2p=0.0, elec=0.0, elecHot=0.0
    double precision  :: pu_s=0.0, pu_o=0.0
  END TYPE temp

  TYPE :: height !h in Steffl Onebox model
    double precision  :: s=0.0, sp=0.0, s2p=0.0, s3p=0.0, s4p=0.0, o=0.0, op=0.0, o2p=0.0, elec=0.0
  END TYPE height



  TYPE :: nu 
    double precision  :: sp_s2p, sp_s3p, sp_s4p, sp_op, sp_o2p, s2p_s3p, s2p_s4p, s2p_op, s2p_o2p
    double precision  :: s3p_s4p, s3p_op, s3p_o2p, s4p_op, s4p_o2p, op_o2p
    double precision  :: sp_elec, s2p_elec, s3p_elec, s4p_elec, op_elec, o2p_elec
    double precision  :: sp_elecHot, s2p_elecHot, s3p_elecHot, s4p_elecHot, op_elecHot, o2p_elecHot, elec_elecHot
  END TYPE nu

  TYPE :: r_ind
    real              :: ish=0.0, isph=0.0, is2ph=0.0, is3ph=0.0, ioh=0.0, ioph=0.0, io2ph=0.0, s_production=0.0, &
                         o_production=0.0, o2s_spike=0.0, o_to_s=0.0
    real              :: cx(CX_SIZE) !(17 different variables in Steffl onebox) cx stands for charge exchange
    real              :: emis_sp(EMIS_SIZE,EMIS_SIZE), emis_s2p(EMIS_SIZE,EMIS_SIZE), &
                         emis_s3p(EMIS_SIZE,EMIS_SIZE), emis_op(EMIS_SIZE,EMIS_SIZE), &
                         emis_o2p(EMIS_SIZE,EMIS_SIZE)
    real              :: emis_temp(EMIS_SIZE), emis_dens(EMIS_SIZE)!EMIS_SIZE and CX_SIZE are declared at the top as parameters
  END TYPE r_ind

  TYPE :: r_dep
    real              :: is=0.0, isp=0.0, is2p=0.0, is3p=0.0, io=0.0, iop=0.0, io2p=0.0, rsp=0.0, rs2p=0.0,&
                         rs3p=0.0, rop=0.0, ro2p=0.0, transport=0.0
  END TYPE r_dep

  TYPE :: lat_dist
    double precision  :: z(LAT_SIZE), sp(LAT_SIZE), s2p(LAT_SIZE), s3p(LAT_SIZE) !LAT_SIZE declared as parameter above
    double precision  :: s4p(LAT_SIZE), op(LAT_SIZE), o2p(LAT_SIZE), elec(LAT_SIZE), elecHot(LAT_SIZE)
  END TYPE lat_dist
 
  TYPE :: ft_int
    double precision  :: cx(CX_SIZE) !charge exchange array
    double precision  :: is=0.0, isp=0.0, is2p=0.0, is3p=0.0, io=0.0, iop=0.0, io2p=0.0, ish=0.0, isph=0.0, is2ph=0.0,&
                         is3ph=0.0, ioh=0.0, ioph=0.0, io2ph=0.0, rsp=0.0, rs2p=0.0, rs3p=0.0, rop=0.0, ro2p=0.0
  END TYPE ft_int

  TYPE:: energy
    real              ::s_ion=0.0, s_cx=0.0, o_ion=0.0, o_cx=0.0, elecHot_eq=0.0, tot_eq=0.0, P_in=0.0, Puv=0.0, Pfast=0.0 &
                      , Ptrans=0.0, Ptrans_elecHot=0.0, P_out=0.0 
  END TYPE

  TYPE :: reac
    type(lat_dist)    ::lat
    type(density)     :: dens
    type(nT)          :: nrg
    type(temp)        :: tmp
    type(height)      :: h
    type(nu)          :: v
    type(r_ind)       :: ind
    type(r_dep)       :: dep
    type(ft_int)      :: ft
  END TYPE reac

  TYPE :: ft_mix
    double precision  ::sp, s2p, s3p, op, o2p
    double precision  ::elec, elecHot, s, o, fc, fh
  END TYPE ft_mix

  TYPE :: recomb
    real              ::O_table(REC_ROWS,REC_O)
    real              ::S_table(REC_ROWS,REC_S)
  END TYPE recomb

!!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

END MODULE varTypes
