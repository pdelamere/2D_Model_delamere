!equivalent to cm3_model_common in Steffl
!holds dybnamic global variables
!constant values are in varTypes.f90 with structure declarations

MODULE GLOBAL

USE varTypes

real                 :: rates,trans,net_source,trans_exp,iondens0,net_source0,tau0,longitude
integer              :: info,runt,dt
real                 :: zoff,rdist, v_r0
real                 :: sys3variations,lag_amp,lag_phase,lag_const,fehot_amp,fehot_phase,fehot_const, lon3
real                 :: Fs, Fsp, Fs2p, Fs3p, Fo, Fop, Fo2p !used in fluxCorrect in timestep module
real                 :: EFsp, EFs2p, EFs3p, EFop, EFo2p, EFelec !Used in updateNT in timestep module

type(recomb)         ::rec_tables !Holds recombinaion tables read in dielectronics.f90

logical              ::trans_type=.false.

real                 ::protons=0.0  !<can be added to inputs.dat

END MODULE GLOBAL
