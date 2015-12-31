MODULE FTMIX

USE varTypes

IMPLICIT NONE

CONTAINS

  function ftint_mix(n, h)
    type(density)     ::n
    type(height)      ::h
    type(ft_mix)      ::ftint_mix
    real              ::tot, sp, s2p, s3p, op, o2p

    ftint_mix%elec = n%elec
    ftint_mix%elecHot = n%elecHot
    ftint_mix%s = n%s
    ftint_mix%o = n%o
    ftint_mix%fc = n%fc
    ftint_mix%fh = n%fh

    sp  = n%sp  * ROOTPI * h%sp
    s2p = n%s2p * ROOTPI * h%s2p
    s3p = n%s3p * ROOTPI * h%s3p
    op  = n%op  * ROOTPI * h%op
    o2p = n%o2p * ROOTPI * h%o2p

    tot = (sp + op + 2.0 * (s2p + o2p) + 3.0 * s3p)/(1.0-n%protons)

    ftint_mix%sp  = sp/tot
    ftint_mix%s2p = s2p/tot
    ftint_mix%s3p = s3p/tot
    ftint_mix%op  = op/tot
    ftint_mix%o2p = o2p/tot

  end function ftint_mix

END MODULE

