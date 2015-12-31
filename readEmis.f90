MODULE ReadEmis

  USE varTypes

  IMPLICIT NONE

  CONTAINS

  subroutine ReadIndices(temp, dens)
    real              ::temp(EMIS_SIZE), dens(EMIS_SIZE)
 
    call ReadIndex('emisTemp.dat', temp)  !make sure to change character lengh if filenames change    
    call ReadIndex('emisDens.dat', dens)      

  end subroutine ReadIndices
  
  subroutine ReadIndex(loc, array)
    real              ::array(EMIS_SIZE)
    character(len=12) ::loc !location or filename of data
    integer           ::i

    
    open(unit=10, file=loc, status="old")
    do i=1, EMIS_SIZE
      read(10,*) array(i)
    end do
    close(10)

  end subroutine ReadIndex

  subroutine ReadEmisTable(loc, emis)
    real              ::emis(EMIS_SIZE,EMIS_SIZE)
    character(len=20) ::loc  !location of data (filename)
    integer           ::i,j

    open(unit=10, file=loc, status="old")

    do i=1, EMIS_SIZE
      do j=1, EMIS_SIZE
        read(10,*), emis(j, i)
      end do
    end do
   
    do i=1, EMIS_SIZE
      do j=1, EMIS_SIZE
        emis(j, i)=emis(j,i)/(1.60217646e-12)
      end do
    end do
!    print *, emis(1,1), emis(1,101) , emis(101,1), emis(101, 101) 

  end subroutine ReadEmisTable

subroutine InitIndependentRates(ind)
  type(r_ind)         ::ind
  character(len=20)   ::loc

  call ReadIndices(ind%emis_temp, ind%emis_dens) !reads in the temp and density tables for inerpolate the emmission tables
! The above function should be unneccesary since values are now calculated rather than searched in a table. Calculation is faster thant searching.

!Read in all emission tables. These tables are used to determine power radiated by the torus
  loc='emisSp.dat'
  call ReadEmisTable(loc, ind%emis_sp)
  loc='emisS2p.dat'
  call ReadEmisTable(loc, ind%emis_s2p)
  loc='emisS3p.dat'
  call ReadEmisTable(loc, ind%emis_s3p)
  loc='emisOp.dat'
  call ReadEmisTable(loc, ind%emis_op)
  loc='emisO2p.dat'
  call ReadEmisTable(loc, ind%emis_o2p)

end subroutine InitIndependentRates

END MODULE ReadEmis

