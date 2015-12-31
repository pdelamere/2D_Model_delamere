MODULE dielectronic

USE GLOBAL
USE VARTYPES
USE PARALLELVARIABLES

IMPLICIT NONE

CONTAINS
  subroutine read_rec_tables()    
    integer           ::i,j
    real              ::temp, O(9), S(3)
    OPEN(unit=10, FILE="recombination_O.dat", status="old")
    OPEN(unit=11, FILE="recombination_S.dat", status="old")
    do i=1, REC_ROWS
      read(10,*) O(1), O(2), O(3), O(4), O(5), O(6), O(7), O(8), O(9)
      read(11,*) S(2), S(3)
      S(1) = O(1)

      do j=1, REC_O
        rec_tables%O_table(i,j)=O(j)
      end do

      do j=1, REC_S
        rec_tables%S_table(i,j)=S(j)
      end do
    end do
    close(10)
    close(11)
  end subroutine read_rec_tables

  real function dielectronic_rate(iz, in, T)
    integer          ::iz, in
    real             ::c
    double precision ::T

    select case(iz)
      case(16)
        select case(in)
          case(16)
            c=alpha(0.137e-8,14.95,T)
          case(15)
            c=alpha(0.80729e-8,17.56,T)
            c=c+alpha(0.11012e-9,7.07,T) 
          case(14)
            c=alpha(0.18172e-7,16.62,T)
            c=c+alpha(0.59195e-10,2.4,T)
          case(13)
            c=alpha(0.17710e-7,13.46,T)
        end select
      case(8)
        select case(in)
          case(8)
            c=alpha(0.969e-9, 15.6,T)
          case(7)
            c=alpha(0.382e-8,18.23,T)
          case(6)
            c=alpha(.55713e-8,38.17,T)
            c=c+alpha(0.35028e-10,1.88,T)
            c=c+alpha(0.54359e-8,17.81,T)
        end select
    end select
  
  dielectronic_rate= c*T**(-1.5)

  END function dielectronic_rate

  real function alpha(a, b, T)
    real             ::a, b
    double precision :: T
   
    alpha = a*exp(-b/T)

  end function alpha

  real function O_chart_interpolate(t, col)
    integer    ::col, numCol=2, i !(column of the recombinaion_S.dat file)
    real       ::t, slope   !t is the log base 10 of temp in kelvin
    real       ::lesser, greater, lvalues(2), gvalues(2), minTemp=1.0, maxTemp=9.0

    if (t<minTemp .or. t>maxTemp) then
      print *, "temp is out of bounds of the Oxygen recombination chart"
      stop
    endif

    if (col<1 .or. col>numCol) then
      print *, "bad column for O chart"
      stop
    endif

    i = 1
    greater = rec_tables%O_table(i,1) 
    gvalues(1) = rec_tables%O_table(i,2)
    gvalues(2) = rec_tables%O_table(i,3)

    do while (t>greater .and. greater<maxTemp)
      i = i+1
      lesser = greater
      lvalues = gvalues
      greater = rec_tables%O_table(i,1) 
      gvalues(1) = rec_tables%O_table(i,2)
      gvalues(2) = rec_tables%O_table(i,3)
    enddo

    slope = (gvalues(col)-lvalues(col))/(greater-lesser)
    O_chart_interpolate= slope * (t-lesser) + lvalues(col)
  end function O_chart_interpolate

  real function S_chart_interpolate(t, col)
    integer    ::col, numCol=2, i !(column of the recombinaion_S.dat file)
    real       ::t, slope  !the log base 10 of temp in kelvin
    real       ::lesser, greater, lvalues(2), gvalues(2), minTemp=1.0, maxTemp=9.0

    if(t<minTemp .or. t>maxTemp) then
      print *, "temp is out of bounds of the Sulphur recombination chart"
      print *, "t = ", t
      print *, "processor #", mype
      stop
    endif

    if(col<1 .or. col>numCol) then
      print *, "bad column for S chart"
      stop
    endif

    i=1
    greater = rec_tables%S_table(i,1) 
    gvalues(1) = rec_tables%S_table(i,2)
    gvalues(2) = rec_tables%S_table(i,3)

    do while (t>greater .and. greater<maxTemp)
      i=i+1
      lesser = greater
      lvalues = gvalues
      greater = rec_tables%S_table(i,1) 
      gvalues(1) = rec_tables%S_table(i,2)
      gvalues(2) = rec_tables%S_table(i,3)
    enddo
    slope = (gvalues(col)-lvalues(col))/(greater-lesser)
    S_chart_interpolate= slope * (t-lesser) + lvalues(col)
  end function S_chart_interpolate

END MODULE dielectronic
