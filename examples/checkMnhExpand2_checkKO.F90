!#PYFT transfo: --checkEmptyParensInMnhExpand Err

SUBROUTINE EXPAND6
IMPLICIT NONE

INTEGER :: JI, JJ
INTEGER :: IOPT
INTEGER, DIMENSION(5) :: ICASE
REAL, DIMENSION(5) :: ZZ
REAL, DIMENSION(5, 6) :: ZZZ

!$mnh_expand_array(JI=1:5)
ZZ(:)=1.
!$mnh_end_expand_array(JI=1:5)

END SUBROUTINE EXPAND6
