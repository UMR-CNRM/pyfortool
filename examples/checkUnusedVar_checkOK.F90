!#PYFT transfo: --checkUnusedLocalVar Err

SUBROUTINE CHECK2(ARG)
   
IMPLICIT NONE   

REAL, INTENT(IN) :: ARG
REAL :: Z1, Z2

Z1=0. !Z1 is used
Z2=1. !Z2 is used

END SUBROUTINE CHECK2
