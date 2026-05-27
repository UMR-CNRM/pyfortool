!#PYFT transfo: --checkONLY Err

SUBROUTINE CHECK_KO

!The following USE instruction has no ONLY clause, the tool should generate an error
USE MOD_VAR

IMPLICIT NONE

REAL :: VAR1, VAR2

END SUBROUTINE CHECK_KO
