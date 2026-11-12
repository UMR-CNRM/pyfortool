!#PYFT transfo: --checkContinuationInArraySection Err

SUBROUTINE FOO(A, B)
  REAL, DIMENSION(:,:,:) :: A, B
  A(1, &
   &2, 3)=1.
END SUBROUTINE FOO
