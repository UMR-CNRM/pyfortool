!#PYFT transfo: --checkContinuationInArraySection Err

SUBROUTINE FOO(A, B)
  REAL, DIMENSION(:,:,:) :: A, B
  A(:,:,:)=B(:, &
  &:,:)
END SUBROUTINE FOO
