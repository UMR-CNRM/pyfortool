!#PYFT transfo: --checkContinuationInArraySection Err

SUBROUTINE FOO(A, B)
  REAL, DIMENSION(:,:,:) :: A, B
  A(:,:,:)=B(:,:,:) + &
          &B(:,:,:)
  CALL F(A, &
        &B)
  A(1,1,1)=F(1, &
            &2)
END SUBROUTINE FOO
