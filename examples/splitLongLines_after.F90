!#PYFT transfo: --maxLineLength 60 --splitLongLines

MODULE MODTEST
  REAL, DIMENSION(10, 10) :: XA
CONTAINS
  SUBROUTINE FOO(PA, PB, PC, PD, PE, PF, PG, PH, PI, PJ, &
    & PK, PL)
    REAL, DIMENSION(:, :), INTENT(INOUT) :: PA, PB, PC, &
      & PD, PE, PF
    REAL, INTENT(IN) :: PG, PH, PI, PJ, PK, PL
    ! This is a very long comment line that must be split on
    ! several lines é à
    PA(:, :) = PB(:, :) + PC(:, :) + PD(:, :) + PE(:, :) &
      & + PF(:, :) * 2.
    PA(1, 1) = PB(1, 1) + XA(1, 1) * FUNC(PG, PH, PI, PJ, &
      & PK, PL) + 3.
    CALL SUB(PA(:, :), PB(:, :), &
      & 'a string, with commas, that is long', PG)
    PA(1, 1) = PG ! a trailing comment that is far too long
                  ! to fit on the line
    PA(1, 1) = PG + PH + PI + PJ ! a trailing comment that
                                 ! must be moved
!$acc kernels present(PA, PB, PC, PD, PE, PF, PG, PH, PI, PJ, PK, PL)
    PA(1, 1) = PG
    PA(:, :) = &
      & PB(PAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA:, :)
!$acc end kernels
  END SUBROUTINE FOO
END MODULE MODTEST
