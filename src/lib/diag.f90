      SUBROUTINE get_dsygvd(H, S, C, E, matrix_size)

      IMPLICIT NONE

      INTEGER :: i, j, matrix_size
      REAL(KIND=8), DIMENSION(matrix_size, matrix_size) :: H, S, C
      REAL(KIND=8), DIMENSION(matrix_size) :: E
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: W
      INTEGER,DIMENSION(:), ALLOCATABLE :: IW

!f2py intent(in, out) :: C
!f2py intent(in, out) :: E

      ALLOCATE (W(100000000), IW(100000000) )

      W = 0.0d0
      IW = 0.0d0

      CALL DSYGVD(1, 'V', 'L', matrix_size, H, matrix_size, S,           &
     &  matrix_size, E, W, 100000000, IW, 100000000, i)

      C = H

!      OPEN(UNIT=22,FILE="eigen_vectors")
!      OPEN(UNIT=23,FILE="eigen_values")

!      DO i=1, matrix_size
!         WRITE(22,*) ((C(i,j)),j= 1, matrix_size)
!      END DO
!      WRITE(23,*) ((E(i)),i = 1, matrix_size)

!      CLOSE(UNIT=22)
!      CLOSE(UNIT=23)

      DEALLOCATE (W, IW)

      END SUBROUTINE get_dsygvd
