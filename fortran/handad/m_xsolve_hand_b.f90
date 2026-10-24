!***********************************************************************
!   XFOILAD - hand-written adjoints of the linear solvers
!
!   These routines replace the Tapenade versions of the same name in
!   generated/m_xsolve_b.f90 (tools/insert_hand_adjoints.py, run by
!   build.sh --tapenade).  Tapenade differentiates the eliminations
!   statement by statement and stores every overwritten matrix entry on
!   its stack (about 2e6 values per BLSOLV call); the adjoint of a linear
!   solve  A x = b  only needs a solve with the transposed matrix:
!
!       b_bar = A^-T x_bar,      A_bar = A_bar - b_bar x^T
!
!   BLSOLV (split mode, see tapenade/run_tapenade.sh): BLSOLV_FWD solves
!   and keeps the factors and the solution (module xad_hand_store),
!   BLSOLV_BWD replays the row operations of the
!   elimination (with the same dropped entries, VACCEL) in reverse order,
!   with the multipliers left in VA / VM by the factorization.  The matrix
!   adjoint uses the formula of the exact solve; it differs from the exact
!   derivative of the approximate solver by terms proportional to the
!   Newton step, which vanish at the converged state the adjoint is
!   evaluated at.
!
!   LUDCMP / BAKSUB: AIJ is only used through BAKSUB after LUDCMP, so
!   BAKSUB_BWD (split mode) gives directly the adjoint of the matrix before
!   factorization, accumulated over all the solves by LUDCMP_B.
!***********************************************************************

  SUBROUTINE LUDCMP_B(nsiz, n, a, ab, indx)
!   adjoint of the unfactored matrix: A_bar = A_bar - Y X^T, with the
!   columns stored by BAKSUB_BWD (one matrix product instead of one rank-1
!   update per solve)
    USE XAD_HAND_STORE, ONLY : bak_y, bak_x, bak_m
    IMPLICIT NONE
    INTEGER :: n, nsiz
    REAL, DIMENSION(nsiz, nsiz) :: a
    REAL, DIMENSION(nsiz, nsiz) :: ab
    INTEGER, DIMENSION(nsiz) :: indx
    IF (bak_m .GT. 0) THEN
      ab(1:n, 1:n) = ab(1:n, 1:n) - MATMUL(bak_y(1:n, 1:bak_m), TRANSPOSE(bak_x(1:n, 1:bak_m)))
      bak_m = 0
    END IF
  END SUBROUTINE LUDCMP_B


  SUBROUTINE BAKSUB_FWD(nsiz, n, a, indx, b)
!   forward sweep: the solve itself (the caller restores b to the solution
!   before BAKSUB_BWD)
    IMPLICIT NONE
    INTEGER :: n, nsiz
    REAL, DIMENSION(nsiz, nsiz) :: a
    REAL, DIMENSION(nsiz) :: b
    INTEGER, DIMENSION(nsiz) :: indx
    CALL BAKSUB(nsiz, n, a, indx, b)
  END SUBROUTINE BAKSUB_FWD


  SUBROUTINE BAKSUB_BWD(nsiz, n, a, ab, indx, b, bb)
!   b holds the solution x of the forward solve, bb its adjoint on entry
!   and the adjoint of the right-hand side on exit
    USE XAD_HAND_STORE, ONLY : bak_y, bak_x, bak_push
    IMPLICIT NONE
    INTEGER :: n, nsiz
    REAL, DIMENSION(nsiz, nsiz) :: a
    REAL, DIMENSION(nsiz, nsiz) :: ab
    REAL, DIMENSION(nsiz) :: b
    REAL, DIMENSION(nsiz) :: bb
    INTEGER, DIMENSION(nsiz) :: indx
    INTEGER :: i, j, ll, k
    REAL :: sum, tmp
!
!---- U^T z = x_bar (forward), then L^T w = z (backward, unit diagonal)
    DO i=1,n
      sum = bb(i)
      DO j=1,i-1
        sum = sum - a(j, i)*bb(j)
      END DO
      bb(i) = sum/a(i, i)
    END DO
    DO i=n,1,-1
      sum = bb(i)
      DO j=i+1,n
        sum = sum - a(j, i)*bb(j)
      END DO
      bb(i) = sum
    END DO
!---- transposed row interchanges, in reverse order
    DO i=n,1,-1
      ll = indx(i)
      tmp = bb(ll)
      bb(ll) = bb(i)
      bb(i) = tmp
    END DO
!
!---- adjoint of the unfactored matrix, applied by LUDCMP_B
    k = bak_push(n)
    bak_y(1:n, k) = bb(1:n)
    bak_x(1:n, k) = b(1:n)
  END SUBROUTINE BAKSUB_BWD


  SUBROUTINE BLSOLV_FWD()
!   forward sweep: primal solve, factors and solution kept for BLSOLV_BWD
    USE I_XFOIL_DIFF
    USE XAD_HAND_STORE, ONLY : bls_stack, bls_push
    IMPLICIT NONE
    INTEGER :: iv, l, nsy, k
!
    nsy = nsys
    k = bls_push(nsy)
    ASSOCIATE (r => bls_stack(k))
      r%ivte1 = isys(iblte(1), 1)
      r%ivz = isys(iblte(2)+1, 2)
      r%vacc(1) = vaccel
      r%vacc(2) = vaccel*2.0/(s(n)-s(1))
      r%vacc(3) = r%vacc(2)
      r%vz = vz
      DO iv=1,nsy
        r%va0(:, :, iv) = va(:, :, iv)
        r%vdel0(:, :, iv) = vdel(:, :, iv)
        r%vb(:, :, iv) = vb(:, :, iv)
        DO l=1,nsy
          r%vm0(:, l, iv) = vm(:, l, iv)
        END DO
      END DO
      CALL BLSOLV()
      DO iv=1,nsy
        r%va(:, :, iv) = va(:, :, iv)
        r%d(:, :, iv) = vdel(:, :, iv)
        DO l=1,nsy
          r%vm(:, l, iv) = vm(:, l, iv)
        END DO
      END DO
    END ASSOCIATE
  END SUBROUTINE BLSOLV_FWD


  SUBROUTINE BLSOLV_BWD()
!   backward sweep: transposed solve with the factors of BLSOLV_FWD
    USE I_XFOIL_DIFF
    USE XAD_HAND_STORE, ONLY : bls_stack, bls_top
    IMPLICIT NONE
    INTEGER :: iv, ivp, ivte1, ivz, k, kv, l, c, nsy
    REAL :: s3, r1, r2, r3, d1, d2, d3
    REAL, DIMENSION(3, 2, izx), SAVE :: rb
!
    ASSOCIATE (r => bls_stack(bls_top))
      nsy = r%nsys
      ivte1 = r%ivte1
      ivz = r%ivz
      DO iv=1,nsy
        rb(:, :, iv) = vdelb(:, :, iv)
      END DO
!
!---- transpose of the back substitution
      DO iv=2,nsy
        DO c=1,2
          s3 = 0.0
          DO kv=1,iv-1
            s3 = s3 + r%vm(1, iv, kv)*rb(1, c, kv) + r%vm(2, iv, kv)*rb(2, c, kv) &
&             + r%vm(3, iv, kv)*rb(3, c, kv)
          END DO
          rb(3, c, iv) = rb(3, c, iv) - s3
        END DO
      END DO
!
!---- transpose of the forward elimination, stations in reverse order
      DO iv=nsy,1,-1
        ivp = iv + 1
        DO c=1,2
          r1 = rb(1, c, iv)
          r2 = rb(2, c, iv)
          r3 = rb(3, c, iv)
          IF (iv .NE. nsy) THEN
            IF (ivp .NE. nsy) THEN
!             lower VM column (same dropped entries as the forward elimination)
              s3 = 0.0
              DO kv=iv+2,nsy
                DO k=1,3
                  IF (ABS(r%vm(k, iv, kv)) .GT. r%vacc(k)) s3 = s3 + r%vm(k, iv, kv)*rb(k, c, kv)
                END DO
              END DO
              r3 = r3 - s3
            END IF
            IF (iv .EQ. ivte1) THEN
!             VZ block
              DO k=1,3
                r1 = r1 - r%vz(k, 1)*rb(k, c, ivz)
                r2 = r2 - r%vz(k, 2)*rb(k, c, ivz)
              END DO
            END IF
!           VB(IV+1) block
            DO k=1,3
              r1 = r1 - r%vb(k, 1, ivp)*rb(k, c, ivp)
              r2 = r2 - r%vb(k, 2, ivp)*rb(k, c, ivp)
              r3 = r3 - r%vm(k, iv, ivp)*rb(k, c, ivp)
            END DO
          END IF
!         inversion of the VA(IV) block
          r2 = r2 - r%va(1, 2, iv)*r1
          r3 = r3 - r%vm(1, iv, iv)*r1 - r%vm(2, iv, iv)*r2
          r3 = r3/r%vm(3, iv, iv)
          r2 = r2 - r%va(3, 2, iv)*r3
          r2 = r2/r%va(2, 2, iv)
          r1 = r1 - r%va(2, 1, iv)*r2 - r%va(3, 1, iv)*r3
          r1 = r1/r%va(1, 1, iv)
          rb(1, c, iv) = r1
          rb(2, c, iv) = r2
          rb(3, c, iv) = r3
        END DO
      END DO
!
!---- adjoints of the matrix blocks:  A_bar = A_bar - sum_c rb_c d_c^T
      DO iv=1,nsy
        DO c=1,2
          d1 = r%d(1, c, iv)
          d2 = r%d(2, c, iv)
          DO k=1,3
            vab(k, 1, iv) = vab(k, 1, iv) - rb(k, c, iv)*d1
            vab(k, 2, iv) = vab(k, 2, iv) - rb(k, c, iv)*d2
          END DO
          IF (iv .GE. 2) THEN
            d1 = r%d(1, c, iv-1)
            d2 = r%d(2, c, iv-1)
            DO k=1,3
              vbb(k, 1, iv) = vbb(k, 1, iv) - rb(k, c, iv)*d1
              vbb(k, 2, iv) = vbb(k, 2, iv) - rb(k, c, iv)*d2
            END DO
          END IF
        END DO
        DO l=1,nsy
          d3 = r%d(3, 1, l)
          d1 = r%d(3, 2, l)
          DO k=1,3
            vmb(k, l, iv) = vmb(k, l, iv) - rb(k, 1, iv)*d3 - rb(k, 2, iv)*d1
          END DO
        END DO
      END DO
      DO c=1,2
        DO k=1,3
          vzb(k, 1) = vzb(k, 1) - rb(k, c, ivz)*r%d(1, c, ivte1)
          vzb(k, 2) = vzb(k, 2) - rb(k, c, ivz)*r%d(2, c, ivte1)
        END DO
      END DO
!
!---- adjoint of the right-hand sides; inputs restored
      DO iv=1,nsy
        vdelb(:, :, iv) = rb(:, :, iv)
        va(:, :, iv) = r%va0(:, :, iv)
        vdel(:, :, iv) = r%vdel0(:, :, iv)
        DO l=1,nsy
          vm(:, l, iv) = r%vm0(:, l, iv)
        END DO
      END DO
    END ASSOCIATE
    bls_top = bls_top - 1
  END SUBROUTINE BLSOLV_BWD
