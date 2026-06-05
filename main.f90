PROGRAM main

  IMPLICIT NONE

  ! -----------------------------------------------------------------------
  ! Parameters
  ! -----------------------------------------------------------------------
  INTEGER,  PARAMETER :: L       = 100           ! square lattice LxL
  REAL*8,   PARAMETER :: h       = 0.005d0      ! RK4 time step
  REAL*8,   PARAMETER :: t_final = 200.d0       ! total integration time
  INTEGER,  PARAMETER :: NSTEPS  = INT(t_final / h)

  ! ---- Initial-condition seed (change this to get a different IC) --------
  INTEGER,  PARAMETER :: BASE_SEED = 42

  ! -----------------------------------------------------------------------
  ! Sweep ranges  <-- edit these to change the parameter scans
  ! -----------------------------------------------------------------------
  REAL*8,  PARAMETER :: PHI_MIN  = 0.5d0
  REAL*8,  PARAMETER :: PHI_MAX  = 1.0d0
  INTEGER, PARAMETER :: N_PHI    = 60           ! number of phi values

  REAL*8,  PARAMETER :: WA_MIN   = 0.0d0
  REAL*8,  PARAMETER :: WA_MAX   = 100.0d0
  INTEGER, PARAMETER :: N_WA     = 120          ! number of wa values

  ! -----------------------------------------------------------------------
  ! Run parameters
  ! -----------------------------------------------------------------------
  REAL*8    :: wt  = 0.0d0
  REAL*8    :: e   = 1.d0
  REAL*8    :: phi, wa

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: eps = 0.1_dp           ! perturbation amplitude

  ! -----------------------------------------------------------------------
  ! Working arrays
  ! -----------------------------------------------------------------------
  REAL*8 :: p(L,L,4)
  REAL*8 :: csi(L,L,4)
  REAL*8 :: k1(L,L,4), k2(L,L,4), k3(L,L,4), k4(L,L,4)

  ! -----------------------------------------------------------------------
  ! Scalars / loop indices
  ! -----------------------------------------------------------------------
  REAL*8   :: entropy_val
  INTEGER  :: k, wc, wc3

  ! Seed array (size determined at run-time to be portable)
  INTEGER              :: seed_size
  INTEGER, ALLOCATABLE :: seed_array(:)

  CHARACTER*100 :: filename_entropy

  ! -----------------------------------------------------------------------
  ! Outer loop over phi
  ! -----------------------------------------------------------------------
  DO wc3 = 0, N_PHI - 1
    IF (N_PHI > 1) THEN
      phi = PHI_MIN + wc3 * (PHI_MAX - PHI_MIN) / REAL(N_PHI - 1, KIND=8)
    ELSE
      phi = PHI_MIN
    END IF

    ! ---- Inner loop over wa ---------------------------------------------
    DO wc = 0, N_WA - 1
      IF (N_WA > 1) THEN
        wa = WA_MIN + wc * (WA_MAX - WA_MIN) / REAL(N_WA - 1, KIND=8)
      ELSE
        wa = WA_MIN
      END IF

      ! ------------------------------------------------------------------
      ! Set reproducible random seed for this (phi, wa) combination
      ! ------------------------------------------------------------------
      CALL random_seed(size = seed_size)
      ALLOCATE(seed_array(seed_size))
      seed_array = BASE_SEED + wc3 * 1000 + wc
      CALL random_seed(put = seed_array)
      DEALLOCATE(seed_array)

      ! ------------------------------------------------------------------
      ! Initial condition: uniform phi/4 + small random perturbation
      ! ------------------------------------------------------------------
      CALL random_number(csi)
      csi = csi - SUM(csi) / REAL(L*L*4, KIND=8)
      csi = SQRT(eps**2 / SUM(csi**2)) * csi
      p   = phi / 4.d0 + csi

      ! ------------------------------------------------------------------
      ! Open entropy output file for this (wa, phi) run
      ! ------------------------------------------------------------------
      WRITE(filename_entropy, '(A,f0.4,A,f0.4,A,f0.4,A)') &
            'ENTROPY_wa_', wa, '_phi_', phi, '_wt_', wt, '.txt'
      OPEN(99, FILE = filename_entropy, STATUS = 'REPLACE')

      ! ------------------------------------------------------------------
      ! RK4 time integration
      ! ------------------------------------------------------------------
      DO k = 1, NSTEPS
        k1 = RHS(p) * h
        k2 = RHS(p + 0.5d0*k1) * h
        k3 = RHS(p + 0.5d0*k2) * h
        k4 = RHS(p + k3) * h
        p  = p + (k1 + 2.d0*k2 + 2.d0*k3 + k4) / 6.d0

        entropy_val = Entropy(p)
        WRITE(99, *) k, entropy_val
      ENDDO

      CLOSE(99)

    ENDDO   ! end wa loop
  ENDDO     ! end phi loop

! ==========================================================================
CONTAINS
! ==========================================================================

  FUNCTION Entropy(p)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: p(L,L,4)
    REAL*8 :: Entropy
    INTEGER :: ii, jj, ss
    Entropy = 0.d0
    DO ii = 1, L
      DO jj = 1, L
        DO ss = 1, 4
          Entropy = Entropy - p(ii,jj,ss) * LOG(p(ii,jj,ss) * 4.d0 / phi)
        ENDDO
      ENDDO
    ENDDO
    Entropy = Entropy / (phi * REAL(L*L, KIND=8))
  END FUNCTION Entropy

  FUNCTION RHS(p)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: p(L,L,4)
    REAL*8 :: RHS(L,L,4)
    REAL*8 :: aux, aaux
    INTEGER :: ii, jj, im, ip, jm, jp

    RHS = 0.d0
    DO ii = 1, L
      im = ii - 1;  IF (ii .EQ. 1) im = L
      ip = ii + 1;  IF (ii .EQ. L) ip = 1
      DO jj = 1, L
        jm = jj - 1;  IF (jj .EQ. 1) jm = L
        jp = jj + 1;  IF (jj .EQ. L) jp = 1

        aux  = SUM(p(ii,jj,:))
        aaux = SUM(p(im,jj,:)) + SUM(p(ip,jj,:)) + SUM(p(ii,jm,:)) + SUM(p(ii,jp,:))

        RHS(ii,jj,1) = wa*p(ii,jm,1)*(1.d0 - e*aux) - wa*p(ii,jj,1)*(1.d0 - e*SUM(p(ii,jp,:))) &
                     + p(ii,jj,2) + p(ii,jj,4) - 2.d0*p(ii,jj,1) &
                     + wt*(p(im,jj,1)+p(ip,jj,1)+p(ii,jm,1)+p(ii,jp,1))*(1.d0 - e*aux) &
                     - wt*p(ii,jj,1)*(4.d0 - e*aaux)

        RHS(ii,jj,2) = wa*p(im,jj,2)*(1.d0 - e*aux) - wa*p(ii,jj,2)*(1.d0 - e*SUM(p(ip,jj,:))) &
                     + p(ii,jj,1) + p(ii,jj,3) - 2.d0*p(ii,jj,2) &
                     + wt*(p(im,jj,2)+p(ip,jj,2)+p(ii,jm,2)+p(ii,jp,2))*(1.d0 - e*aux) &
                     - wt*p(ii,jj,2)*(4.d0 - e*aaux)

        RHS(ii,jj,3) = wa*p(ii,jp,3)*(1.d0 - e*aux) - wa*p(ii,jj,3)*(1.d0 - e*SUM(p(ii,jm,:))) &
                     + p(ii,jj,2) + p(ii,jj,4) - 2.d0*p(ii,jj,3) &
                     + wt*(p(im,jj,3)+p(ip,jj,3)+p(ii,jm,3)+p(ii,jp,3))*(1.d0 - e*aux) &
                     - wt*p(ii,jj,3)*(4.d0 - e*aaux)

        RHS(ii,jj,4) = wa*p(ip,jj,4)*(1.d0 - e*aux) - wa*p(ii,jj,4)*(1.d0 - e*SUM(p(im,jj,:))) &
                     + p(ii,jj,1) + p(ii,jj,3) - 2.d0*p(ii,jj,4) &
                     + wt*(p(im,jj,4)+p(ip,jj,4)+p(ii,jm,4)+p(ii,jp,4))*(1.d0 - e*aux) &
                     - wt*p(ii,jj,4)*(4.d0 - e*aaux)

      ENDDO
    ENDDO
  END FUNCTION RHS

END PROGRAM main
