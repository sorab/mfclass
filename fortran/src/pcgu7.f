      SUBROUTINE PCGU7U1AR(IN, MXITER, HICLOSE, ITER1, IFDPARAM, IPCGUM)
C     ******************************************************************
C     ALLOCATE STORAGE FOR PCG ARRAYS AND READ PCGU DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE SimModule
      USE SOLUTIONMODULE, ONLY: neq,nja,ia,ja,amat,rhs,x
      USE PCGUMODULE
      IMPLICIT NONE
C     + + + DUMMY VARIABLES + + +
      INTEGER, INTENT(IN) :: IN
      INTEGER, INTENT(IN) :: MXITER
      DOUBLE PRECISION, INTENT(IN) :: HICLOSE
      INTEGER, INTENT(IN) :: ITER1
      INTEGER, INTENT(IN) :: IFDPARAM
      INTEGER, INTENT(INOUT) :: IPCGUM
C     + + + LOCAL VARIABLES + + +
      CHARACTER*200 LINE
      CHARACTER (LEN= 10) :: clin(1:2)
      CHARACTER (LEN= 20) :: cipc(0:3)
      CHARACTER (LEN= 20) :: cscale(0:2)
      CHARACTER (LEN= 25) :: corder(0:2)
      CHARACTER (LEN=16), DIMENSION(2) :: ccnvgopt
      CHARACTER (LEN=  4) :: cval
      INTEGER :: lloc
      INTEGER :: istart, istop
      INTEGER :: i, n
      INTEGER :: i0
      INTEGER :: iscllen, iolen
      REAL :: r
C     + + + PARAMETERS + + + 
      INTEGER, PARAMETER :: IZERO = 0
      REAL, PARAMETER :: RZERO = 0.0
      DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER :: DONE  = 1.0D0
C       DATA
      DATA clin  /'CG        ',
     2            'BCGS      '/
      DATA cipc  /'NONE                ',
     2            'JACOBI              ',
     3            'INCOMPLETE LU       ',
     4            'MOD. INCOMPLETE LU  '/
      DATA cscale/'NO SCALING          ',
     2            'SYMMETRIC SCALING   ',
     3            'L2 NORM SCALING     '/
      DATA corder/'ORIGINAL ORDERING        ',
     2            'RCM ORDERING             ',
     3            'MINIMUM DEGREE ORDERING  '/
      DATA ccnvgopt  /'INFINITY NORM   ',
     2                'L2 NORM         '/
C       OUTPUT FORMATS
02010 FORMAT (1X,/,14X,'SOLUTION BY THE CONJUGATE-GRADIENT METHOD',
     &        /,1X,66('-'),/,
     &        ' MAXIMUM OF ',I6,' CALLS OF SOLUTION ROUTINE',/,
     &        ' MAXIMUM OF ',I6,
     &        ' INTERNAL ITERATIONS PER CALL TO SOLUTION ROUTINE',/,
     &        ' LINEAR ACCELERATION METHOD            =',1X,A,/,
     &        ' MATRIX PRECONDITIONING TYPE           =',1X,A,/,
     &        ' MATRIX SCALING APPROACH               =',1X,A,/,
     &        ' MATRIX REORDERING APPROACH            =',1X,A,/,
     &        ' HEAD CHANGE CRITERION FOR CLOSURE     =',E15.5,/,
     &        ' RESIDUAL CHANGE CRITERION FOR CLOSURE =',E15.5,/,
     &        ' RESIDUAL CONVERGENCE OPTION           =',I9,/,
     &        ' RESIDUAL CONVERGENCE NORM             =',1X,A,/,
     &        ' RELAXATION FACTOR                     =',E15.5,/,
     &        '  ONLY USED WITH MILU0 PRECONDITIONER',//)
02020 FORMAT (///,1X,'PCGU DATA INPUT ERROR:',
     &          /,2X,'SCALING MUST BE USED (ISCL.GT.0) IF USING',
     &          /,2X,'THE ILU0 OR MILU0 PRECONDITIONERS (IPC.EQ.2 OR',
     &          /,2X,'IPC.EQ.3) WITH MATRIX REORDERING (IORD.GT.0)')
2030  FORMAT(1X,A20,1X,6(I6,1X))
2040  FORMAT(1X,20('-'),1X,6(6('-'),1X))
2050  FORMAT(1X,62('-'),/)
C     ------------------------------------------------------------------
C
C-------ALLOCATE VARIABLES
        ALLOCATE(ILINMETH)
        ALLOCATE(ITER1C)
        ALLOCATE(IPC)
        ALLOCATE(ISCL)
        ALLOCATE(IORD)
        ALLOCATE(ICNVGOPT)
        ALLOCATE(IACPC)
        ALLOCATE(NITERC)
        ALLOCATE(NIABCGS)
        ALLOCATE(NIAPC)
        ALLOCATE(NJAPC)
        ALLOCATE(NNZAPC)
        ALLOCATE(HCLOSEPCGU)
        ALLOCATE(RCLOSEPCGU)
        ALLOCATE(RELAXPCGU)
        ALLOCATE(EPFACT)
        ALLOCATE(L2NORM0)
      
C
C-------TRANSFER COMMON VARIABLES FROM SMS TO UPCG
      ILINMETH = 0
      
      HCLOSEPCGU = HICLOSE
      ITER1C = ITER1
      IACPC = 1
      
      ICNVGOPT = 0
C
C-------PRINT A MESSAGE IDENTIFYING UPCG PACKAGE
      WRITE (IOUT,2000)
02000 FORMAT (1X,/1X,'PCGU -- UNSTRUCTURED CONJUGATE-GRADIENT SOLUTION',
     &        ' PACKAGE, VERSION 7.02, 08/13/2013')
C
C-------READ AND PRINT COMMENTS
      CALL URDCOM(IN,IOUT,LINE)
      IF ( IFDPARAM.EQ.0 ) THEN
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IPC,R,IOUT,IN)
        cval = LINE(ISTART:ISTOP)
        SELECT CASE (cval)
          CASE ( 'CG' )
            ILINMETH = 1
          CASE ( 'BCGS' )
            ILINMETH = 2
          CASE DEFAULT
            ILINMETH = 1
            READ (CVAL,*) IPC
        END SELECT
        IF ( cval.EQ.'CG  ' .OR. cval.EQ.'BCGS' ) THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPC,R,IOUT,IN)
        END IF
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISCL,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IORD,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,RCLOSEPCGU,IOUT,IN)
        IF ( RCLOSEPCGU.LT.RZERO ) THEN
          RCLOSEPCGU = ABS( RCLOSEPCGU )
          ICNVGOPT = 1
        END IF
        IF ( IPC.EQ.3 ) THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,RELAXPCGU,-IOUT,IN)
          IF ( ISTART.EQ.200 ) THEN
            IF ( ISTOP.EQ.200 ) RELAXPCGU = 0.97
          END IF
          IF ( RELAXPCGU.EQ.RZERO ) THEN
            IF ( LINE(200:200).EQ.'E' ) RELAXPCGU = 0.97
          END IF
        ELSE
          RELAXPCGU = 0.0
        END IF
      ELSE
        CALL SET_PCGUINPUT(IFDPARAM)
      END IF
      IPCGUM = ILINMETH
C
C-------ERROR CHECKING FOR OPTIONS
      IF ( IPC.LT.0  ) IPC  = 0
      IF ( IPC.GT.3  ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: IPC  MUST BE .LE. 3'
        CALL USTOP('PCGU7AR: IPC  MUST BE .LE. 3')
      END IF
      IF ( ISCL.LT.0 ) ISCL = 0
      IF ( ISCL.GT.2  ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: ISCL MUST BE .LE. 2'
        CALL USTOP('PCGU7AR: ISCL MUST BE .LE. 2')
      END IF
      IF ( IORD.LT.0 ) IORD = 0
      IF ( IORD.GT.2  ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: IORD MUST BE .LE. 2'
        CALL USTOP('PCGU7AR: IORD MUST BE .LE. 2')
      END IF
      IF ( RCLOSEPCGU.EQ.0.0 ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: RCLOSEPCGU MUST .NE. 0.0'
        CALL USTOP('PCGU7AR: RCLOSEPCGU MUST .NE. 0.0')
      END IF
      IF ( RELAXPCGU.LT.0.0 ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: RELAXPCGU MUST BE .GE. 0.0'
        CALL USTOP('PCGU7AR: RELAXPCGU MUST BE .GE. 0.0')
      END IF
      IF ( RELAXPCGU.GT.1.0 ) THEN
        WRITE( IOUT,'(A)' ) 'PCGU7AR: RELAXPCGU MUST BE .LE. 1.0'
        CALL USTOP('PCGU7AR: RELAXPCGU MUST BE .LE. 1.0')
      END IF
C
C-------PRINT MXITER,ITER1C,IPC,ISCL,IORD,HCLOSEPCGU,RCLOSEPCGU
      WRITE (IOUT,2010) MXITER, ITER1C, 
     2                   clin(ILINMETH), cipc(IPC), 
     3                   cscale(ISCL), corder(IORD), 
     4                   HCLOSEPCGU, RCLOSEPCGU, 
     5                   ICNVGOPT,ccnvgopt(ICNVGOPT+1),
     6                   RELAXPCGU
C
C-------ENSURE THAT SCALING IS USED WITH THE ILU0 AND MILU0
C       PRECONDITIONERS IF RCM OR MINIMUM DEGREE ORDERING IS USED
      IF ( IPC.EQ.2 .OR. IPC.EQ.3 ) THEN
        IF ( IORD.NE.0 ) THEN
          IF ( ISCL.EQ.0 ) THEN
            WRITE ( IOUT,2020 )
            CALL USTOP('SCALING MUST BE USED FOR ILU0 AND MILU0 '//
     2                 'WITH REORDERING')
          END IF
        END IF
      END IF
C
C-------INITIALIZE PCGU VARIABLES
      NITERC = 0
C
C-------ALLOCATE AND INITIALIZE MEMORY FOR PCGU
      iscllen  = 1
      IF ( ISCL.NE.0 ) iscllen  = NEQ
      ALLOCATE ( DSCALE(iscllen), DSCALE2(iscllen) )
C       ALLOCATE MEMORY FOR PRECONDITIONING MATRIX
      NIAPC  = NEQ
      NJAPC  = NJA
      NNZAPC = NJA
      IF ( IPC.EQ.0 ) THEN
        NIAPC  = 1
        NJAPC  = 1
        NNZAPC = 1
      ELSE IF ( IPC.EQ.1 ) THEN
        NIAPC  = 1
        NJAPC  = 1
        NNZAPC = NJA
      END IF
      ALLOCATE( IAPC(NIAPC+1) )
      ALLOCATE( JAPC(NJAPC) )
      ALLOCATE( APC(NNZAPC) )
C       ALLOCATE MEMORY FOR ILU0 AND MILU0 NON-ZERO ROW ENTRY VECTOR
      ALLOCATE( IW(NIAPC) )
      ALLOCATE( W(NIAPC) )
C       GENERATE IAPC AND JAPC
      IF ( IPC.EQ.2 .OR. IPC.EQ.3 ) THEN
        CALL SPCGU_PCCRS(NEQ,NJA,IA,JA,
     2                   IAPC,JAPC)
      END IF
C       ALLOCATE SPACE FOR PERMUTATION VECTOR
      i0     = 1
      iolen  = 1
      IF ( IORD.NE.0 ) THEN
        i0     = NEQ
        iolen  = NJA
      END IF
      ALLOCATE( LORDER(i0)  )
      ALLOCATE( IORDER(i0)  )
      ALLOCATE( IARO(i0+1)  )
      ALLOCATE( JARO(iolen) )
      ALLOCATE( ARO(iolen)  )
C       ALLOCATE WORKING VECTORS FOR PCGU SOLVER      
      ALLOCATE( ID(NEQ) )
      ALLOCATE( D(NEQ) )
      ALLOCATE( P(NEQ) )
      ALLOCATE( Q(NEQ) )
      ALLOCATE( Z(NEQ))
C       ALLOCATE MEMORY FOR BCGS WORKING ARRAYS
      NIABCGS = 1
      IF ( ILINMETH.EQ.2 ) THEN
        NIABCGS = NEQ
      END IF
      ALLOCATE( T(NIABCGS) )
      ALLOCATE( V(NIABCGS) )
      ALLOCATE( DHAT(NIABCGS) )
      ALLOCATE( PHAT(NIABCGS) )
      ALLOCATE( QHAT(NIABCGS) )
C       INITALIZE PCGU VECTORS
      DO n = 1, iscllen
        DSCALE(n)  = DONE
        DSCALE2(n) = DONE
      END DO
      DO n = 1, NNZAPC
        APC(n)  = DZERO
      END DO
C       WORKING VECTORS
      DO n = 1, NEQ
        ID(n)    = IZERO
        X(n)    = DZERO
        D(n)    = DZERO
        P(n)    = DZERO
        Q(n)    = DZERO
        Z(n)    = DZERO
      END DO
      DO n = 1, NIAPC
        IW(n)   = IZERO
        W(n)    = DZERO
      END DO
C       BCGS WORKING VECTORS
      DO n = 1, NIABCGS
        T(n)    = DZERO
        V(n)    = DZERO
        DHAT(n) = DZERO
        PHAT(n) = DZERO
        QHAT(n) = DZERO
      END DO
C-------REORDERING VECTORS
      DO n = 1, i0 + 1
        IARO(n) = IZERO
      END DO
      DO n = 1, iolen
        JARO(n) = IZERO
        ARO(n)  = DZERO
      END DO
C
C-------REVERSE CUTHILL MCKEE ORDERING
C       NOTE - USING GNRCM AND ODRV SUBROUTINES IN THE XMD SOLVER SOURCE CODE
C              SPECIFICALLY IN xmblib.f
      IF ( IORD.NE.0 ) THEN
        CALL SSPCGU_CALC_ORDER( IOUT,IORD,NEQ,NJA,IA,JA,
     2                          LORDER,IORDER )
      END IF
C
C-------RETURN
      RETURN
      END SUBROUTINE PCGU7U1AR
C
      SUBROUTINE PCGU7U1AP(ICNVG,KSTP,KITER,IN_ITER)
C
C     ******************************************************************
C     SOLUTION BY THE CONJUGATE GRADIENT METHOD -
C                                          UP TO ITER1 ITERATIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE SimModule
      USE SOLUTIONMODULE, ONLY: neq,nja,ia,ja,amat,rhs,x
      USE PCGUMODULE
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      INTEGER, INTENT(INOUT)                          :: ICNVG
      INTEGER, INTENT(IN)                             :: KSTP
      INTEGER, INTENT(IN)                             :: KITER
      INTEGER, INTENT(INOUT)                          :: IN_ITER
C     + + + LOCAL DEFINITIONS + + +
      INTEGER :: n
      INTEGER :: innerit
      INTEGER :: irc
      INTEGER :: itmax
      DOUBLEPRECISION :: tv
      DOUBLEPRECISION :: rmax
C     + + + PARAMETERS + + + 
      DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
      DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C
C     + + + CODE + + +
C
C-------SET EPFACT BASED ON MFUSG TIMESTEP
      IF ( KSTP.EQ.1 ) THEN
        EPFACT = 0.01
      ELSE
        EPFACT = 0.10
      END IF

C-------SCALE PROBLEM
      IF ( ISCL.NE.0 ) THEN
        CALL SPCGU_SCALE(0,ISCL,NEQ,NJA,IA,JA,AMAT,X,RHS,
     2                   DSCALE,DSCALE2)
      END IF
C
C-------PERMUTE ROWS, COLUMNS, AND RHS
      IF ( IORD.NE.0 ) THEN
        CALL DPERM(NEQ,AMAT,JA,IA,ARO,JARO,IARO,LORDER,ID,1)
        CALL VPERM(NEQ,  X, LORDER)  
        CALL VPERM(NEQ, RHS, LORDER)
        IA0 => IARO
        JA0 => JARO
        A0  => ARO
      ELSE
        IA0 => IA
        JA0 => JA
        A0  => AMAT
      END IF
C
C-------UPDATE PRECONDITIONER
      CALL SPCGU_PCU(IOUT,NJA,NEQ,NIAPC,NJAPC,NNZAPC,IPC,RELAXPCGU,
     2               A0,IA0,JA0,APC,IAPC,JAPC,IW,W)
C-------INITILIZE SOLUTION VARIABLE AND ARRAYS
      IF ( KITER.EQ.1 ) NITERC = 0
      irc    = 1
      ICNVG  = 0
      DO n = 1, NEQ
        D(n) = DZERO
        P(n) = DZERO
        Q(n) = DZERO
        Z(n) = DZERO
      END DO
C-------CALCULATE INITIAL RESIDUAL
      CALL SPCGU_MV(NJA,NEQ,A0,X,D,IA0,JA0)
      rmax = DZERO
      L2NORM0 = DZERO
      DO n = 1, NEQ
        tv     = D(n)
        D(n) = RHS(n) - tv
        IF ( ABS( D(n) ).GT.rmax ) rmax = ABS( D(n) )
        L2NORM0 = L2NORM0 + D(n) * D(n)
      END DO
      L2NORM0 = SQRT(L2NORM0)
C-------CHECK FOR EXACT SOLUTION
      itmax = ITER1C
      IF ( rmax.EQ.DZERO ) itmax = 0
C-------SOLUTION BY THE CONJUGATE GRADIENT METHOD      
      IF ( ILINMETH.EQ.1 ) THEN
        CALL SPCGU_CG(ICNVG,itmax,innerit)
C-------SOLUTION BY THE BICONJUGATE GRADIENT STABILIZED METHOD      
      ELSE IF ( ILINMETH.EQ.2 ) THEN
        CALL SPCGU_BCGS(ICNVG,itmax,innerit)
      END IF
C
C-------BACK PERMUTE AMAT, SOLUTION, AND RHS
      IF ( IORD.NE.0 ) THEN
        CALL DPERM(NEQ,A0,JA0,IA0,AMAT,JA,IA,IORDER,ID,1)
        CALL VPERM(NEQ,  X, IORDER)  
        CALL VPERM(NEQ, RHS, IORDER)  
      END IF
C
C-------UNSCALE PROBLEM
      IF ( ISCL.NE.0 ) THEN
        CALL SPCGU_SCALE(1,ISCL,NEQ,NJA,IA,JA,AMAT,X,RHS,
     2                   DSCALE,DSCALE2)
      END IF
C
C-------SET SMS INNER ITERATION NUMBER (IN_ITER) TO NUMBER OF 
C       PCGU INNER ITERATIONS (innerit)
      IN_ITER = innerit
C
C-------RETURN
      RETURN
C
      END SUBROUTINE PCGU7U1AP
C
C-------ROUTINE TO CALCULATE LORDER AND IORDER FOR REORDERING
C       NOTE - USING GNRCM AND ODRV SUBROUTINES IN THE XMD SOLVER SOURCE CODE
C              SPECIFICALLY IN xmblib.f
      SUBROUTINE SSPCGU_CALC_ORDER( IOUT,IORD,NEQ,NJA,IA,JA,
     2                              LORDER,IORDER )
      USE SMSMODULE,   ONLY: IPRSMS
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IOUT
        INTEGER, INTENT(IN) :: IORD
        INTEGER, INTENT(IN) :: NEQ
        INTEGER, INTENT(IN) :: NJA
        INTEGER, DIMENSION(NEQ+1), INTENT(IN)  :: IA
        INTEGER, DIMENSION(NJA),   INTENT(IN)  :: JA
        INTEGER, DIMENSION(NEQ), INTENT(INOUT) :: LORDER
        INTEGER, DIMENSION(NEQ), INTENT(INOUT) :: IORDER
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        INTEGER :: nsp
        INTEGER, DIMENSION(:), ALLOCATABLE :: iwork0, iwork1
        INTEGER :: iflag
        INTEGER :: i,j
C     + + + PARAMETERS + + + 
        INTEGER, PARAMETER :: IZERO = 0
C     + + + FUNCTIONS + + +
C     + + + FORMATS + + +
2030  FORMAT(1X,A20,1X,6(I6,1X))
2040  FORMAT(1X,20('-'),1X,6(6('-'),1X))
2050  FORMAT(1X,62('-'),/)
C     + + + CODE + + +
        DO n = 1, NEQ
          LORDER(n) = IZERO
          IORDER(n) = IZERO
        END DO
        ALLOCATE ( iwork0(NEQ)  )
        SELECT CASE ( IORD )
          CASE ( 1 )
            ALLOCATE ( iwork1(NEQ) )
            CALL GENRCM(NEQ, NJA, IA, JA,
     2                  LORDER, iwork0, iwork1 )
          CASE ( 2 )
            nsp = 3 * NEQ + 4 * NJA
            ALLOCATE ( iwork1(nsp)  )
            CALL ODRV( IA, JA, LORDER, iwork0, iwork1,
     2                 NEQ, NJA, nsp, iflag )
            IF ( iflag.NE.0 ) THEN
              CALL USTOP('ERROR CREATING MINIMUM DEGREE ORDER'//
     2                   'PERMUTATION ') 
            END IF
        END SELECT
C
C         GENERATE INVERSE OF LORDER
        DO n = 1, NEQ
          IORDER( LORDER(n) ) = n
        END DO
C
C         WRITE SUMMARY OF REORDERING INFORMATION
C         TO LIST FILE
        IF ( IPRSMS.EQ.2 ) THEN
          DO i = 1, NEQ, 6
            WRITE (IOUT,2030) 'ORIGINAL NODE      :',
     2                        (j,j=i,MIN(i+5,NEQ))
            WRITE (IOUT,2040)
            WRITE (IOUT,2030) 'REORDERED INDEX    :',
     2                        (LORDER(j),j=i,MIN(i+5,NEQ))
            WRITE (IOUT,2030) 'REORDERED NODE     :',
     2                        (IORDER(j),j=i,MIN(i+5,NEQ))
            WRITE (IOUT,2050)
         END DO
         END IF
C         DEALLOCATE TEMPORARY STORAGE
        DEALLOCATE ( iwork0, iwork1 )
C---------RETURN
        RETURN
      END SUBROUTINE SSPCGU_CALC_ORDER
C
C-------ROUTINE TO SCALE THE COEFFICIENT MATRIX (AMAT), 
C       THE RHS (B), AND THE ESTIMATE OF X (X)
      SUBROUTINE SPCGU_SCALE(IOPT,ISCL,NEQ,NJA,IA,JA,AMAT,X,B,
     2                       DSCALE,DSCALE2)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IOPT
        INTEGER, INTENT(IN) :: ISCL
        INTEGER, INTENT(IN) :: NEQ
        INTEGER, INTENT(IN) :: NJA
        INTEGER, DIMENSION(NEQ+1), INTENT(IN) :: IA
        INTEGER, DIMENSION(NJA),   INTENT(IN) :: JA
        DOUBLEPRECISION, DIMENSION(NJA),  INTENT(INOUT) :: AMAT
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT) :: X
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT) :: B
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT) :: DSCALE
        DOUBLEPRECISION, DIMENSION(NEQ), INTENT(INOUT)  :: DSCALE2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, n
        INTEGER :: id, jc
        INTEGER :: i0, i1
        DOUBLEPRECISION :: v, c1, c2
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C
C---------SCALE SCALE AMAT, X, AND B
        IF ( IOPT.EQ.0 ) THEN
C-----------SYMMETRIC SCALING
          SELECT CASE ( ISCL )
            CASE ( 1 )
              DO n = 1, NEQ
                id   = IA(n)
                v    = AMAT(id)
                c1   = 1.0D0 / SQRT( ABS( v ) )
                DSCALE(n)  = c1
                DSCALE2(n) = c1
              END DO
C               SCALE AMAT -- AMAT = DSCALE(row) * AMAT(i) * DSCALE2(col)
              DO n = 1, NEQ
                c1 = DSCALE(n)
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  jc = JA(i)
                  c2 = DSCALE2(jc)
                  AMAT(i) = c1 * AMAT(i) * c2 
                END DO
              END DO
C-----------L-2 NORM SCALING
            CASE ( 2 )
C               SCALE EACH ROW SO THAT THE L-2 NORM IS 1
              DO n = 1, NEQ
                c1 = 0.0D0
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  c1 = c1 + AMAT(i) * AMAT(i)
                END DO
                c1 = SQRT( c1 )
                IF ( c1.EQ.0.0D0 ) THEN
                  c1 = 1.0D0
                ELSE
                  c1 = 1.0D0 / c1
                END IF
                DSCALE(n) = c1 
C                 INITIAL SCALING OF AMAT -- AMAT = DSCALE(row) * AMAT(i)              
                DO i = i0, i1
                  AMAT(i) = c1 * AMAT(i)
                END DO
              END DO
C               SCALE EACH COLUMN SO THAT THE L-2 NORM IS 1
              DO n = 1, NEQ
                DSCALE2(n) = 0.0D0
              END DO
              c2 = 0.0D0
              DO n = 1, NEQ
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  jc = JA(i)
                  c2 = AMAT(i)
                  DSCALE2(jc) = DSCALE2(jc) + c2 * c2
                END DO
              END DO
              DO n = 1, NEQ
                c2 = DSCALE2(n)
                IF ( c2.EQ.0.0D0 ) THEN
                  c2 = 1.0D0
                ELSE
                  c2 = 1.0D0 / SQRT( c2 )
                END IF
                DSCALE2(n) = c2
              END DO
C               FINAL SCALING OF AMAT -- AMAT = DSCALE2(col) * AMAT(i)              
              DO n = 1, NEQ
                i0 = IA(n)
                i1 = IA(n+1) - 1
                DO i = i0, i1
                  jc = JA(i)
                  c2 = DSCALE2(jc)
                  AMAT(i) = c2 * AMAT(i)
                END DO
              END DO
          END SELECT
C-----------SCALE X AND B
          DO n = 1, NEQ
            c1    = DSCALE(n)
            c2    = DSCALE2(n)
            X(n)  = X(n) / c2
            B(n)  = B(n) * c1
          END DO
C---------UNSCALE SCALE AMAT, X, AND B
        ELSE
          DO n = 1, NEQ
            c1 = DSCALE(n)
            i0 = IA(n)
            i1 = IA(n+1) - 1
C             UNSCALE AMAT
            DO i = i0, i1
              jc = JA(i)
              c2 = DSCALE2(jc)
              AMAT(i) = ( 1.0D0 / c1 ) * AMAT(i) * ( 1.0D0 / c2 ) 
            END DO
C             UNSCALE X AND B
            c2   = DSCALE2(n)
            X(n) = X(n) * c2
            B(n) = B(n) / c1
          END DO     
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_SCALE
C
C-------ROUTINE TO UPDATE THE PRECONDITIONER
      SUBROUTINE SPCGU_PCU(IOUT,NJA,NEQ,NIAPC,NJAPC,NNZAPC,IPC,RELAX,
     2                     AMAT,IA,JA,APC,IAPC,JAPC,IW,W)
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IOUT
        INTEGER, INTENT(IN) :: NJA
        INTEGER, INTENT(IN) :: NEQ
        INTEGER, INTENT(IN) :: NIAPC
        INTEGER, INTENT(IN) :: NJAPC
        INTEGER, INTENT(IN) :: NNZAPC
        INTEGER, INTENT(IN) :: IPC
        REAL, INTENT(IN) :: RELAX
        DOUBLEPRECISION, DIMENSION(NJA),  INTENT(IN)     :: AMAT
        INTEGER, DIMENSION(NEQ+1), INTENT(IN)    :: IA
        INTEGER, DIMENSION(NJA), INTENT(IN)      :: JA
        DOUBLEPRECISION, DIMENSION(NNZAPC), INTENT(INOUT) :: APC
        INTEGER, DIMENSION(NIAPC+1), INTENT(INOUT) :: IAPC
        INTEGER, DIMENSION(NJAPC), INTENT(INOUT)   :: JAPC
        INTEGER, DIMENSION(NIAPC), INTENT(INOUT)   :: IW
        DOUBLEPRECISION, DIMENSION(NIAPC), INTENT(INOUT) :: W
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: izero
        DOUBLEPRECISION :: delta
C     + + + FUNCTIONS + + +
C     + + + FORMATS + + +
!2000    FORMAT (/,' MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT.  CHECK',
!     &          ' INPUT FILES.',/,' -- STOP EXECUTION (SPCGU_PCU)')
2000    FORMAT (/,' MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT.',
     &          /,' ADDING SMALL VALUE TO PIVOT (SPCGU_PCU)')
C     + + + CODE + + +
        izero = 0
        delta = 0.0D0
        SELECT CASE(IPC)
C           NO PRE-CONDITIONER
          CASE (0)
C           JACOBI PRE-CONDITIONER
          CASE (1)
            CALL SPCGU_PCJ(NJA,NEQ,AMAT,APC,IA,JA)
C           ILU0 AND MILU0
          CASE (2,3)
            LUPC: DO
              CALL SPCGU_PCILU0(NJA,NEQ,AMAT,IA,JA,
     2                          APC,IAPC,JAPC,IW,W,
     3                          RELAX,izero,delta)
              IF ( izero.LT.1 ) THEN
                EXIT LUPC
              END IF
              delta = 1.5D0 * delta + 0.001
              izero = 0
              IF ( delta.GT.0.5D0 ) THEN
                WRITE(IOUT,2000)
                delta = 0.5D0
                izero = 2
              END IF
            END DO LUPC
C           ADDITIONAL PRECONDITIONERS - ILUT, etc.
        END SELECT
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCU
C
C-------JACOBI PRECONDITIONER - INVERSE OF DIAGONAL 
      SUBROUTINE SPCGU_PCJ(NJA,NEQ,AMAT,APC,IA,JA)
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NJA
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NJA),  INTENT(IN)      :: AMAT
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT)   :: APC
        INTEGER, DIMENSION(NEQ+1), INTENT(IN) :: IA
        INTEGER, DIMENSION(NJA),   INTENT(IN) :: JA
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, n
        INTEGER :: ic0, ic1
        INTEGER :: id
        DOUBLEPRECISION :: tv
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NEQ
            ic0 = IA(n)
            ic1 = IA(n+1) - 1
            id = IA(n)
            DO i = ic0, ic1
              IF ( JA(i).EQ.n ) THEN
                id = i
                EXIT
              END IF
            END DO
            tv  = AMAT(id)
            IF ( ABS( tv ).GT.DZERO ) tv = DONE / tv
            APC(n) = tv
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCJ

      SUBROUTINE SPCGU_JACA(NEQ,A,D1,D2)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT) :: D2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NEQ
          tv     = A(n) * D1(n)
          D2(n) = tv
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_JACA

      SUBROUTINE SPCGU_PCILU0(NJA,NEQ,AMAT,IA,JA,
     2                        APC,IAPC,JAPC,IW,W,
     3                        RELAX,IZERO,DELTA)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NJA
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NJA),  INTENT(IN)     :: AMAT
        INTEGER, DIMENSION(NEQ+1), INTENT(IN)    :: IA
        INTEGER, DIMENSION(NJA), INTENT(IN)      :: JA
        DOUBLEPRECISION, DIMENSION(NJA), INTENT(INOUT)   :: APC
        INTEGER, DIMENSION(NEQ+1), INTENT(INOUT) :: IAPC
        INTEGER, DIMENSION(NJA), INTENT(INOUT)   :: JAPC
        INTEGER, DIMENSION(NEQ), INTENT(INOUT)   :: IW
        DOUBLEPRECISION, DIMENSION(NEQ), INTENT(INOUT)   :: W
        REAL, INTENT(IN) :: RELAX
        INTEGER, INTENT(INOUT) :: IZERO
        DOUBLEPRECISION, INTENT(IN) :: DELTA
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1
        INTEGER :: iic0, iic1
        INTEGER :: iu, iiu
        INTEGER :: j, n
        INTEGER :: jj
        INTEGER :: jcol, jw
        INTEGER :: jjcol
        DOUBLEPRECISION :: drelax
        DOUBLEPRECISION :: sd1
        DOUBLEPRECISION :: tl
        DOUBLEPRECISION :: rs
        DOUBLEPRECISION :: d
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
        DOUBLEPRECISION, PARAMETER :: DTINY = 1.0D-06
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        drelax = DBLE( RELAX )
        DO n = 1, NEQ
          IW(n)  = 0
          W(n)   = DZERO
        END DO
        MAIN: DO n = 1, NEQ
          ic0 = IA(n)
          ic1 = IA(n+1) - 1
          DO j = ic0, ic1
            jcol      = JA(j)
            IW(jcol) = 1
            W(jcol) = W(jcol) + AMAT(j)
          END DO
          ic0 = IAPC(n)
          ic1 = IAPC(n+1) - 1
          iu  = JAPC(n)
          rs   = DZERO
          LOWER: DO j = ic0, iu-1
            jcol     = JAPC(j)
            iic0     = IAPC(jcol) 
            iic1     = IAPC(jcol+1) - 1
            iiu      = JAPC(jcol)
            tl       = W(jcol) * APC(jcol)
            W(jcol) = tl
            DO jj = iiu, iic1
              jjcol = JAPC(jj)
              jw    = IW(jjcol)
              IF ( jw.NE.0 ) THEN
                W(jjcol) = W(jjcol) - tl * APC(jj)
              ELSE
                rs = rs + tl * APC(jj)
              END IF
            END DO
          END DO LOWER
C           DIAGONAL - CALCULATE INVERSE OF DIAGONAL FOR SOLUTION
          d   = W(n)
          tl  = ( DONE + DELTA ) * d - ( drelax * rs )
C-----------ENSURE THAT THE SIGN OF THE DIAGONAL HAS NOT CHANGED AND IS NOT ZERO
          sd1 = SIGN(d,tl)
          IF ( sd1.NE.d ) THEN
C             USE SMALL VALUE IF DIAGONAL SCALING IS NOT EFFECTIVE FOR ELIMINATING
C             PIVOTS THAT CHANGE THE SIGN OF THE DIAGONAL
            IF ( IZERO.GT.1 ) THEN
              tl = SIGN(DTINY,d)
C             DIAGONAL SCALING CONTINUES TO BE EFFECTIVE
            ELSE
              IZERO = 1
              EXIT MAIN
            END IF
          END IF
          IF ( ABS(tl).EQ.DZERO ) THEN
C             USE SMALL VALUE IF DIAGONAL SCALING IS NOT EFFECTIVE FOR ELIMINATING
C             ZERO PIVOTS
            IF ( IZERO.GT.1 ) THEN
              tl = SIGN(DTINY,d)
C             DIAGONAL SCALING CONTINUES TO BE EFFECTIVE FOR ELIMIATING ZERO PIVOTS
            ELSE
              IZERO = 1
              EXIT MAIN
            END IF
          END IF
          APC(n) = DONE / tl
C           RESET POINTER FOR IW TO ZERO
          IW(n) = 0
          W(n)  = DZERO
          DO j = ic0, ic1
            jcol = JAPC(j)
            APC(j) = W(jcol)
            IW(jcol) = 0
            W(jcol) = DZERO
          END DO
      END DO MAIN
!      
!        do n = 1, NEQ
!          ic0 = IAPC(n)
!          ic1 = IAPC(n+1) - 1
!          iu  = JAPC(n)
!          do j = ic0,iu-1
!            jcol = JAPC(j)
!            write (10991,'(i10,i10,g15.7)' ) n, jcol, 
!     2                                       APC(j)
!          end do
!          write (10991,'(i10,i10,g15.7)' ) n, n, 
!     2                                     APC(n)
!          do j = iu,ic1
!            jcol = JAPC(j)
!            write (10991,'(i10,i10,g15.7)' ) n, jcol, 
!     2                                       APC(j)
!          end do
!        end do
C
C---------RESET IZERO IF SUCCESSFUL COMPLETION OF MAIN
        IZERO = 0
C
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCILU0

      SUBROUTINE SPCGU_ILU0A(NJA,NEQ,APC,IAPC,JAPC,R,D)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NJA
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NJA),  INTENT(INOUT)  :: APC
        INTEGER, DIMENSION(NEQ+1), INTENT(IN) :: IAPC
        INTEGER, DIMENSION(NJA), INTENT(IN)   :: JAPC
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(IN)     :: R
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT)  :: D
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1
        INTEGER :: iu
        INTEGER :: jcol
        INTEGER :: j, n
        DOUBLEPRECISION :: tv
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C         FORWARD SOLVE - APC * D = R
        FORWARD: DO n = 1, NEQ
          tv   = R(n)
          ic0 = IAPC(n)
          ic1 = IAPC(n+1) - 1
          iu  = JAPC(n) - 1
          LOWER: DO j = ic0, iu
            jcol = JAPC(j)
            tv    = tv - APC(j) * D(jcol)
          END DO LOWER
          D(n) = tv
        END DO FORWARD
C         BACKWARD SOLVE - D = D / U
        BACKWARD: DO n = NEQ, 1, -1
          ic0 = IAPC(n)
          ic1 = IAPC(n+1) - 1
          iu  = JAPC(n)
          tv   = D(n)
          UPPER: DO j = iu, ic1
            jcol = JAPC(j)
            tv    = tv - APC(j) * D(jcol)
          END DO UPPER
C           COMPUTE D FOR DIAGONAL - D = D / U
          D(n) =  tv * APC(n)
        END DO BACKWARD
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_ILU0A

      SUBROUTINE SPCGU_CG(ICNVG,ITMAX,INNERIT)
        USE SOLUTIONMODULE, ONLY: neq,nja,ia,ja,amat,rhs,x
        USE PCGUMODULE
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(INOUT) :: ICNVG
        INTEGER, INTENT(IN)    :: ITMAX
        INTEGER, INTENT(INOUT) :: INNERIT
C       + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        INTEGER :: iiter
        DOUBLEPRECISION :: dhclose, drclose
        DOUBLEPRECISION :: tv
        DOUBLEPRECISION :: deltax
        DOUBLEPRECISION :: rmax
        DOUBLEPRECISION :: l2norm
        DOUBLEPRECISION :: rcnvg
        DOUBLEPRECISION :: alpha, beta
        DOUBLEPRECISION :: rho, rho0
        DOUBLEPRECISION :: machprec
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C         + + + FUNCTIONS + + +
        DOUBLEPRECISION :: SPCGU_DP
C
C         + + + CODE + + +
        INNERIT  = 0
        machprec = EPSILON( DZERO )
        dhclose  = DBLE( HCLOSEPCGU )
        drclose  = DBLE( RCLOSEPCGU )
C
C-------INNER ITERATION          
        INNER: DO iiter = 1, itmax
           INNERIT = INNERIT + 1 
           NITERC  = NITERC  + 1 
C----------APPLY PRECONDITIONER
          SELECT CASE (IPC)
C             NO PRECONDITIONER
            CASE (0)
              DO n = 1, NEQ
                Z(n) = D(n)
              END DO
C             JACOBI PRECONDITIONER
            CASE (1)
              CALL SPCGU_JACA(NEQ,APC,D,Z)
            CASE (2,3)
              CALL SPCGU_ILU0A(NJA,NEQ,APC,IAPC,JAPC,D,Z)
          END SELECT
          rho = SPCGU_DP( NEQ, D, Z )
C-----------COMPUTE DIRECTIONAL VECTORS
          IF (IITER.EQ.1) THEN
            DO n = 1, NEQ
              P(n) = Z(n)
            END DO
          ELSE
            beta = rho / rho0
            DO n = 1, NEQ
              P(n) = Z(n) + beta * P(n)
            END DO
          END IF
C-----------COMPUTE ITERATES
C           UPDATE qc
          CALL SPCGU_MV(NJA,NEQ,A0,P,Q,IA0,JA0)

          alpha = rho / SPCGU_DP( NEQ, P, Q )
C-----------UPDATE X AND RESIDUAL
          deltax = DZERO
          rmax   = DZERO
          l2norm = DZERO
          DO n = 1, NEQ
            tv      = alpha * P(n)
            X(n)  = X(n) + tv
            deltax = MAX( ABS(tv), deltax )
            tv      = D(n)
            tv      = tv - alpha * Q(n)
            D(n)  = tv
            rmax   = MAX( ABS(tv), rmax )
            l2norm = l2norm + tv * tv
          END DO
          l2norm = SQRT(l2norm)
C-----------TEST FOR SOLVER CONVERGENCE
          IF ( ICNVGOPT.EQ.1 ) THEN
            rcnvg = l2norm
          ELSE
            rcnvg = rmax
          END IF
          CALL SPCGU_TESTCNVG( ICNVGOPT,ICNVG,
     2                         deltax,rcnvg,
     3                         L2NORM0,EPFACT,dhclose,drclose )
C           CHECK FOR EXACT SOLUTION
          IF ( rcnvg.EQ.DZERO ) ICNVG = 1
          IF ( ICNVG.EQ.1 ) EXIT INNER
C-----------CHECK THAT CURRENT AND PREVIOUS rho ARE DIFFERENT
          IF ( ABS( rho - rho0 ).LT.machprec ) THEN
            rho0 = rho
            EXIT INNER
          END IF
C-----------SAVE CURRENT INNER ITERATES
          rho0 = rho
        END DO INNER
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_CG

      SUBROUTINE SPCGU_BCGS(ICNVG,ITMAX,INNERIT)
        USE SOLUTIONMODULE, ONLY: neq,nja,ia,ja,amat,rhs,x
        USE PCGUMODULE
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(INOUT) :: ICNVG
        INTEGER, INTENT(IN)    :: ITMAX
        INTEGER, INTENT(INOUT) :: INNERIT
C       + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        INTEGER :: iiter
        DOUBLEPRECISION :: dhclose, drclose
        DOUBLEPRECISION :: tv
        DOUBLEPRECISION :: deltax
        DOUBLEPRECISION :: rmax
        DOUBLEPRECISION :: l2norm
        DOUBLEPRECISION :: rcnvg
        DOUBLEPRECISION :: alpha, alpha0 
        DOUBLEPRECISION :: beta
        DOUBLEPRECISION :: rho, rho0
        DOUBLEPRECISION :: omega, omega0
        DOUBLEPRECISION :: machprec
        DOUBLEPRECISION :: numer, denom
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
        DOUBLEPRECISION, PARAMETER :: DTWO  = 2.0D0
C         + + + FUNCTIONS + + +
        DOUBLEPRECISION :: SPCGU_DP
!        DOUBLEPRECISION :: SPCGU_RES
!        DOUBLEPRECISION :: SPCGU_L2NORM
C
C         + + + CODE + + +
        INNERIT  = 0
        machprec = EPSILON( DZERO )
        dhclose  = DBLE( HCLOSEPCGU )
        drclose  = DBLE( RCLOSEPCGU )
        
        alpha = DZERO
        beta  = DZERO
        rho   = DZERO
        rho0  = DZERO
C        
C-------SAVE INITIAL RESIDUAL
        DO n = 1, NEQ
          DHAT(n) = D(n)
        END DO
C
C-------INNER ITERATION          
        INNER: DO iiter = 1, itmax
           INNERIT = INNERIT + 1 
           NITERC = NITERC + 1 
C----------CALCULATE rho
          rho = SPCGU_DP( NEQ, DHAT, D )
C-----------COMPUTE DIRECTIONAL VECTORS
          IF (IITER.EQ.1) THEN
            DO n = 1, NEQ
              P(n) = D(n)
            END DO
          ELSE
            beta = ( rho / rho0 ) * ( alpha0 / omega0 )
            DO n = 1, NEQ
              P(n) = D(n) + beta * ( P(N) - omega * V(n) )
            END DO
          END IF
C----------APPLY PRECONDITIONER TO UPDATE PHAT
          SELECT CASE (IPC)
C             NO PRECONDITIONER
            CASE (0)
              DO n = 1, NEQ
                PHAT(n) = P(n)
              END DO
C             JACOBI PRECONDITIONER
            CASE (1)
              CALL SPCGU_JACA(NEQ,APC,P,PHAT)
            CASE (2,3)
              CALL SPCGU_ILU0A(NJA,NEQ,APC,IAPC,JAPC,P,PHAT)
          END SELECT
C-----------COMPUTE ITERATES
C           UPDATE V WITH A AND PHAT
          CALL SPCGU_MV(NJA,NEQ,A0,PHAT,V,IA0,JA0)
C           UPDATE alpha WITH DHAT AND V
          denom = SPCGU_DP( NEQ, DHAT, V )
          !denom = denom + SIGN(dsmall,denom)
          alpha = rho / denom
C-----------UPDATE Q
          DO n = 1, NEQ
            Q(n) = D(n) - alpha * V(n)
          END DO
C-----------CALCULATE INFINITY NORM OF Q - TEST FOR TERMINATION
C           TERMINATE IF rmax IS LESS THAN MACHINE PRECISION (machprec)
          rmax = DZERO
          DO n = 1, NEQ
              tv = Q(n)
              IF ( ISCL.NE.0 ) tv = tv / DSCALE(n)
              IF ( ABS(tv).GT.ABS(rmax) ) rmax = tv
          END DO
          IF ( ABS(rmax).LE.machprec ) THEN
            deltax = DZERO
            DO n = 1, NEQ
              tv      = alpha * PHAT(n)
              IF ( ISCL.NE.0 ) THEN
                tv = tv * DSCALE(n)
              END IF
              X(n)  = X(n) + tv
              IF ( ABS(tv).GT.ABS(deltax) ) deltax = tv
            END DO
            CALL SPCGU_TESTCNVG( ICNVGOPT,ICNVG,
     2                           deltax,rmax,
     3                           rmax,EPFACT,dhclose,drclose )
            IF ( ICNVG.EQ.1 ) EXIT INNER
          END IF
C-----------APPLY PRECONDITIONER TO UPDATE QHAT
          SELECT CASE (IPC)
C             NO PRECONDITIONER
            CASE (0)
              DO n = 1, NEQ
                QHAT(n) = Q(n)
              END DO
C             JACOBI PRECONDITIONER
            CASE (1)
              CALL SPCGU_JACA(NEQ,APC,Q,QHAT)
            CASE (2,3)
              CALL SPCGU_ILU0A(NJA,NEQ,APC,IAPC,JAPC,Q,QHAT)
          END SELECT
C           UPDATE T WITH A AND QHAT
          CALL SPCGU_MV(NJA,NEQ,A0,QHAT,T,IA0,JA0)
C-----------UPDATE omega
          numer = SPCGU_DP( NEQ, T, Q )
          denom = SPCGU_DP( NEQ, T, T )
          denom = denom + SIGN(machprec,denom)
          omega = numer / denom
C-----------UPDATE X AND RESIDUAL
          deltax = DZERO
          rmax   = DZERO
          l2norm = DZERO
          DO n = 1, NEQ
C-------------X AND DX            
            tv      = alpha * PHAT(n) + omega * QHAT(n)
            X(n)  = X(n) + tv
            IF ( ISCL.NE.0 ) THEN
              tv = tv * DSCALE(n)
            END IF
            IF ( ABS(tv).GT.ABS(deltax) ) deltax = tv
C-------------RESIDUAL
            tv      = Q(n) - omega * T(n)
            D(n)  = tv
            IF ( ISCL.NE.0 ) THEN
              tv = tv / DSCALE(n)
            END IF
            IF ( ABS(tv).GT.ABS(rmax) ) rmax = tv
            l2norm = l2norm + tv * tv
          END DO
C-----------TEST FOR SOLVER CONVERGENCE
          IF ( ICNVGOPT.EQ.1 ) THEN
            rcnvg = l2norm
          ELSE
            rcnvg = rmax
          END IF
          CALL SPCGU_TESTCNVG( ICNVGOPT,ICNVG,
     2                         deltax,rcnvg,
     3                         L2NORM0,EPFACT,dhclose,drclose )
C           CHECK FOR EXACT SOLUTION
          IF ( rcnvg.EQ.DZERO ) ICNVG = 1
          IF ( ICNVG.EQ.1 ) EXIT INNER
C-----------CHECK THAT CURRENT AND PREVIOUS rho, alpha, AND omega ARE DIFFERENT
          IF ( ABS( rho - rho0 ).LT.machprec ) THEN
            rho0 = rho
            EXIT INNER
          END IF
          IF ( ABS( alpha - alpha0 ).LT.machprec ) THEN
            alpha0 = alpha
            EXIT INNER
          END IF
          IF ( ABS( omega - omega0 ).LT.machprec ) THEN
            omega0 = omega
            EXIT INNER
      END IF
C-----------SAVE CURRENT INNER ITERATES
          rho0   = rho
          alpha0 = alpha
          omega0 = omega
        END DO INNER
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_BCGS
C
C---------TEST FOR SOLVER CONVERGENCE
        SUBROUTINE SPCGU_TESTCNVG( Icnvgopt,Icnvg,
     2                             Hmax,Rmax,
     3                             Rmax0,Epfact,Hclose,Rclose )
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN)         :: Icnvgopt
        INTEGER, INTENT(INOUT)      :: Icnvg
        DOUBLEPRECISION, INTENT(IN) :: Hmax
        DOUBLEPRECISION, INTENT(IN) :: Rmax
        DOUBLEPRECISION, INTENT(IN) :: Rmax0
        DOUBLEPRECISION, INTENT(IN) :: Epfact
        DOUBLEPRECISION, INTENT(IN) :: Hclose
        DOUBLEPRECISION, INTENT(IN) :: Rclose
C     + + + LOCAL DEFINITIONS + + +
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        IF ( Icnvgopt.EQ.0 ) THEN
          IF ( ABS(Hmax).LE.Hclose .AND. ABS(Rmax).LE.Rclose ) THEN
            Icnvg = 1
          END IF
        ELSE
          IF ( Rmax.LE.Rclose ) THEN
            Icnvg = 1
          ELSE IF ( Rmax.LE.Rmax0*Epfact ) THEN
            Icnvg = 1
          END IF
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_TESTCNVG
C
C---------GENERATE IAPC AND JAPC FROM IA AND JA
C         JAPC(1:NEQ) HAS THE POSITION OF THE UPPER ENTRY FOR A ROW
C         JAPC(NEQ+1:NJA) IS THE COLUMN POSITION FOR ENTRY
C         APC(1:NEQ) PRECONDITIONED INVERSE OF THE DIAGONAL
C         APC(NEQ+1:NJA) PRECONDITIONED ENTRIES FOR OFF DIAGONALS
        SUBROUTINE SPCGU_PCCRS(NEQ,NJA,IA,JA,
     2                         IAPC,JAPC)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN)         :: NEQ
        INTEGER, INTENT(IN)         :: NJA
        INTEGER, DIMENSION(NEQ+1), INTENT(IN)    :: IA
        INTEGER, DIMENSION(NJA), INTENT(IN)      :: JA
        INTEGER, DIMENSION(NEQ+1), INTENT(INOUT) :: IAPC
        INTEGER, DIMENSION(NJA), INTENT(INOUT)   :: JAPC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n, j
        INTEGER :: i0, i1
        INTEGER :: nlen
        INTEGER :: ic,ip
        INTEGER :: jcol
        INTEGER, DIMENSION(:), ALLOCATABLE :: iarr
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        ip = NEQ + 1
        DO n = 1, NEQ
          i0 = IA(n)
          i1 = IA(n+1) - 1
          nlen = i1 - i0
          ALLOCATE( iarr(nlen) )
          ic = 0
          DO j = i0, i1
            jcol = JA(j)
            IF ( jcol.EQ.n ) CYCLE
            ic = ic + 1
            iarr(ic) = jcol
          END DO
          CALL SPCGU_ISORT(nlen,iarr)
          IAPC(n) = ip
          DO j = 1, nlen
            jcol = iarr(j)
            JAPC(ip) = jcol
            ip = ip + 1
          END DO
          DEALLOCATE(iarr)
        END DO
        IAPC(NEQ+1) = NJA + 1
C---------POSITION OF THE FIRST UPPER ENTRY FOR ROW         
        DO n = 1, NEQ
          i0 = IAPC(n)
          i1 = IAPC(n+1) - 1
          JAPC(n) = IAPC(n+1)
          DO j = i0, i1
            jcol = JAPC(j)
            IF ( jcol.GT.n ) THEN
              JAPC(n) = j
              EXIT
            END IF
          END DO
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_PCCRS
C
C-------SIMPLE INPLACE SORTING ROUTINE FOR AN INTEGER ARRAY      
      SUBROUTINE SPCGU_ISORT(NVAL,IARRAY)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER,INTENT(IN) :: NVAL
        INTEGER,DIMENSION(NVAL),INTENT(INOUT) :: IARRAY
C     + + + LOCAL DEFINITIONS + + +
        integer :: i, j, itemp
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO i = 1, NVAL-1
            DO j = i+1, NVAL
                if(IARRAY(i).GT.IARRAY(j)) then
                    itemp = IARRAY(j)
                    IARRAY(j) = IARRAY(i)
                    IARRAY(i) = itemp
                END IF
            END DO
        END DO
      END SUBROUTINE SPCGU_ISORT
      
C       COPY ONE DOUBLEPRECISION VECTOR TO ANOTHER
      SUBROUTINE SSPCGU_DCOPY(NR, V, R)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NR
        DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)    :: V
        DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT) :: R
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NR
          R(n) = V(n)
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SSPCGU_DCOPY

C       COPY ONE INTEGER VECTOR TO ANOTHER
      SUBROUTINE SSPCGU_ICOPY(NR, V, R)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NR
        INTEGER, DIMENSION(NR), INTENT(IN)    :: V
        INTEGER, DIMENSION(NR), INTENT(INOUT) :: R
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NR
          R(n) = V(n)
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SSPCGU_ICOPY
      
      SUBROUTINE SPCGU_MV(NJA,NEQ,A,D1,D2,IA,JA)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NJA
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NJA),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(INOUT) :: D2
        INTEGER, DIMENSION(NEQ+1), INTENT(IN) :: IA
        INTEGER, DIMENSION(NJA), INTENT(IN)   :: JA
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1
        INTEGER :: icol
        INTEGER :: m, n
        DOUBLEPRECISION :: tv
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NEQ
C           ADD DIAGONAL AND OFF-DIAGONAL TERMS
          tv     = DZERO
          ic0   = IA(n)
          ic1   = IA(n+1)-1
          DO m = ic0, ic1
            icol = JA(m)
            tv  = tv + A(m) * D1(icol)
          END DO
          D2(n) = tv
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SPCGU_MV

      DOUBLEPRECISION FUNCTION SPCGU_DP(NEQ,A,B) RESULT(C)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NEQ),  INTENT(IN)    :: B
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
C     + + + PARAMETERS + + + 
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        C = DZERO
        DO n = 1, NEQ
          C = C + A(n) * B(n)
        END DO
C---------RETURN
        RETURN
      END FUNCTION SPCGU_DP
C
C-------CALCULATE THE L2 NORM OF A VECTOR
      DOUBLEPRECISION FUNCTION SPCGU_L2NORM(NEQ,V) RESULT(value)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NEQ
        DOUBLEPRECISION, DIMENSION(NEQ), INTENT(IN) :: V
C     + + + LOCAL DEFINITIONS + + +
        DOUBLEPRECISION :: dotp
C     + + + FUNCTIONS + + +
        DOUBLEPRECISION :: SPCGU_DP
C     + + + CODE + + +
        dotp  = SPCGU_DP(NEQ,V,V)
        value = SQRT(dotp)
C---------RETURN
        RETURN
      END FUNCTION SPCGU_L2NORM
C
C-------CALCULATE THE RESIDUAL FOR A NODE USING CURRENT VALUES OF
C       AMAT, X, AND THE RHS FOR THE NODE PASSED TO FUNCTION     
      DOUBLEPRECISION FUNCTION SPCGU_RES(N,IA,JA,AMAT,X,B) 
     2  RESULT(value)
        USE SOLUTIONMODULE, ONLY: neq,nja
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: N
        INTEGER, DIMENSION(NEQ+1), INTENT(IN) :: IA
        INTEGER, DIMENSION(NJA), INTENT(IN)   :: JA
        DOUBLEPRECISION, DIMENSION(NJA), INTENT(IN)   :: AMAT
        DOUBLEPRECISION, DIMENSION(NEQ), INTENT(IN)   :: X
        DOUBLEPRECISION, INTENT(IN)   :: B
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: j
        INTEGER :: i0, i1
        INTEGER :: jcol
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        value = 0.0
        i0 = IA(N)
        i1 = IA(N+1) - 1
        DO j = i0, i1
          jcol = JA(j)
          value = value + AMAT(j) * X(jcol)
        END DO
        value = B - value
C---------RETURN
        RETURN
      END FUNCTION SPCGU_RES
C
C-------ROUTINES FROM SPARSKIT TO PERMUTATE A LINEAR SYSTEM OF EQUATIONS
C       IN ORDER TO REORDER THE MATRIX TO MINIMIZE THE BANDWIDTH USING
C       THE REVERSE CUTHILL MCKEE ALGORITHM
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! This routine permutes the rows and columns of a matrix stored in CSR
! format. i.e., it computes P A Q, where P, Q are permutation matrices.
! P maps row i into row perm(i) and Q maps column j into column qperm(j):
!      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
! In the particular case where Q is the transpose of P (symmetric
! permutation of A) then qperm is not needed.
! note that qperm should be of length ncol (number of columns) but this
! is not checked.
!-----------------------------------------------------------------------
! Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n       = dimension of the matrix
! a, ja,
!    ia = input matrix in a, ja, ia format
! perm       = integer array of length n containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix -- also the destination of column i in case
!         permutation is symmetric (job .le. 2)
!
! qperm      = same thing for the columns. This should be provided only
!         if job=3 or job=4, i.e., only in the case of a nonsymmetric
!        permutation of rows and columns. Otherwise qperm is a dummy
!
! job      = integer indicating the work to be done:
! * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
!             job = 1      permute a, ja, ia into ao, jao, iao
!             job = 2 permute matrix ignoring real values.
! * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
!             job = 3      permute a, ja, ia into ao, jao, iao
!             job = 4 permute matrix ignoring real values.
!
! on return:
!-----------
! ao, jao, iao = input matrix in a, ja, ia format
!
! in case job .eq. 2 or job .eq. 4, a and ao are never referred to
! and can be dummy arguments.
! Notes:
!-------
!  1) algorithm is in place
!  2) column indices may not be sorted on return even  though they may be
!     on entry.
!----------------------------------------------------------------------c
! local variables
      integer locjob, mod
!
!     locjob indicates whether or not real values must be copied.
!
      locjob = mod(job,2)
!
! permute rows first
!
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
!
! then permute columns
!
      locjob = 0
!
      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob)
      else
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob)
      endif
!
      return
!-------end-of-dperm----------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.
! rperm  computes B = P A  where P is a permutation matrix.
! the permutation P is defined through the array perm: for each j,
! perm(j) represents the destination row number of row number j.
! Youcef Saad -- recoded Jan 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n       = dimension of the matrix
! a, ja, ia = input matrix in csr format
! perm       = integer array of length nrow containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job .ne. 1 :  ignore real values.
!                     (in which case arrays a and ao are not needed nor
!                      used).
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
! note :
!        if (job.ne.1)  then the arrays a and ao are not used.
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1)
!
!     determine pointers for output matix.
!
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
!
! get pointers from lengths
!
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
!
! copying
!
      do 100 ii=1,nrow
!
! old row = ii  -- new row = iperm(ii) -- ko = new pointer
!
         ko = iao(perm(ii))
         do 60 k=ia(ii), ia(ii+1)-1
            jao(ko) = ja(k)
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
!
      return
!---------end-of-rperm -------------------------------------------------
!-----------------------------------------------------------------------
      end

!-----------------------------------------------------------------------
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*8 a(*), ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the columns of a matrix a, ja, ia.
! the result is written in the output matrix  ao, jao, iao.
! cperm computes B = A P, where  P is a permutation matrix
! that maps column j into column perm(j), i.e., on return
!      a(i,j) becomes a(i,perm(j)) in new matrix
! Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow       = row dimension of the matrix
!
! a, ja, ia = input matrix in csr format.
!
! perm      = integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job .ne. 1 :  ignore real values ao and ignore iao.
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
!
! Notes:
!-------
! 1. if job=1 then ao, iao are not used.
! 2. This routine is in place: ja, jao can be the same.
! 3. If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return.
!
!----------------------------------------------------------------------c
! local parameters:
      integer k, i, NJA
!
      NJA = ia(nrow+1)-1
      do 100 k=1,NJA
         jao(k) = perm(ja(k))
 100  continue
!
!     done with ja array. return if no need to touch values.
!
      if (job .ne. 1) return
!
! else get new pointers -- and copy values too.
!
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
!
      do 2 k=1, NJA
         ao(k) = a(k)
 2    continue
!
      return
!---------end-of-cperm--------------------------------------------------
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine vperm (n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of a real vector x 
! according to the permutation array perm(*), i.e., on return, 
! the vector x satisfies,
!
!	x(perm(j)) :== x(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! x	= input vector
!
! on return:
!---------- 
! x	= vector x permuted according to x(perm(*)) :=  x(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables 
      real*8 tmp, tmp1
!
      init      = 1
      tmp       = x(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
!     
! loop
! 
 6    k = k+1
!
! save the chased element --
! 
      tmp1      = x(ii) 
      x(ii)     = tmp
      next      = perm(ii) 
      if (next .lt. 0 ) goto 65
!     
! test for end 
!
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
!
! end loop 
!
      goto 6
!
! reinitilaize cycle --
!
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp       = x(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
!     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
!     
      return
!-------------------end-of-vperm--------------------------------------- 
!-----------------------------------------------------------------------
      end
!
      SUBROUTINE SET_PCGUINPUT(IFDPARAM)
      USE PCGUMODULE, ONLY:  IPC,ISCL,IORD,RCLOSEPCGU,RELAXPCGU,ILINMETH
      INTEGER IFDPARAM
C Simple option
      SELECT CASE ( IFDPARAM )
      CASE(1)
        ILINMETH=1
        IPC = 2
        ISCL = 0
        IORD = 0
        RCLOSEPCGU = 1.0e-1
        RELAXPCGU = 0.0
C Moderate
      CASE(2)
        ILINMETH=2
        IPC = 3
        ISCL = 0
        IORD = 0
        RCLOSEPCGU = 1.0e-1
        RELAXPCGU = 0.97
C Complex
      CASE(3)
        ILINMETH=2
        IPC = 3
        ISCL = 0
        IORD = 0
        RCLOSEPCGU = 1.0e-1
        RELAXPCGU = 0.97
      END SELECT
      RETURN
      END
      