      SUBROUTINE SMS7U1AR(IN)

      USE SOLUTION, ONLY: neq,nja,ia,ja,amat,rhs,x
      USE SimModule
      USE SMSMODULE
      USE XMDMODULE, ONLY: IACL
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
      INTRINSIC INT
      EXTERNAL URDCOM, URWORD, UPARLSTAL
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: IN
!     ------------------------------------------------------------------
!     LOCAL VARIABLES
!     ------------------------------------------------------------------
      INTEGER :: LLOC, ISTART, ISTOP
      INTEGER :: I, K, N
      INTEGER :: IFDPARAM, MXVL, NPP
      INTEGER :: IPCGUM
      CHARACTER(LEN=200) :: LINE
      REAL :: R, HCLOSEDUM, HICLOSEDUM, THETADUM, AMOMENTDUM
      REAL :: AKAPPADUM, GAMMADUM, BREDUCDUM,BTOLDUM,RESLIMDUM
      INTEGER :: ISYMFLG=1

!     ------------------------------------------------------------------
!
C1------IDENTIFY PACKAGE AND INITIALIZE.
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'SMS -- SPARSE MATRIX SOLVER PACKAGE, VERSION 7',
     &', 5/2/2005',/,9X,'INPUT READ FROM UNIT',I3)
      ALLOCATE (HCLOSE, HICLOSE,BIGCHOLD,BIGCH)
      ALLOCATE (RELAXOLD)
      ALLOCATE (RES_PREV)
      ALLOCATE (RES_NEW)
      ALLOCATE (RES_IN)
      ALLOCATE (IBCOUNT)
      ALLOCATE (ICNVG)
      ALLOCATE (ITER1,THETA,MXITER,LINMETH,NONMETH,IPRSMS)
      ALLOCATE (AKAPPA,GAMMA,AMOMENTUM,BREDUC,BTOL,RES_LIM,
     &  NUMTRACK,IBFLAG)
C
      I = 1
      NUMTRACK = 0
      BTOL = 0
      BREDUC = 0.
      RES_LIM = 0.
      IBFLAG = 0
! CHECK IF DEFAULT SOLVER VALUES WILL BE USED
      LLOC = 1
      IFDPARAM = 0
      CALL URDCOM(IN, IOUT, LINE)
      NPP = 0
      MXVL = 0
      CALL UPARLSTAL(IN,IOUT,LINE,NPP,MXVL)
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SIMPLE') THEN
        IFDPARAM = 1
         WRITE(IOUT,21)
   21    FORMAT(1X,'SIMPLE OPTION:',/,
     1     1X,'DEFAULT SOLVER INPUT VALUES FOR FAST SOLUTIONS')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'MODERATE') THEN
         IFDPARAM=2
         WRITE(IOUT,23)
   23    FORMAT(1X,'MODERATE OPTION:',/,1X,'DEFAULT SOLVER',
     1         ' INPUT VALUES REFLECT MODERETELY NONLINEAR MODEL')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'COMPLEX') THEN
         IFDPARAM=3
         WRITE(IOUT,25)
   25    FORMAT(1X,'COMPLEX OPTION:',/,1X,'DEFAULT SOLVER',
     1 ' INPUT VALUES REFLECT STRONGLY NONLINEAR MODEL')
      ELSE
        BACKSPACE IN
        WRITE(IOUT,'(A)') ' ALL SOLVER INPUT DATA WILL BE READ ',
     +                     'FROM THE SOLVER INPUT FILE. '
      END IF
C2------READ NONLINEAR ITERATION PARAMETERS AND LINEAR SOLVER SELECTION INDEX
      LLOC = 1
      CALL URDCOM(IN, IOUT, LINE)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, HCLOSEDUM, IOUT, IN)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, HICLOSEDUM, IOUT, IN)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, MXITER, R, IOUT, IN)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, ITER1, R, IOUT, IN)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, IPRSMS, R, IOUT, IN)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NONMETH, R, IOUT, IN)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, LINMETH, R, IOUT, IN)
      IF(NONMETH.NE.0)THEN
        IF ( IFDPARAM.EQ.0 ) THEN
        LLOC = 1
        CALL URDCOM(IN, IOUT, LINE)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, THETADUM, IOUT, IN)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I,AKAPPADUM, IOUT, IN)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, GAMMADUM, IOUT, IN)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3,I,AMOMENTDUM, IOUT, IN)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NUMTRACK, R, IOUT, IN)
        THETA = THETADUM
        AKAPPA = AKAPPADUM
        GAMMA = GAMMADUM
        AMOMENTUM = AMOMENTDUM
        IF( NUMTRACK.GT.0 ) THEN
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I,  BTOLDUM, IOUT, IN)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I,BREDUCDUM, IOUT, IN)
        CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I,RESLIMDUM, IOUT, IN)
        BTOL = BTOLDUM
        BREDUC = BREDUCDUM
        RES_LIM = RESLIMDUM
        ENDIF
        ELSE
        CALL SET_RELAX(IFDPARAM)
        END IF
      END IF
C
      HCLOSE = HCLOSEDUM
      HICLOSE = HICLOSEDUM
      IF ( THETA.LT.CLOSEZERO ) THETA = 1.0E-3
C
C3------ECHO INPUT OF NONLINEAR ITERATIN PARAMETERS AND LINEAR SOLVER INDEX
      WRITE(IOUT,9002) HCLOSE,HICLOSE,MXITER,ITER1,IPRSMS,
     * NONMETH,LINMETH
C
 9002 FORMAT(1X,'OUTER ITERATION CONVERGENCE CRITERION (HCLOSE) = ',
     &  E15.6,
     &      /1X,'INNER ITERATION CONVERGENCE CRITERION (HICLOSE) = ',
     &  E15.6,
     &      /1X,'MAXIMUM NUMBER OF OUTER ITERATIONS (MXITER)     = ',I9,
     &      /1X,'MAXIMUM NUMBER OF INNER ITERATIONS (ITER1)      = ',I9,
     &      /1X,'SOLVER PRINTOUT INDEX             (IPRSMS)      = ',I9,
     &      /1X,'NONLINEAR ITERATION METHOD    (NONLINMETH)      = ',I9,
     &      /1X,'LINEAR SOLUTION METHOD           (LINMETH)      = ',I9)
C
      IF(NONMETH.NE.0)THEN
        WRITE(IOUT,9003)THETA,AKAPPA,GAMMA,AMOMENTUM,NUMTRACK
        IF(NUMTRACK.NE.0) WRITE(IOUT,9004) BTOL,BREDUC,RES_LIM
      ENDIF
9003  FORMAT(1X,'D-B-D WEIGHT REDUCTION FACTOR      (THETA)      = ',
     &  E15.6,
     &      /1X,'D-B-D WEIGHT INCREASE INCREMENT    (KAPPA)      = ',
     &  E15.6,
     &      /1X,'D-B-D PREVIOUS HISTORY FACTOR      (GAMMA)      = ',
     &  E15.6,
     &      /1X,'MOMENTUM TERM                  (AMOMENTUM)      = ',
     &  E15.6,
     &      /1X,'MAXIMUM NUMBER OF BACKTRACKS    (NUMTRACK)      = ',I9)
9004  FORMAT(1X,'BACKTRACKING TOLERANCE FACTOR       (BTOL)      = ',
     &  E15.6,
     &      /1X,'BACKTRACKING REDUCTION FACTOR     (BREDUC)      = ',
     &  E15.6,
     &      /1X,'BACKTRACKING REIDUAL LIMIT       (RES_LIM)      = ',
     &  E15.6)
      IF(MXITER.LE.0) THEN
        WRITE(*,5)
        CALL USTOP(' ')
      ELSEIF(ITER1.LE.0) THEN
        WRITE(*,7)
        CALL USTOP(' ')
      ENDIF
    5 FORMAT(/1X,'ERROR: OUTER ITERATION NUMBER MUST BE > 0.')
    7 FORMAT(/1X,'ERROR: INNER ITERATION NUMBER MUST BE > 0.')
C
      ISYMFLG = 1
      IF ( NONMETH.GT.0 )THEN
        WRITE(IOUT,*) '***NEWTON LINEARIZATION WILL BE USED***'
        WRITE(IOUT,*)
        ISYMFLG = 0
      ELSEIF ( NONMETH.EQ.0 )THEN
        WRITE(IOUT,*) '***PICARD LINEARIZATION WILL BE USED***'
        WRITE(IOUT,*)
      ELSEIF ( NONMETH.LT.0 )THEN
        WRITE(IOUT,*) '***PICARD LINEARIZATION WILL BE USED WITH RELAXAT
     *ION***'
        WRITE(IOUT,*)
      ELSE
        WRITE(IOUT,*) '***INCORRECT VALUE FOR VARIABLE NONMETH WAS ',
     +                'SPECIFIED. CHECK INPUT.***'
        WRITE(IOUT,*)
        CALL USTOP('  ')
      END IF
C4------CALL SECONDARY SUBROUTINE TO INITIALIZE AND READ LINEAR SOLVER PARAMETERS
      IF ( LINMETH==1 )THEN
C4A-------FOR XMD SOLVER
        WRITE(IOUT,*) '***XMD LINEAR SOLVER WILL BE USED***'
        CALL XMD7U1AR(IN,IFDPARAM)
        WRITE(IOUT,*)
        ISYMFLG = 0
        IF ( IACL.EQ.0 ) ISYMFLG = 1
      ELSEIF ( LINMETH==2 )THEN
C4B-------FOR PCGU SOLVER
        WRITE(IOUT,*) '***PCGU LINEAR SOLVER WILL BE USED***'
        CALL PCGU7U1AR(IN, MXITER, HICLOSE, ITER1, IFDPARAM, IPCGUM)
        WRITE(IOUT,*)
        ISYMFLG = 0
        IF ( IPCGUM.EQ.1 ) ISYMFLG = 1
      ELSE
C4C-----INCORRECT LINEAR SOLVER FLAG
        WRITE(IOUT,*) '***INCORRECT VALUE FOR LINEAR SOLUTION METHOD ',
     +                'SPECIFIED. CHECK INPUT.***'
        WRITE(IOUT,*)
        CALL USTOP('  ')
      END IF
C
C---------------------------------------------------------------------------------
C5-----ALLOCATE SPACE FOR NONLINEAR ARRAYS AND INITIALIZE
      ALLOCATE (HTEMP(NEQ))
      ALLOCATE (HNCG(MXITER),LRCH(3,MXITER))
      IF(NONMETH.GT.0)THEN
        ALLOCATE (AMATFL(NJA))
      ELSE
        AMATFL => AMAT
      ENDIF
      IF(IABS(NONMETH).EQ.1)THEN
        ALLOCATE (WSAVE(NEQ),HCHOLD(NEQ),DEOLD(NEQ))
        WSAVE = 0.
        HCHOLD = 0.
        DEOLD = 0.
      ENDIF
      HNCG = 0.0D0
      LRCH = 0
C6------RETURN
      RETURN
      END
C
      SUBROUTINE SET_RELAX(IFDPARAM)
        USE SMSMODULE, ONLY: AKAPPA,GAMMA,AMOMENTUM,
     &                       BREDUC,BTOL,NUMTRACK,THETA,
     &                       RES_LIM
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IFDPARAM
C SIMPLE OPTION
      SELECT CASE ( IFDPARAM )
      CASE ( 1 )
        THETA = 1.0
        AKAPPA = 0.0
        GAMMA = 0.0
        AMOMENTUM = 0.0
        NUMTRACK = 0
        BTOL = 0.0
        BREDUC = 0.0
        RES_LIM = 0.0
C MODERATE
       CASE ( 2 )
        THETA = 0.9
        AKAPPA = 0.0001
        GAMMA = 0.0
        AMOMENTUM = 0.0
        NUMTRACK = 0
        BTOL = 0.0
        BREDUC = 0.0
        RES_LIM = 0.0
C COMPLEX
       CASE ( 3 )
        THETA = 0.8
        AKAPPA = 0.0001
        GAMMA = 0.0
        AMOMENTUM = 0.0
        NUMTRACK = 20
        BTOL = 1.05
        BREDUC = 0.1
        RES_LIM = 0.002
      END SELECT
      RETURN
      END SUBROUTINE SET_RELAX

C-----------------------------------------------------------------------------
C
      SUBROUTINE GLO2SMS1AP(KITER,KSTP,KPER)
C******************************************************************
C PERFORM RESIDUAL REDUCTION AND NEWTON LINEARIZATION AND
C PREPARE FOR SPARSE SOLVER, AND CHECK CONVERGENCE OF NONLINEARITIES
C******************************************************************
      USE SimModule
      USE SOLUTION, ONLY: neq,nja,ia,ja,amat,rhs,x,IBOUND=>active
      USE SMSMODULE
      USE XMDMODULE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KITER
      INTEGER, INTENT(IN) :: KSTP
      INTEGER, INTENT(IN) :: KPER
      !LOCAL VARIABLES
      INTEGER :: J, N
      INTEGER :: IN_ITER
      INTEGER :: NB
      DOUBLEPRECISION :: ABIGCH,HDIF,AHDIF,ADIAG,BIG
C---------------------------------------------------------------------
C
C--------------------------------------------------------------------
C2------PERFORM RESIDUAL REDUCTION CYCLES IF REQUIRED
      IF(NUMTRACK.GT.0)THEN
        IF(KITER.EQ.1.AND.IBFLAG.EQ.0)THEN
C2A-------WRITE HEADER FOR SOLVER OUTPUT SUMMARY WITH BACKTRACKING
          WRITE(IOUT,11)
11        FORMAT(/' Outer-Iteration  Inner-Iteration  Backtracking  ',
     1    'Number of        Incoming       Outgoing  Maximum Head',
     1    'Change      Maximum Head Change'/
     1    '     Number           Count           Flag       Backtracks',
     1    7X,'Residual       Residual           Value              ',
     1    'Location')
        ENDIF
C
C2B-------CALL SUBROUTINE TO DETERMINE IF BACKTRACKING IS NEEDED
        CALL SGLO2SMS1RR(KITER,KSTP,KPER)
      ELSE
C2D-------WRITE HEADER FOR SOLVER OUTPUT SUMMARY WITHOUT BACKTRACKING
        IF(KITER.EQ.1)THEN
          WRITE(IOUT,12)
12        FORMAT(/' Outer-Iteration  Inner-Iteration    Maximum Head ',
     1    'Change  Maximum Head Change'/
     1    '     Number           Count               Value',
     1    14X,'Location')
        ENDIF
      ENDIF
C-----------------------------------------------------------------------------
C2-------CALCULATE AND FILL DERIVATIVE TERMS IN JACOBIAN
C2-------FOR NEWTON METHOD FOR LAYERS WITH LAYCON = 4
C2-------THE DERIVATIVE TERMS INCLUDE DKR FOR UNCONFINED AND DWADI FOR VERTICAL FLOW CORRECTION
      IF(NONMETH.GT.0) THEN
C----------------------------------------------------------------------------
C2D-------TRANSFER FLOW TERMS INTO AMATFL FOR USE IN MASS BALANCE CALCULATION
        DO J=1,NJA
          AMATFL(J) = AMAT(J)
        ENDDO
      ENDIF
C
C3-----TAKE CARE OF LOOSE ENDS FOR ALL NODES BEFORE CALL TO SOLVER
      BIG = 1.0E20
      DO N=1,NEQ
C3a-------STORE X IN TEMPORARY LOCATION
        HTEMP(N) = X(N)
C3b-------SET DIRICHLET BOUNDARY AND NO-FLOW CONDITION
        IF(IBOUND(N).LE.0)THEN
          AMAT(IA(N)) = 1.0*BIG
          RHS(N) = X(N)*BIG
        ELSE
C3c---------TAKE CARE OF ZERO ROW DIAGONAL
          ADIAG = ABS(AMAT(IA(N)))
          IF(ADIAG.LT.1.0E-15)THEN
            AMAT(IA(N)) = 1.0E06
            RHS(N) = RHS(N) + X(N)*1.0E06
          ENDIF
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
C4------call linear solver
      CALL SOLVERS(KITER,KSTP,KPER,IN_ITER)
C------------------------------------------------------------
C5------CHECK OUTER ITERATION CONVERGENCE
      NB=1
      ICNVG=0
      BIGCH=0.0
      ABIGCH=0.0
      DO N=1,NEQ
        IF(IBOUND(N).EQ.0) CYCLE
        HDIF=X(N)-HTEMP(N)
        AHDIF=ABS(HDIF)
        IF(AHDIF.GE.ABIGCH)THEN
          BIGCH= HDIF
          ABIGCH= AHDIF
          NB = N
        ENDIF
      ENDDO
C
      IF(ABIGCH.LE.HCLOSE) ICNVG=1
C
C5a------STORE MAXIMUM CHANGE VALUE AND LOCATION
      HNCG(KITER) = BIGCH
C
      LRCH(1,KITER) = NB
      IF(NUMTRACK.GT.0)THEN
        write(iout,22)kiter,IN_ITER,bigch,nb
22      format(i9,I17,69x,g15.6,9x,i9,11x,'GWF-node number')
      ELSE
        write(iout,23)kiter,IN_ITER,bigch,nb
23      format(I9,I17,10X,G16.5,9X,I9,11x,'GWF-node number')
      ENDIF
C
C-----NOT CONVERGED, IF EITHER IS NOT CONVERGED
      IF(ICNVG.EQ.0) ICNVG = 0
204   CONTINUE
!C
!C7------USE CONVERGE OPTION TO FORCE CONVERGENCE
!      IF(IFRCNVG.EQ.1.AND.KITER.EQ.MXITER) THEN
!        ICNVG=1
!      ENDIF
C
C-----------------------------------------------------------
C8-------PERFORM UNDERRELAXATION WITH DELTA-BAR-DELTA
      IF(NONMETH.NE.0.AND.ICNVG.EQ.0) CALL GLO2SMS1UR(KITER)
C
C9------WRITE ITERATION SUMMARY FOR CONVERGED SOLUTION
      IF(ICNVG.EQ.0 .AND. KITER.NE.MXITER) GOTO 600
      IF(KSTP.EQ.1) WRITE(IOUT,1000)
 1000 FORMAT(/1X)
      WRITE(IOUT,1010) KITER,KSTP,KPER
 1010 FORMAT(1X,I5,' CALLS TO SPARSE MATRIX SOLVER PACKAGE ',
     & ' IN FLOW TIME STEP',I4,' STRESS PERIOD',I4)
C
C9A------FOR SOLUTION NODES
       CALL SSMS2BCFU1P(HNCG,LRCH,KITER,MXITER)
C
  600 CONTINUE
C10-----RETURN
      RETURN
      END
C-----------------------------------------------------------------------------
       SUBROUTINE SOLVERS(KITER,KSTP,KPER,IN_ITER)
C******************************************************************
C PREPARE FOR SPARSE SOLVER, AND CHECK CONVERGENCE OF NONLINEARITIES
C******************************************************************
      USE SimModule
      USE SOLUTION, ONLY: neq,nja,ia,ja,amat,rhs,x,IBOUND=>active
      USE SMSMODULE
      USE XMDMODULE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KITER
      INTEGER, INTENT(IN) :: KSTP
      INTEGER, INTENT(IN) :: KPER
      INTEGER, INTENT(INOUT) :: IN_ITER
      !LOCAL VARIABLES
      INTEGER :: I, N
      INTEGER :: I1, I2
      INTEGER :: ITESTMAT
      INTEGER :: ITER
      DOUBLEPRECISION :: STOL
C---------------------------------------------------------------------
C1-------PRINT MATRIX AND RHS FOR CHECKING
      itestmat = 0
      if(itestmat.eq.1)then
        open(99,file='mat_USGs.TXT')
        WRITE(99,*)'NODE, RHS, AMAT FOLLOW'
        DO N=1,NEQ
          I1 = IA(N)
          I2 = IA(N+1)-1
          WRITE(99,66)N,RHS(N),(AMAT(I),I=I1,I2)
        ENDDO
66      FORMAT(I9,1X,G15.6,2X,10G15.6)
        CLOSE (99)
        stop
      endif
C
C-----------------------------------------------------------
C2-------CALL LINEAR SOLVER
C-----------------------------------------------------------
        IF(LINMETH.EQ.1)THEN
C
C2A---------CALL XMD SOLVER
C2A1-------- ILU FACTORIZATION
          IF (IDROPTOL.EQ.0) THEN
C2A2--------numerical factorization only for level based scheme
            call xmdnfctr(amat, rhs, ia, ja, nja, NEQ, ierr)
          ELSE
C2A3--------level/drop tolerance preconditioning
            call xmdprecd(amat, rhs, epsrn, ia, ja, nja, NEQ, level,
     &            ierr)
          ENDIF
C
C2A4---------solve matrix
          iter = iter1
          call xmdsolv(amat, rhs, x, stol, rrctol, ia, ja, nja, NEQ,
     &                 north, iter, iacl, ierr)
          IN_ITER = ITER
C--------------------------------------------------------------
        ELSEIF(LINMETH.EQ.2)THEN
C2B---------CALL PCGC SOLVER
C
          CALL PCGU7U1AP(ICNVG,KSTP,KITER,IN_ITER)
C          
        ENDIF
C
C ----------------------------------------------------------------------
C3------PRINT SOLUTION FOR CHECKING
      if(itestmat.eq.1)then
        WRITE(IOUT,*)'MATRIX SOLUTION FOLLOWS'
        WRITE(IOUT,67)(n,X(N),N=1,NEQ)
67    FORMAT(10(I8,G15.6))
        stop
      endif
C
C4-------Return
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE SSMS2BCFU1P(CNCG,LRCH,ITP,MXITER)
C******************************************************************
C PRINT MAXIMUM HEAD CHANGES FOR ALL ITERATIONS FOR POROUS MEDIUM NODES
C******************************************************************
C
      USE SimModule
      IMPLICIT  NONE
      DOUBLEPRECISION, DIMENSION(MXITER), INTENT(IN) :: CNCG
      INTEGER, DIMENSION(3,MXITER), INTENT(IN) :: LRCH
      INTEGER, INTENT(IN) :: ITP
      INTEGER, INTENT(IN) :: MXITER
      !LOCAL VARIABLES
      INTEGER :: J
C---------------------------------------------------------------------------
C1-----FOR POROUS MATRIX NODES
      WRITE(IOUT,2)ITP
2     FORMAT(/1X,'TOTAL OF ',I7,'OUTER ITERATIONS')
        WRITE(IOUT,15)
   15   FORMAT(1X,' MAXIMUM CHANGE FOR EACH ITERATION:'
     &      /1X, 5('  MAX. CHANGE        NODE')/1X,132('-'))
        WRITE(IOUT,20) (CNCG(J),LRCH(1,J),J=1,ITP)
   20   FORMAT((1X,5(G13.5,', ',I10)))
C---------------------------------------------------------------------------
C2------RETURN
      RETURN
      END
C-------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      SUBROUTINE GLO2SMS1UR(KITER)
C     ******************************************************************
C     UNDER RELAX AS PER DELTA-BAR-DELTA OR COOLEY FORMULA
C     ******************************************************************
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
      USE SOLUTION, ONLY: neq,nja,ia,ja,amat,rhs,x,IBOUND=>active
      USE SMSMODULE
      IMPLICIT NONE
!     -----------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      INTEGER :: KITER
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      DOUBLEPRECISION :: WW, HSAVE, DELH, RELAX, ES,
     &                   AES, AMOM, CLOSEBOT, FELEV
      INTEGER :: K,N,NNDLAY,NSTRT,I,kk
!     -----------------------------------------------------------------
      closebot = 0.9
      IF(ABS(NONMETH).EQ.1) THEN
C1-------OPTION FOR USING DELTA-BAR-DELTA SCHEME TO UNDER-RELAX SOLUTION FOR ALL EQUATIONS
        DO N=1,NEQ
C
C2---------COMPUTE NEWTON STEP-SIZE (DELTA H) AND INITIALIZE D-B-D PARAMETERS
          DELH = X(N) - HTEMP(N)
          IF ( kiter.EQ.1 )THEN
            Wsave(N) = 1.0D0
            Hchold(N) = 1.0E-20
            DEold(N) = 0.0D0
          END IF
C
C3---------COMPUTE NEW RELAXATION TERM AS PER DELTA-BAR-DELTA
          ww = Wsave(N)
          IF ( DEold(N)*DELH.LT.0.0D0 ) THEN
C3A-----------FOR FLIP-FLOP CONDITION, DECREASE FACTOR
c             ww = -0.5* DEold(N) / DELH
c             if(ww.gt.1.0) ww = 1.0
c             if(ww.lt.1.0e-08) ww = 1.0e-08
             ww = Theta*Wsave(N)
          ELSE
C3B---------WHEN CHANGE IS OF SAME SIGN, INCREASE FACTOR
c             ww = 1.0
            ww = Wsave(N) + akappa
          END IF
          IF ( ww.GT.1.0d0 ) ww = 1.0d0
          Wsave(N) = ww
C4----------COMPUTE EXPONTENTIAL AVERAGE OF PAST CHANGES IN Hchold
          If(kiter.eq.1)then !this method does it right after Newton - need to do it after underrelaxation and backtracking.
            Hchold(N) = DELH
          Else
            Hchold(N) = (1-gamma) * DELH + gamma * Hchold(N)
          Endif
C5--------STORE SLOPE (CHANGE) TERM FOR NEXT ITERATION
         DEold(N) = DELH
C
C6----------COMPUTE ACCEPTED STEP-SIZE AND NEW HEAD
          amom = 0.0
          if(kiter.gt.4) amom = amomentum
          DELH = DELH * ww + amom * Hchold(N)
          X(N) = HTEMP(N) + DELH
!C7----------ACCOUNT FOR ICONCV=0 CONDITION FOR LAYCON=4 CASE
!          IF ( ICONCV.EQ.0.AND.IUNSAT.EQ.0) THEN
!            IF(N.LE.NODES)THEN  !-------FOR POROUS MEDIUM NODES
!            DO K=1,NLAY
!              NNDLAY = NODLAY(K)
!              NSTRT = NODLAY(K-1)+1
!              IF(N.GE.NSTRT.AND.N.LE.NNDLAY)THEN
!               KK = K
!               GO TO 11
!              ENDIF
!            ENDDO
!11          CONTINUE
!              IF ( LAYCON(KK).EQ.4 ) THEN
!                IF ( X(N).LT.Bot(N) ) THEN
!                  hsave = X(N)
!                  X(N) = HTEMP(N)*(1.0-closebot) + Bot(N)*closebot
!                  DELH = X(N) - hsave
!                END IF
!              ENDIF
!            ENDIF
!          ENDIF
          ENDDO
C---------------------------------------------------------------------------
      ELSEIF(ABS(NONMETH).EQ.2) THEN
C8-------DO COOLEY UNDERRELAXATION
        IF(KITER.EQ.1)THEN
          RELAX = 1.0
          RELAXOLD = 1.0
          BIGCHOLD = BIGCH
        ELSE
C9---------COMPUTE RELAXATION FACTOR
          ES = BIGCH / (BIGCHOLD*RELAXOLD)
          AES = ABS(ES)
          IF(ES.LT.-1.0E0)THEN
            RELAX = 0.5/AES
          ELSE
            RELAX = (3.0+ES) / (3.0+AES)
          ENDIF
        ENDIF
        RELAXOLD = RELAX
C10---------MODIFY COOLEY TO USE EXPONENTIAL AVERAGE OF PAST CHANGES AS PER LINE BELOW
        BIGCHOLD = (1-gamma)*BIGCH  + gamma*BIGCHOLD  !this method does it right after Newton - need to do it after underrelaxation and backtracking.
C        if(numtrack.eq.0) BIGCHOLD = (1-gamma)*BIGCH  + gamma*BIGCHOLD
        IF(RELAX.LT.1.0)THEN
C11---------COMPUTE NEW HEAD AFTER UNDER-RELAXATION
          DO N = 1, NEQ
            DELH = X(N) - HTEMP(N)
            X(N) = HTEMP(N) + RELAX * DELH
          ENDDO
        ENDIF
!C12-------ACCOUNT FOR ICONCV=0 CONDITION APPROPRIATELY
!         IF ( ICONCV.EQ.0.AND.IUNSAT.EQ.0) THEN
!          DO K=1,NLAY
!            NNDLAY = NODLAY(K)
!            NSTRT = NODLAY(K-1)+1
!            IF ( LAYCON(K).EQ.4 ) THEN
!              DO N=NSTRT,NNDLAY
!                IF ( X(N).LT.Bot(N) ) THEN
!                  X(N) = HTEMP(N)*(1.0-closebot) + Bot(N)*closebot
!                END IF
!              ENDDO
!            ENDIF
!          ENDDO
!101       CONTINUE
!C
!        ENDIF
C---------------------------------------------------------------------------
      ENDIF
C13-----RETURN
      RETURN
      END SUBROUTINE GLO2SMS1UR
C
C-----------------------------------------------------------------------
      SUBROUTINE SGLO2SMS1RR(KITER,KSTP,KPER)
C     ******************************************************************
C     COMPUTE RESIDUAL AND EVALUATE IF REDUCTION IS REQUIRED
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE SimModule
      USE SOLUTION, ONLY: neq,nja,ia,ja,amat,rhs,x,IBOUND=>active
      USE SMSMODULE, ONLY: IBFLAG,BTOL,BREDUC,NUMTRACK,HCLOSE,HTEMP,
     &                     GAMMA,HCHOLD,BIGCHOLD,BIGCH,RES_LIM,NONMETH,
     &                     RES_PREV,RES_NEW,RES_IN,IBCOUNT 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KITER
      INTEGER, INTENT(IN) :: KSTP
      INTEGER, INTENT(IN) :: KPER
C     LOCAL VARIABLES
      INTEGER :: N
      DOUBLEPRECISION :: CHMAX,DELH,ABSDELH,HDIF,AHDIF,ABIGCH
C     ------------------------------------------------------------------
C
C1------AT FIRST ITERATION
      IF(KITER.EQ.1)THEN
c1A------RESET RESIDUAL REDUCTION COUNT AND FLAG
        IBCOUNT = 0
        IBFLAG = 0
C
C1B------INITIALIZE PREVIOUS RESIDUAL AND RETURN
        CALL RES_FUNC(RES_PREV)
        WRITE(IOUT,66) KITER,IBFLAG,IBCOUNT,RES_PREV,RES_PREV
        RETURN
      ENDIF
66    FORMAT(I9,17X,I13,I14,7X,G15.6,2X,G15.6)
C--------------------------------------------------------------
C2------COMPUTE CURRENT RESIDUAL
      CALL RES_FUNC(RES_NEW)
      IF(IBCOUNT.EQ.0) RES_IN = RES_NEW
C
C3-----SET COUNT AND FLAG IF DESIRED RESIDUAL REDUCTION DID NOT OCCUR
      IF(RES_NEW.GT.RES_PREV * BTOL)THEN
C
C3A-------BUT NO BACKTRACKING IF MAXIMUM TRACKS ARE EXCEEDED SO RETURN
        IF(IBCOUNT.GE.NUMTRACK)THEN
          IBFLAG = 2
          WRITE(IOUT,66) KITER,IBFLAG,IBCOUNT,RES_IN,RES_PREV
          IBFLAG = 0
          IBCOUNT = 0
          RES_PREV = RES_NEW
          RETURN
        ENDIF
C
C3B-------BUT NO BACKTRACKING IF RESIDUAL IS SMALLER THAN LIMIT SO RETURN
        IF(RES_NEW.LT.RES_LIM)THEN
          IBFLAG = 3
          WRITE(IOUT,66) KITER,IBFLAG,IBCOUNT,RES_IN,RES_NEW
          IBFLAG = 0
          IBCOUNT = 0
          RES_PREV = RES_NEW
          RETURN
        ENDIF
C
C3C-------ALSO NO BACKTRACKING IF MAXIMUM CHANGE IS LESS THAN CLOSURE SO RETURN
        CHMAX = 0.0
        DO N=1,NEQ
          DELH = BREDUC*(X(N) - HTEMP(N))
          ABSDELH = ABS(DELH)
          IF(ABSDELH.GT.CHMAX) CHMAX = ABSDELH
        ENDDO
        IF(CHMAX.LT.HCLOSE)THEN
          IBFLAG = 4
          WRITE(IOUT,66) KITER,IBFLAG,IBCOUNT,RES_IN,RES_NEW
          IBFLAG = 0
          IBCOUNT = 0
          RES_PREV = RES_NEW
          RETURN
        ENDIF
C
C4-------PERFORM BACKTRACKING IF FREE OF CONSTRAINTS AND SET COUNTER AND FLAG
        DO N=1,NEQ
          DELH = BREDUC*(X(N) - HTEMP(N))
          X(N) = HTEMP(N) + DELH
        ENDDO
        IBCOUNT = IBCOUNT + 1
        IBFLAG = 1
C
C5------RESET COUNT AND FLAG IF DESIRED RESIDUAL REDUCTION DID OCCUR
      ELSE
        WRITE(IOUT,66) KITER,IBFLAG,IBCOUNT,RES_IN,RES_NEW
        IBFLAG = 0
        IBCOUNT = 0
        RES_PREV = RES_NEW
      ENDIF
C
C7------RETURN.
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE RES_FUNC(RES)
C     ******************************************************************
C     COMPUTE RESIDUAL - USING MEAN SQUARE RESIDUAL
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE SOLUTION, ONLY: neq,nja,ia,ja,amat,rhs,x,IBOUND=>active
      USE SMSMODULE, ONLY: IBFLAG,BTOL,BREDUC,NUMTRACK,HCLOSE,HTEMP
      IMPLICIT NONE
      DOUBLEPRECISION, INTENT(INOUT) :: RES
C     LOCAL VARIABLES      
      INTEGER :: J, N
      INTEGER :: JJ
      DOUBLEPRECISION :: ROWSUM,RESIDUAL
C     ------------------------------------------------------------------
C
C1------COMPUTE Q FOR ALL NODES
      RESIDUAL = 0.0
      DO N=1,NEQ
        IF(IBOUND(N).GT.0)THEN
          ROWSUM = 0.0
          DO J = IA(N),IA(N+1)-1
            JJ = JA(J)
            ROWSUM = ROWSUM + AMAT(J) * X(JJ)
          ENDDO
C2----------COMPUTE MEAN SQUARE RESIDUAL FROM Q OF EACH NODE
          RESIDUAL = RESIDUAL +  (ROWSUM - RHS(N))**2
        ENDIF
      ENDDO
      RES = RESIDUAL
C3------RETURN
      RETURN
      END SUBROUTINE RES_FUNC
