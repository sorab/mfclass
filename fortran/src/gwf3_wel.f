      SUBROUTINE GWF2WEL7U1AR(IN)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR WELL PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,NODES,IUNSTR,NEQS
      USE GWFWELMODULE, ONLY:NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL,NPWEL,
     1                       IWELPB,NNPWEL,WELAUX,WELL,IWELQV,NNPWCLN,
     2                       IAFR
C
      CHARACTER*200 LINE
C     ------------------------------------------------------------------
      ALLOCATE(NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL,IAFR)
      ALLOCATE(NPWEL,IWELPB,NNPWEL,IWELQV,NNPWCLN)
C
C1------IDENTIFY PACKAGE AND INITIALIZE NWELLS.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'WEL -- WELL PACKAGE, VERSION 7, 5/2/2005',
     1' INPUT READ FROM UNIT ',I4)
      NWELLS=0
      NNPWEL=0
      NNPWCLN=0
      IWELQV=0
      IAFR=0
C
C2------READ MAXIMUM NUMBER OF WELLS AND UNIT OR FLAG FOR
C2------CELL-BY-CELL FLOW TERMS.
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPWEL,MXPW)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(3I10)') MXACTW,IWELCB
         LLOC=21
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTW,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IWELCB,R,IOUT,IN)
      END IF
      WRITE(IOUT,3) MXACTW
    3 FORMAT(1X,'MAXIMUM OF ',I6,' ACTIVE WELLS AT ONE TIME')
      IF(IWELCB.LT.0) WRITE(IOUT,7)
    7 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF(IWELCB.GT.0) WRITE(IOUT,8) IWELCB
    8 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
      WRITE(IOUT,9) MXACTW,IWELCB
    9 FORMAT(1X,'MAXIMUM NUMBER OF ACTIVE WELLS (MXACTW) =',I7
     *  /1X,'C-B-C FLUX FLAG OR UNIT NUMBER (IWELCB) =',I3)
C
C3------READ AUXILIARY VARIABLES AND PRINT FLAG.
      ALLOCATE(WELAUX(20))
      NAUX=0
      IPRWEL=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.20) THEN
            NAUX=NAUX+1
            WELAUX(NAUX)=LINE(ISTART:ISTOP)
            WRITE(IOUT,12) WELAUX(NAUX)
   12       FORMAT(1X,'AUXILIARY WELL VARIABLE: ',A)
         END IF
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'AUTOFLOWREDUCE') THEN
         WRITE(IOUT,16)
   16    FORMAT(1X,'WELL FLUX WILL BE REDUCED WHEN SATURATED ',
     1       'THICKNESS IS LESS THAN 1 PERCENT OF CELL THICKNESS')
         IWELQV = 1
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'IUNITAFR') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IAFR,R,IOUT,IN)
         WRITE(IOUT,25) IAFR
   25    FORMAT(1X,'WELL REDUCTION INFO WILL BE WRITTEN TO UNIT: ',
     1       I5)
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,'LISTS OF WELL CELLS WILL NOT BE PRINTED')
         IPRWEL = 0
         GO TO 10
      END IF
C3A-----THERE ARE FOUR INPUT VALUES PLUS ONE LOCATION FOR
C3A-----CELL-BY-CELL FLOW.
      NWELVL=5+NAUX
C
C4------ALLOCATE SPACE FOR THE WELL DATA.
      IWELPB=MXACTW+1
      MXWELL=MXACTW+MXPW
      ALLOCATE (WELL(NWELVL,MXWELL))
C
C5------READ NAMED PARAMETERS.
      WRITE(IOUT,18) NPWEL
   18 FORMAT(1X,//1X,I5,' Well parameters')
      IF(NPWEL.GT.0) THEN
        LSTSUM=IWELPB
        DO 120 K=1,NPWEL
          LSTBEG=LSTSUM
          CALL UPARLSTRP(LSTSUM,MXWELL,IN,IOUT,IP,'WEL','Q',1,
     &                   NUMINST)
          NLST=LSTSUM-LSTBEG
          IF(NUMINST.EQ.0) THEN
C5A-----READ PARAMETER WITHOUT INSTANCES.
            IF(IUNSTR.EQ.0)THEN
              CALL ULSTRD(NLST,WELL,LSTBEG,NWELVL,MXWELL,1,IN,
     &        IOUT,'WELL NO.  LAYER   ROW   COL   STRESS FACTOR',
     &        WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
            ELSE
             CALL ULSTRDU(NLST,WELL,LSTBEG,NWELVL,MXWELL,1,IN,
     &        IOUT,'WELL NO.      NODE       STRESS FACTOR',
     &        WELAUX,20,NAUX,IFREFM,NEQS,4,4,IPRWEL)
            ENDIF
          ELSE
C5B-----READ INSTANCES.
            NINLST=NLST/NUMINST
            DO 110 I=1,NUMINST
            CALL UINSRP(I,IN,IOUT,IP,IPRWEL)
            IF(IUNSTR.EQ.0)THEN
              CALL ULSTRD(NINLST,WELL,LSTBEG,NWELVL,MXWELL,1,IN,
     &        IOUT,'WELL NO.  LAYER   ROW   COL   STRESS FACTOR',
     &        WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
            ELSE
             CALL ULSTRDU(NINLST,WELL,LSTBEG,NWELVL,MXWELL,1,IN,
     &        IOUT,'WELL NO.      NODE       STRESS FACTOR',
     &        WELAUX,20,NAUX,IFREFM,NEQS,4,4,IPRWEL)
            ENDIF
            LSTBEG=LSTBEG+NINLST
  110       CONTINUE
          END IF
  120   CONTINUE
      END IF
C
C6------RETURN
      RETURN
      END
      SUBROUTINE GWF2WEL7U1RP(IN)
C     ******************************************************************
C     READ WELL DATA FOR A STRESS PERIOD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,NODES,IUNSTR,NEQS,INCLN
      USE GWFWELMODULE, ONLY:NWELLS,MXWELL,NWELVL,IPRWEL,NPWEL,
     1                       IWELPB,NNPWEL,WELAUX,WELL,NNPWCLN
C
      CHARACTER*6 CWELL
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'WEL -- WELL PACKAGE, VERSION 7, 5/2/2005',
     1' INPUT READ FROM UNIT ',I4)
C
C1----READ NUMBER OF WELLS (OR FLAG SAYING REUSE WELL DATA).
C1----AND NUMBER OF PARAMETERS
      NNPWCLN=0
      ITMPCLN=0
      IF(INCLN.GT.0)THEN
        IF(IFREFM.EQ.0) THEN
          READ(IN,'(3I10)') ITMP,NP,ITMPCLN
        ELSE
          READ(IN,*) ITMP,NP,ITMPCLN
        END IF
      ELSE
        IF(IFREFM.EQ.0) THEN
          READ(IN,'(2I10)') ITMP,NP
        ELSE
          READ(IN,*) ITMP,NP
        END IF
      ENDIF
C
C------Calculate some constants.
      NAUX=NWELVL-5
      IOUTU = IOUT
      IF (IPRWEL.EQ.0) IOUTU=-IOUTU
C
C1A-----IF ITMP LESS THAN ZERO REUSE NON-PARAMETER DATA. PRINT MESSAGE.
C1A-----IF ITMP=>0, SET NUMBER OF NON-PARAMETER WELLS EQUAL TO ITMP.
      IF(ITMP.LT.0) THEN
         WRITE(IOUT,6)
    6    FORMAT(1X,/
     1    1X,'REUSING NON-PARAMETER WELLS FROM LAST STRESS PERIOD')
      ELSE
         NNPWEL=ITMP
      END IF
C
      IF(INCLN.GT.0)THEN
C       TREAT STRUCTURED CASE WITH CLN
        IF(ITMPCLN.LT.0) THEN
           WRITE(IOUT,7)
    7      FORMAT(1X,/
     1    1X,'REUSING NON-PARAMETER CLN WELLS FROM LAST STRESS PERIOD')
        ELSE
         NNPWCLN=ITMPCLN
        END IF
      ENDIF


C
C1B-----IF THERE ARE NEW NON-PARAMETER WELLS, READ THEM.
      MXACTW=IWELPB-1
      IF(ITMP.GT.0.OR.ITMPCLN.GT.0) THEN
        IF(NNPWEL.GT.MXACTW) THEN
          WRITE(IOUT,99) NNPWEL,MXACTW
   99     FORMAT(1X,/1X,'THE NUMBER OF ACTIVE WELLS (',I6,
     1                   ') IS GREATER THAN MXACTW(',I6,')')
          CALL USTOP(' ')
        END IF
C
        IF(ITMP.GT.0) THEN
          IF(IUNSTR.EQ.0)THEN
            CALL ULSTRD(NNPWEL,WELL,1,NWELVL,MXWELL,1,IN,IOUT,
     1           'WELL NO.  LAYER   ROW   COL   STRESS RATE',
     2            WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
          ELSE
            CALL ULSTRDU(NNPWEL,WELL,1,NWELVL,MXWELL,1,IN,IOUT,
     &      'WELL NO.      NODE       STRESS FACTOR',
     &       WELAUX,20,NAUX,IFREFM,NEQS,4,4,IPRWEL)
          ENDIF
        ENDIF
C
        IF(ITMPCLN.GT.0) THEN
          CALL ULSTRDU(NNPWCLN,WELL,NNPWEL+1,NWELVL,MXWELL,1,IN,
     &    IOUT, 'WELL NO.  CLN-NODE       STRESS FACTOR',
     &    WELAUX,20,NAUX,IFREFM,NEQS,4,4,IPRWEL)
        ENDIF
      END IF
      NWELLS=NNPWEL+NNPWCLN
C
C1C-----IF THERE ARE ACTIVE WELL PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('Q')
      NREAD=NWELVL-1
      IF(NP.GT.0) THEN
         DO 30 N=1,NP
         CALL UPARLSTSUB(IN,'WEL',IOUTU,'Q',WELL,NWELVL,MXWELL,NREAD,
     1                MXACTW,NWELLS,4,4,
     2            'WELL NO.  LAYER   ROW   COL   STRESS RATE',
     3            WELAUX,20,NAUX)
   30    CONTINUE
      END IF
C
C3------PRINT NUMBER OF WELLS IN CURRENT STRESS PERIOD.
      CWELL=' WELLS'
      IF(NWELLS.EQ.1) CWELL=' WELL '
      WRITE(IOUT,101) NWELLS,CWELL
  101 FORMAT(1X,/1X,I6,A)
C
C-------FOR STRUCTURED GRID, CALCULATE NODE NUMBER AND PLACE IN LAYER LOCATION
      IF(ITMP.GT.0.AND.IUNSTR.EQ.0)THEN
        DO L=1,NWELLS
          IR=WELL(2,L)
          IC=WELL(3,L)
          IL=WELL(1,L)
          N = IC + NCOL*(IR-1) + (IL-1)* NROW*NCOL
          WELL(1,L) = N
        ENDDO
      ENDIF
C
C-------FOR CONDUIT-NODES, CALCULATE GLOBAL NODE NUMBER
      IF(ITMPCLN.GT.0)THEN
        DO L=NNPWEL+1,NNPWEL+NNPWCLN
          WELL(1,L) = WELL(1,L) + NODES
        ENDDO
      ENDIF
C
C6------RETURN
      RETURN
      END
      SUBROUTINE GWF2WEL7U1FM
C     ******************************************************************
C     SUBTRACT Q FROM RHS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,AMAT,IA,TOP,BOT,HNEW,NODES
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFWELMODULE, ONLY:NWELLS,WELL,IWELQV
      DOUBLE PRECISION QTHIK,X,Y,Q,QA,QEPS,DQ,EPS,BOTT,THCK,HD
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF WELLS <= 0 THEN RETURN.
      IF(NWELLS.LE.0) RETURN
C
C2------PROCESS EACH WELL IN THE WELL LIST.
      DO 100 L=1,NWELLS
      N=WELL(1,L)
      Q=WELL(4,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(IBOUND(N).LE.0) GO TO 100
C
C2B-----IF THE CELL IS VARIABLE HEAD THEN SUBTRACT Q FROM
C       THE RHS ACCUMULATOR.
      IF(IWELQV.EQ.1.AND.Q.LT.0)THEN
        IPIV = IA(N)
C-------HONOR SUPPLY/DEMAND CONDITIONS FOR EXTRACTION WELLS (NEWTON METHOD)
        HD = HNEW(N)
        IF(N.GT.NODES)THEN
          ICLN = N-NODES
          !!langevin mf2015 todo: remove thisCALL CLNV(ICLN,THCK)
          QTHIK = THCK * 0.01 ! OVER 1 PERCENT OF CELL THICKNESS
          BOTT = ACLNNDS(ICLN,5)
        ELSE
          QTHIK = (TOP(N) - BOT(N)) * 0.01 ! OVER 1 PERCENT OF CELL THICKNESS
          BOTT = BOT(N)
        ENDIF
        X = (HD - BOTT) /QTHIK
        CALL SMOOTH(X,Y)
        QA = Q * Y
C-------CALCULATE DQ/DH
        EPS = 0.01 * QTHIK
        X = (HD+EPS - BOTT) /QTHIK
        CALL SMOOTH(X,Y)
        QEPS = Q*Y
        DQ = (QEPS - QA) / EPS
        AMAT(IPIV) = AMAT(IPIV) + DQ
        RHS(N) = RHS(N) - QA + DQ*HD
      ELSE
        RHS(N)=RHS(N)-Q
      ENDIF
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
c ---------------------------------------------------
      subroutine smooth(x,y)
C     COMPUTES THE S CURVE FOR SMOOTH DERIVATIVES BETWEEN X=0 AND X=1
      double precision x,y
      if (x.le.0.0) then
        y = 0.0
      elseif(x.lt.1.0)then
        y = - 2.0 * x**3 + 3*x*x
      else
        y = 1.0
      endif
      return
      end
c ---------------------------------------------------
      SUBROUTINE GWF2WEL7U1BD(KSTP,KPER)
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR WELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,
     1                  IUNSTR,TOP,BOT,HNEW,NEQS,INCLN
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNCB 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFWELMODULE,ONLY:NWELLS,IWELCB,WELL,NWELVL,WELAUX,IWELQV,IAFR
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,QTHIK,X,Y,HD,THCK,BOTT
      DATA TEXT(1) /'           WELLS'/
      DATA TEXT(2) /'       CLN WELLS'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      NWELLAFR=0
      IF(IWELCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IWELCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF(IBD.EQ.2) THEN
         NAUX=NWELVL-5
         IF(IAUXSV.EQ.0) NAUX=0
         IF(IUNSTR.EQ.0)THEN 
           CALL UBDSV4(KSTP,KPER,TEXT(1),NAUX,WELAUX,IWELCB,NCOL,NROW,
     1          NLAY,NWELLS,IOUT,DELT,PERTIM,TOTIM,IBOUND)
         ELSE 
           CALL UBDSV4U(KSTP,KPER,TEXT(1),NAUX,WELAUX,IWELCB,NODES,
     1          NWELLS,IOUT,DELT,PERTIM,TOTIM,IBOUND)
         ENDIF
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NEQS
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO WELLS, DO NOT ACCUMULATE FLOW.
      IF(NWELLS.EQ.0) GO TO 200
C
C5------LOOP THROUGH EACH WELL CALCULATING FLOW.
      DO 100 L=1,NWELLS
C
C5A-----GET NODE NUMBER OF CELL CONTAINING WELL.
      N=WELL(1,L)
      Q=ZERO
C
C5B-----IF THE CELL IS NO-FLOW OR CONSTANT_HEAD, IGNORE IT.
      IF(IBOUND(N).LE.0)GO TO 99
C
C5C-----GET FLOW RATE FROM WELL LIST.
      Q=WELL(4,L)
      IF(IWELQV.EQ.1.AND.Q.LT.0)THEN
C-------HONOR SUPPLY/DEMAND CONDITIONS FOR EXTRACTION WELLS
        HD = HNEW(N)
        IF(N.GT.NODES)THEN
          ICLN = N-NODES
          !!langevin mf2015 todo: remove thisCALL CLNV(ICLN,THCK)
          QTHIK = THCK * 0.01 ! OVER 1 PERCENT OF CELL THICKNESS
          BOTT = ACLNNDS(ICLN,5)
        ELSE
          QTHIK = (TOP(N) - BOT(N)) * 0.01 ! OVER 1 PERCENT OF CELL THICKNESS
          BOTT = BOT(N)
        ENDIF
        X = (HD - BOTT) /QTHIK
        CALL SMOOTH(X,Y)
        Q = Q * Y
      ENDIF
      QQ=Q
C
C5D-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT(1),KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0)THEN
          IF(N.GT.NODES)THEN
            WRITE(IOUT,64)N-NODES,Q
   64       FORMAT(1X,'CLN NODE',I6,'   RATE',1PG15.6)
          ELSE
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'WELL ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL '
     1       ,I5,'   RATE ',1PG15.6)
          ENDIF
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'WELL ',I6,'    NODE ',I8,'   RATE ',1PG15.6)
        ENDIF
         IBDLBL=1
      END IF
C
C5E-----ADD FLOW RATE TO BUFFER.
      BUFF(N)=BUFF(N)+Q
C
C5F-----SEE IF FLOW IS POSITIVE OR NEGATIVE.
      IF(Q.GE.ZERO) THEN
C
C5G-----FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
        RATIN=RATIN+QQ
      ELSE
C
C5H-----FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
        RATOUT=RATOUT-QQ
      END IF
C
C5I-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C5I-----COPY FLOW TO WELL LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN 
        IF(IUNSTR.EQ.0)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
          CALL UBDSVB(IWELCB,NCOL,NROW,IC,IR,IL,Q,
     1                  WELL(1,L),NWELVL,NAUX,5,IBOUND,NLAY)
        ELSE
          CALL UBDSVBU(IWELCB,NODES,N,Q,
     1                  WELL(1,L),NWELVL,NAUX,5,IBOUND)
        ENDIF
      ENDIF
      WELL(NWELVL,L)=Q
C
C5J-----WRITE FLOW REDUCTION INFO IF REQUESTED
      IF (IWELQV.GT.0. .AND. IAFR.GT.0) THEN
        IF(Q.GT.WELL(4,L)) THEN
            IF (NWELLAFR.EQ.0) THEN
              WRITE(IAFR, *)
              WRITE(IAFR,300) KPER,KSTP
            END IF
            IF(IUNSTR.EQ.0) THEN
              IL = (N-1) / (NCOL*NROW) + 1
              IJ = N - (IL-1)*NCOL*NROW
              IR = (IJ-1)/NCOL + 1
              IC = IJ - (IR-1)*NCOL
              IF(NWELLAFR.EQ.0) WRITE(IAFR,400)
              WRITE(IAFR,500) L,IL,IR,IC,WELL(4,L),Q,HD,BOTT
            ELSE
              IF(NWELLAFR.EQ.0) WRITE(IAFR,401)
              WRITE(IAFR,501) L,N,WELL(4,L),Q,HD,BOTT
            ENDIF
            NWELLAFR = NWELLAFR + 1
        ENDIF
      END IF
  300 FORMAT(' WELLS WITH REDUCED PUMPING FOR STRESS PERIOD ',I5,
     1      ' TIME STEP ',I5)
  400 FORMAT('WELL.NO   LAY   ROW   COL         APPL.Q          ACT.Q',
     1       '        GW_HEAD       CELL_BOT')
  401 FORMAT('WELL.NO    NODE         APPL.Q          ACT.Q',
     1       '        GW_HEAD       CELL_BOT')
  500 FORMAT(I7,3I6,4(1PG15.6))
  501 FORMAT(I7,1X,I9,4(1PG15.6))
C
C-----END OF WELL LOOP
  100 CONTINUE
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT(1),IWELCB,BUFF(1),NCOL,
     1                   NROW,NLAY,IOUT)
        IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
          IF(ICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT(2),ICLNCB,
     1         BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
        ENDIF
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT(1),IWELCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
C
C7------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVL(3,MSUM)=RIN
      VBVL(4,MSUM)=ROUT
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT(1)
C
C8------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUM=MSUM+1
C
C9------RETURN
      RETURN
      END
      SUBROUTINE GWF2WEL7U1DA
C  Deallocate WEL MEMORY
      USE GWFWELMODULE
C
        DEALLOCATE(NWELLS)
        DEALLOCATE(MXWELL)
        DEALLOCATE(NWELVL)
        DEALLOCATE(IWELCB)
        DEALLOCATE(IWELQV)
        DEALLOCATE(IPRWEL)
        DEALLOCATE(IAFR)
        DEALLOCATE(NPWEL)
        DEALLOCATE(IWELPB)
        DEALLOCATE(NNPWEL)
        DEALLOCATE(WELAUX)
        DEALLOCATE(WELL)
C
      RETURN
      END

