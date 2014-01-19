module tdismodule
  integer :: nper,kper,kstp
  integer :: itrss
  real,    save,    dimension(:),    allocatable  ::perlen
  integer, save,    dimension(:),    allocatable  ::nstp
  real,    save,    dimension(:),    allocatable  ::tsmult
  integer, save,    dimension(:),    allocatable  ::issflg

  contains

  subroutine tdis_ar(fname,iout)
! ******************************************************************************
! Initialize variables in the tdismodule using the values in the tdis input
! file called fname.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  implicit none
  character(len=*),intent(in) :: fname
  integer,intent(in) :: iout
  character(len=200) :: line
  integer :: n,i,in,lloc,istart,istop,iss,itr
  real :: r,zero
! ------------------------------------------------------------------------------
!
! -- Get a unit number for tdis and open the file
  call freeunitnumber(in)
  OPEN(UNIT=in,FILE=fname,STATUS='unknown',FORM='FORMATTED',          &
       ACCESS='SEQUENTIAL')
!
! -- Identify package
      WRITE(IOUT,11) IN
   11 FORMAT(1X,/1X,'TDIS -- TEMPORAL DISCRETIZATION PACKAGE,',                &
        ' VERSION 1 : 1/15/2014 - INPUT READ FROM UNIT ',I4)
!
! -- Read comments and the first line following the comments
  CALL URDCOM(IN,IOUT,LINE)
!
! -- Read NPER
  lloc=1
  CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPER,R,IOUT,in)
      WRITE(IOUT,20) NPER
   20 FORMAT(1X,I4,' STRESS PERIOD(S) IN SIMULATION')
!
! -- Allocate arrays for temporal information
  ALLOCATE (PERLEN(NPER),NSTP(NPER),TSMULT(NPER),ISSFLG(NPER))
!
! -- Read and prepare stress period information
      WRITE(IOUT,161)
  161 FORMAT(1X,//1X,'STRESS PERIOD     LENGTH       TIME STEPS',              &
                  '     MULTIPLIER FOR DELT    SS FLAG',/1X,76('-'))
      ISS=0
      ITR=0
      DO 200 N=1,NPER
      READ(IN,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,PERLEN(N),IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSTP(N),R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TSMULT(N),IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF (LINE(ISTART:ISTOP).EQ.'TR') THEN
         ISSFLG(N)=0
         ITR=1
      ELSE IF (LINE(ISTART:ISTOP).EQ.'SS') THEN
         ISSFLG(N)=1
         ISS=1
      ELSE
         WRITE(IOUT,162)
  162    FORMAT(' SSFLAG MUST BE EITHER "SS" OR "TR"',                         &
            ' -- STOP EXECUTION (SGWF2BAS7U1ARDIS)')
         CALL USTOP(' ')
      END IF
      WRITE (IOUT,163) N,PERLEN(N),NSTP(N),TSMULT(N),LINE(ISTART:ISTOP)
  163 FORMAT(1X,I8,1PG21.7,I7,0PF25.3,A11)
!
!-----STOP IF NSTP LE 0, PERLEN EQ 0 FOR TRANSIENT STRESS PERIODS,
!-----TSMULT LE 0, OR PERLEN LT 0..
      IF(NSTP(N).LE.0) THEN
         WRITE(IOUT,164)
  164    FORMAT(1X,/1X,                                                        &
        'THERE MUST BE AT LEAST ONE TIME STEP IN EVERY STRESS PERIOD')
         CALL USTOP(' ')
      END IF
      ZERO=0.
      IF(PERLEN(N).EQ.ZERO .AND. ISSFLG(N).EQ.0) THEN
         WRITE(IOUT,165)
  165    FORMAT(1X,/1X,                                                        &
        'PERLEN MUST NOT BE 0.0 FOR TRANSIENT STRESS PERIODS')
         CALL USTOP(' ')
      END IF
      IF(TSMULT(N).LE.ZERO) THEN
         WRITE(IOUT,170)
  170    FORMAT(1X,/1X,'TSMULT MUST BE GREATER THAN 0.0')
         CALL USTOP(' ')
      END IF
      IF(PERLEN(N).LT.ZERO) THEN
         WRITE(IOUT,175)
  175    FORMAT(1X,/1X,                                                        &
        'PERLEN CANNOT BE LESS THAN 0.0 FOR ANY STRESS PERIOD')
         CALL USTOP(' ')
      END IF
  200 CONTINUE
!
!-----Assign ITRSS.
      IF(ISS.EQ.0 .AND. ITR.NE.0) THEN
         ITRSS=1
         WRITE(IOUT,270)
  270    FORMAT(/,1X,'TRANSIENT SIMULATION')
      ELSE IF(ISS.NE.0 .AND. ITR.EQ.0) THEN
         ITRSS=0
         WRITE(IOUT,275)
  275    FORMAT(/,1X,'STEADY-STATE SIMULATION')
      ELSE
         ITRSS=-1
         WRITE(IOUT,280)
  280    FORMAT(/,1X,'COMBINED STEADY-STATE AND TRANSIENT SIMULATION')
      END IF


  end subroutine tdis_ar

end module tdismodule

