MODULE SimModule
  implicit none
  integer :: iout
END MODULE SimModule



program mf2015fort
! ******************************************************************************
! Main MODFLOW-2015 program.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  use SimModule
  use SolutionModule
  use ModelModule
  use CrossModule
  implicit none
  class(modeltype), pointer :: m
  class(solutiontype), pointer :: s
  class(crosstype), pointer :: c
  integer :: is,im,nsg,ic,ipos
  CHARACTER*40 VERSION
  CHARACTER*10 MFVNAM
  PARAMETER (VERSION='0.1.00 01/14/2014')
  PARAMETER (MFVNAM='-2015')
  INTEGER :: IBDT(8)
  INTEGER :: I
! ------------------------------------------------------------------------------
!
! -- Write banner to screen and define constants
  WRITE (*,1) MFVNAM,VERSION
1 FORMAT (/,34X,'MODFLOW',A,/,                                                 &
  4X,'U.S. GEOLOGICAL SURVEY MODULAR FINITE-DIFFERENCE',                       &
  ' GROUNDWATER FLOW MODEL',/,29X,'Version ',A/)
  nsg=1  !hardwire the number of solution groups

! -- Get current date and time, assign to IBDT, and write to screen
      CALL DATE_AND_TIME(VALUES=IBDT)
      WRITE(*,2) (IBDT(I),I=1,3),(IBDT(I),I=5,7)
    2 FORMAT(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ',             &
      I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2,/)

  ! -- Read the name file and initialize instances
  call simulation_init('simulation.nam')
  !
  ! -- Initialize each solution in solutionlist
  do is=1,solutionlist%nsolutions
    call solutionlist%getsolution(s,is)
    call s%initialize()
  enddo
  !
  ! -- Initialize each cross in crosslist
  do ic=1,crosslist%ncrosses
    call crosslist%getcross(c,ic)
    call c%initialize()
  enddo
  !
  ! -- Assign crosses to correct solutions
  do is=1,solutionlist%nsolutions
    call solutionlist%getsolution(s,is)
    allocate(s%crosslist%crosses(1))
    ipos=1
    do ic=1,crosslist%ncrosses
      call crosslist%getcross(c,ic)
      if(c%m1%solutionid==is) then
        call s%addcross(c,ipos)
        ipos=ipos+1
        cycle
      elseif(c%m2%solutionid==is) then
        call s%addcross(c,ipos)
        ipos=ipos+1
        cycle
      endif
    enddo
  enddo
  !
  ! -- Run the connect routines for each solution.  The connect routines build the
  ! -- solution ia and ja arrays as well as the mapping arrays for each
  ! -- model that is part the solution.
  do is=1,solutionlist%nsolutions
    call solutionlist%getsolution(s,is)
    call s%connect()
    call s%smsinit()
  enddo
  !
  ! -- Read stress period data
  do is=1,solutionlist%nsolutions
    call solutionlist%getsolution(s,is)
    do im=1,s%modellist%nmodels
      call s%modellist%getmodel(m,im)
      call m%read_prepare()
    enddo
  enddo
  !
  ! -- Solve each solution
  do is=1,solutionlist%nsolutions
    call solutionlist%getsolution(s,is)
    call s%solve()
  enddo
  !
  ! -- Save each solution
  do is=1,solutionlist%nsolutions
    call solutionlist%getsolution(s,is)
    call s%save('output.dat')
  enddo

  CALL GLO1BAS6ET(IOUT,IBDT,1)

! -- Deallocate everything here
    
  print *, 'ending...'

end program mf2015fort

subroutine simulation_init(simfile)
! ******************************************************************************
! Read the simulation name file and initialize the models, crosses
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  use SimModule
  use TdisModule
  use SolutionModule
  use ModelModule
  use GWFModule
  use CrossModule
  implicit none
  character(len=*),intent(in) :: simfile
  character(len=300) :: line,fname,tag
  character(len=20) :: filtyp
  integer :: lloc,ityp1,ityp2,n,inunit,id,sid,mid,sgid
  integer :: nsolutions, nmodels, ncrosses, m1, m2,nsol,nmod,im
  real :: r
  class(solutiontype), pointer :: sp
  class(modeltype), pointer :: mp
! ------------------------------------------------------------------------------
!
! -- Set defaults and open the simulation name file on inunit
  inunit=99
  iout=-1
  open(unit=inunit, file=simfile)
  WRITE(*,490)' Using Simulation file: ',simfile
490 FORMAT(A,A)
!
! -- Read each line of the simulation name file, skip blank or comment lines and
! -- then parse the line
  do
    ! -- read line and skip if necessary
    read(inunit,'(a)',end=100) line
    if(line.eq.' ') cycle
    if(line(1:1).eq.'#') cycle
    !
    !find number of models
    lloc=1
    CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
    filtyp=line(ityp1:ityp2)
    !
    !open the simulation list file - temporary solution???
    if (filtyp.eq.'LIST') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
      fname=line(ityp1:ityp2)
      call freeunitnumber(iout)
      open(unit=iout,file=fname,status='unknown')
    !
    !set nmodels
    elseif(filtyp.eq.'NMODELS') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,nmodels,R,IOUT,INUNIT)
      call model_list_init(nmodels)
    !
    !create models
    elseif(filtyp=='MODEL') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
      fname=line(ityp1:ityp2)
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,id,R,IOUT,INUNIT)
      call model_create(fname,id)
    !
    elseif(filtyp=='GWFMODEL') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
      fname=line(ityp1:ityp2)
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,id,R,IOUT,INUNIT)
      call gwfmodel_create(fname,id)
    !
    elseif(filtyp=='TDIS') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
      fname=line(ityp1:ityp2)
      call tdis_ar(fname,iout)
    !
    !set nsolutions
    elseif(filtyp.eq.'NSOLUTIONS') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,nsolutions,R,IOUT,INUNIT)
      call solution_list_init(nsolutions)
    !
    !create solutions
    elseif(filtyp=='SOLUTION') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
      fname=line(ityp1:ityp2)
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,id,R,IOUT,INUNIT)
      call solution_create(fname,id)
    !    
    !set ncross
    elseif(filtyp.eq.'NCROSSES') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,ncrosses,R,IOUT,INUNIT)
      call cross_list_init(ncrosses)
    !
    !create solutions
    id=1
    elseif(filtyp=='XRS') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
      fname=line(ityp1:ityp2)
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,m1,R,IOUT,INUNIT)
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,m2,R,IOUT,INUNIT)
      call cross_create(fname,id,m1,m2)
      id=id+1
    !
    !solution groups
    elseif(filtyp=='SOLUTION_GROUP') then
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,sgid,R,IOUT,INUNIT)
      solutiongrouplist(1)%id=sgid
      read(inunit,*) fname,solutiongrouplist(sgid)%mxiter
      read(inunit,*) fname,nsol
      solutiongrouplist(sgid)%nsolutions=nsol
      allocate(solutiongrouplist(sgid)%solutions(nsol))
      do n=1,nsol
        read(inunit,*) fname,sid
        call solutionlist%getsolution(sp,sid)
        read(inunit,*) fname,nmod
        allocate(sp%modellist%models(nmod))
        do im=1,nmod
          read(inunit,*) fname,mid
          call modellist%getmodel(mp,mid)
          call sp%addmodel(mp,im)
        enddo
      enddo
    endif
    cycle
100 exit        
  enddo
  close(inunit)
end subroutine simulation_init

subroutine freeunitnumber(iu)
! ******************************************************************************
! Assign a free unopened unit number to the iu dummy argument.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  implicit none
  integer :: lastunitnumber
  parameter(lastunitnumber=10000)
  integer,intent(inout) :: iu
  integer, save :: nextunitnumber=1000
  integer i
  logical :: opened
! ------------------------------------------------------------------------------
!
  do i=nextunitnumber,lastunitnumber
    inquire(unit=i,opened=opened)
    if(.not. opened) exit
  enddo
  nextunitnumber=i
  iu=nextunitnumber
end subroutine freeunitnumber

      SUBROUTINE GLO1BAS6ET(IOUT,IBDT,IPRTIM)
!     ******************************************************************
!     Get end time and calculate elapsed time
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------
      INTEGER IBDT(8), IEDT(8), IDPM(12)
      DATA IDPM/31,28,31,30,31,30,31,31,30,31,30,31/ ! Days per month
      DATA NSPD/86400/  ! Seconds per day
!     ------------------------------------------------------------------
!
!     Get current date and time, assign to IEDT, and write.
      CALL DATE_AND_TIME(VALUES=IEDT)
      WRITE(*,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
 1000 FORMAT(1X,'Run end date and time (yyyy/mm/dd hh:mm:ss): ',               &
             I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)
      IF(IPRTIM.GT.0) THEN
        WRITE(IOUT,'(1X)')
        WRITE(IOUT,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
      END IF
!
!     Calculate elapsed time in days and seconds
      NDAYS=0
      LEAP=0
      IF (MOD(IEDT(1),4).EQ.0) LEAP = 1
      IBD = IBDT(3)            ! BEGIN DAY
      IED = IEDT(3)            ! END DAY
!     FIND DAYS
      IF (IBDT(2).NE.IEDT(2)) THEN
!       MONTHS DIFFER
        MB = IBDT(2)             ! BEGIN MONTH
        ME = IEDT(2)             ! END MONTH
        NM = ME-MB+1             ! NUMBER OF MONTHS TO LOOK AT
        IF (MB.GT.ME) NM = NM+12
        MC=MB-1
        DO 10 M=1,NM
          MC=MC+1                ! MC IS CURRENT MONTH
          IF (MC.EQ.13) MC = 1
          IF (MC.EQ.MB) THEN
            NDAYS = NDAYS+IDPM(MC)-IBD
            IF (MC.EQ.2) NDAYS = NDAYS + LEAP
          ELSEIF (MC.EQ.ME) THEN
            NDAYS = NDAYS+IED
          ELSE
            NDAYS = NDAYS+IDPM(MC)
            IF (MC.EQ.2) NDAYS = NDAYS + LEAP
          ENDIF
   10   CONTINUE
      ELSEIF (IBD.LT.IED) THEN
!       START AND END IN SAME MONTH, ONLY ACCOUNT FOR DAYS
        NDAYS = IED-IBD
      ENDIF
      ELSEC=NDAYS*NSPD
!
!     ADD OR SUBTRACT SECONDS
      ELSEC = ELSEC+(IEDT(5)-IBDT(5))*3600.0
      ELSEC = ELSEC+(IEDT(6)-IBDT(6))*60.0
      ELSEC = ELSEC+(IEDT(7)-IBDT(7))
      ELSEC = ELSEC+(IEDT(8)-IBDT(8))*0.001
!
!     CONVERT SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
      NDAYS = ELSEC/NSPD
      RSECS = MOD(ELSEC,86400.0)
      NHOURS = RSECS/3600.0
      RSECS = MOD(RSECS,3600.0)
      NMINS = RSECS/60.0
      RSECS = MOD(RSECS,60.0)
      NSECS = RSECS
      RSECS = MOD(RSECS,1.0)
      MSECS = NINT(RSECS*1000.0)
      NRSECS = NSECS
      IF (RSECS.GE.0.5) NRSECS=NRSECS+1
!
!     Write elapsed time to screen
        IF (NDAYS.GT.0) THEN
          WRITE(*,1010) NDAYS,NHOURS,NMINS,NRSECS
 1010     FORMAT(1X,'Elapsed run time: ',I3,' Days, ',I2,' Hours, ',I2,        &
            ' Minutes, ',I2,' Seconds',/)
        ELSEIF (NHOURS.GT.0) THEN
          WRITE(*,1020) NHOURS,NMINS,NRSECS
 1020     FORMAT(1X,'Elapsed run time: ',I2,' Hours, ',I2,                     &
            ' Minutes, ',I2,' Seconds',/)
        ELSEIF (NMINS.GT.0) THEN
          WRITE(*,1030) NMINS,NSECS,MSECS
 1030     FORMAT(1X,'Elapsed run time: ',I2,' Minutes, ',                      &
            I2,'.',I3.3,' Seconds',/)
        ELSE
          WRITE(*,1040) NSECS,MSECS
 1040     FORMAT(1X,'Elapsed run time: ',I2,'.',I3.3,' Seconds',/)
        ENDIF
!
!     Write times to file if requested
      IF(IPRTIM.GT.0) THEN
        IF (NDAYS.GT.0) THEN
          WRITE(IOUT,1010) NDAYS,NHOURS,NMINS,NRSECS
        ELSEIF (NHOURS.GT.0) THEN
          WRITE(IOUT,1020) NHOURS,NMINS,NRSECS
        ELSEIF (NMINS.GT.0) THEN
          WRITE(IOUT,1030) NMINS,NSECS,MSECS
        ELSE
          WRITE(IOUT,1040) NSECS,MSECS
        ENDIF
      ENDIF
!
      RETURN
      END
