module gwfmodule
  use ModelModule
  use PackageModule
  use global,only:gwfglotype
  use GWFBASMODULE,only:gwfbastype
  use GWFBCFMODULE,only:gwfbcftype
  private
  public :: gwfmodel_create
  type, extends(modeltype) :: gwfmodeltype
    type(gwfglotype) :: gwfglodat
    type(gwfbastype) :: gwfbasdat
    type(gwfbcftype) :: gwfbcfdat
  end type gwfmodeltype

  contains

  subroutine gwfmodel_create(filename,id,ioutsim)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id,ioutsim
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,iout,inunit,npackages,ipak,iunstr
    integer :: nfile,istart,istop,iu,inam1,inam2,iflen,LENVER,INDENT
    logical :: lop
    real :: r
    CHARACTER*40 VERSION,SPACES
    CHARACTER*10 MFVNAM
    PARAMETER (VERSION='GWF 3 01/18/2014')
    PARAMETER (MFVNAM='-2015')    
    type(gwfmodeltype), pointer :: gwfmodel
    
    SPACES=' '
    LENVER=LEN_TRIM(VERSION)
    INDENT=40-(LENVER+8)/2    

    !create a new model and add it to the modellist container
    allocate(gwfmodel)
    call modellist%setmodel(gwfmodel,id)
!
! -- temporarily set iout to ioutsim.  if 'list' is present in name file
! -- then iout will be set to that unit number
    iout=ioutsim
!    
    print *,'Creating gwfmodel: ', id
    call freeunitnumber(inunit)
    print *, 'opening model namefile on unit: ', inunit
    
    open(unit=inunit,file=filename,status='old')
    ipak=1
    nfile=0
    do
      !read line and skip if necessary
      read(inunit,'(a)',end=100) line
      if(line.eq.' ') cycle
      if(line(1:1).eq.'#') cycle
      !
      !decode
      LLOC=1
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
      FILTYP=LINE(ITYP1:ITYP2)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IU,R,IOUT,INUNIT)
      CALL URWORD(LINE,LLOC,INAM1,INAM2,0,N,R,IOUT,INUNIT)
      IFLEN=INAM2-INAM1+1
      FNAME(1:IFLEN)=LINE(INAM1:INAM2)
      INQUIRE(UNIT=IU,OPENED=LOP)
      IF(LOP) THEN
         IF(IOUT.EQ.0) THEN
            WRITE(*,11) FNAME(1:IFLEN),IU
   11       FORMAT(1X,/1X,'CANNOT OPEN ',A,' ON UNIT',I4,                      &
                    ' BECAUSE UNIT IS ALREADY BEING USED')
         ELSE
            WRITE(IOUT,11) FNAME(1:IFLEN),IU
         END IF
         CALL USTOP(' ')
      END IF
!      
      IF(NFILE.EQ.0) THEN
        IF(FILTYP.EQ.'LIST') THEN
          IOUT=IU
          OPEN(UNIT=IU,FILE=FNAME(1:IFLEN),STATUS='REPLACE',                   &
                FORM='FORMATTED',ACCESS='SEQUENTIAL')
          WRITE(IOUT,60) MFVNAM,SPACES(1:INDENT),VERSION(1:LENVER)
60        FORMAT(34X,'MODFLOW',A,/,                                            &
                   6X,'U.S. GEOLOGICAL SURVEY MODULAR',                        &
                   ' FINITE-DIFFERENCE GROUNDWATER FLOW MODEL',/,              &
                   A,'VERSION ',A,/)
          WRITE(IOUT,78) FNAME(1:IFLEN),IOUT
78        FORMAT(1X,'LIST FILE: ',A,/25X,'UNIT ',I4)
        ELSE
          WRITE(*,*) ' FIRST ENTRY IN NAME FILE MUST BE "LIST".'
          CALL USTOP(' ')
        END IF
! -- Get next file name
        NFILE=1
        cycle
      END IF      
      !
      if(filtyp.eq.'DISU') then
        iunstr=1
        open(unit=iu,file=FNAME(1:IFLEN),status='old')
        call GWF3DIS9AR(iu,iout,iunstr)
      !
      elseif(filtyp.eq.'DIS') then
        iunstr=0
        open(unit=iu,file=FNAME(1:IFLEN),status='old')
        call GWF3DIS9AR(iu,iout,iunstr)
      !
      elseif(filtyp.eq.'BAS6') then
        open(unit=iu,file=FNAME(1:IFLEN),status='old')
        call GWF3BAS9AR(iu,VERSION,MFVNAM)
      !
      !create packages
      elseif(filtyp.eq.'NPACKAGES') then
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,npackages,R,IOUT,INUNIT)
        !allocate(model%packages(npackages))
        ipak=1
        cycle
      !
      !create all other package types
      else
        print *,'Not doing anything yet with: ', FNAME(1:IFLEN)
        ipak=ipak+1
        cycle
      endif
    cycle
  100 exit        
    enddo
    close(inunit)
    !
    ! -- save the pointers
    call gwfmodel%gwfglodat%pntsav
    call gwfmodel%gwfbasdat%pntsav
    return
  end subroutine gwfmodel_create

  subroutine package_create(fname,filtyp,ipakid,packages)
    use PackageModule
    use welmodule
    use ghbmodule
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: filtyp
    integer,intent(in) :: ipakid
    type(packagelist),intent(inout) :: packages
    class(packagetype), pointer :: packobj
    if(filtyp=='WEL') then
      call wel_create(packobj,fname,ipakid)
    elseif(filtyp=='GHB') then
      call ghb_create(packobj,fname,ipakid)
    endif
    call packages%setpackage(packobj,ipakid)
  end subroutine package_create

end module gwfmodule