module gwfmodule
  use ModelModule
  use PackageModule
  use global,only:gwfglotype,niunit
  use GWFBASMODULE,only:gwfbastype
  use GWFBCFMODULE,only:gwfbcftype
  private
  public :: gwf3ar
  type, extends(modeltype) :: gwfmodeltype
    type(gwfglotype) :: gwfglodat
    type(gwfbastype) :: gwfbasdat
    type(gwfbcftype) :: gwfbcfdat
  end type gwfmodeltype
  CHARACTER*40 VERSION
  CHARACTER*10 MFVNAM
  integer :: maxunit
  CHARACTER*80 HEADNG(2)
  PARAMETER (VERSION='GWF 3 01/18/2014')
  PARAMETER (MFVNAM='-2015') 
  INTEGER,DIMENSION(:),POINTER :: IUNIT
  CHARACTER*4 CUNIT(NIUNIT)
  DATA CUNIT/   'BCF6', 'WEL ', 'DRN ', 'RIV ', 'EVT ', '    ', 'GHB ',  & !  7
                'RCH ', '    ', '    ', '    ', 'OC  ', 'SMS ', 'PCB ',  & ! 14
                'BCT ', 'FHB ', 'RES ', 'STR ', 'IBS ', 'CHD ', 'HFB6',  & ! 21
                'LAK ', 'LPF ', 'DIS ', 'DISU', 'PVAL', '    ', 'HOB ',  & ! 28
                'CLN ', '    ', 'ZONE', 'MULT', 'DROB', 'RVOB', 'GBOB',  & ! 35
                'GNC ', '    ', 'CHOB', 'ETS ', 'DRT ', '    ', 'GMG ',  & ! 42
                'hyd ', 'SFR ', '    ', 'GAGE', 'LVDA', '    ', 'lmt6',  & ! 49
                'MNW1', '    ', '    ', 'KDEP', 'SUB ', 'UZF ', 'gwm ',  & ! 56
                'SWT ', '    ', '    ', '    ', '    ', '    ', '    ',  & ! 63
                37*'    '/
  contains

  subroutine gwf3ar(filename,id,ioutsim)
! ******************************************************************************
! Main MODFLOW-2015 program.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id,ioutsim
    integer :: inunit
    type(gwfmodeltype), pointer :: gwfmodel
! ------------------------------------------------------------------------------
    
    !create a new model and add it to the modellist container
    allocate(gwfmodel)
    call modellist%setmodel(gwfmodel,id)    
    print *,'Creating gwfmodel: ', id
    !
    ! -- open the gwf name file
    call freeunitnumber(inunit)
    print *, 'opening model namefile on unit: ', inunit
    open(unit=inunit,file=filename,status='old')
    !
    ! -- allocate and read bas and dis and then save to pointers
    CALL GLO2BAS8AR(INUNIT,CUNIT,VERSION,24,31,32,MAXUNIT,12,                  &
                      HEADNG,26,MFVNAM,29,27,30,36)
    call gwfmodel%gwfglodat%pntsav
    call gwfmodel%gwfbasdat%pntsav
    !
    ! -- make a pointer to iunit for easier syntax within these routines
    iunit=>gwfmodel%gwfglodat%iunit
    !
    ! -- allocate and read bcf and then save to pointers
    IF(IUNIT(1).GT.0.OR.IUNIT(23).GT.0)THEN                                     
        CALL GWF2BCFU1AR(IUNIT(1),IUNIT(22),IUNIT(23))
        call gwfmodel%gwfbcfdat%pntsav
    ENDIF
    !
    ! -- set and allocate gwfmodel variables and arrays
    gwfmodel%neq=gwfmodel%gwfglodat%neqs
    gwfmodel%nja=gwfmodel%gwfglodat%nja    
    allocate(gwfmodel%ia(gwfmodel%neq+1))
    allocate(gwfmodel%ja(gwfmodel%nja))
    allocate(gwfmodel%cond(gwfmodel%nja))
    allocate(gwfmodel%rhs(gwfmodel%neq))
    allocate(gwfmodel%idxglo(gwfmodel%nja))
    !
    return
  end subroutine gwf3ar

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