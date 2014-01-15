module gwfmodule
  use ModelModule
  use PackageModule
  use global,only:gwfglotype
  private
  public :: gwfmodel_create
  type, extends(modeltype) :: gwfmodeltype
    type(gwfglotype) :: gwf_data
  end type gwfmodeltype

  contains

  subroutine gwfmodel_create(filename,id)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,iout,inunit,npackages,ipak
    real :: r
    type(gwfmodeltype), pointer :: gwfmodel

    !create a new model and add it to the modellist container
    allocate(gwfmodel)
    call modellist%setmodel(gwfmodel,id)

    print *,'Creating gwfmodel: ', id
    call freeunitnumber(inunit)
    print *, 'opening model namefile on unit: ', inunit
    
    open(unit=inunit,file=filename)
    ipak=1
    do
      !read line and skip if necessary
      read(inunit,'(a)',end=100) line
      if(line.eq.' ') cycle
      if(line(1:1).eq.'#') cycle
      !
      !decode
      lloc=1
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
      filtyp=line(ityp1:ityp2)
      !
      if(filtyp.eq.'DIS') then
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
        fname=line(ityp1:ityp2)
        call gwfmodel%disread(fname,id)
      !
      !create packages
      elseif(filtyp.eq.'NPACKAGES') then
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,npackages,R,IOUT,INUNIT)
        !allocate(model%packages(npackages))
        ipak=1
      !
      !create all other package types
      else
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
        fname=line(ityp1:ityp2)
        call package_create(fname,filtyp,ipak,gwfmodel%packages)
        ipak=ipak+1
      endif
    cycle
  100 exit        
    enddo
    close(inunit)
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