module ModelModule
use PackageModule

implicit none

private
public :: model_list_init
public :: modeltype
public :: modellisttype
public :: getmodel
public :: model_create
public :: modellist

type :: modelcontainer
    class(modeltype), pointer :: obj
end type modelcontainer

type modellisttype
    integer :: nmodels = 0
    type(modelcontainer), allocatable, dimension(:) :: models
    contains
    procedure :: setmodel
    procedure :: getmodel
end type modellisttype
type(modellisttype) :: modellist

type modeltype
    character(len=20) :: name
    integer :: id
    integer :: neq
    integer :: nja
    integer :: offset
    integer :: solutionid=0
    integer, allocatable, dimension(:) :: ia
    integer, allocatable, dimension(:) :: ja
    double precision, pointer, dimension(:) :: x
    double precision, pointer, dimension(:) :: rhs
    double precision, allocatable, dimension(:) :: cond
    integer, allocatable, dimension(:) :: idxglo
    type(packagelist) :: packages
    integer, dimension(10) :: models
    integer, dimension(10) :: crosses
    contains
    procedure :: disread
    procedure :: modelst
    procedure :: modelrp
    procedure :: modelad
    procedure :: modelfmcalc
    procedure :: modelfmfill
    procedure :: modelbd
    procedure :: printname
end type modeltype

contains

subroutine model_list_init(nmodels)
    implicit none
    integer,intent(in) :: nmodels
    allocate(modellist%models(nmodels))
end subroutine model_list_init

subroutine setmodel(this, newmodel, ipos)
    implicit none
    class(modellisttype) :: this
    class(modeltype), target, intent(in) :: newmodel
    integer, intent(in) :: ipos
    this%models(ipos)%obj => newmodel
    this%nmodels=max(ipos,this%nmodels)
end subroutine setmodel

subroutine getmodel(this, themodel, ipos)
    implicit none
    class(modellisttype) :: this
    class(modeltype), pointer :: themodel
    integer, intent(in) :: ipos
    themodel => this%models(ipos)%obj
end subroutine getmodel

subroutine model_create(filename,id)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,iout,inunit,npackages,ipak
    real :: r
    type(modeltype), pointer :: model

    !create a new model and add it to the modellist container
    allocate(model)
    call modellist%setmodel(model,id)

    print *,'Creating model: ', id
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
            call model%disread(fname,id)
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
            call package_create(fname,filtyp,ipak,model%packages)
            ipak=ipak+1
        endif
    cycle
100 exit        
    enddo
    close(inunit)

end subroutine model_create

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


subroutine disread(this,filename,id)
    implicit none
    class(modeltype) :: this
    character(len=*), intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line
    character(len=20) :: name
    integer :: lloc,ityp1,ityp2,n,iout,inunit
    real :: r
    this%id=id
    call freeunitnumber(inunit)
    print *, 'opening model dis file on unit: ', inunit
    open(unit=inunit,file=filename)
    !
    !problem size: neq, nja
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo    
    lloc=1
    CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,this%neq,R,IOUT,INUNIT)
    CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,this%nja,R,IOUT,INUNIT)
    !
    !allocate arrays
    allocate(this%ia(this%neq+1))
    allocate(this%ja(this%nja))
    allocate(this%cond(this%nja))
    allocate(this%rhs(this%neq))
    allocate(this%idxglo(this%nja))
    !
    !ia
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%ia
    !
    !ja
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%ja
    !
    !cond
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%cond
    !
    !rhs
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%rhs
    close(inunit)

end subroutine disread

subroutine modelst(this)
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
    !
    !stress timing
    print *,'modelst'
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        call p%packagest()
    enddo
end subroutine modelst

subroutine modelrp(this)
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
    !
    !read and prepare
    print *,'modelrp'
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        call p%packagerp()
    enddo
end subroutine modelrp

subroutine modelad(this)
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
    !
    !advance
    print *,'modelad'
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        call p%packagead()
    enddo
end subroutine modelad

subroutine modelfmcalc(this)
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
    !
    !this is where conductance could be calculated
    !
    !calculate the rhs and hcof terms for each package
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        call p%fmcalc()
    enddo
end subroutine modelfmcalc

subroutine modelfmfill(this,amatsln,njasln)
    !fill amatsln and rhssln with amat and rhs from this model
    implicit none
    class(modeltype) :: this
    double precision,dimension(njasln),intent(inout) :: amatsln
    integer,intent(in) :: njasln
    class(packagetype), pointer :: p
    integer :: ip,n,i,ipos
    print *,'modelfmfill'
    !
    !copy the model conductance into the solution amat
    do ipos=1,this%nja
        amatsln(this%idxglo(ipos))=this%cond(ipos)
    enddo
    !
    !copy the package rhs and hcof into the rhs and amat
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        do i=1,p%nbound
            n=p%nodelist(i)
            this%rhs(n)=this%rhs(n)+p%rhs(i)
            ipos=this%ia(n)
            amatsln(this%idxglo(ipos))=amatsln(this%idxglo(ipos))+p%hcof(i)
        enddo
    enddo
end subroutine modelfmfill

subroutine modelbd(this)
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
    !
    !advance
    print *,'modelbd'
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        call p%packagebd(this%x)
    enddo
end subroutine modelbd

subroutine printname(this)
  implicit none
  class(modeltype), intent(in) :: this
  print *, 'Model Name: ', this%name
end subroutine printname

end module ModelModule
