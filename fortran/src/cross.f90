module CrossModule
use ModelModule
use PackageModule
use ghbmodule

implicit none

private
public :: crosstype
public :: crosslisttype
public :: getcross
public :: cross_create
public :: cross_list_init
public :: crosslist

type :: crosscontainer
    class(crosstype), pointer :: obj
end type crosscontainer

type crosslisttype
    integer :: ncrosses = 0
    type(crosscontainer), allocatable, dimension(:) :: crosses
    contains
    procedure :: setcross
    procedure :: getcross
end type crosslisttype
type(crosslisttype) :: crosslist

type crosstype
    character(len=20) :: name
    integer :: id
    class(modeltype), pointer :: m1
    class(modeltype), pointer :: m2
    logical :: implicit=.TRUE.
    integer :: ncross
    integer, allocatable, dimension(:) :: nodem1
    integer, allocatable, dimension(:) :: nodem2
    double precision, allocatable, dimension(:) :: cond
    integer, allocatable, dimension(:) :: idxglo
    integer, allocatable, dimension(:) :: idxsymglo
    integer :: m1pakid
    integer :: m2pakid
    type(ghbtype),pointer :: p1
    type(ghbtype),pointer :: p2
    contains
    procedure :: crossinit
    procedure :: crossfmfill
    procedure :: crossbd
    !procedure :: fill_explicit
    !procedure :: fill_implicit
end type crosstype

contains

subroutine setcross(this, newcross, ipos)
    implicit none
    class(crosslisttype) :: this
    class(crosstype), target, intent(in) :: newcross;
    integer, intent(in) :: ipos
    this%crosses(ipos)%obj => newcross
    this%ncrosses=max(ipos,this%ncrosses)
end subroutine setcross

subroutine getcross(this, thecross, ipos)
    implicit none
    class(crosslisttype) :: this
    class(crosstype), pointer :: thecross
    integer, intent(in) :: ipos
    thecross => this%crosses(ipos)%obj
end subroutine getcross

subroutine cross_list_init(ncrosses)
    implicit none
    integer,intent(in) :: ncrosses
    allocate(crosslist%crosses(ncrosses))
end subroutine cross_list_init

subroutine cross_create(filename,id,m1id,m2id)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id,m1id,m2id
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,iout,inunit,npackages,ipak
    real :: r
    type(crosstype), pointer :: cross
    class(modeltype), pointer :: mp

    !create a new cross and add it to the crosslist container
    allocate(cross)
    call crosslist%setcross(cross,id)

    print *,'Creating cross.'
    call modellist%getmodel(mp, m1id)
    cross%m1=>mp
    call modellist%getmodel(mp, m2id)
    cross%m2=>mp
    !
    !open the file
    call freeunitnumber(inunit)
    print *, 'opening cross file on unit: ', inunit
    open(unit=inunit,file=filename)
    !
    !problem size: ncross
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo    
    lloc=1
    CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,cross%ncross,R,IOUT,INUNIT)
    !
    !allocate arrays
    allocate(cross%nodem1(cross%ncross))
    allocate(cross%nodem2(cross%ncross))
    allocate(cross%cond(cross%ncross))
    allocate(cross%idxglo(cross%ncross))
    allocate(cross%idxsymglo(cross%ncross))
    !
    !nodem1
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) n,cross%nodem1
    !
    !nodem2
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) n,cross%nodem2
    !
    !cond
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) cross%cond
    !
    close(inunit)

end subroutine cross_create

subroutine crossinit(this)
    implicit none
    class(crosstype) :: this
    class(packagetype), pointer :: packobj1,packobj2
    character(len=16) :: packname
    integer :: i,ipakid
    print *,'initialize'
    !
    !check if m1 and m2 are part of the same solution.  if not
    !then they must be coupled explicitly
    if(this%m1%solutionid/=this%m2%solutionid) then
        this%implicit=.FALSE.
        !
        !create a ghb-like package for model 1
        ipakid=this%m1%packages%npackages+1
        packname='XRS '//this%m2%name
        packname=adjustr(packname)
        call ghb_new(this%p1,packname,0,this%ncross,this%ncross, &
                     this%nodem1,this%cond)
        packobj1=>this%p1
        call this%m1%packages%setpackage(packobj1,ipakid)
        !
        !create a ghb-like package for model 2
        ipakid=this%m2%packages%npackages+1
        packname='XRS '//this%m1%name
        packname=adjustr(packname)
        call ghb_new(this%p2,packname,0,this%ncross,this%ncross, &
                     this%nodem2,this%cond)
        packobj2=>this%p2
        call this%m2%packages%setpackage(packobj2,ipakid)
    endif
end subroutine crossinit

subroutine crossfmfill(this)
    implicit none
    class(crosstype) :: this
    !class(packagetype), pointer :: p1,p2
    integer :: i,n1,n2
    if(.not. this%implicit) then
        do i=1,this%ncross
            n1=this%nodem1(i)
            n2=this%nodem2(i)
            this%p1%stage(i)=this%m2%x(n2)
            this%p2%stage(i)=this%m1%x(n1)
        enddo
    endif
end subroutine crossfmfill

subroutine crossbd(this)
! -- langevin mf2015 todo -- need to write budget terms to cbc
! -- langevin mf2015 todo: need to put ibound check in here
  implicit none
  class(crosstype) :: this
  character(len=16) :: packname
  integer :: i,n1,n2
  real :: rin,rout,rate,zero
  double precision :: ratin,ratout,rrate
  zero=0.
  ratin=zero
  ratout=zero
  if(this%implicit) then
    do i=1,this%ncross
      n1=this%nodem1(i)
      n2=this%nodem2(i)
      rrate=this%cond(i)*this%m2%x(n2)-this%cond(i)*this%m1%x(n1)
      rate=rrate
      if(rate<zero) then
        ratout=ratout-rrate
      else
        ratin=ratin+rrate
      endif
    enddo
    rin=ratin
    rout=ratout
!
! -- Add the budget terms to model 1
    packname='XRS '//this%m2%name
    packname=adjustr(packname)
    call this%m1%modelbdentry(packname,rin,rout)
!
! -- Add the budget terms to model 2
    packname='XRS '//this%m1%name
    packname=adjustr(packname)
    call this%m2%modelbdentry(packname,rout,rin)
  endif
end subroutine crossbd

subroutine printname(this)
  implicit none
  class(crosstype), intent(in) :: this
  print *, 'cross Name: ', this%name
end subroutine printname

subroutine ghb_new(packobj,name,inunit,maxbound,nbound,nodelist,cond)
    !create a new ghb package for an explicit cross connection
    implicit none
    type(ghbtype), pointer :: packobj
    character(len=*),intent(in) :: name
    integer,intent(in) :: inunit
    integer,intent(in) :: maxbound
    integer,intent(in) :: nbound
    integer,dimension(maxbound),target :: nodelist
    double precision, dimension(maxbound),target :: cond
    type(ghbtype), pointer :: ghbobj
    allocate(ghbobj)
    packobj => ghbobj
    ghbobj%name=name
    ghbobj%inunit=inunit
    call ghbobj%ghb_allocate(maxbound)
    ghbobj%nbound=nbound
    ghbobj%nodelist=>nodelist
    ghbobj%cond=>cond
end subroutine ghb_new

end module CrossModule

