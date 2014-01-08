module CrossModule
use ModelModule
use PackageModule

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
    type(modeltype), pointer :: m1
    type(modeltype), pointer :: m2
    logical :: implicit=.TRUE.
    integer :: ncross
    integer, allocatable, dimension(:) :: nodem1
    integer, allocatable, dimension(:) :: nodem2
    double precision, allocatable, dimension(:) :: cond
    integer, allocatable, dimension(:) :: idxglo
    integer, allocatable, dimension(:) :: idxsymglo
    contains
    !procedure :: initialize
    !procedure :: fill
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
    read(inunit,*) cross%nodem1
    !
    !nodem2
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) cross%nodem2
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



subroutine printname(self)
  implicit none
  class(crosstype), intent(in) :: self
  print *, 'cross Name: ', self%name
end subroutine printname

end module CrossModule

