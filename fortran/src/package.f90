module PackageModule
private
public :: packagetype
public :: packagelist

type :: packagecontainer
    class(packagetype), pointer :: obj
end type PackageContainer

!langevin mf2015 todo: note maxpackages is hardwired for new-style
!packages.  Need to fix this so it is based on entries in namefile
integer, parameter :: maxpackages=20
type packagelist
    integer :: npackages = 0
    type(PackageContainer), dimension(maxpackages) :: packages
    contains
    procedure :: setpackage
    procedure :: getpackage
end type packagelist

type packagetype
    character(len=20) :: name
    integer :: inunit
    integer :: maxbound
    integer :: nbound
    real :: rin
    real :: rout 
    integer,pointer,dimension(:) :: nodelist
    double precision, allocatable, dimension(:) :: hcof
    double precision, allocatable, dimension(:) :: rhs
contains
    procedure :: packagest
    procedure :: packagerp
    procedure :: packagead
    procedure :: fmcalc => packagefmcalc
    procedure :: packagebd
    procedure :: pack_allocate
end type packagetype

contains

subroutine setpackage(this, newpackage, ipos)
    implicit none
    class(packagelist) :: this
    class(packagetype), target, intent(in) :: newpackage
    integer, intent(in) :: ipos
    this%packages(ipos)%obj => newpackage
    this%npackages=max(ipos,this%npackages)
end subroutine setpackage

subroutine getpackage(this, thepackage, ipos)
    implicit none
    class(packagelist) :: this
    class(packagetype), pointer :: thepackage
    integer, intent(in) :: ipos
    thepackage => this%packages(ipos)%obj
end subroutine getpackage

subroutine pack_allocate(this,maxbound)
    implicit none
    class(packagetype) :: this
    integer,intent(in) :: maxbound
    this%maxbound=maxbound
    allocate(this%nodelist(maxbound))
    allocate(this%hcof(maxbound))
    allocate(this%rhs(maxbound))
end subroutine pack_allocate

subroutine packagest(this)
    implicit none
    class(packagetype) :: this
    !this package has no ST routine
end subroutine packagest

subroutine packagerp(this)
    implicit none
    class(packagetype) :: this
    !this package has no RP routine
end subroutine packagerp

subroutine packagead(this)
    implicit none
    class(packagetype) :: this
    !this package has no AD routine
end subroutine packagead

subroutine packagefmcalc(this)
    implicit none
    class(packagetype) :: this
    print *, 'you should never see this.  this should be overridden.'
end subroutine packagefmcalc

subroutine packagebd(this,x)
! -- langevin mf2015 todo: need to put ibound check in here
    implicit none
    class(packagetype) :: this
    double precision,dimension(*),intent(in) :: x
    integer :: i,node
    real :: rate,zero
    double precision ratin,ratout,rrate
    zero=0.
    ratin=zero
    ratout=zero
    do i=1,this%nbound
        node=this%nodelist(i)
        rrate=this%hcof(i)*x(node)-this%rhs(i)
        rate=rrate
        if(rate<zero) then
          ratout=ratout-rrate
        else
          ratin=ratin+rrate
        endif
    enddo
    this%rin=ratin
    this%rout=ratout
end subroutine packagebd

end module PackageModule

