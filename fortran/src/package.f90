module PackageModule
private
public :: packagetype
public :: packagelist

type :: packagecontainer
    class(packagetype), pointer :: obj
end type PackageContainer

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
    integer,pointer,dimension(:) :: nodelist
    double precision, allocatable, dimension(:) :: hcof
    double precision, allocatable, dimension(:) :: rhs
contains
    procedure :: packagest
    procedure :: packagerp
    procedure :: fmcalc => packagefmcalc
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
    print *, 'packagest.  if you see this, this package has no ST routine.'
end subroutine packagest

subroutine packagerp(this)
    implicit none
    class(packagetype) :: this
    print *, 'packagerp.  if you see this, this package has no RP routine'
end subroutine packagerp

subroutine packagefmcalc(this)
    implicit none
    class(packagetype) :: this
    print *, 'you should never see this.  this should be overridden.'
end subroutine packagefmcalc

end module PackageModule

