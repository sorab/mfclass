! -- langevin: this is a start on a new-style well package
! -- that can be added as a package to any model.  It builds
! -- on the information in wel8subs.f, which contains the
! -- modflow-usg well package and some edits to make it
! -- more general.  Next steps are to try and remove all
! -- use statements in the well package.
!
module wel8module
use PackageModule
use GWFWELMODULE,only:gwfweltype
private
public :: wel8_create
type, extends(packagetype) :: weltype
    double precision, allocatable, dimension(:) :: q
    contains
    procedure :: wel_ar
    procedure :: wel_allocate
    procedure :: packagerp=>welrp
    procedure :: fmcalc => welfmcalc
end type weltype
type(gwfweltype) :: gwfweldat

contains

subroutine wel8_create(packobj,fname,id,inunit,iout)
    implicit none
    class(packagetype), pointer :: packobj
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    integer,intent(in) :: inunit
    integer,intent(in) :: iout
    type(weltype), pointer :: welobj
    integer :: maxbound
    allocate(welobj)
    packobj => welobj
    welobj%inunit=inunit
    welobj%iout=iout
    call welobj%wel_ar(fname,id)
end subroutine wel8_create

subroutine wel_ar(this,fname,id)
    implicit none
    class(weltype) :: this
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    write(this%iout,*) 'Creating wel8 package: ', id
    this%filtyp='WEL'
    write(this%name,'(a,i1)') 'WEL8_',id
    write(this%iout,*) 'opening well input file on unit: ', this%inunit
    open(unit=this%inunit,file=fname,status='old')
    call GWF3WEL8U1AR(this%inunit,this%iout)
    call gwfweldat%pntsav
    this%maxbound=gwfweldat%mxwell
    !
    !allocate arrays in package superclass and in weltype
    call this%wel_allocate(this%maxbound)
end subroutine wel_ar

subroutine wel_allocate(this,maxbound)
    !allocate wel package of size maxbound
    implicit none
    class(weltype) :: this
    integer,intent(in) :: maxbound
    call this%pack_allocate(maxbound)
    allocate(this%q(maxbound))
end subroutine wel_allocate

subroutine welrp(this)
    implicit none
    class(weltype) :: this
    integer :: i
    call gwfweldat%pntset
    call GWF3WEL8U1RP(this%inunit,this%iout)
    this%nbound=gwfweldat%nwells
    this%nodelist=gwfweldat%well(1,:)
end subroutine welrp

subroutine welfmcalc(this)
    implicit none
    class(weltype) :: this
    integer :: i
    call gwfweldat%pntset
    this%hcof(:)=0.
    this%rhs(:)=0.
    call GWF3WEL8U1FM(this%rhs)
end subroutine welfmcalc

end module wel8module

