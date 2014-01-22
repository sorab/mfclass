MODULE SimModule
  implicit none
  integer :: iout
  integer :: iverbose=1 !0: print nothing
                        !1: print first level subroutine headings
END MODULE SimModule

subroutine sim_message(iv,message)
  use simmodule
  implicit none
  integer,intent(in) :: iv
  character(len=*),intent(in) :: message
  if(iv>iverbose) then
    write(iout,'(a)') message
  endif
end subroutine sim_message

  
