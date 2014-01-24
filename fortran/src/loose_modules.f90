MODULE SimModule
  implicit none
  integer :: iout
  integer :: iverbose=0 !0: print nothing
                        !1: print first level subroutine information
END MODULE SimModule

subroutine sim_message(iv,message)
! -- iv is the verbosity level of this message
! --  (1) means primary subroutine for simulation, cross, model, 
! --      solution, package, etc.
! -- message is a character string message to write
  use simmodule
  implicit none
  integer,intent(in) :: iv
  character(len=*),intent(in) :: message
  if(iv<=iverbose) then
    write(iout,'(a)') message
  endif
end subroutine sim_message

  
