! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

                                                                   
module log_tools

  !! This module reuses/modifies code that was originally developed !!
  !!  as part of the Glimmer ice-sheet modelling project.           !!
  !!                                                                !!
  !! It provides file logging and error/message handling.Six levels !!
  !!  of message/error are defined:                                 !!
  !!                                                                !!
  !!  msg_diagnostic : Diagnostic messages                          !!
  !!  msg_timestep   : Timestep enumeration and related information !!
  !!  msg_info       : Information messages                         !!
  !!  msg_warning    : Warning messages                             !!
  !!  msg_error      : Error messages                               !!
  !!  msg_fatal      : Fatal error messages                         !!
  !!                                                                !!
  !! These are numbered 1--6, with increasing severity, and the     !!
  !!  level of message output may be set to output all messages,    !!
  !!  only those above a particular severity, or none at all. It    !!
  !!  should be noted that even if all messages are turned off, if  !!
  !!  a fatal error is encountered this will still cause a halt!    !!
  !!                                                                !!
  !! The other point to note is that when calling the messaging     !!
  !!  routines, the numerical identifier of a message level should  !!
  !!  be replaced by the appropriate parameter (msg_<see above>).   !!

  implicit none

  character,parameter :: dirsep = '/'

  integer,parameter :: msg_diagnostic = 1 ! Numerical identifier for diagnostic messages.
  integer,parameter :: msg_timestep   = 2 ! Numerical identifier for timestep messages.
  integer,parameter :: msg_info       = 3 ! Numerical identifier for information messages.
  integer,parameter :: msg_warning    = 4 ! Numerical identifier for warning messages.
  integer,parameter :: msg_error      = 5 ! Numerical identifier for (non-fatal) error messages.
  integer,parameter :: msg_fatal      = 6 ! Numerical identifier for fatal error messages.

  integer, parameter            :: msg_levels = 6
  logical, private, dimension(msg_levels) :: msg_show = .false.

  character(len=*), parameter, dimension(0:msg_levels), private :: msg_prefix = (/ &
       '* UNKNOWN      ', &
       '*              ', &
       '*              ', &
       '               ', &
       '* WARNING:     ', &
       '* ERROR:       ', &
       '* FATAL ERROR :' /)

  character(len=100),private :: logfile_name    ! name of log file
  integer,private            :: logfile_unit=6  ! log unit

  character(len=400),public :: message

contains
  !
  !----------------------------------------------------------------------------
  !
  function logname(fname)
    ! derives name of log file from file name by stripping directories and appending .log
    implicit none
    character(len=*), intent(in) :: fname
    character(len=100)           :: logname

    character(len=*), parameter :: suffix='.log'
    integer i
    i = scan(fname,dirsep,.True.)
    if (i.ne.0) then
       logname = trim(fname(i+1:))//suffix
    else
       logname = trim(fname)//suffix
    end if
  end function logname
  !
  !----------------------------------------------------------------------------
  ! 
  subroutine open_log(unit,fname)
    ! opens log file
    implicit none
    integer, optional          :: unit   ! file unit to use
    character(len=*), optional :: fname  ! name of log file

    ! local variables..
    character(len=8) :: date
    character(len=10) :: time

    if (present(unit)) then
       logfile_unit = unit
    end if
    if (present(fname)) then
       logfile_name = adjustl(trim(fname))
    else
       logfile_name = 'log.log'
    end if

    if (logfile_unit.ne.6) then
       open(unit=logfile_unit,file=logfile_name,status='unknown')
    end if

    call date_and_time(date,time)
    call write_log_div
    write(unit=logfile_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Started logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div

    write(*,*)" Program is running.  All log output will go to: ",trim(logfile_name)

  end subroutine open_log
  !
  !----------------------------------------------------------------------------
  !
  subroutine write_log( message , type , file , line )
    ! write to log
    implicit none
    integer,intent(in),optional          :: type    ! Type of error to be generated (see list above).
    character(len=*),intent(in)          :: message ! message to be written
    character(len=*),intent(in),optional :: file    ! the name of the file which triggered the message
    integer,intent(in),optional          :: line    ! the line number at the which the message was triggered

    ! local variables..
    character(len=len_trim(message)+50) :: msg
    integer :: local_type 
    character(len=6) :: line_num

    local_type = 0
    if (present(type)) then
       if (type.ge.1 .or. type.le.msg_levels) then
          local_type = type
       end if
    else
       local_type = msg_info
    end if

    ! constructing message
    if (present(file) .and. present(line)) then
       write(line_num,'(I6)')line
       write(msg,*) trim(msg_prefix(local_type))//' (',trim(file),':',trim(adjustl(line_num)),') '//trim(message)
    else
       write(msg,*) trim(msg_prefix(local_type))//' '//trim(message)
    end if
    ! messages are always written to file log
    write(logfile_unit,*) trim(msg)
    ! and maybe to std out
    if (local_type.ne.0) then
       if (msg_show(local_type)) write(*,*) trim(msg)
    end if
    ! stop logging if we encountered a fatal error
    if (local_type.eq.msg_fatal) then
       call close_log( local_type ) 
    end if
  end subroutine write_log
  !
  !----------------------------------------------------------------------------
  !
  subroutine write_log_div
    ! start a new section
    implicit none
    write(logfile_unit,*) '*******************************************************************************'
  end subroutine write_log_div
  !
  !----------------------------------------------------------------------------
  !
  subroutine close_log( local_type )
    ! close log file
    implicit none
    ! arguments..
    integer,intent(in) :: local_type
    ! local variables..
    character(len=8) :: date
    character(len=10) :: time
    call date_and_time(date,time)
    call write_log_div
    write(unit=logfile_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Finished logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div

    close(logfile_unit)

    if ( local_type .eq. msg_fatal ) then
      write(*,*)"Fatal message encountered.  See logfile for details.  Stopping."
      stop
    else
      write(*,*)"Log closed."
    endif
  end subroutine close_log
  !
  !----------------------------------------------------------------------------
  !
  subroutine sync_log
    ! synchronise log to disk
    implicit none
    close(logfile_unit)
    open(unit=logfile_unit,file=logfile_name, position="append", status='old')
  end subroutine sync_log
  !
  !----------------------------------------------------------------------------
  !
  subroutine set_msg_level(level)
    ! Sets the output message level.
    integer, intent(in) :: level ! The message level (6 is all messages; 0 is no messages). 
    integer :: i

    do i = 1,msg_levels
       if ( i > ( msg_levels - level ) ) then
          msg_show(i) = .true.
       else
          msg_show(i) = .false.
       endif
    enddo

  end subroutine set_msg_level
  !
  !----------------------------------------------------------------------------
  !
  integer function get_logunit()
    ! return log unit
    implicit none

    get_logunit = logfile_unit
  end function get_logunit
  !
  !----------------------------------------------------------------------------
  !  
end module log_tools
!
!----------------------------------------------------------------------------
!
