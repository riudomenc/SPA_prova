! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_cmd_line

  !! Routines to handle command line input. !!

  implicit none

  integer,private :: cmd_line_len = 5000

  public :: parse_spa_cmd_line

contains
  !
  !----------------------------------------------
  !
  subroutine parse_spa_cmd_line( config_filename )

    ! Looks through the command line, and does appropriate   !
    ! things with arguments etc. Returns name of config-file !

    use log_tools
    use scale_declarations, only: fname_length

    implicit none

    ! arguments..
    character(len=fname_length),intent(out) :: config_filename 

    ! local variables..
    character(len=cmd_line_len) :: entire_cmd_line
    integer                     :: i, nos_of_args
    character(len=fname_length) :: argument, this_program

    ! Get the program name..
    call get_command_argument( 0 , this_program )

    ! Get the number of arguments..
    nos_of_args = command_argument_count()

    ! if user hasn't supplied any arguments at all, then ask for name..
    if ( nos_of_args .gt. 0 ) then

       ! loop over command line arguments
       i = 0
       do
          i = i + 1
          if ( i .gt. nos_of_args ) exit

          call get_command_argument( i , argument )

          ! check if it is an option
          if ( argument(1:1) .eq. '-' ) then
             select case ( trim( argument ) )
             case ('-h') ! print help and exit..

                call print_cmd_line_help( this_program )

             case default ! print help and exit..

                write(*,*) 'Unknown option ',trim( argument )
                call print_cmd_line_help( this_program )

             end select
          else
             ! it's not an option - assume it's our config file...
             call get_command_argument( 1 , config_filename )
             call write_log( "Will use runtime parameters from : "//trim(config_filename) )
          end if
       end do

    else ! user didn't supply any arguments, so...

       write(*,*) 'Enter name of configuration file to be read:'
       read(*,'(a)') config_filename
       write(*,*) 'Thank you'

    end if

  end subroutine parse_spa_cmd_line
  !
  !-------------------------------------------------------------------
  !
  subroutine print_cmd_line_help( this_program )

    ! Print a help statement !

    implicit none

    character(len=*),intent(in) :: this_program

    write(*,*) 'Usage: ',trim(this_program),' [options] configname'
    write(*,*) 'where [options] are'
    write(*,*) '  -h:          this message'

    stop "Stopped."

  end subroutine print_cmd_Line_help
  !
  !----------------------------------------------
  !
end module spa_cmd_line
!
!----------------------------------------------
!
