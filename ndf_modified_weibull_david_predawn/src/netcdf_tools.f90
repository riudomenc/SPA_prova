! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module netcdf_tools

  !! This module declares generic derived types and procedures    !!
  !! that allow us to more easily manipulate the data associated  !!
  !! with netcdf files.                                           !!
  !!                                                              !!
  !! The derived types allow us to declare a general 'type', such !!
  !! as nc_dimension, with a bunch of properties, such as 'name', !!
  !! and then build more complex things, such as the 'nc_met_file'!!
  !! which has a nc_header (containing several dimensions) and a  !!
  !! series of 'nc_variables..                                    !!

  use scale_declarations, only: fname_length

  implicit none

  ! A generic dimension...
  type nc_dimension
     character(fname_length)   :: name     ! name of dimension
     integer                   :: dim_id = -999    ! dim id
     integer                   :: var_id = -999    ! var id
     integer                   :: length =  0      ! size of dimension
     real,dimension(:),pointer :: values => null() ! values of dimension
  end type nc_dimension

  ! A generic variable..
  type nc_variable
     character(fname_length)    :: name = ''          ! name of variable
     integer                    :: id   = -999        ! var id
     real                       :: mdi  = -1e32       ! missing data indicator
     real,dimension(:,:,:),pointer :: data  => null() ! contains the variable's values
     type(nc_variable),pointer  :: next     => null() ! points to next variable in list!
     type(nc_variable),pointer  :: previous => null() ! points to next variable in list!
  end type nc_variable

  ! A typical file header (comprising file info and dimensions)..
  type nc_header
     character(fname_length)    :: name    =  ''        ! file name
     integer                    :: id      =  -999      ! file handle
     logical                    :: define_mode =.False.   ! indicates read/write (false) or define (true)
     logical                    :: exists = .False.     ! did the file exist before opening?
     type(nc_dimension),pointer :: time    => null()    ! time
     type(nc_dimension),pointer :: lat     => null()    ! latitude
     type(nc_dimension),pointer :: lon     => null()    ! longitude
  end type nc_header

  ! A generic output file - note this is a linked list (using the %next pointer)
  !  so may contain information about multiple output files.
  type nc_output
     type(nc_header),pointer :: header    => null() ! information about the file, such as name & dims
     integer                 :: frequency =   999   ! (days) between writes.
     type(nc_output),pointer :: next      => null() ! link to next output file.
     type(nc_output),pointer :: previous  => null() ! link to next output file.
  end type nc_output

  ! --- private to this module ---
  private :: create_var

  ! --- interface to handle coord-dims/var writing ---
  interface write_to_nc
    module procedure put_dim_nc, put_var_nc
  end interface

contains
  !
  !-------------------------------------------------------------------
  !
  subroutine open_nc_file( file_info , for_output )

    ! opens a netcdf file for reading/writing. !

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    type(nc_header),pointer     :: file_info    ! handle to file information.
    logical,intent(in),optional :: for_output   ! flag to indicate whether to open file for writing

    ! local variables..
    integer :: status
    logical :: write_mode,file_exists

    ! if not given, assume open read-only...
    if (.not. present( for_output ) ) then
       write_mode = .False.
    else
       write_mode = for_output
    endif

    ! open appropriately..
    if ( write_mode ) then

       ! check if file already exists...
       inquire( file=trim(file_info%name) , exist=file_info%exists ,iostat=status )
       if ( status .ne. 0 ) call write_log("Could not get information about nc file",msg_fatal,__FILE__,__LINE__)

       if ( file_info%exists ) then
          ! open in write-mode..
          status = nf90_open( trim(file_info%name) , nf90_write , file_info%id )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          file_info%define_mode = .False.
          write(message,*)"Opened existing file "//trim(file_info%name)//" for writing to."
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
       else
          ! create new file
          status = nf90_create( trim(file_info%name) , nf90_noclobber , file_info%id )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          file_info%define_mode = .True.
          write(message,*)"Created new file "//trim(file_info%name)//" for writing to."
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
       endif
       call handle_nc_error( __FILE__ , status , __LINE__ )
    else
       ! opening file for read-only..
       status = nf90_open( trim(file_info%name) , nf90_nowrite , file_info%id )
       call handle_nc_error( __FILE__ , status , __LINE__ )
       file_info%define_mode = .False.
       write(message,*)"Opened existing file "//trim(file_info%name)//" for reading from"
       call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
       file_info%exists = .True.
    endif

  endsubroutine open_nc_file
  !
  !-------------------------------------------------------------------
  !
  subroutine change_mode( nc_handle , define )

    ! Take a given ncdf handle and puts the file into define !
    ! (True) or write (False) mode.  As part of this, it     !
    ! first check what mode the file is currently in.        !

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    type(nc_header) :: nc_handle
    logical         :: define    ! True => want file in define mode
    ! False => want file in read/write mode

    ! local variables..
    integer :: status

    ! Only do something if the file's current mode
    !  differs from what is wanted...
    if ( nc_handle%define_mode .neqv. define ) then

       if ( define ) then ! assume file is in read/write mode..

          status = nf90_redef( nc_handle%id )
          call handle_nc_error( __FILE__ , status , __LINE__ )
          nc_handle%define_mode = .True.
          call write_log( "Changed file to define mode" , msg_info , __FILE__ , __LINE__ )

       else ! assume file is in define mode..

          status = nf90_enddef( nc_handle%id )
          if ( status .ne. 0) then
            call handle_nc_error( __FILE__ , status , __LINE__ )
          endif
          nc_handle%define_mode = .False.
          call write_log( "Changed file to read/write mode" , msg_info , __FILE__ , __LINE__ )

       endif

    endif

  end subroutine change_mode
  !
  !-------------------------------------------------------------------
  !
  subroutine get_dim_info( file_id , dim )

    ! Returns the dim-id,var-id and size of !
    ! a given dimension in a given file.    !

    use netcdf

    ! arguments..
    integer,intent(in) :: file_id
    type(nc_dimension),pointer :: dim

    ! local variables..
    integer :: status, lala
    character(len=100) :: blah

    ! get the dimension id..
    status = nf90_inq_dimid( file_id, dim%name, dim%dim_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! ...the variable id..
    status = nf90_inq_varid( file_id, dim%name , dim%var_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! ...the dimension size...
    status = nf90_inquire_dimension( file_id, dim%dim_id, len=dim%length )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! and finally the dimension data...
    allocate( dim%values(dim%length) )
    status = nf90_get_var( file_id, dim%var_id, dim%values )
    call handle_nc_error( __FILE__, status , __LINE__ )

  end subroutine get_dim_info
  !
  !-------------------------------------------------------------------
  !
  subroutine create_dim( file_id , dim )

    ! Creates a new dimension in a file. !

    use netcdf

    ! arguments..
    integer,intent(in) :: file_id
    type(nc_dimension) :: dim

    ! local variables..
    integer :: status

    ! create dimension
    status = nf90_def_dim( file_id , dim%name , dim%length , dim%dim_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine create_dim
  !
  !-------------------------------------------------------------------
  !
  subroutine create_coord_dim( file_id , dim )

    ! calls appropriate other routines to make !
    ! a coordinate dimension.                  !

    implicit none

    ! arguments..
    integer,intent(in) :: file_id
    type(nc_dimension) :: dim

    ! local variables..

    call create_dim( file_id , dim )
    call create_var( file_id , dim%name , dim%var_id ,  (/ dim%dim_id /) )

  end subroutine create_coord_dim
  !
  !-------------------------------------------------------------------
  !
  subroutine create_nc_var( file_id , var , dim_ids )

    ! calls appropriate other routines to make !
    ! a coordinate dimension.                  !

    implicit none

    ! arguments..
    integer,intent(in) :: file_id, dim_ids(3)
    type(nc_variable) :: var

    call create_var( file_id , var%name , var%id ,dim_ids )

  end subroutine create_nc_var
  !
  !-------------------------------------------------------------------
  !
  subroutine create_var( file_id , var_name , var_id , dim_ids , from_file_id, from_var_id )

    ! Create the structure for a variable in a netcdf file. !
    ! If given another file's id and variable id it will    !
    ! copy the "units" and "long-name" attributes of that   !
    ! variable to this one.                                 !

    use netcdf

    implicit none

    ! arguments..
    integer,intent(in)          :: file_id, dim_ids(:) 
    character(len=*),intent(in) :: var_name
    integer,intent(in),optional :: from_file_id, from_var_id
    integer,intent(out)         :: var_id

    ! local variables..
    integer :: status

    ! create variable..
    status = nf90_def_var( file_id , var_name , nf90_float , dim_ids , var_id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    if ( present(from_file_id) .and. present(from_var_id) ) then
       ! copy units attribute..
       status = nf90_copy_att( from_file_id, from_var_id, 'units', file_id, var_id)
       call handle_nc_error( __FILE__ , status , __LINE__ )

       ! copy long-name attribute..
       status = nf90_copy_att( from_file_id, from_var_id, 'long_name', file_id, var_id)
       call handle_nc_error( __FILE__ , status , __LINE__ )
    endif

  end subroutine create_var
  !
  !-------------------------------------------------------------------
  !
  subroutine copy_global_atts( from_nc , to_nc )

    ! Copy the global attributes from one file to another. !

    use netcdf

    implicit none

    ! arguments..
    type(nc_header),intent(in) :: from_nc, to_nc

    ! local variables..
    integer :: i, nos_gl_atts, status
    character(len=nf90_max_name) :: name

    ! get number of global attributes..
    status = nf90_inquire( from_nc%id, nAttributes=nos_gl_atts )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    do i = 1 , nos_gl_atts

       ! get the name associated with each attribute..
       status = nf90_inq_attname( from_nc%id, nf90_global, i, name)
       call handle_nc_error( __FILE__ , status , __LINE__ )

       ! specify a copy of this attribute by name..
       status = nf90_copy_att( from_nc%id, nf90_global, trim(name), to_nc%id, nf90_global )
       call handle_nc_error( __FILE__ , status , __LINE__ )

    enddo

  end subroutine copy_global_atts
  !
  !-------------------------------------------------------------------
  !
  subroutine set_dim_pointer( dim, dimname, dimlength, dimvalues )

    ! Put the given data into the given dimension !

    use log_tools

    implicit none

    ! arguments..
    type(nc_dimension),pointer :: dim
    character(len=*),intent(in) :: dimname
    integer,intent(in) :: dimlength
    real,dimension(:),intent(in) :: dimvalues

    dim%name   = dimname
    dim%length = dimlength
    if ( associated( dim%values ) ) then
       call write_log( trim("Re-allocating "//dimname) , msg_warning , __FILE__ , __LINE__ )
       deallocate( dim%values )
    endif
    allocate( dim%values( dimlength ) )
    dim%values = -1e30
    dim%values = dimvalues

  end subroutine set_dim_pointer
  !
  !-------------------------------------------------------------------
  !
  subroutine set_var_pointer( var, varname, dim1_len , dim2_len , dim3_len )

    ! Setup the variable pointer !

    use log_tools

    implicit none

    ! arguments..
    type(nc_variable),pointer :: var
    character(len=*),intent(in) :: varname
    integer,intent(in) :: dim1_len, dim2_len, dim3_len

    var%name   = varname
    if ( associated( var%data ) ) then
       call write_log( trim("Re-allocating "//varname) , msg_warning , __FILE__ , __LINE__ )
       deallocate( var%data )
    endif
    allocate( var%data( dim1_len , dim2_len , dim3_len ) )
    var%data = -1e30

  end subroutine set_var_pointer
  !
  !-------------------------------------------------------------------
  !
  subroutine get_nc_var( file_id, var_handle )

    ! for the given var (assumes %name is already put in)
    ! it (1) checks to see if name is in file,
    ! (2) gets dim-ids for var from file
    ! (3) allocates var%data to be correct size..

    use log_tools
    use netcdf

    implicit none

    ! arguments..
    integer :: file_id
    type(nc_variable) :: var_handle

    ! local variables..
    integer                      :: i,n_dims,status,datatype
    integer,allocatable          :: dim_id_list(:), dim_length(:)
    character(len=nf90_max_name),allocatable :: dim_name_list(:)
    real,allocatable             :: data_4d(:,:,:,:)
    real(kind(0.d0)),allocatable :: dble_data_3d(:,:,:), dble_data_4d(:,:,:,:)

    ! get the variable's id..
    status = nf90_inq_varid( file_id , var_handle%name , var_handle%id )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! get the number of dim's that variable has..
    status = nf90_inquire_variable( file_id , var_handle%id , xtype=datatype , ndims=n_dims )
    call handle_nc_error(  __FILE__ , status , __LINE__ )

    ! get the id of the dimensions associated with a particular variable..
    allocate( dim_id_list( n_dims ) )
    status = nf90_inquire_variable( file_id , var_handle%id , dimids=dim_id_list )
    call handle_nc_error( __FILE__ , status , __LINE__ )

    ! and the sizes of those dims..
    allocate( dim_length( n_dims ) , dim_name_list( n_dims ) )
    do i=1,n_dims
       status = nf90_inquire_dimension( file_id , dim_id_list(i) , len=dim_length(i) , name=dim_name_list(i) )
       call handle_nc_error( __FILE__ , status , __LINE__ )
    enddo

    ! now get the data, and trim a dimension if it happens to be a 4-d object..
    select case (n_dims)
    case (3)
       call write_log("Data variable is 3-d" , msg_info , __FILE__ , __LINE__ )
       if ( associated(var_handle%data ) ) deallocate(var_handle%data)
       allocate( var_handle%data( dim_length(1) , dim_length(2) , dim_length(3) ) )
       ! check in case input-data is a double (we will Not pass along!)..
       if ( datatype .eq. nf90_double ) then
          allocate( dble_data_3d( dim_length(1) , dim_length(2) , dim_length(3) ) )
          status = nf90_get_var( file_id , var_handle%id , dble_data_3d )
          call handle_nc_error(  __FILE__ , status , __LINE__ )
          var_handle%data = real( dble_data_3d )
       else
          status = nf90_get_var( file_id , var_handle%id , var_handle%data )
          call handle_nc_error(  __FILE__ , status , __LINE__ )
       endif
    case (4)
       write(message,*)"Data variable is 4-d - need to remove a dimension!"
       call write_log( message , msg_info , __FILE__ , __LINE__ )
       allocate( data_4d( dim_length(1) , dim_length(2) , dim_length(3) , dim_length(4) ) )
       ! check in case input-data is a double (we will Not pass along!)..
       if ( datatype .eq. nf90_double ) then
          allocate( dble_data_4d( dim_length(1) , dim_length(2) , dim_length(3) , dim_length(4) ) )
          status = nf90_get_var( file_id , var_handle%id , dble_data_4d )
          call handle_nc_error( __FILE__ ,  status , __LINE__ )
          data_4d = real( dble_data_4d , kind(data_4d) )
       else
          status = nf90_get_var( file_id , var_handle%id , data_4d )
          call handle_nc_error( __FILE__ ,  status , __LINE__ )
       endif
       ! now check which, if any, of the dimensions can be concatenated (ie are of size 1)
       if ( dim_length(1) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(1))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
          if ( associated(var_handle%data) ) deallocate(var_handle%data)
          allocate( var_handle%data( dim_length(2), dim_length(3), dim_length(4) ) )
          var_handle%data = data_4d( 1 , 1:dim_length(2) , 1:dim_length(3) , 1:dim_length(4) )
       else if ( dim_length(2) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(2))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
          if ( associated(var_handle%data) ) deallocate(var_handle%data)
          allocate( var_handle%data( dim_length(1), dim_length(3), dim_length(4) ) )
          var_handle%data = data_4d( 1:dim_length(1) , 1 , 1:dim_length(3) , 1:dim_length(4) )
       else if ( dim_length(3) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(3))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
          if ( associated(var_handle%data) ) deallocate(var_handle%data)
          allocate( var_handle%data( dim_length(1), dim_length(2), dim_length(4) ) )
          var_handle%data = data_4d( 1:dim_length(1) , 1:dim_length(2) , 1 , 1:dim_length(4) )
       else if ( dim_length(4) .eq. 1 ) then
          write(message,*)"dimension ",trim(dim_name_list(4))," is of length 1 and will be removed!"
          call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
          if ( associated(var_handle%data) ) deallocate(var_handle%data)
          allocate( var_handle%data( dim_length(1), dim_length(2), dim_length(3) ) )
          var_handle%data = data_4d( 1:dim_length(1) , 1:dim_length(2) , 1:dim_length(3) , 1 )
       endif
    case default ! n-dims case
       write(message,*) "I don't know what to do with a variable with this many dimensions:",n_dims
       call write_log( trim(message) , msg_fatal, __FILE__ , __LINE__ )
    end select

    if ( allocated(  dim_id_list  ) )  deallocate(  dim_id_list  )
    if ( allocated(   dim_length  ) )  deallocate(   dim_length  )
    if ( allocated( dim_name_list ) )  deallocate( dim_name_list )
    if ( allocated(    data_4d    ) )  deallocate(    data_4d    )

  end subroutine get_nc_var
  !
  !-------------------------------------------------------------------
  !
  subroutine put_dim_nc( file_id , dim )

    ! Write a dimension's values to a netcdf file. !

    use netcdf

    implicit none

    ! arguments..
    integer,intent(in)         :: file_id   ! handle to nc-file to write to
    type(nc_dimension),pointer :: dim

    ! local variables..
    integer :: status

    status = nf90_put_var( file_id , dim%var_id , dim%values )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  end subroutine put_dim_nc
  !
  !-------------------------------------------------------------------
  !
  subroutine put_var_nc( file_id , var )

    ! Write a variable's data to a netcdf file. !

    use netcdf

    implicit none

    ! arguments..
    integer,intent(in)         :: file_id   ! handle to nc-file to write to
    type(nc_variable),pointer :: var

    ! local variables..
    integer :: status

    status = nf90_put_var( file_id , var%id , var%data )
    call handle_nc_error( __FILE__ , status , __LINE__ )

  endsubroutine put_var_nc
  !
  !-------------------------------------------------------------------
  !
  subroutine handle_nc_error( file , status , line )

    ! handle netCDF error !

    use netcdf
    use log_tools

    implicit none

    ! arguments..
    character(len=*), intent(in) :: file    ! name of f90 file error occured in
    integer,          intent(in) :: status  ! netCDF return value
    integer,optional             :: line    ! line number error occured at

    ! local variables..
    integer :: linenos

    if ( .not. present(line) ) then
       linenos = -1
    else
       linenos = line
    endif

    if ( status .ne. nf90_noerr ) then
       call write_log( nf90_strerror(status) , type=msg_fatal , file=file , line=linenos )
    end if

  end subroutine handle_nc_error
  !
  !-------------------------------------------------------------------
  !
  function close_nc_output( oc )

    ! close output file, and remove element from linked list !

    use netcdf
    use log_tools

    implicit none

    ! arguments..
    type(nc_output), pointer :: close_nc_output
    type(nc_output), pointer :: oc

    ! local variables..
    integer :: status

    if ( associated(oc) ) then
       ! If there was a previous file, change that file's 'next' link
       ! to point below me..
       if (associated(oc%previous)) then
          oc%previous%next => oc%next
       end if
       ! If there is a next file..
       !    change that file's 'previous' link to point above me.
       !    Also change 'me' to point to that next file.
       ! Else there are no more files to examine.
       if (associated(oc%next)) then
          oc%next%previous => oc%previous
          close_nc_output => oc%next
       else
          close_nc_output => NULL()
       end if
       ! Now close the current file..
       status = nf90_close(oc%header%id)
       call write_log_div
       call write_log( 'Closing output file '//trim(oc%header%name) )
       ! And remove its pointers..
       deallocate(oc)
    end if

  end function close_nc_output
  !
  !----------------------------------------------------------------------------
  !
end module netcdf_tools
