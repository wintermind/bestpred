!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME:         bestpred_log.f90
! VERSION:      2.0 beta
! RELEASED:     25 JANUARY 2008
! AUTHORS:      John B. Cole (john.cole@ars.usda.gov)
! DESCRIPTION:  Logging tools for the best prediction programs.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------
subroutine log_check_file_exists(LogFileName,FileExists)
  integer, intent(out) :: FileExists
  character*76, intent(in) :: LogFileName
  ! The STATUS='OLD' keyword specifies that the file LogFileName must
  ! exist or IOSTAT will return a code other than 0 (2 for Absoft).
  OPEN(FILE=LogFileName, UNIT=39, STATUS='OLD', ACTION='READ', IOSTAT=FileExists)
  if ( FileExists == 0 ) then
    ! The logfile already exists
    CLOSE(39)
  end if
end subroutine log_check_file_exists
!--------------------------------------------------------------------
subroutine log_create_file(LogFileName,CreatedFile)
  integer, intent(out) :: CreatedFile
  character*76, intent(in) :: LogFileName
  character*8               :: curr_date, curr_time
  character*10              :: long_time, long_date

  OPEN(FILE=LogFileName, UNIT=39, STATUS='NEW', ACTION='WRITE', IOSTAT=CreatedFile)
  if ( CreatedFile == 0 ) then
    call DATE_AND_TIME(curr_date,long_time)
    curr_time = long_time(1:2)//':'//long_time(3:4)//':'//long_time(5:6)
    long_date = curr_date(5:6)//'/'//curr_date(7:8)//'/'//curr_date(1:4)
    write (39,*) '# ', trim(LogFileName), ' created on ', long_date, ' at ', trim(curr_time)
    CLOSE(39)
  else
    print *,'[bestpred_log]: Unable to create the log file ', LogFileName
  end if
end subroutine log_create_file
!--------------------------------------------------------------------
! 1. Check to see if file exists
! 2. If file exists:
!    2a. Open the fule
!    2b. Write the event
! 3. If the file does not exist:
!    3a. Create the file
!    3b. Open the fule
!    3c. Write the event
! 4. Close the file
subroutine log_message(LogFile,LogMessage)
  character*80, intent(in)  :: LogMessage
  character*64, intent(in)  :: LogFile
  integer                   :: OpenStatus, FileExists, CreatedFile
  character*80              :: LogFileName
  character*8               :: curr_date, curr_time
  character*10              :: long_time

  ! First, form the complete file name.
  call DATE_AND_TIME(curr_date,long_time)
  curr_time = long_time(1:2)//':'//long_time(3:4)//':'//long_time(5:6)
  LogFileName = trim(LogFile)//'.'//trim(curr_date)//'.log'
  ! 1. Check to see if file exists
  call log_check_file_exists(LogFileName,FileExists)
  ! 2. If file exists:
  if ( FileExists == 0 ) then
    !    2a. Open the fule
    open(FILE=LogFileName, UNIT=39, POSITION='APPEND', STATUS='OLD', ACTION='WRITE', IOSTAT=OpenStatus)
    if ( OpenStatus == 0 ) then
      !    2b. Write the event
      write (39,*) trim('['//curr_time//']: '//trim(LogMessage))
      ! 4. Close the file
      close(39)
    end if
  ! 3. If the file does not exist:
  else
    !    3a. Create the file
    call log_create_file(LogFileName,CreatedFile)
    if ( CreatedFile == 0 ) then
      !    3b. Open the fule
      open(FILE=LogFileName, UNIT=39, POSITION='APPEND', STATUS='OLD', ACTION='WRITE', IOSTAT=OpenStatus)
      if ( OpenStatus == 0 ) then
        !    3c. Write the event
        write (39,*) trim('['//curr_time//']: '//trim(LogMessage))
        ! 4. Close the file
        close(39)
      end if
    end if
  end if
end subroutine log_message
