!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME:         bestpred_parm.f90
! VERSION:      2.0 beta
! RELEASED:     01 AUGUST 2007
! AUTHORS:      John B. Cole (john.cole@ars.usda.gov)
! DESCRIPTION:  Reads program parameters from the file bestpred.par.
!               This program is part of the bestpred package from AIPL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_parms(laclen, maxlen, dailyfreq, plotfreq             &
      , use3X, mtrait, GLOBALmtrait, GRAFplot, PERSfloor, PERSceiling &
      , source, WRITEcurve, CURVEfile, WRITEdata, DATAfile, INfile    &
      , OUTfile, maxprnt, ONscreen, obs, maxshow, maxtd, INTmethod    &
      , INTmethodSCS, DEBUGmsgs, UNITSin, UNITSout, breed11, breedUNK &
      , dim0, dim0flag, LOGon, LOGfile, LOGfreq, CURVEsmall           &
      , CURVEsingle, region, season)
    ! Dummy variables
    integer :: GRAFplot(4),dim0(8),source,DEBUGmsgs,maxprnt   &
      ,ONscreen,mtrait,use3X,GLOBALmtrait                     &
      ,obs,maxshow,laclen,dailyfreq,maxlen,WRITEcurve         &
      ,WRITEdata,plotfreq,maxtd, breed11, breedUNK, dim0flag  &
      ,LOGon, LOGfreq, CURVEsmall, CURVEsingle, region, season
    character :: INTmethod, INTmethodSCS, UNITSin, UNITSout
    real*8 :: PERSfloor, PERSceiling
    character*64 :: INfile, OUTfile, CURVEfile, DATAfile, LOGfile

    ! Subroutine Variables
    integer :: DEBUGparms, j
    namelist /bestpred/ laclen, maxlen, dailyfreq, plotfreq           &
      , use3X, mtrait, GLOBALmtrait, GRAFplot, PERSfloor, PERSceiling &
      , source, WRITEcurve, CURVEfile, WRITEdata, DATAfile, INfile    &
      , OUTfile, maxprnt, ONscreen, obs, maxshow, maxtd, INTmethod    &
      , INTmethodSCS, DEBUGparms, DEBUGmsgs, UNITSin, UNITSout        &
      , breed11, breedUNK, dim0, dim0flag, LOGon, LOGfile, LOGfreq    &
      , CURVEsmall, CURVEsingle, region, season

    !
    ! IMPORTANT NOTE: If you are getting strange error messages
    ! when you call read_parm() check bestpred.par to make sure
    ! that it has valid <eol> characters in it. Some editors do
    ! not set them correctly, particularly on the last line of
    ! the file.
    !
    OPEN(50, file='bestpred.par', access='sequential', form='formatted')
    READ(50, NML=bestpred)
    CLOSE(50)

    ! Check the values of laclan and dailyfreq and set them to
    ! default values if they're outside the range of valid values,
    if ( laclen < 0 .or. laclen > 999 ) laclen = 305
    if ( maxlen < 0 .or. maxlen > 999 ) maxlen = 365
    if ( maxlen < laclen ) maxlen = laclen
    if ( dailyfreq < 0 .or. dailyfreq > 999 ) dailyfreq = 6
    if ( plotfreq < 0 .or. plotfreq > 999 ) plotfreq = 6
    if ( use3X < 0 .or. use3X > 3 ) use3X = 3
    if ( mtrait < 1 .or. mtrait > 4 ) mtrait = 3
    if ( GLOBALmtrait < 1 .or. GLOBALmtrait > 4 ) GLOBALmtrait = 3
    do j = 1, 4
      if ( GRAFplot(j) < 0 .or. GRAFplot(j) > 2 ) GRAFplot(j) = 2
    end do
    if ( PERSfloor > 0 ) PERSfloor = -9.99
    if ( PERSceiling < 0 ) PERSceiling = 9.99
    if ( source /= 10 .and. source /= 11 .and. source /= 12 .and. source /= 14 .and. source /= 15 .and. source /= 24 ) source = 11
    if ( WRITEcurve .ne. 0 .and. WRITEcurve .ne. 1 ) WRITEcurve = 0
    if ( len(CURVEfile) <= 0 ) CURVEfile = 'cowcurve'
    if ( WRITEdata .ne. 0 .and. WRITEdata .ne. 1 ) WRITEdata = 0
    if ( len(DATAfile) <= 0 ) DATAfile = 'cowdata'
    if ( len(INfile) <= 0 ) INfile = 'pcdart.bpi'
    if ( len(OUTfile) <= 0 ) OUTfile = 'pcdart.bpo'
    if ( maxprnt < 1 ) maxprnt = 1
    if ( ONscreen .ne. 0 .and. ONscreen .ne. 1 ) ONscreen = 1
    if ( obs < 0 ) obs = 99999999
    if ( maxshow < 0 .or. maxshow > obs ) maxshow = 5
    if ( maxtd < 1 ) maxtd = 50
    if ( INTmethod .ne. 'L' .and. INTmethod .ne. 'W' .and. INTmethod .ne. 'R' ) INTmethod = 'L'
    if ( INTmethodSCS .ne. 'L' .and. INTmethodSCS .ne. 'G' .and. INTmethodSCS .ne. 'S' ) INTmethodSCS = 'L'
    if ( DEBUGmsgs < 0 .or. DEBUGmsgs > 2 ) DEBUGmsgs = 0
    if ( UNITSin .ne. 'P' .and. UNITSin .ne. 'K' ) UNITSin = 'P'
    if ( UNITSout .ne. 'P' .and. UNITSout .ne. 'K' ) UNITSout = 'P'
    if ( breed11 < 1 .or. breed11 > 6 ) breed11 = 4
    if ( breedUNK < 1 .or. breedUNK > 6 ) breedUNK = 4
    if ( dim0(1) < 1 .or. dim0(1) > laclen ) dim0(1) = 115
    if ( dim0(2) < 1 .or. dim0(2) > laclen ) dim0(2) = 115
    if ( dim0(3) < 1 .or. dim0(3) > laclen ) dim0(3) = 150
    if ( dim0(4) < 1 .or. dim0(4) > laclen ) dim0(4) = 155
    if ( dim0(5) < 1 .or. dim0(5) > laclen ) dim0(5) = 161
    if ( dim0(6) < 1 .or. dim0(6) > laclen ) dim0(6) = 152
    if ( dim0(7) < 1 .or. dim0(7) > laclen ) dim0(7) = 159
    if ( dim0(8) < 1 .or. dim0(8) > laclen ) dim0(8) = 148
    if ( dim0flag < 0 .or. dim0flag > 1 ) dim0flag = 0
    if ( LOGon < 0 .or. LOGon > 1 ) LOGon = 1
    if ( len(LOGfile) <= 0 ) LOGfile = 'example'
    if ( LOGfreq <= 0 ) LOGfreq = 0
    if ( CURVEsmall .ne. 0 .and. CURVEsmall .ne. 1 ) CURVEsmall = 0
    if ( CURVEsingle .ne. 0 .and. CURVEsingle .ne. 1 ) CURVEsingle = 0
    if ( region < 1 .or. region > 7 ) region = 2       ! Midwest
    if ( season < 1 .or. season > 4 ) season = 1       ! Spring

    ! Finished with checking parms.
    if ( DEBUGparms == 1 ) then
      print *,'--------------------------------------------------', &
              '--------------------'
      print *,'Parameters read from bestpred.par:'
      print *,'================================='
      write(*,NML=bestpred)
      print *,'--------------------------------------------------', &
              '--------------------'
    end if

  return
end subroutine read_parms
