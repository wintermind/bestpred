!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME:         bestpred_main.f90
! VERSION:      2.0 beta
! RELEASED:     01 AUGUST 2007
! AUTHORS:      Paul M. VanRaden (paul@aipl.arsusda.gov)
!               John B. Cole (john.cole@ars.usda.gov)
! DESCRIPTION:  Prepares data for subroutine fmt4DCR. This program is part
!               of the bestpred package from AIPL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program bestpred_main
!
!      For the love of humanity...if you declare a variable as a parameter and
!      later try to assign it a value (like when reading from the parameter
!      file NAMELIST) gfortran (the GNU Fortan compiler) will build an
!      executable that segfaults. Other compilers, including G95 (Windows)
!      and Absoft (Linux), do NOT behave this way. The wonderful thing about
!      standards is that they can be interpreted in so many ways...
!
!      integer    obs,          maxshow,   maxtd, READparms
!      parameter (obs=99999999, maxtd=50, READparms=0)
       integer :: obs = 99999999
       integer :: maxshow
       integer :: maxtd = 50
       integer :: READparms = 0
!
!     PARAMETER DEFINITIONS:
!                           source = 10 inputs test day data in format 4
!                                       as defined on aipl.arsusda.gov
!                                  = 11 uses testing plans coded by user
!                                  = 12 inputs USDA master file records
!                                  = 13 inputs test file of RIP records
!				   = 14 inputs TD data in DRMS format
!                                  = 15 as 10, but also reads 305-d means
!                                       from the file format4.means.
!                                  = 24 reads filenames for format 14
!                                       files from a text file.
!                           obs     is maximum records to read
!                           maxshow is maximum cow records to display
!                           maxtd   is maximum test days in a lactation
!                           READparms is a flag passed downstream from
!                                     bestpred to fmt4dcr and then to
!                                     aipldcr that indicates which parms
!                                     should be used by that particular
!                                     subroutine -- those passed in from
!                                     bestpred (0) or those stored in the
!                                     file bestpred.par (1).
!-----------------------------------------------------------------------
!     NOTES ON USE OF AIPLDCR:
!                           1) Change parameters above if desired.
!                           2) Set data set names in OPEN stmts below.
!                           3) To combine data from multiple records,
!                              sort format 4T or master file records by
!                              cow and fresh date. Just one record will
!                              be output for each lactation including
!                              all test days and all herds with data.
!                           4) Edit file 11 if you want to examine DCR
!                              for particular testing plans.
!----------------------------------------------------------------------
      real*8 DCRm,DCRc,DCRs,YLDvec(2,16),PERSvec(4),RELyld(4),RELpers(4) &
            ,Yvec(2,4),herd305(2,4),bump(4),ageadj,PROJact(16)

      ! What's coming out of fmt4dcr?
      ! DCRx: Data collection ratings
      ! YLDvec: BP of lactation yields
      ! PERSvec: BP of persistencies
      ! RELyld: Yield reliabilities
      ! RELpers: Persistency reliabilities
      ! Yvec:
      ! herd305: Vector of herd averages for m,f,p,s

      ! DAILYbp and DAILYherd will contain daily BP of yield for cows
      ! and herds, respectively, unless dailyfreq=0. If dailyfreq is
      ! 0 then they will contain only zeroes.
      real*8 DAILYbp(2,4,999)
      real*8 DAILYherd(4,999)

      !character*710 format4
      character*1400 format4
      character*248 f248,n248
      character*100 BLUPout(10)
      !character*18, dimension(:,:), allocatable ::  GRAFout
      character*18 ::  GRAFout(500,4)
!      character*23  segment(maxtd),seg(maxtd)
      character*23, dimension(:), allocatable :: segment, seg
      character*80  a80
      character*17  cowid,cow,cowidmeans
!      character*14  shrtseg(maxtd)
      character*14, dimension(:), allocatable :: shrtseg
      character*8   herd,nherd,fresh,nfresh,birth,nbirth,freshmeans
      character*4   fmoday
      character(len=4) practicec
      integer i,j,k,m,maxprnt,practice,nomore,nread                    &
             ,ntd,nseg,multhrd,lacno,nlacno                            &
             ,ncow,length                                              &
             ,BLUPn(10),herdys,skip,plan(9),dims,fyr                   &
             ,herd14(3),avgage,lnwt,me305(3),dev305(3)                 &
             ,dimm,milk,fatpc,propc,scs,freq, maxlen, breed11          &
             ,breedUNK
!      integer iyld(4,maxtd), dim(maxtd),Xmilk(maxtd),weigh(maxtd),     &
!             , sample(maxtd), MRD(maxtd), super(maxtd), status(maxtd)  &
!             , pctship(maxtd)
      integer, dimension(:,:), allocatable :: iyld
      integer, dimension(:), allocatable   :: dim, Xmilk, weigh, sample &
             , MRD, super, status, pctship
      integer :: use3X = 3, mtrait = 3
      integer*2 doprev,ndoprev
      integer :: dim0(8), dim0flag

!!!   real, dimension(:,:,:), allocatable, save :: dyield

!                                                            Added by JBC
!      character*23  gSEG(maxtd)
      character*23, dimension(:), allocatable :: gSEG
      character*7   cowid7
      character*2   breed2
      integer       BPfreq, nreadalbert, MTswitch
      integer       gDIM,gMRD,gMILKSHIPPED, gNDOPREV, gNLACNO,ios,     &
                    gHERD305(4), iosseg, grfstub,                      &
                    laclen, dailyfreq, plotfreq
!      integer :: gDIMVEC(maxtd)
      integer, dimension(:), allocatable :: gDIMVEC
!       integer :: grdim, grfmilk, grffat, grfprot, grfscs
!       real*8        gTDDIM(365,4)
      character*4   gMILK
      character*2   gFAT, gPROT, gSCS
      character*1   gSUPCODE,gTDSTATUS,gMILKFREQ,gNMILKWEIGHED,        &
                    gNMILKSAMPLED
      logical ::    TDstatus(305)
!                                                  Print plots for MFPS
      integer :: GRAFplot(4) = (/1,0,0,0/) &
                 ,source = 11 ,DEBUGmsgs = 0 ,ONscreen = 0             &
                 ,GLOBALmtrait = 3                                     &
                 ,WRITEcurve = 0, WRITEdata = 0, LOGon = 0             &
                 ,LOGfreq = 0, CURVEsmall = 0, CURVEsingle = 0
      character*64 :: INfile = 'pcdart.bpi'                            &
                      ,OUTfile = 'pcdart.bpo'                          &
                      ,CURVEfile = 'cowcurve'                          &
                      ,DATAfile = 'cowdata'                            &
                      ,LOGfile = 'logfile'                             &
                      ,MEANSfile = 'format4.means'
      real*8  :: TEMPestyld(4), TEMPesthrd(4), PERSfloor, PERSceiling
      !!!integer :: TEMPestdim(4)
      integer ::  TDflag
      character(7) :: GRAFname(8)
      character :: INTmethod, INTmethodSCS, UNITSin, UNITSout
      character(len=17) cowidfmt11

      ! Used for processing source 24 files.
      integer :: ios24 = 0
      integer :: n24 = 0
      integer :: iostest = 0
      integer :: eoftest = 0
      integer :: oldsource = 0
      integer :: saved_use3X = 0
      character*64 :: fmt14file

      integer :: region=2, season=1, mnth, ndup, year

      ! breed11breed is used to map the breed11 parameter to a 2-character
      ! breed code when forming Format 4 records based on the testing plans
      ! coded in DCRexample.txt. It's used only when source == 11.
      character(2) :: breed11breed(6)

!      ! This interface is necessary so that interpolate() can accept optional arguments.
!      INTERFACE
!        SUBROUTINE interpolate(trait, month, dyield, lacn, dyld, sum, meanyld, meanp &
!          , sd, dsd, method, DEBUGmsgs, maxlen, breed, region, season)
!          character :: method
!          integer :: month, lacn, DEBUGmsgs, maxlen, trait, breed
!          real :: dyield(maxlen,4,2), dyld(12,4,2), dsd(12,4,2)
!          real*8 :: meanyld(4,maxlen,2), sd(4,maxlen,2), meanp(2,4), sum
!          integer :: region, season
!        END SUBROUTINE
!      END INTERFACE

!      character*80 :: LogMessage

      data practice       /0/
      data nomore         /0/
      data ntd            /0/
      data nseg           /1/
      data multhrd        /0/
      data lacno          /1/
      data length         /0/
      data skip           /0/
      data dev305         /3*0/
      data breed11breed   /'AY','BS','GU','HO','JE','MS'/
      data GRAFname       /'MT Milk','MT  Fat','MT Prot'                    &
                          ,'MT  SCS','ST Milk','ST  Fat','ST Prot','ST  SCS'/

!      ios = 0
!      iosseg = 0
!      grfstub = -1
!      TDflag = 0
!      Xmilk = 0

! ---------------------------------------------- Load the parameter file
      !
      ! IMPORTANT NOTE: If you are getting strange error messages
      ! when you call read_parm() check bestpred.par to make sure
      ! that it has valid <eol> characters in it. Some editors do
      ! not set them correctly, particularly on the last line of
      ! the file.
      !
      call read_parms(laclen, maxlen, dailyfreq, plotfreq                   &
        ,use3X, mtrait, GLOBALmtrait, GRAFplot, PERSfloor, PERSceiling      &
        ,source, WRITEcurve, CURVEfile, WRITEdata, DATAfile, INfile         &
        ,OUTfile, maxprnt, ONscreen, obs, maxshow, maxtd, INTmethod         &
        ,INTmethodSCS, DEBUGmsgs, UNITSin, UNITSout, breed11, breedUNK      &
        ,dim0, dim0flag, LOGon, LOGfile, LOGfreq, CURVEsmall, CURVEsingle   &
        ,region, season)

!    print *, '[bestpred_main]: maxprnt: ', maxprnt
!    print *, '[bestpred_main]: maxshow: ', maxshow
!     print *, '[bestpred_main]: source: ', source
!     print *, '[bestpred_main]: CURVEsingle = ', CURVEsingle

! Okay, now we need to allocate all of the "new" allocatable structures that
! are dimensioned using maxtd. That's done here to make sure that values from
! the parameter file are used when calling from the command line.

      allocate(segment(maxtd))
      allocate(seg(maxtd))
      allocate(shrtseg(maxtd))
      allocate(iyld(4,maxtd))
      allocate(dim(maxtd))
      allocate(Xmilk(maxtd))
      allocate(weigh(maxtd))
      allocate(sample(maxtd))
      allocate(MRD(maxtd))
      allocate(super(maxtd))
      allocate(status(maxtd))
      allocate(pctship(maxtd))
      allocate(gSEG(maxtd))
      allocate(gDIMVEC(maxtd))

      ios = 0
      iosseg = 0
      grfstub = -1
      TDflag = 0
      Xmilk = 0

!     ...................................................... INPUT FILES
      if(source == 10 .or. source == 15) OPEN (10, file='format4.dat')
      if(source == 15) OPEN (15, file='format4.means')
      if(source == 11) OPEN (11, file='DCRexample.txt')
      if(source == 12) OPEN (12, file='input.dcr')
      if(source == 13) OPEN (13, status='old', file='/data1/lori/AIPL.research.only')
! Source 14 was added for Albert de Vries.
      if ( source == 14 ) then
        if ( len(trim(INfile)) > 0 ) then
          if ( ONscreen == 1 ) print *,'Opening ', trim(INfile), ' for input'
          OPEN (14, file=trim(INfile))
        else
          if ( ONscreen == 1 ) print *,'Opening pcdart.bpi for input'
          OPEN (14, file='pcdart.bpi')
        end if
        if ( len(trim(OUTfile)) > 0 ) then
          if ( ONscreen == 1 ) print *,'Opening ', trim(OUTfile), ' for output'
          OPEN (45, file=trim(OUTfile))
        else
          if ( ONscreen == 1 ) print *,'Opening pcdart.bpo for output'
          OPEN (45, file='pcdart.bpo')
        end if
        if ( ONscreen == 1 ) print *,'Opening pcdart.fmt48 for ' &
          ,'intermediate fmt4 processing'
        OPEN (44, file='pcdart.fmt48')
!        print *,'Files opened'
      end if
! Source 24 was added so that you can give BESTPRED a file that contains a list
! of Format 14 files to process. The file of filenames MUST be named "pcdart_files.txt".
      if ( source == 24 ) then
        !if ( ONscreen == 1 )
        print *,'Opening pcdart_files.txt for input'
        OPEN (64, file='pcdart_files.txt')
        if ( len(trim(OUTfile)) > 0 ) then
          if ( ONscreen == 1 ) print *,'Opening ', trim(OUTfile), ' for output'
          OPEN (45, file=trim(OUTfile))
        else
          if ( ONscreen == 1 ) print *,'Opening pcdart.bpo for output'
          OPEN (45, file='pcdart.bpo')
        end if
        if ( ONscreen == 1 ) print *,'Opening pcdart.fmt48 for ' &
          ,'intermediate fmt4 processing'
        OPEN (44, file='pcdart.fmt48')
      end if
! Back to your regulary scheduled programming,,,
      if ( source .ne. 12 ) then
!       ................................................... OUTPUT FILES
        OPEN (20, file='results_v2.dcr')
!                      (Data Collection Rating, 305 yield,
!                      persistency, reliability)
        OPEN (21, file='lctcurve.dat')
!                      (Lactation curve for cow, contemps)
      else
!                                                         File for REMLD
        OPEN (24, file='yldpers.mfp')
!                                                     Files for MTDFREML
        OPEN (25, file='yldpers.mlk')
        OPEN (26, file='yldpers.fat')
        OPEN (27, file='yldpers.pro')
        OPEN (28, form='UNFORMATTED', file='Rinvers.mlk')
        OPEN (29, form='UNFORMATTED', file='Rinvers.fat')
        OPEN (30, form='UNFORMATTED', file='Rinvers.pro')
        OPEN (31, form='UNFORMATTED', file='Rinvers.m')
        OPEN (32, form='UNFORMATTED', file='Rinvers.f')
        OPEN (33, form='UNFORMATTED', file='Rinvers.p')
!       Note that unit 39 is used by the logging system.
      end if
      !maxprnt = maxshow

      ! This is a kludge. The bestpred subroutine gets called twice for
      ! each record when we run from the command line, once to get the
      ! actuals and once for the MEs. So, we need to double maxprnt and
      ! maxshow in order to get the desired number of plots and summaries.
      maxprnt = maxprnt * 2
      maxshow = maxshow * 2

! ------------------------------------- BEGIN READING FORMAT 4T RECORDS
      nread = 0
!      ndoprev = 0
      nreadalbert = 0

      if ( source == 11 ) go to 100
!
! The oldsource variable is used to keep track of when we're reading from a
! source 24 file.
 10   if ( oldsource == 24 ) then
        source = 24
        if ( nomore == 1 ) go to 200
        if ( nomore /= 1 ) go to 2400
      end if
      if ( nomore == 1 ) go to 200
!                                                         Fixed portion
 11   if ( source .ne. 14 .and. source .ne. 24 ) then
        if ( nread > 0 ) then
          herd305(:,:) = 0.d0
          f248 = n248
!          read(f248,12) cowid,herd,fresh
! 12       format(t3,a17,t107,a8,t128,a8,t159,i2)
! Added to support previous days open (thanks to BJH).
          read(f248,12) cowid,herd,fresh,ndoprev
 12       format(t3,a17,t107,a8,t128,a8,t159,i2,t246,i3)
          if ( source == 15 ) then
            read(15,13) cowidmeans,freshmeans,herd305(1,1),herd305(1,2),herd305(1,3),herd305(1,4)
 13         format(a17,1x,a8,1x,f5.0,1x,f4.0,1x,f4.0,2x,f3.0)
            if ( cowid .ne. cowidmeans ) then
              print *,'[ERROR]: Cow IDs do not match for Source 15 records: ', cowid, ' versus ', cowidmeans, '!'
              if ( fresh .ne. freshmeans ) print *,'[ERROR]: Fresh dates do not match for Source 15 records: ', &
                fresh, ' versus ', freshmeans, '!'
              print *, '[WARNING]: Setting means to 0!'
              herd305(1,:)= 0.d0
            end if
            !print *, '[bestpred_main]: herd305 after reading format4.means', herd305(1,:)
          end if
          doprev = ndoprev
!                                                1 segment per test day
          do 15 i=1,nseg
            ntd = ntd + 1
            if ( ntd > maxtd ) go to 15
            segment(ntd) = seg(i)
            read(seg(i),'(i3)') dims
            length = max(length,dims)
!			if ( ntd > 20 ) print *, 'segment ', ntd, ' ', segment(ntd)
 15       continue
        end if
!        print *, '[bestpred_main]: length read from format 4T:', length
        write(f248(136:138),'(i3)') length
      end if
!                                                     Read next record
 18   if ( source == 10 .or. source == 15 ) then
        read(10,20,end=50) n248,nseg,(seg(i),i=1,nseg)
 20     format(a248,i2,50a23)
!        read(n248,12) cow,nherd,nfresh
        read(n248,181) cow,nherd,nfresh
 181    format(t3,a17,t107,a8,t128,a8)
!        print *, '[bestpred_main]: n248   : ', n248
        ! This needs more thorough investigation when I have time -- I think
        ! that f248 and n248 are redundant and a single variable can be used
        ! consistently in bestpred_main.f90. This fixes Miels's problem with
        ! incomplete Format information being passed to bestpred_fmt4.f90. The
        ! values written to the format4 variable were being taken from f248 in
        ! all cases, never from n248, so the data from the input file were never
        ! being passed downstream.
        f248 = n248
      end if

!***
!if(cow /= 'HOUSA00093WPJ3927') go to 18
      if(source == 12) then
        read(12,25,end=50) cow,nherd,nfresh,nlacno, &
                            nseg,(seg(i),i=1,nseg)
 25     format(a17,a8,a8,i1,i2,50a23)
        write(n248,12) cow,nherd,nfresh,nlacno
      end if

      if(source == 13) then
        read(13,30,end=50) cow,nherd,nfresh,nbirth,nlacno, &
                           lnwt,ndoprev,herd14,avgage,     &
                           nseg,(shrtseg(i),i=1,nseg)
 30     format(a17,3a8,i1,i3,i3,i5,i4,i4,i2,i2,50a14)
!                                                  Adjust herd average
        ageadj = 1.0 + .0056*(72 - avgage)
        if(avgage > 72) ageadj = 1.0
        do 35 j=1,3
 35       herd14(j) = herd14(j)*ageadj
!                                                       End adjustment
        do 37 j=1,nseg
          read(shrtseg(j),36) dimm,milk,fatpc,propc,scs,freq
 36       format(i3,i4,3i2,i1)
 37       write(seg(j),38) dimm,'11',freq,freq,freq,'01100', &
                           milk,fatpc,propc,scs
 38       format(i3,a,3i1,a,i4,3i2)
        write(n248,40) cow,nbirth,nherd,nfresh,nlacno,herd14,dev305
 40     format(t3,a17,t71,a8,t107,a8,t128,a8,t159,i2, &
               t188,i5,2i4,i6,2i5)
      end if

! If the source is 24 then we need to loop over the contents of the file and
! handle each row as a separate source 14 file. IF we've read all of the
! records from a source 14 file then it's time to read the next filename.
         !print *, 'Handling source == 24...'

! I need to handle at least three distinctly different cases:
!   1. This is the first time through the loop and we need to get the first filename.
!   2. This is not the first time through the loop, but we haven't finished reading
!      all of the records from the Format 14 file.
!   3. This is not the first time through the loop, but we've finished reading the
!      Format 14 file and we need to read the next file name.
!
!   1. First time through the loop and we need to get the first filename.
!      I'm using an input unit of 64 because 24 would collide with some
!      old code of PVR's, and we don't want to leave that kind of thing
!      hanging around to blow up on us later.
2400  if ( source == 24 .and. oldsource == 0 ) then
        read(64, 2490, IOSTAT=ios24) fmt14file
2490    format(a64)
        n24 = n24 + 1
        if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]: Read filename number ', n24, ' (', trim(fmt14file), ') ', &
          ' from the file pcdart_files.txt (IOSTAT=', ios24, ')'
        ! Hopefully, we can open the format 14 file here and let the handle fall
        ! through to the code beginning at 1400.
        if ( len(trim(fmt14file)) > 0 ) then
          if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]: Opening the file ', trim(fmt14file), &
            ' for input after reading a format 24 file'
          OPEN (14, file=trim(fmt14file))
        else
          if ( DEBUGmsgs > 0 ) print *, '[bestpred_main]: Empty line read from the format 24!'
        end if
        oldsource = 24
        source = 14
        nreadalbert = 0
        goto 1400
!   2. This is not the first time through the loop, but we haven't finished reading
!      all of the records from the Format 14 file.
      else if ( source == 24 .and. oldsource == 24 .and. ios == 0 ) then
        ! Check to see if we're at the end of the file (UNIT 14 is correct, we're checking
        ! the testday data file)
        read(14, "(A)", IOSTAT=iostest) eoftest
        ! Move the pointer back to the start of the line
        if ( iostest /= 0 ) then
          ios = iostest
          !nomore = 1
          go to 10
        else
          backspace(14)
          if ( DEBUGmsgs > 1 ) then
            print *, "[bestpred_main]: This is not the first time through the loop, but ", &
              "we haven't finished reading all of the records from the Format 14 file."
            print *, '                 nread     : ', nread
            print *, '                 source    : ', source
            print *, '                 oldsource : ', oldsource
            print *, '                 ios       : ', ios
            print *, '                 ios24     : ', ios24
          end if
          source = 14
          goto 1400
        end if
!   3. This is not the first time through the loop, but we've finished reading the
!      Format 14 file and we need to read the next file name.
      else if ( source == 24 .and. oldsource == 24 .and. ios /= 0 ) then
        ! Try and get the name of the next file to read.
        read(64, 2490, IOSTAT=ios24) fmt14file
        ! If there are no more names then we're done!
        if ( ios24 /= 0 ) then
           if ( DEBUGmsgs == 1 ) print *, '[bestpred_main]: There are no more source 14 filenames to read!'
           !source = 14
           !nreadalbert = 0
           !go to 1400
           nomore = 1
           if ( DEBUGmsgs > 1 ) then
               print *, "[bestpred_main]: This is not the first time through the loop, but ", &
                      "we have finished reading all of the records from the Format 14 file."
               print *, '                 nread     : ', nread
               print *, '                 source    : ', source
               print *, '                 oldsource : ', oldsource
               print *, '                 ios       : ', ios
               print *, '                 ios24     : ', ios24
           end if
           go to 10
        else
          n24 = n24 + 1
          if ( DEBUGmsgs == 1 ) print *,'[bestpred_main]: Read filename number ', n24, ' (', trim(fmt14file), ') ', &
            ' from the file pcdart_files.txt (IOSTAT=', ios24, ')'
          ! Open the next file of source 14 records using the filename we just read
          ! from unit 64.
          if ( len(trim(fmt14file)) > 0 ) then
            if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]: Opening the file ', trim(fmt14file), &
              ' for input after reading a format 24 file'
            OPEN (14, file=trim(fmt14file))
          else
            if ( DEBUGmsgs > 0 ) print *, '[bestpred_main]: Empty line read from the format 24 file!'
          end if
          ! Clean up our flags.
          source = 14
          nreadalbert = 0
          goto 1400
        end if
!   4. I don't think we'll ever get here, but if we do this tells me where to start looking
!      the source of the problem.
!      else
!        if ( DEBUGmsgs == 1 ) print *, '[bestpred_main]: An unexpected case was found ', &
!          'while processing source 24 records!'
!        stop
      end if

! The first line of the format 14 file contains herd average milk, fat, protein
! and SCS. Following lines contain individual cow records.
!                                                      Header Portion
 1400     if ( source == 14 .and. nreadalbert == 0 ) then
            if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]: Reading Albert header'
            read(14,1401,IOSTAT=ios) herd,gHERD305(1),gHERD305(2),  &
                                     gHERD305(3),gHERD305(4)
 1401       format(a8,i5,i4,i4,i2)
            do 1497 i=1,4
 1497         herd305(1,i) = gHERD305(i)
            nreadalbert = 1
          end if
 !                                                Detail & TD Segments
          !print *, '[bestpred_main]: source = ', source
          if ( source == 14 .and. nreadalbert == 1 ) then
            if ( ios /= 0 ) then
              nomore = 1
              go to 10
            end if

            ntd = 0
            if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]: Reading Albert detail and TD segments'
            read(14,1402,IOSTAT=ios) cowid7,breed2,nbirth,gNDOPREV,   &
                                     gNLACNO,nfresh,MTswitch,         &
                                     BPfreq,nseg,(gSEG(i),i=1,nseg)
 1402       format(t9,a7,a2,a8,i3,i2,a8,i1,i3,i2,22a21)

            if(gNDOPREV > 99) gNDOPREV = 99
            if(gNLACNO > 9) gNLACNO = 9
            if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]:   Read record for cow ', cowid7
            if ( DEBUGmsgs > 0 ) print *,'[bestpred_main]:     Reading ', nseg, ' segments'
 !                                     Form TD segments that conform to Fmt 4.
            gMILKSHIPPED = 0
            gDIMVEC = 0
            TDstatus = .false.
            do 1498 i=1,nseg
              read(gSEG(i),1403,IOSTAT=iosseg) gDIM,gSUPCODE,gTDSTATUS, &
                                 gMILKFREQ,gNMILKWEIGHED,gNMILKSAMPLED, &
                                 gMRD,gMILK,gFAT,gPROT,gSCS
 1403         format(i4,t5,a1,t6,a1,t7,a1,t8,a1,t9,a1,t10,i2,t12,a4,    &
                t16,a2,t18,a2,t20,a2)
!               print *,'      ios: ', ios, ' iosseg: ', iosseg
              if ( ios /= 0 ) go to 1406
!               if ( iosseg /= 0 ) go to 1499
              if ( iosseg /= 0 ) go to 1406
              if ( gDIM > 999 ) gDIM = 999
              gDIMVEC(i) = gDIM
              if ( gDIM <= 305 ) TDstatus(gDIM) = .true.
                write(seg(i),1404) gDIM,gSUPCODE,gTDSTATUS,gMILKFREQ,     &
                                 gNMILKWEIGHED,gNMILKSAMPLED,gMRD,        &
                                 gMILKSHIPPED,gMILK,gFAT,gPROT,gSCS
 1404         format(i3,a1,a1,a1,a1,a1,i2,i3,a4,a2,a2,a2)
 !             print *,'      Segment ', i, ' : ', seg(i)
              ntd = ntd + 1
              if(ntd > maxtd) go to 1499
              segment(ntd) = seg(i)
              read(seg(i),'(i3)') dims
 1498         length = max(length,dims)
 1499       continue
!              Now we need to form the 248-byte "header" for the fmt4 record.
 1406       write(f248,1405) breed2,cowid7,nbirth,herd,nfresh,gNLACNO
 1405       format(t3,a2,t8,a7,t71,a8,t107,a8,t128,a8,t159,i2)
            write(f248(136:138),'(i3)') length
            if ( DEBUGmsgs > 0 ) print *, '[bestpred_main]: Wrote fmt4 ', &
              'header for cow ', cowid7
      end if

      !!! 12/20/2007 JBC Commented-out so that records with
      !!!                no TD will be processed, e.g. for
      !!!                the bias study.
      nread = nread + 1
      if ( nread > obs ) go to 50
      go to 52
 50   nomore = 1
      if ( nread == 0 ) go to 200
      !!! This was used to skip the loop at 52, but we actually need
      !!! to enter that code in order to catch the last record in a
      !!! file with no TD data.
      !cow = '00000000000000000'
!                                      Check for > 1 record / lactation
!  52   if ( cowid == cow .and. fresh == nfresh ) then
 52    if ( nomore == 1 ) go to 55
       if ( cowid == cow .and. fresh == nfresh ) then
!          print *, '[bestpred_main]: >1 record/lactation'
         if ( herd == nherd ) then
           ndup = ndup + 1
         else
           multhrd = multhrd + 1
         end if
         go to 10
       end if

!                                                         Fill format 4
 55   if ( nseg == 0 ) then
!         print *, '[bestpred_main]: Filling format 4 w/no segments ', cow
        write(format4,57) f248
        57 format(a248)
      else
        write(format4,20) f248,ntd,(segment(i),i=1,ntd)
      end if
!      print *, '[bestpred_main]: format4: ', format4

!     ..................................... CALCULATE DCR AND ME YIELDS
!
!     If the source is 14, individual animal records contain the
!     MTswitch variable, which takes the same values as mtrait, and
!     serves the same function. If the GLOBALmtrait parameter is
!     set to a valid value of mtrait (1, 3, or 4) then GLOBALmtrait's
!     value is used for all records read from a file, overriding the
!     individual switches.
      if ( source == 14 .or. source == 24 ) then
        if ( GLOBALmtrait > 0 ) then
          if ( DEBUGmsgs > 0 ) print *,'Setting mtrait = GLOBALmtrait'
          mtrait = GLOBALmtrait
        else
          if ( DEBUGmsgs > 0 ) print *,'Setting mtrait = MTswitch'
          mtrait = MTswitch
        end if
        ! We need to call fmt4dcr twice -- once to get the
        ! projected/actual yields and store them in PROJact
        ! and once to get the ME yields. Note that the
        ! maxprnt and ONscreen parameters are hard-coded to
        ! 0 for the call to get projected/actual yields.
        saved_use3X = use3X
        use3X = 0
!        if ( DEBUGmsgs > 1 ) print *, 'herd305', herd305
!        if ( DEBUGmsgs > 0 ) print *, 'Calling fmt4dcr() to get projected/actual 305-d yields'
        !call bestpred_fmt4(format4,doprev,0                    &
        !        ,DCRm,DCRc,DCRs,YLDvec,PERSvec,RELyld,RELpers  &
        !        ,Yvec,herd305,bump,BLUPout,BLUPn               &
        !        ,BPfreq,GRAFplot,DEBUGmsgs,0                   &
        !        ,mtrait,use3X,laclen,dailyfreq                 &
        !        ,INTmethod,maxlen,DAILYbp,DAILYherd            &
        !        ,WRITEcurve,CURVEfile,WRITEdata,DATAfile       &
        !        ,plotfreq, INTmethodSCS, READparms             &
        !        ,UNITSin, UNITSout, breedUNK, dim0, dim0flag   &
        !        ,LOGon,LOGfile,LOGfreq,maxshow,CURVEsmall      &
        !        ,CURVEsingle, region, season                   &
        !)

        ! I think that some other things also need to be hard-coded to 0
        ! in order to avoid records being written to output files more
        ! once, such as WRITEcurve and WRITEdata.
        call bestpred_fmt4(format4,doprev,0                    &
                ,DCRm,DCRc,DCRs,YLDvec,PERSvec,RELyld,RELpers  &
                ,Yvec,herd305,bump,BLUPout,BLUPn               &
                ,BPfreq,GRAFplot,DEBUGmsgs,0                   &
                ,mtrait,use3X,laclen,dailyfreq                 &
                ,INTmethod,maxlen,DAILYbp,DAILYherd            &
                ,0,CURVEfile,0,DATAfile                        &
                ,plotfreq, INTmethodSCS, READparms             &
                ,UNITSin, UNITSout, breedUNK, dim0, dim0flag   &
                ,LOGon,LOGfile,LOGfreq,0,0                     &
                ,0, region, season                             &
        )

        PROJact = YLDvec(2,:)
        !print *, 'PROJact: ', PROJact
        use3X = saved_use3X
      end if

!      if ( DEBUGmsgs > 1 ) print *, '================================================================================'
      if ( DEBUGmsgs > 1 ) print *, '[bestpred_main]: herd305 before bestpred_fmt4()', herd305
!      print *, '[bestpred_fmt4]: Entering bestpred_fmt4()'
!      print *, 'Calling bestpred_fmt4 for cow ', cowid7
      call bestpred_fmt4(format4,doprev,maxprnt                 &
                  ,DCRm,DCRc,DCRs,YLDvec,PERSvec,RELyld,RELpers &
                  ,Yvec,herd305,bump,BLUPout,BLUPn              &
                  ,BPfreq,GRAFplot,DEBUGmsgs,ONscreen           &
                  ,mtrait,use3X,laclen,dailyfreq                &
                  ,INTmethod,maxlen,DAILYbp,DAILYherd           &
                  ,WRITEcurve,CURVEfile,WRITEdata,DATAfile      &
                  ,plotfreq, INTmethodSCS, READparms            &
                  ,UNITSin, UNITSout, breedUNK, dim0, dim0flag  &
                  ,LOGon,LOGfile,LOGfreq,maxshow, CURVEsmall    &
                  ,CURVEsingle, region, season                  &
                  )
!	  print *, '[bestpred_fmt4]: Leaving bestpred_fmt4()'
      if ( DEBUGmsgs > 1 ) print *, '[bestpred_main]: herd305 after bestpred_fmt4()', herd305
!      print *, '[bestpred_main]: herd305 after bestpred_fmt4()', herd305
      skip = 1
!
!     ........................................ OUTPUT DCR AND ME YIELDS
!                                       DCRc is mean of fat,protein DCR
!      print *, '[bestpred_main]: I am not 12'
      if(source == 12) go to 65
!***  if(source < 12) then
!                             Write Albert-formatted records to unit 45
!      print *, '[bestpred_main]: I am not 14'
      if ( source == 14 .or. source == 24 ) then
        ! Persistencies are ~N(0,1), so we can safely set any values
        ! smaller than -9.99 or larger than +9.99 to those values as
        ! a floor and ceiling, respectively, stored in the PERSfloor
        ! and PERSceiling parameters.
        !do i = 1,8
        do i = 1,4
          if ( PERSvec(i) < PERSfloor ) then
            if ( ONscreen == 1 ) print *,GRAFname(i),' persistency for '  &
              & ,'cow ',cowid7,' was changed from ',PERSvec(i),' to a '   &
              & ,'floor of ',PERSfloor
            PERSvec(i) = PERSfloor
          end if
          if ( PERSvec(i) > PERSceiling ) then
            if ( ONscreen == 1 ) print *,GRAFname(i),' persistency for '  &
              & ,'cow ',cowid7,' was changed from ',PERSvec(i),' to a '   &
              & ,'ceiling of ',PERSceiling
            PERSvec(i) = PERSceiling
          end if
        end do
!        print *,'cowid   DIM LTDm    LTDf   LTDp   LTDs'
        do i = 1,305
          if ( i <= 305 ) then
            !do j = 1,4
            !    read(GRAFout(i,j),'(i3,2f6.1)') TEMPestdim(j), &
            !        TEMPestyld(j), TEMPesthrd(j)
            !end do
            ! Set the test day status flag
            TDflag = 0
            if ( TDstatus(i) ) TDflag = 1
            ! Write the record to the OUTfile.
            write(45,"(a7,1x,i3,1x,i1,1x,f5.1,1x,f3.1,1x,f3.1,1x,f3.1,1x,      &
                  3(i5,1x,i4,1x,i4,1x,i4,1x),                                  &
                  3(i3,1x),4(f5.2,1x),8(i2,1x),f5.1,1x,f3.1,1x,f3.1,1x,f3.1)") &
              cowid7, i, TDflag                                                &
              , TEMPestyld(1), TEMPestyld(2), TEMPestyld(3), TEMPestyld(4)                         &
              , int(YLDvec(1,9)), int(YLDvec(1,10)), int(YLDvec(1,11)), int(YLDvec(1,12))          &
              , int(PROJact(1)), int(PROJact(2)), int(PROJact(3)), int(PROJact(4))                 &
              , int(YLDvec(1,1)), int(YLDvec(1,2)), int(YLDvec(1,3)), int(YLDvec(1,4))             &
              , int(DCRm), int(DCRc), int(DCRs)                                                    &
              , PERSvec(1), PERSvec(2), PERSvec(3), PERSvec(4)                                     &
              , int(RELyld(1)*100), int(RELyld(2)*100), int(RELyld(3)*100), int(RELyld(4)*100)     &
              , int(RELpers(1)*100), int(RELpers(2)*100), int(RELpers(3)*100), int(RELpers(4)*100) &
              , TEMPesthrd(1), TEMPesthrd(2), TEMPesthrd(3), TEMPesthrd(4)
            if ( TDstatus(i) ) then
              write (*,'(1x,a7,1x,i3,1x,f7.1,1x,f6.1,1x,f6.1,1x,f5.2)')  &
                cowid7, i ,YLDvec(1,9), YLDvec(1,10), YLDvec(1,11),      &
                YLDvec(1,12)
            end if
          end if
        end do
        ! If we're reading source 24 and there are still records in the current
        ! format 14 file then we need to go back up.
        if ( DEBUGmsgs > 1 ) then
            print *, '[bestpred_main]: Finished writing results, looping back up to process more records.'
            print *, '                 nread     : ', nread
            print *, '                 source    : ', source
            print *, '                 oldsource : ', oldsource
            print *, '                 ios       : ', ios
            print *, '                 ios24     : ', ios24
        end if
        if ( source == 14 .and. oldsource == 24 .and. nomore == 0 ) goto 10
      end if

      if ( source < 99 ) then
!		print *, '[bestpred_main]: I am < 99'
!        if ( DEBUGmsgs > 0 ) print *, 'herd305', herd305
        write(20,60) cowid,fresh,length &
          ,DCRm,DCRc,DCRs,YLDvec(1,:),PERSvec,RELyld,RELpers,Yvec(1,:),herd305(1,:),bump
 60     format(a17,1x,a8,i4,1x,3f5.0,4(3f7.0,f7.2),4f6.2,8f5.2   &
          ,2(3f7.0,f7.2),1x,3f7.4,f9.2)
      else
        write(20,61) cowid,fresh,lacno,length,DCRm,DCRc,DCRs,YLDvec(1,:)
 61     format(a17,a8,i1,i3,3f4.0,4(f6.0,2f5.0,f4.2))
      end if
!
!     ........................................... OUTPUT FILES FOR BLUP
!      print *, '[bestpred_main]: Entering OUTPUT FILES FOR BLUP section'
 65   if(skip == 1) go to 91
      if(source == 10) ncow = nread - multhrd - ndup - 1
      if(source == 11) ncow = practice
      if(source == 12) read(cowid,'(9x,i8)') ncow
      if(source .ne. 12) go to 99
!      print *, '[bestpred_main]: Reading fresh date'
      read(fresh,'(2x,i2,i2)') year,mnth
!	  print *, '[bestpred_main]: Reading herd'
      read(herd,'(i8)') herdys
      herdys = herdys*10 + (year - 90)
!     if(nread<maxshow) print *,'herd = ',herd,'  herdys = ',herdys
!                                                  Data for sire model
      write(BLUPout(1),70) ncow,herdys,mnth,lacno,Yvec(1,:)
 70   format(i8,i12,2i3,4f7.2)
!                                                Data for animal model
!!!   This may not work as intended -- Yvec has been resized so that it
!!!   is only 4 elements long since only ST or MT calculations are done
!!!   in a single call, instead of the original code in which both were
!!!   done.
      write(BLUPout(2),72) ncow,herdys,mnth,lacno,Yvec(1,1),Yvec(1,4)
      write(BLUPout(3),72) ncow,herdys,mnth,lacno,Yvec(1,2),Yvec(1,1)
      write(BLUPout(4),72) ncow,herdys,mnth,lacno,Yvec(1,3),Yvec(1,2)
 72   format(i8,i12,2i3,3(f7.0,f7.2))
!                                                        Record length
      do 80 k=1,4
        if(k == 1) BLUPn(k) = 82
 80     if(k >= 2) BLUPn(k) = 40
!                                                     Write BLUP files
       do 90 k=1,10
         if(k <= 4) write(k+23,'(a)') BLUPout(k)(1:BLUPn(k))
 90      if(k >= 5) write(k+23) BLUPout(k)(1:BLUPn(k))
!      ............................................ Output lact curves
 91    if ( source == 14 ) then
!       if ( source == 14 ) then
         write(21,921) 0,cowid7,nfresh,nlacno
 921     format(i3,1x,a7,1x,a8,i3)
       else
         write(21,922) 0,cowid,fresh,lacno
 922     format(i3,1x,a17,1x,a8,i3)
       end if

 99    ntd = 0
       length = 0
       if(nread > 0) then
         !print *, '[bestpred_main]: I am jumping to 10 after processing cow ', cowid
         go to 10
       end if
! ................................................ PROCESS EXAMPLE DATA
 100     i = 0
!                                                          Herd average
      me305(1) = 20000
      me305(2) = 700
      me305(3) = 600
!                                                     Read testing plan
      if ( practice == 0 ) print *,'Reading example data from DCRexample.txt'
 110  read(11,111,end=200) a80
 111  format(a80)
      if ( a80(1:1) == '_' ) go to 110
      do 120 m=1,80
        if ( a80(m:m) > ' ' ) go to 130
 120  continue
      practice = practice + 1
      go to 55
 130  read(a80,*) (plan(j),j=1,9)
!                                                       Construct segments
      !!!dim = 0
      do 150 k=plan(7),plan(8),plan(6)
        i = i + 1
        dim(i) = k
        super(i) = plan(1)
        Xmilk(i) = plan(2)
        weigh(i) = plan(3)
        sample(i) = plan(4)
        MRD(i) = plan(5)
        pctship(i) = 100
        status(i) = 0
!                                                        Simulate data
        iyld(1,i) = 15.d0 + (me305(1)*.003 - .12*(dim(i) - 150) + &
                    15.d0*(.7 - (.7/MRD(i))*dsin(k/11.d0)))*10
        iyld(2,i) = (3.6 + .4*dsin((11 - k)/11.d0))*10
        iyld(3,i) = (3.2 + .3*dsin((11 - k)/11.d0))*10
        iyld(4,i) = (3.3 + 2.*dsin((60 - k)/60.d0))*10
        do 140 m=2,4
 140      if(sample(i) == 0) iyld(m,i) = 0
!                                              Write test day segments
!
          write(segment(i),145) dim(i),super(i),status(i),Xmilk(i)   &
            ,weigh(i),sample(i),MRD(i),pctship(i) &
            ,(iyld(m,i),m=1,4)
 145      format(i3,5i1,i2,i3,i4,3i2)
        if ( DEBUGmsgs > 1 ) then
          print *, segment(i)
        end if
        length = max(length,dim(i))
!        print *, '[bestpred_main]: length: ', length
 150    continue
      ntd = i
!                                              Construct fixed portion
!      cowid = 'HOUSA.EX.COW.'
      cowid = breed11breed(breed11) // 'USA.EX.COW.'
      herd = '12345678'
      lacno = plan(9)
      fresh = '19990401'
!     To convert from character to integer for uniquely numbering
!     output files in aipldcr.f90.
      cowidfmt11 = cowid
      practicec = cowidfmt11(14:17)
      write(practicec,'(i4.4)') practice

!                                                  Calculate birth day
      read(fresh,'(i4,a4)') fyr,fmoday
      write(birth,'(i4,a4)') fyr - lacno - 1,fmoday
!                                                  Write fixed portion
      !
      write(f248,40) cowidfmt11,birth,herd,fresh,lacno,me305,dev305
      write(f248(136:138),'(i3)') length
!      print *, '[bestpred_main]: f248(136:138): ', f248(136:138)
      doprev = 200
      go to 110

 200  if ( source == 10 ) then
        if ( ONscreen == 1 ) print *,nread,' records read from format 4T ' &
          ,'input file, unit 10'
      else if (source == 12) then
        if ( ONscreen == 1 ) print *,nread,' records read from USDA ' &
          ,'master data file, unit 12'
      else if(source == 13) then
        if ( ONscreen == 1 ) print *,nread,' records read from USDA ' &
          ,'master data file, unit 13'
      else if(source == 11) then
        if ( ONscreen == 1 ) print *,practice,' examples input from ' &
          ,'unit 11'
      else if(source == 14) then
        if ( ONscreen == 1 ) print *,nread-1,' lactation records read ' &
          ,'from pcdart.bpi'
      else if(source == 24) then
        if ( ONscreen == 1 ) print *, nread-1, &
          ' lactation records read from the ', n24, ' files listed in pcdart_files.txt'
      else
        if ( ONscreen == 1 ) print *,multhrd,' records for cows that ' &
          ,'changed herds'
        if ( ONscreen == 1 ) print *,ndup,' duplicate records for same ' &
          ,'herd-cow-lactation'
      end if

stop
end program bestpred_main
