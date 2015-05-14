!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME:         bestpred_fmt4.f90
! VERSION:      2.0 beta
! RELEASED:     01 AUGUST 2007
! AUTHORS:      Paul M. VanRaden (paul@aipl.arsusda.gov)
!               John B. Cole (john.cole@ars.usda.gov)
! DESCRIPTION:  Prepares format 4 data for subroutine aiplDCR. This program
!               is part of the bestpred package from AIPL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine bestpred_fmt4( &
!                                                                 Input
          format4,doprev,maxprnt &
!                                                                Output
          ,DCRm,DCRc,DCRs,YLDvec,PERSvec,RELyld,RELpers &
!                                                              Optional
          ,Yvec,herd305,bump,BLUPout,BLUPn              &
          ,bpfreqin,GRAFplot,DEBUGmsgs,ONscreen         &
          ,mtrait,use3X                                 &
          ,laclen,dailyfreq,INTmethod,maxlen            &
          ,DAILYbp,DAILYherd,WRITEcurve,CURVEfile       &
          ,WRITEdata,DATAfile,plotfreq, INTmethodSCS    &
          ,READparms,UNITSin, UNITSout, breedUNK, dim0  &
          ,dim0flag, LOGon, LOGfile, LOGfreq,maxshow    &
          ,CURVEsmall, CURVEsingle, region, season      &
          )
! .....................................................................
!   Definitions (MT=multi-trait,ST=single-trait,m=milk,f=fat,p=protein)
!                                                                 Input
!       format4 = AIPL format 4 record (710 bytes) as defined
!                 on aipl.arsusda.gov web site
!       doprev  = days not pregnant in previous lactation
!       maxprnt = number of graphs to display (must not be negative)
!               = -1 displays only the input variables
!                                                                Output
!       DCRm   = MT m DCR, DCRc = MT f,p DCR,  DCRs = MT scs DCR
!       YLDvec = MT m,f,p,scs, ST m,f,p,scs, non-305 m,f,p,scs
!       PERSvec= MT m,f,p,scs persistency, ST m,f,p,scs persistency
!       RELyld = MT m,f,p,scs Reliability, ST m,f,p,scs Reliability
!       RELpers= MT m,f,p,scs Rel(persist),ST m,f,p,scs Rel(persist)
!       length = days from fresh date to last test or cull date
!                                                              Optional
!       GRAFplot = Array of switches toggling lactation curve graphs
!                  on (1) or off (0); the order of traits is: M,F,P,S.
!       DEBUGmsgs = Used to display debugging messages.
!       ONscreen = toggle output to the screen on (1) and off (0).
!       laclen specifies default lactation length
!       dailyfreq specifies how often actual dailies are calculated
!       INTmethod specifies the inerpolation method to be used for
!           calculating lactation curves.
!       READparms is a flag passed downstream from maindcr to fmt4dcr
!           and then to aipldcr that indicates which parms should be
!           used by that particular subroutine -- those passed in from
!           maindcr (0) or those stored in the file maindcr.par (1).
! .....................................................................
      integer    mtrait,  use3X,  maxtd, READparms
      parameter (maxtd=50)
!                                            aiplDCR control variables:
!                         mtrait = 1 does only ST to save CPU time
!                                = 3 does MT m,f,p; 4 does MT m,f,p,scs
!                          use3X = 0 doesn't adjust for 3X milking,
!                                  1 uses old factors, 2 new factors,
!                                  3 uses phased-in factors over time
      integer*2 agebase
!                        agebase = 0 uses mature equivalent age base
!                                =36 uses 36 month age base, etc.
      integer   showage
!                        showage = 1 displays age factors at startup
!                                = 0 doesn't
! ....................................................................
      !character*710 :: format4
      character*1400 :: format4
      character*100 :: BLUPout(10)
      character*23  :: segment(maxtd)
      character*17  :: cowid
      character*8   :: fresh,herd,birth
      character     ::  brd
      character(len=*)  :: CURVEfile, DATAfile, LOGfile
      real*8 DCRm,DCRc,DCRs,YLDvec(2,16),PERSvec(4),RELyld(4),RELpers(4), &
             yield(2,4,maxtd),Yvec(2,4),bump(4),agefac(4),avgfac(4),        &
             herd305(2,4),me305(3),dev305(3),zero
      real*8 oldfac(4)
      integer ntd,parity,maxprnt,length,BLUPn(10),                      &
              dim(maxtd),Xmilk(maxtd),weigh(maxtd),sample(maxtd),       &
              MRD(maxtd),super(maxtd),status(maxtd),                    &
              adjscs4,i,WRITEcurve,WRITEdata,LOGon,LOGfreq,maxshow
      integer   count, CURVEsmall, CURVEsingle
      integer*2 pday,age,doprev,docur,lacno,frmo,fryr,state         &
               ,fill2,yld2(3),adjyld2(3),i305,scs305
!                                                 Simple test day means
      real :: ybar(4)
!                               User-specified multiple trait parameter
      integer bpfreqin, GRAFplot(4), DEBUGmsgs, ONscreen, laclen        &
        , dailyfreq, maxlen, plotfreq, breedUNK
      character INTmethod, INTmethodSCS, UNITSin, UNITSout

      ! DAILYbp and DAILYherd will contain daily BP of yield for cows
      ! and herds, respectively, unless dailyfreq=0. If dailyfreq is
      ! 0 then they will contain only zeroes.
      real*8 DAILYbp(4,999)
      real*8 DAILYherd(4,999)
      character*2 :: herdstate

      ! integer :: breed11
      integer :: dim0(8), dim0flag, j, k
      integer :: region, season

      external pday

      data ybar    /70.,2.5,2.2,3.2/
      data count   /0/
      data docur   /90/
      data fill2   /0/
      data i305    /305/
      data agebase /0/
      data showage /0/
!      data showage /1/

!      print *, '[bestpred_fmt4]: CURVEsingle = ', CURVEsingle

! ....................................................................
      ! Make sure that mtrait has a valid value. 1 is used as a defauly
      ! because that's its value in PVR's original code, where mtrait
      ! had a hard-coded value.
      if (mtrait == 1 .or. mtrait == 3 .or. mtrait == 4) then
        continue
      else
        mtrait = 1
      end if
!      print *,'Frequency of BP: ', bpfreqin

      zero = 0.

      !herd305(:,:) = 0.d0

      do 5 j=1,3
 5      yld2(j) = 10000
!                                                 Display age factors
      count = count + 1
      if(count == 1 .and. showage == 1) then
        do i=1,6
          if (i == 1) brd = 'H'
          if (i == 2) brd = 'J'
          if (i == 3) brd = 'B'
          if (i == 4) brd = 'G'
          if (i == 5) brd = 'A'
          if (i == 6) brd = 'M'
          if (i == 1 ) print '(26x,a)','Milk    Fat     Prot    SCS   prevDO'
          do fryr = 1990,1966,-6
            !print *,'Brd = ',brd,', Fresh year = ',fryr
            do age=22,99,14
              lacno = (age - 6)/14
              parity = lacno
              if(doprev < 1) doprev = 140
              if(lacno == 1) doprev = 0
              avgfac = 0.
              do frmo = 1,12
                do k=1,5
                  if (k == 1) state = 21
                  if (k == 2) state = 35
                  if (k == 3) state = 42
                  if (k == 4) state = 74
                  if (k == 5) state = 93
                  scs305 = 329
                  call aiplage(brd,age,fryr,frmo,lacno,state,doprev, &
                       agebase,agefac(1),agefac(2),agefac(3))
                  !print *, "Calling adjscs4()"
                  agefac(4) = adjscs4(brd,lacno,frmo,state, &
                       age,i305,scs305)/(scs305*1.d0)
                  do j=1,4
                    if(j < 4) oldfac(j) = adjyld2(j) / (yld2(j) * 1.d0)
                    avgfac(j) = avgfac(j) + agefac(j) / 60.
                    end do
                  end do
               print '(a,2i3,a,4f8.2,i8)','Age-par',age,lacno,' factors ',agefac,doprev
               print '(a,2i3,a,4f8.2,i8)','Age-par',age,lacno,' factors ',oldfac,doprev
!               print *
                end do
              print '(a,2i3,a,4f8.2,i8)','Age-par',age,lacno,' factors ',avgfac,doprev
              end do
            end do
          end do
        end if
!.......................................................... Read record
      read(format4,11) cowid,birth,herd,fresh,length,lacno,doprev,ntd, &
                      (segment(i),i=1,ntd)
! 11   format(t3,a17,t71,a8,t107,a8,t128,a8 &
!            ,t136,i3,t159,i2,t249,i2,50a23)
! Added to support previous days open (thanks to BJH).
 11   format(t3,a17,t71,a8,t107,a8,t128,a8 &
        ,t136,i3,t159,i2,t246,i3,t249,i2,50a23)
      !print *, cowid, " ", birth, " ", herd, " ", fresh, " ", length, " ", lacno, " ", ntd, " ", doprev
!                                                Read test day segments
      do 16 i=1,ntd
        read(segment(i),13) dim(i),super(i),status(i), &
                   Xmilk(i),weigh(i),sample(i),MRD(i), &
                   (yield(1,k,i),k=1,4)
 13     format(i3,5i1,i2,3x,f4.1,3f2.1)
        !print *, "    ", dim(i), super(i), status(i), Xmilk(i), weigh(i), sample(i), MRD(i), yield(1,1,i) & 
        !  , yield(1,2,i), yield(1,3,i), yield(1,4,i)
        if ( status(i) == 3 ) then
          do 14 k=1,4
 14         yield(1,k,i) = zero
        end if
        if(status(i) == 2) then
          yield(1,2,i) = zero
          yield(1,3,i) = zero
          ! JBC 08/11/2009 Commented this out so that records with SCS but no
          ! fat or protein test make it through to bestpred().
          !yield(1,4,i) = zero
        end if
!        if ( yield(1,4,i) == 0.d0 ) then
!            yield(1,4,i) = -0.01
!        end if
 16     continue

! Read the state code stored in bytes 107-108 of Format 4
      herdstate = herd(1:2)
      scs305 = 329
      if ( DEBUGmsgs > 0 ) print *, '[bestpred_fmt4] herdstate: ', herdstate
!                                       Read cow's deviation and me305
!      print *, "Format 4: ", format4
      read(format4,17) me305(1),me305(2),me305(3), &
                       dev305(1),dev305(2),dev305(3)
 17   format(t188,f5.0,2f4.0,f6.0,2f5.0)
      if ( DEBUGmsgs > 0 ) print *, '[bestpred_fmt4]: herd305 before bestpred_fmt4()', herd305
!                                    Backsolve to get mature herd mean
!      print *, '[bestpred_fmt4]: herd305 before back-calculation in fmt4', herd305
!      print *, '[bestpred_fmt4]: herd305 before zero check in bestpred_fmt4()', herd305(1,:)
      do j = 1,4
!        herd305(j) = me305(j) - dev305(j)
! I added this IF to help with Albert's work on using some sort of PPA by
! increasing the herd average for a later-lactation cow to reflect the
! repeatability of her earlier performance. The herd average for a trait
! is only back-calculated from me305 and dev305 if the herd average in
! herd305 is <= 0.
        if ( herd305(1,j) <= 0 ) then
          if ( j < 4 ) then
            herd305(1,j) = me305(j) - dev305(j)
            if ( me305(j) == 0.d0 ) herd305(1,j) = 0.d0
          else
            herd305(1,j) = 0.d0
            !herd305(1,4) = scs305
          end if
        end if
      end do
      herd305(2,:) = herd305(1,:)
!      print *, '[bestpred_fmt4]: herd305 after back-calculation in fmt4', herd305

!      print *, '[bestpred_fmt4]: herd305 after zero check in bestpred_fmt4()', herd305(1,:)
!      herd305(4) = 0.d0
      if(birth == '        ') birth = '19500101'
      if(birth  < '19500101') birth = '19500101'
      age = (pday(fresh) - pday(birth))/30.5
      if(fresh == '        ') age = 84
!                                         Assign missing lactation num
      if(lacno <= 0) lacno = max(1,(age-6)/13)
      if(birth == '19500101') age = 12 + 15*lacno
      read(fresh(1:4),'(i4)') fryr
      read(fresh(5:6),'(i2)') frmo
      read(herd(1:2),'(i2)') state
      if(frmo < 1) frmo = 1
      brd = cowid(1:1)
      if(index('ABGHJM',brd) == 0) brd = 'H'
      parity = lacno
      if(doprev < 1) doprev = 140
      if(lacno == 1) doprev = 0
      !scs305 = 329
!                                    Get factor to adjust 305-d records
!          ageadj4                               is official C function
!          ageadjF                       is substitute Fortran function
!          aiplage             also adjusts pre-1990 data by C function
!     call ageadj4(brd,lacno,frmo,state,age, &
!                  fill2,doprev,docur,yld2,adjyld2)

      if ( DEBUGmsgs > 1 ) then
          print *, 'Age factors before call to aiplage(): ', agefac
          print *, 'aiplage() parameters'
          print *, '--------------------'
          print *, '    brd     : ', brd
          print *, '    age     : ', age
          print *, '    fryr    : ', fryr
          print *, '    frmo    : ', frmo
          print *, '    lacno   : ', lacno
          print *, '    state   : ', state
          print *, '    doprev  : ', doprev
          print *, '    agebase : ', agebase
      end if

      ! If the state code is empty/zero then default to 35, which is the
      ! code for Wisconsin. This prevents a nasty little side-effect with
      ! source 10 records (partial Format 4 records).
      if ( state == 0 ) then
          print *, '[bestpred_fmt4]: Changing state from ', &
              state, ' to 35 (Wisconsin).'
          state = 35
      end if

      call aiplage(brd,age,fryr,frmo,lacno,state,doprev, &
                   agebase,agefac(1),agefac(2),agefac(3))
!                                        Get factor to adjust 305-d SCS
!                 adjscs4                        is official C function
!                 adjscsF                is substitute Fortran function
      !agefac(4) = adjscs4(brd,lacno,frmo,state, &
      !            age,i305,scs305)/(scs305*1.d0)
      !print *, "Calling adjscs4()"
      agefac(4) = adjscs4(brd,lacno,frmo,state, &
                  age,i305,scs305)/(scs305*1.d0)
!      print *, 'Age factors after call to aiplage() :', agefac

      ! JBC 08/11/2008 -- We need to divide herd average SCS by 100 so that it is on the proper
      ! scale before we call bestpred(). The edits system in the routine run processing sends
      ! in SCS in X.XX format, not (X.XX*100). If we do the dividion in bestpred() we end
      ! up with herd average SCS near 0, which results in a systematic underestimation of
      ! ME SCS, and is particularly pronounced for early RIPs.
      herd305(:,4) = herd305(:,4) / 100

!      print *, '[bestpred_fmt4]: length in fmt4 before entry to bestpred: ', length
!      print *, '[bestpred_fmt4: maxprnt: ', maxprnt
!      print *, '[bestpred_fmt4: maxshow: ', maxshow
      !print *, "ntd before calling bestpred(): ", ntd
!      print *, '[bestpred_fmt4]: yield going into bestpred(): ', yield

      call bestpred(mtrait,use3X,maxprnt                        &
                  ,ntd,dim,super,Xmilk,weigh,sample,MRD,yield   &
                  ,herd305,agefac,cowid,fresh,parity,length     &
                  ,DCRm,DCRc,DCRs,YLDvec,PERSvec,RELyld,RELpers &
!                                                              Optional
                  ,Yvec,bump,BLUPout,BLUPn,GRAFplot             &
                  ,DEBUGmsgs,ONscreen,laclen,dailyfreq          &
                  ,INTmethod,maxlen,DAILYbp,DAILYherd           &
                  ,WRITEcurve,CURVEfile,WRITEdata,DATAfile      &
                  ,plotfreq,INTmethodSCS,READparms,UNITSin      &
                  ,UNITSout, breedUNK, dim0, dim0flag           &
                  ,LOGon,LOGfile,LOGfreq,maxshow                &
                  ,herdstate,CURVEsmall, CURVEsingle            &
                  ,region, season                               &
      )
!	  print *, '[bestpred_fmt4]: Leaving bestpred()'
      return
    end subroutine bestpred_fmt4
!----------------------------------------------------------------------
!                                 Fortran routine to standardize yield
!                                !!! Example factors for age,season !!!
!                                !!! Poor substitute for C function !!
!
    subroutine ageadjF(brd,lacno,frmo,state,age, &
                        fill2,doprev,docur,yld2,adjyld2)
!
      character brd*1
      integer*2 lacno,frmo,state,age,fill2,doprev,docur,yld2(3) &
               ,adjyld2(3)
      real agefac
      integer begin, j
      data begin   /0/

      if(begin == 0) then
        print *,'Age factors for example only, use C function instead'
        begin = 1
        end if
!                                                     Age adjustment
      if(lacno == 0) lacno = max(1,(age-6)/13)
      agefac = 1.0 + .0056*(72 - age)
      if(age > 72) agefac = 1.0
!                                                  Season adjustment
      agefac = agefac*(1.0 + .02*sin((frmo - 3)*3.1416/6.0))
!                                             Convert actual to mature
      do 20 j=1,3
 20     adjyld2(j) = yld2(j)*agefac
      return
    end subroutine ageadjF
!----------------------------------------------------------------------
!                              !!Example factors to standardize SCS !!!
!                                !!! Poor substitute for C function !!
!
    integer function adjscsF(brd,lacno,frmo,state,age,dim,scs305)
!
      character brd*1
      integer*2 lacno,frmo,state,age,dim,scs305
      real agefac
      integer begin
      data begin   /0/
      if(begin == 0) then
        print *,'SCS factors for example only, use C function instead'
        begin = 1
        end if
!                                     Adjustment to average age (45 mo)
      if(lacno == 0) lacno = max(1,(age)/13)
      agefac = 1.0 + .0094*(45 - age)
      if(age > 72) agefac = .75
!                                                     Season adjustment
      agefac = agefac*(1.0 + .04*sin((frmo - 9)*3.1416/6.0))
!                                                     Adjust input data
      adjscsF = agefac*scs305
      return
    end function adjscsF
! -------------------------- Compute perpetual days --------------------
!                                                  by Chris Schanefelter
!                                     Year 2000 changes by Paul VanRaden
    integer*2 function pday(date8)
      implicit none
      character  date8*8
      integer modays(12,2) &
        ,i,year,month,day,leap,balance_years        &
        ,number_of_years,year_groups

      data modays /31,28,31,30,31,30,31,31,30,31,30,31  &
                  ,31,29,31,30,31,30,31,31,30,31,30,31/
      data year   /0/
      data month  /0/
!                                                    Accepts 8-byte date
      read(date8,'(i4,2i2)') year,month,day
      year = year - 1900
!                                           0 represents an unknown date
!                                               set to December 31, 1889
      if (year <= 0) then
        pday = -21916
        return
      endif
      leap = year/4                  !     leap years are divisible by 4
      number_of_years = year - 60    !       number years from base year
      year_groups = abs(number_of_years) / 4    !          4 year groups
!                        Calculate number of years not in a 4 year group
      balance_years = abs(number_of_years) - (year_groups * 4)
      if(year .eq. (leap*4)) then
        leap = 2                                            !  leap year
      else
        leap = 1                                       ! not a leap year
      end if
!                                        Establish minimums and maximums
      if (month <= 0) month = 1
      if (month > 12) month = 12
      if (day .eq. 0) day = 1
      if (day .gt. modays(month,leap)) day = modays(month,leap)
      do i = 1, month-1
        day = day + modays(i,leap)
      end do
      if (number_of_years .le. -1 ) then !  before 1/1/1960
!         Adjust days when going away from zero in a negative direction
        day = (365 - day)
!                                              Calculate perpetual days
        pday = -((year_groups * 1461) +  day &
                + (balance_years - 1) * 365)
      else                               ! after and including 1/1/1960
!                For years not in a year group, add 1 leap and 365 days
!                                               for each year completed
        if (balance_years > 0) &
          balance_years = (((balance_years - 1)* 365)+ 366)

!                                              Calculate perpetual days
        pday = (year_groups * 1461) + day + balance_years
      end if
      pday = pday - 1
      return
    end function pday
