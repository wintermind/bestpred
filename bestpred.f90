!e!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME:         bestpred.f90
! VERSION:      2.0 beta
! RELEASED:     01 AUGUST 2007
! AUTHORS:      Paul M. VanRaden (paul@aipl.arsusda.gov)
!               John B. Cole (john.cole@ars.usda.gov)
! DESCRIPTION:  Reads program parameters from the file bestpred.par.
!               This program is part of the bestpred package from AIPL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         aiplDCR:    Fortran / written by Paul VanRaden 11-99
!                      Animal Improvement Programs Lab, USDA
!                      paul@aipl.arsusda.gov
!
!         Computes:    1) Data Collection Rating from test day data.
!                      2) yields (projected 305, partial, 305, 365)
!                      3) milk, fat, protein, and somatic cell score
!                      4) persistency and reliability of persistency
!                      5) factors to adjust for 3X milking
!
!          Methods:    Multi-trait or single-trait best prediction using
!                      a correlation matrix among test day yields.
!
!       References:    1996 J. Dairy Science 79(Suppl. 1):143 (abstr.)
!                      1997 J. Dairy Science 80:3015
!                      1998 Proc. 6th World Conference on Genetics
!                           Applied to Livestock Production 23: 347
! .....................................................................
!   Definitions (ST=single-trait,MT=multi-trait,m=milk,f=fat,p=protein)
!                                                     Control variables
!                         mtrait = 1 does only ST to save CPU time
!                                = 3 does MT m,f,p; 4 does MT m,f,p,scs
!                         mature = 1 adjusts for age etc.,  0 doesn't,
!                                 -1 is for pre-adjusted deviations
!                                    and inputs f,p lbs instead of %
!                          use3X = 0 doesn't adjust for 3X milking,
!                                  1 uses old factors, 2 new factors,
!                                  3 uses phased-in factors over time
!                        maxprnt = number of graphs to display, use
!                                 -1 to display the input variables
subroutine bestpred &
!                                                               Control
      (mtrait,use3X,maxprnt &
!                                                             Input
      ,tstdays,dim,super,Xmilk,weigh,sample,MRD, yield &
      ,herdavg,agefac,cowid,fresh,parity,length &
!                                                            Output
      ,DCRm,DCRc,DCRs,YLDvec,PERSvec,RELyld,RELpers &
!                                                          Optional
      ,Yvec,bump,BLUPout,BLUPn                            &
!                                                      Added by JBC
      ,GRAFplot,DEBUGmsgs,ONscreen                        &
      ,laclen,dailyfreq, INTmethod, maxlen                &
      ,DAILYbp,DAILYherd,WRITEcurve,CURVEfile             &
      ,WRITEdata,DATAfile,plotfreq,INTmethodSCS           &
      ,READparms,UNITSin, UNITSout, breedUNK, dim0        &
      ,dim0flag, LOGon, LOGfile, LOGfreq, maxshow         &
      ,herdstate, CURVEsmall, CURVEsingle                 &
      ,region, season                                     &
      )
! .....................................................................
!                                                                 Input
!       tstdays is number of test days
!       dim     contains days in milk at each test day
!         next: supervision codes, milking and sampling frequencies
!       MRD     is number of milk recorded days averaged together (LER)
!       yield   contains m,f,p, and scs data for each test day
!       herdavg 305-d herd means for m,f,p,scs adjusted for age, not 3X;
!!!             called herd305 in bestpref_main and bestpred_fmt.
!       agefac  has adjustment factors for m,f,p,scs
!       cowid   is 17-byte American ID of cow
!       fresh   is 8-byte YYYYMMDD fresh date
!       parity  is lactation number (assumed / if parity = 0)
!       length  is days from fresh date to termination date
!                                                                Output
!       DCRm   = MT  m  DCR
!       DCRc   = MT f,p DCR
!       DCRs   = MT scs DCR
!       YLDvec = ME    : 305 m,f,p,scs, 365 m,f,p,scs, maxlen m,f,p,scs, laclen m,f,p,scs
!                Actual: 305 m,f,p,scs, 365 m,f,p,scs, maxlen m,f,p,scs, laclen m,f,p,scs
!       PERSvec= MT m,f,p,scs persistency, ST m,f,p,scs persistency
!       RELyld = MT m,f,p,scs Reliability, ST m,f,p,scs Reliability
!       RELpers= MT m,f,p,scs Rel(persist),ST m,f,p,scs Rel(persist)
!                ST values are substituted if no MT values
!                                                              Optional
!       GRAFplot indicates which MFPS plots should be drawn
!       DEBUGmsgs is used to display debugging messages (0/1)
!       ONscreen toggles output to the screen on (1) and off (0)
!       laclen specifies default lactation length that is used to
!          figure out what actuals to return.
!       dailyfreq specifies how often actual dailies are calculated
!       INTmethod specifies the inerpolation method to be used for
!           calculating lactation curves.
!       maxlen specifies the longest possible lactation length
!           to be used in calculating lactation curves.
!       READparms is a flag passed downstream from maindcr to fmt4dcr
!           and then to aipldcr that indicates which parms should be
!           used by that particular subroutine -- those passed in from
!           maindcr (0) or those stored in the file maindcr.par (1).
!       maxshow specified how many plots should be drawn on the screen.
! ....................................................................
!                                                         OTHER OPTIONS
!  integer,save :: mbreed
  integer,parameter :: last=2, part=1
  integer,parameter :: nbump=0
  integer,save :: nbmp
  integer,save :: lacn
  integer,save :: month
  integer,save :: ncall=0
  integer,save :: precise=1
  integer :: READparms, maxshow
  real,save :: maxbump
  real, dimension(:,:,:), allocatable, save :: dyield
  real*8,save :: lplus
  real*8, dimension(:,:,:), allocatable, save :: meanyld
  real*8, dimension(:,:,:,:,:), allocatable, save :: covari
  real*8,save :: corr305(4,4)
  real*8, dimension(:,:,:,:,:), allocatable, save :: covd
  real*8,save :: meanp(2,4)
  real*8,save :: varp(2,4)
  real*8,save :: lacdif
  real*8,save :: sddif
  real*8,save :: stdvar(4,last)
  real*8, dimension(:,:,:), allocatable, save :: sd
  real*8,save :: summ  ! I renamed this because it was colliding with the
                       ! Fortran 90 array function SUM().
  real*8,save :: dV1(4,4)
  real*8,save :: dVd(4,4)
  !real*8,save :: d0(4)
  real*8 :: d0(2,4)
  real*8,save :: vijkl
  real*8,save :: qpq=0.d0
  real*8,save :: tdhi
  real*8,save :: tdlo
  real*8,save :: zero
  integer :: dim0(8), dim0flag
  real :: dim0wgts(2) = (/ 0.285, 0.715 /) ! Should sum to 1.

  ! DAILYbp and DAILYherd will contain daily BP of yield for cows
  ! and herds, respectively, unless dailyfreq=0. If dailyfreq is
  ! 0 then they will contain only zeroes.
  real*8 :: DAILYbp(2,4,maxlen)
  real*8 :: DAILYherd(4,maxlen)

  !character, parameter :: mybreed='H'
  character :: INTmethod, INTmethodSCS, UNITSin, UNITSout
  character*2 :: herdstate
  integer :: region, season
  integer, save :: oldregion, oldseason
  integer, parameter :: maxy=200, maxtd=50
!
!                         mybreed= breed letter for stats display
!                         last   = last lactation group (2 is 2nd&later)
!                         nbump  reduces REL of worst 1/10**nbump lacts
!                                nbump ranges from 1-5 (see Fval table)
!                                = 0 will not reduce REL of bad data
!                         maxtd  is maximum test days in a lactation
!                         maxy   is maximum number of milk, fat, and
!                                protein obs in any lact (<= 4*maxtd)
  integer mtrait,use3X,maxprnt, &
     dim(maxtd),Xmilk(maxtd),weigh(maxtd),sample(maxtd), &
     MRD(maxtd),super(maxtd),tstdays,length, &
     parity,BLUPn(10), DEBUGmsgs             &
     ,ONscreen, laclen, maxlen
  character(100) BLUPout(10)
  real*8 DCRm,DCRc,DCRs,YLDvec(2,16),PERSvec(4),RELyld(4),RELpers(4), &
    var(maxy,maxy),ymean,yield(2,4,maxtd),                            &
    vary,covary,vari(4,4),dev(maxy,1),                              &
    multi(4,1), &
    herdavg(2,4), &
    herd305(2,4),hratio(4),herd365(2,4),herd999(2,4),herdpart(2,4), &
    sum305(4),sum365(4),sum999(4),sumpart(4), &
    dvari(4,4),dcov(4,maxy),dcovp(maxy,4), &
    persist(4,1),pers(4,1), &
    partial(4,maxy),partrec(4,1), &
    xsingl(4),xmulti(4),covsum(4,maxy), &
    agefac(4),varfac(4),test3X(4,maxtd),fact3X(4),part3X(4), &
    Yvec(2,4),DCRvec(4),lacwt(4),milk3X, &
    qCVC1(4,4),Rvar(2,2),Rvec(12), &
    Rvecs,Rvecm,Rvecf,Rvecp,Rcorr, &
    Zdev(maxtd),Zdif,regrel(4), &
    tdregr,bump(4),bumpsd(4),bigbump(4),varb(4)
  real*8, dimension(:,:), allocatable :: curve, curvem, graph1
  real*8, dimension(:,:,:), allocatable :: grafall
  real*8, dimension(:,:,:), allocatable :: curves
  real*8, dimension(:), allocatable :: freq3X

  real*8, dimension(:,:,:), allocatable :: std, tempSTD!, tempDEV
  real*8, dimension(:,:), allocatable :: stdslice, tempSTDslice, tempTD

! These matrices are used to catch covariance between TD and a given DIM.
! '305' is used for 305 d lactations, '365' for 365 d lactations, and
! '999' for user-defined lactation lengths. The 'part's are used for
! calculating actual yields based on DIM of most recent test for RIP.
  real*8  cov305(4,maxy),covp305(maxy,4),                  &
        covar305(4,maxy),vari305(4,4),multi305(4,1),     &
        cov365(4,maxy),covp365(maxy,4),                  &
        covar365(4,maxy),vari365(4,4),multi365(4,1),     &
        cov999(4,maxy),covp999(maxy,4),                  &
        covar999(4,maxy),vari999(4,4),                   &
        multi999(4,1),single305(4,1),single365(4,1),     &
        single999(4,1)

  integer dailyfreq, st_start, plotfreq, breedUNK

  character*8 :: BESTPREDversion = '2.0rc7'
  character*10 :: BESTPREDdate = '11/26/2014'
  character*64 :: BESTPREDname = 'John B. Cole'
  character*64 :: BESTPREDemail = 'john.cole@ars.usda.gov'
  character*64 :: BESTPREDwebsite = 'http://www.aipl.arusda.gov/software/bestpred/'
!
  real, save :: zero0=0.d0, lb=2.205
!
!                               Daily yields of milk, fat, protein, SCS
!                              tested monthly    (first lactation curve)
  real, save :: dyld(12,4,last)=reshape((/ &
      53.3,63.,66.3,65.5,63.9,62.,60.1,58.1,55.6,52.3,49.,47. &
     ,2.3,2.32,2.31,2.28,2.25,2.22,2.19,2.16,2.11,2.03,1.95,1.9 &
     ,1.81,1.9,1.99,2.01,2.01,1.99,1.96,1.92,1.87,1.8,1.74,1.67 &
     ,3.38,2.63,2.41,2.44,2.48,2.52,2.57,2.6,2.66,2.74,2.82,2.9 &
!                                           later lactation curve
     ,75.1,86.,86.9,82.2,77.2,72.2,67.2,62.1,56.2,49.8,46.,42. &
     ,3.33,3.2,3.0,2.85,2.71,2.58,2.45,2.29,2.12,1.92,1.81,1.7 &
     ,2.6,2.55,2.55,2.49,2.4,2.3,2.17,2.04,1.89,1.72,1.62,1.55 &
     ,3.35,2.84,2.8,2.9,3.02,3.15,3.28,3.43,3.58,3.76,3.88,4.0/) &
     ,(/12,4,last/)) &
!                           Daily s.d. of milk, fat, protein, SCS
     ,dsd(12,4,last)=reshape((/ &
      10.7,10.6,9.9,9.6,9.5,9.5,9.5,9.6,9.8,10.3,10.7,11.2 &
     ,.51,.45,.41,.39,.38,.38,.37,.37,.37,.39,.40,.41 &
     ,.34,.31,.28,.28,.28,.29,.29,.29,.30,.32,.33,.34 &
     ,1.51,1.45,1.42,1.41,1.4,1.41,1.41,1.41,1.41,1.42,1.42,1.42 &
!                                            later lactation s.d.
     ,15.,15.3,14.4,13.8,13.3,13.,12.9,13.,13.6,14.7,15.,15.5 &
     ,.75,.69,.62,.58,.55,.52,.51,.51,.53,.57,.58,.59 &
     ,.50,.43,.40,.38,.38,.38,.38,.39,.42,.46,.47,.49 &
     ,1.63,1.69,1.71,1.69,1.66,1.61,1.56,1.51,1.47,1.45,1.44,1.4/) &
     ,(/12,4,last/))
!                                          Mature breed means for 1995
  character breed6*6
  character, save :: breed(6)=(/'A', 'B',  'G',  'H',  'J',  'M'/)
  real, save :: brdyld(6,4,2)=reshape((/ &
      15080., 16616., 13864., 20845., 14120., 14472. &
     ,  589.,   668.,   625.,   763.,   662.,   518. &
     ,  505.,   586.,   481.,   656.,   531.,   476. &
     ,  3.16,   3.22,   3.35,   3.20,   3.31,   2.87 &
     ,18532., 21577., 16850., 25658., 18161., 17475. &  ! and for 2000
     ,  713.,   867.,   745.,   935.,   833.,   624. &
     ,  582.,   714.,   551.,   771.,   644.,   544. &
     ,  2.95,   2.93,   3.29,   3.08,   3.32,   3.06/),(/6,4,2/))
  real, save :: brdsd(6,4)=reshape((/ &
       2300.,  2530.,  2324.,  2946.,  2128.,  2332. &
     ,  115.,   115.,   112.,   119.,   109.,   112. &
     ,   92.,    92.,    92.,    97.,    89.,    92. &
     ,  1.28,   1.21,   1.35,   1.34,   1.18,   1.30/),(/6,4/))
!                        Squared correlations of monthly, daily testing
  real, save :: rmonth(4,last)=reshape((/ &
      .962,.960,.962,.943, &
      .958,.956,.958,.960/),(/4,last/))
!                                       F values for test day bumpiness
   real Fvalue
   real, save :: Fval(5,10)=reshape((/ &
!       F1        F01        F001      F0001      F00001
     2.70555,   6.63495,   10.8277,   15.1369,   19.5118, &
     2.30260,   4.60521,    6.9079,    9.2105,   11.5132, &
     2.08381,   3.78166,    5.4222,    7.0360,    8.6341, &
     1.94487,   3.31921,    4.6168,    5.8783,    7.1185, &
     1.84728,   3.01729,    4.1031,    5.1491,    6.1714, &
     1.77412,   2.80202,    3.7430,    4.6428,    5.5180, &
     1.71673,   2.63937,    3.4746,    4.2683,    5.0371, &
     1.67021,   2.51131,    3.2656,    3.9786,    4.6666, &
     1.63153,   2.40737,    3.0975,    3.7468,    4.3713, &
     1.59873,   2.32096,    2.9589,    3.5565,    4.1298/),(/5,10/))
!                                Expansion amount: full=1 or reduced=.6
!     real, save :: expamt=1.0
  real, save :: expamt=0.6 &
!                                   Repeatabilities and correlations of
!                               random effects in yield and persistency
     ,repty(4)=(/.55, .55, .55, .30/) &
     ,reptp(4)=(/.44, .44, .44, .20/) &
     ,reptyp(4)=(/.00, .00, .00, .00/)
!                              If no repeated records, substitute
!                          heritabilities and genetic correlation
!    ,repty(4)=(/.25, .25, .25, .10/) &
!    ,reptp(4)=(/.10, .10, .10, .05/) &
!    ,reptyp(4)=(/.00, .00, .00, .00/)
!                                         DIM midpoints for persistency
!                                     If value is 0, DIM will be chosen
!                                   to make yield and pers uncorrelated
!                        Now passed in from parm file (JBC, 11/14/2007)
!     integer, save :: dim0(4)=(/161,159,166,155/)
!     integer, save :: dim0(4)=(/  0,  0,  0,  0/)
!                                              LSCS is sum / 305, * 100
  real, save :: X305(4)=(/1., 1., 1., 305./)
  real, save :: X100(4)=(/1., 1., 1., 100./)
  integer q(6),qq(6,maxy),size,i,j,k,l,m,n,i6,brd &
     ,ntests(4),maxyld,minyld,middle(maxtd) &
     ,dimmilk(maxtd),dimsort(maxtd) &
     ,mstart,mstop,kk,MRDi &
     ,testday,yearfr
  integer, save :: oldbrd
  integer usetd
  integer :: maxyield(4),minyield(4)
  integer, dimension(:), allocatable :: dimvec
  integer :: LOGon, LOGfreq
  character cowid*17,fresh*8,brd1*1
  character, save :: plus*1='+'
  character(4) percent
  character(4), dimension(:), allocatable :: fatpct, propct, scs
  character plot(40),plt !,tests(maxlen/plotfreq)
  character, dimension(:), allocatable :: tests
!                        tests(maxlen/6)
  character, save :: capital*13='LUSBACOQPRVXD',small*13='lusbacoqprvxd'
  character(2) MTorST(4)
  character(7), save :: trt(-3:8)=(/'devMilk','dev.Fat','devProt','dev.SCS' &
      ,            '   Milk','    Fat','   Prot','    SCS' &
      ,            'ME.Milk','ME..Fat','ME.Prot','ME..SCS'/)
!
  integer :: GRAFplot(4), GRAFplotMIXED, GRAFplotSUM
  character(4), save :: GRAFname(4)=(/'Milk','Fat ','Prot','SCS '/)
  character(1), save :: SHORTtrait(4) = (/'M','F','P','S'/)
  character(2), save :: STorMTlabel(2) = (/'ST','MT'/)
  integer :: WRITEcurve, WRITEdata, CURVEsmall
  character(len=*) :: CURVEfile, DATAfile, LOGfile
  integer :: YIELDunits(4) = (/5,1,1,1/)
  integer :: YIELDsteps(4) = (/1,1,1,1/)
  character*80 :: LogMessage
  character*4 :: int2str
  integer :: overlaps, ntd, ndup, CURVEsingle

  external vary, covary
  equivalence (breed,breed6)

!     .................................................. STARTUP MESSAGE
  if ( ONscreen == 1 .and. ncall == 0 ) then
    print *,'--------------------------------------------------------------------------------'
    print *,'BESTPRED version ', trim(BESTPREDversion)
    print *,'    Released: ', trim(BESTPREDdate)
    print *,'    Support : ', trim(BESTPREDname), ' (',trim(BESTPREDemail),')'
    print *,'    Website : ', trim(BESTPREDwebsite)
    print *,'--------------------------------------------------------------------------------'
  end if

  ! Log program start-up
  ! You should assign the message to LogMessage rather than passing it as a string literal
  ! in the subroutine call. If you pass the string without assignment it will be papped to
  ! 80 characters using NULLs, which will be written to the log file.
  !if ( LOGon == 1 ) then
  !    LogMessage = 'Starting BESTPRED '//trim(BESTPREDversion)
  !    call log_message(LOGfile,LogMessage)
  !end if

! ---------------------------------------------- Allocate the arrays declared
!                                                as allocatable above.
!  print *, "tstdays on bestpred() entry: ", tstdays
  if ( DEBUGmsgs > 1 ) print *,'Allocating memory'
  if ( ncall == 0 ) then
    if ( LOGon == 1 .and. LOGfreq > 0 ) then
      LogMessage = 'Allocating matrices for means, SD, and covariances'
      call log_message(LOGfile,LogMessage)
    end if
    allocate(dyield(maxlen,4,last))
    allocate(meanyld(4,maxlen,last))
    allocate(sd(4,maxlen,last))
    allocate(covari(4,4,maxlen,maxlen,last))
    allocate(covd(4,4,maxlen,maxlen,last))
  end if
  if ( LOGon == 1 .and. LOGfreq > 0 .and. mod(ncall,LOGfreq) == 0 ) then
    LogMessage = 'Allocating arrays in call '//trim(char(ncall+ICHAR('0')))
    call log_message(LOGfile,LogMessage)
  end if
  allocate(curve(maxlen,maxy))
  allocate(curves(maxlen,maxy,4))
  allocate(curvem(maxlen,maxy))
  allocate(graph1(maxlen,1))
  allocate(grafall(2,maxlen,4))
  allocate(freq3X(maxlen))
  allocate(std(2,4,maxlen))
  allocate(stdslice(1,maxlen))
  allocate(tempTD(4,maxlen))
  allocate(tempSTD(2,4,maxlen))
  allocate(tempSTDslice(1,maxlen))
!  allocate(tempDEV(2,4,maxlen))
  allocate(fatpct(maxlen))
  allocate(propct(maxlen))
  allocate(scs(maxlen))
  allocate(tests(maxlen/plotfreq))
  allocate(dimvec(maxlen))

  ! tempTD is used for storing TD values to be written to an output file.
  ! If a given day is a test day then tempTD stores the milk/components
  ! test value. Otherwise, it holds -999.0.
  std = -999.0
  tempTD = -999.0
  tempSTD = -999.0
!  tempDEV = -999.0
  DAILYbp = -999.0
  DAILYherd = 0.0

  sum305 = 0.
  sum365 = 0.
  sum999 = 0.
  sumpart = 0.

  !herdstate = '00'

!!!  print *, dim(1:10)
!!!  print *, yield(1,1:10)

!  print *, '[bestpred]: length upon entry to bestpred: ', length
!  print *, yield(1,1,1:tstdays)
!  print *, yield(1,2,1:tstdays)
!  print *, yield(1,3,1:tstdays)
!  print *, yield(1,4,1:tstdays)
!  print *, '[bestpred]: herdavg upon entry to bestpred: ', herdavg
!  print *, '[bestpred]: dev: ', dev

  if ( DEBUGmsgs > 1 ) then
    print *, 'UNITSin: ', UNITSin
    print *, 'UNITSout: ', UNITSout
  end if

!  if ( DEBUGmsgs > 1 ) print *, '[bestpred]: herd305 passed from bestpred_fmt4()', herdavg

  ! If the user is sending in pounds we need to convert them to kg.
  if ( UNITSin == 'P' ) then
    if ( LOGon == 1 .and. LOGfreq > 0 .and. mod(ncall,LOGfreq) == 0 ) then
      LogMessage = 'Converting from lbs to kg for interal calculations'
      call log_message(LOGfile,LogMessage)
    end if
    yield(:,1,:) = yield(:,1,:) / lb
    herdavg(:,1:3) = herdavg(:,1:3) / lb
    ! Note that dyld, dsd, brdyld, and brdsd are SAVEd so we only
    ! need to convert them once.
    if ( ncall == 0 ) then
      dyld(:,1:3,:) = dyld(:,1:3,:) / lb
      dsd(:,1:3,:) = dsd(:,1:3,:) / lb
      brdyld(:,1:3,:) = brdyld(:,1:3,:) / lb
      brdsd(:,1:3) = brdsd(:,1:3) / lb
    end if
  end if

! The herd average SCS DOES NOT need to be divided by 100. For Format 4 records that
! is done on bestpred_fmt4() before bestpred is called(). For the edits system, that
! is sone before herd averages are passed to bestpred().
!  herdavg(:,4) = herdavg(:,4) / 100

  if ( dim0flag == 1 ) then
    dim0 = 0
    d0 = 0
  else
    d0 = 0
    do j=1,4
      d0(1,j) = dim0(j)
      d0(2,j) = dim0(4+j)
    end do
  end if

  ! Loop over the elements of GRAFplot to determine the appropriate
  ! labels for the output. This is a little tricky b/c all, some,
  ! or no plats can be requested, each with a separate flag.
  GRAFplotMIXED = 0
  GRAFplotSUM = sum(GRAFplot)
  if ( GRAFplotSUM > 0 ) then
    do j = 1, 4
      if ( j == 1 .and. GRAFplot(j) == 1 ) GRAFplotMIXED = 1
      if ( j == 1 .and. GRAFplot(j) == 2 ) GRAFplotMIXED = 2
      if ( j > 1 .and. GRAFplot(j) == 1 .and. GRAFplotMIXED == 2 ) GRAFplotMIXED = 3
      if ( j > 1 .and. GRAFplot(j) == 2 .and. GRAFplotMIXED == 1 ) GRAFplotMIXED = 3
    end do
  end if

 !                               Print input parameters
 if ( ncall == 0 .and. DEBUGmsgs > 0 ) then
     print *,'[bestpred]: Tipping points (1): ', dim0(1:4)
     print *,'[bestpred]: Tipping points (2): ', dim0(5:8)
 end if
 !if ( maxprnt < 0 ) then
 if ( DEBUGmsgs > 0 ) then
    print *,'======================================================'
    print *,'Inputs to bestpred for cow=',cowid,' fresh=',fresh
    print *,'       tstdays=',tstdays,' parity=',parity, &
                     ' length=',length, ' maxprnt=',maxprnt
    print *,' dim sup frq way sam mrd milk   fat  prot   scs'
    do j = 1,min(tstdays,maxtd)
      if ( UNITSin == 'P' ) then
        print 1,dim(j),super(j),Xmilk(j),weigh(j),sample(j),mrd(j),yield(1,1,j)*lb,yield(1,2:4,j)
      else
        print 1,dim(j),super(j),Xmilk(j),weigh(j),sample(j),mrd(j),yield(1,:,j)
      end if
 1        format(6i4,4f6.1)
      if(tstdays > maxtd) print *,tstdays,' tests but',maxtd,' used'
    enddo
    if ( UNITSin == 'P' ) then
      print 2,'Herd average mlk fat pro scs =', herdavg(:,1:3)*lb, herdavg(:,4)*100
    else
      print 2,'Herd average mlk fat pro scs =', herdavg
    end if
 2      format(a,8f10.2)
    print *,'======================================================'
end if
if (  DEBUGmsgs > 1 ) then
    print *,'======================================================'
    print *,"John's arguments"
    print *,"    GRAFplot     :", GRAFplot
    print *,"    DEBUGmsgs    :", DEBUGmsgs
    print *,"    ONscreen     :", ONscreen
    print *,"    laclen       :", laclen
    print *,"    dailyfreq    :", dailyfreq
    print *,"    INTmethod    :", INTmethod
    print *,"    maxlen       :", maxlen
    print *,"    WRITEcurve   :", WRITEcurve
    print *,"    CURVEfile    :", CURVEfile
    print *,"    WRITEdata    :", WRITEdata
    print *,"    DATAfile     :", DATAfile
    print *,"    plotfreq     :", plotfreq
    print *,"    INTmethodSCS :", INTmethodSCS
    print *,"    READparms    :", READparms
    print *,"    UNITSin      :", UNITSin
    print *,"    UNITSout     :", UNITSout
    print *,"    breedUNK     :", breedUNK
    print *,"    dim0         :", dim0
    print *,"    dim0         :", dim0
    print *,"    dim0flag     :", dim0flag
    print *,"    LOGon        :", LOGon
    print *,"    LOGfile      :", LOGfile
    print *,"    LOGfreq      :", LOGfreq
    print *,"    agefac       :", agefac
    print *,'======================================================'
  end if
!  print *,"    agefac       :", agefac

  ! Moved up from line ~566 so that breed can be passed to interpolate
  read(cowid,'(a1)') brd1
  brd = index(breed6,brd1)
  ! If brd == 0 then then brd1 wasn't found in breed6. Set the default to
  ! breedUNK and print a warning.
  if ( brd == 0 ) then
    if ( DEBUGmsgs == 1 ) print *, "[WARNING]: Unknown breed ", brd1, &
      "provided; defaulting to ", breedUNK, " (", breed(breedUNK), ")."
    if ( LOGon == 1 .and. LOGfreq > 0 .and. mod(ncall,LOGfreq) == 0 ) then
      LogMessage = 'Unknown breed '//trim(char(brd+ICHAR('0')))//' provided; defaulting to ' &
        //trim(char(breedUNK+ICHAR('0')))//' ('//trim(breed(breedUNK))//')'
      call log_message(LOGfile,LogMessage)
    end if
    brd = breedUNK
  end if

  ! If herdstate is 00 then bestpred() is being called by the edits system. We don't want to call interpolate()
  ! every tine the season of calving changes since we're not using regional or seasonal effects, so we're fixing
  ! the region and season at the defaults.
  if ( ( INTmethod == 'L' .or. INTmethod == 'W' ) .and. ( INTmethodSCS == 'L' .or. INTmethodSCS == 'G' ) ) then
    herdstate = '00' ! Avoid flip-flopping on interpolate() unless were using seasons and/or regions
  end if
  if ( herdstate .eq. '00' ) then
    region = 2 ! Midwest
    season = 1 ! Spring
  else
    call state_to_region(herdstate, region, DEBUGmsgs)
    call date_to_season(fresh, season, DEBUGmsgs)
  end if

  if ( ncall == 0 ) then
    oldbrd = brd
    oldregion = region
    oldseason = season
  end if
! -------------------------------------- SECTION DONE DURING FIRST CALL
  if ( ncall > 0 .and. brd .eq. oldbrd .and. dim0flag == 0 .and. season .eq. oldseason .and. region .eq. oldregion ) go to 20
  if ( maxprnt > 0 .and. ncall == 0 ) then
    print *,'Lactation yields are MATURE (adjusted for age, etc.) and ACTUAL (unadjusted) values'
    if ( GRAFplotMIXED == 0 ) print *, 'No lactation curves are being plotted.'
    if ( GRAFplotMIXED == 1 ) print *, 'Plotted lactation curves are MATURE (adjusted for age, etc.)'
    if ( GRAFplotMIXED == 2 ) print *, 'Plotted lactation curves are ACTUAL (unadjusted for age, etc.)'
    if ( GRAFplotMIXED == 3 ) &
      print *, 'Plotted lactation curves are a mix of MATURE (adjusted for age, etc.) and ACTUAL (unadjusted) values'
    if ( mtrait > 1 ) then
      print *,'Multi-trait methods include data for',trt(1:mtrait)
    else
      print *,'!! Only single-trait methods were used !!'
    end if
    if ( use3X == 0 ) print *,'3X adjustments were not used'
    if ( use3X == 1 ) print *,'Old 3X factors were used'
    if ( use3X == 2 ) print *,'New 3X factors were used'
    if ( use3X == 3 ) print *,'New 3X factors phased in over time'
  end if
  if ( nbump > 0 ) then
    if ( maxprnt > 0 ) print *, &
      'DCR was reduced for bumpiest 1 /',10**nbump,'lacts'
    nbmp = nbump
    maxbump = 1.0
  else
    maxbump = 99.
    if ( maxprnt > 0 .and. ncall == 0 ) print *,'DCRs were not adjusted for bumpiness'
  end if
!                                       Compute means, s.d. of test days
!                                           within each lactation group
  do 19 lacn = last,1,-1
!                               Interpolate between monthly means and sd
!                 Note that interpolate() has to be called once for each
!                 trait under evaluation.
    do j = 1,4
      if ( j <= 3 ) then
        if ( DEBUGmsgs > 0 ) then
          print *,'[bestpred]: Calling interpolate for trait ', j, ' with method ', INTmethod, ' and breed ', brd &
            , ' and region ', region, ' and season ', season
        end if
        if ( LOGon == 1 .and. LOGfreq > 0 .and. mod(ncall,LOGfreq) == 0 ) then
          LogMessage = '[bestpred]: Calling interpolate(): trait '//trim(char(j+ICHAR('0')))//', method ' &
            //trim(INTmethod)//', breed '//trim(char(brd+ICHAR('0')))//', region '            &
            //trim(char(region+ICHAR('0')))//', and season '//trim(char(season+ICHAR('0')))
          call log_message(LOGfile,LogMessage)
        end if
        call interpolate(j, month, dyield, lacn, dyld, summ, meanyld, meanp &
              , sd, dsd, INTmethod, DEBUGmsgs, maxlen, brd, region, season)
      else
        if ( LOGon == 1 .and. LOGfreq > 0 .and. mod(ncall,LOGfreq) == 0 ) then
          !LogMessage = 'Calling interpolate with method '//trim(INTmethod)//' and breed ' &
          !  //trim(char(brd+ICHAR('0')))
          LogMessage = '[bestpred]: Calling interpolate(): trait '//trim(char(j+ICHAR('0')))//', method ' &
            //trim(INTmethodSCS)//', breed '//trim(char(brd+ICHAR('0')))//', region '         &
            //trim(char(region+ICHAR('0')))//', and season '//trim(char(season+ICHAR('0')))
          call log_message(LOGfile,LogMessage)
        end if
        if ( DEBUGmsgs > 0 ) then
          print *,'[bestpred]: Calling interpolate for trait ', j, ' with method ', INTmethodSCS, ' and breed ', brd &
            , ' and region ', region, ' and season ', season
        end if
        call interpolate(j, month, dyield, lacn, dyld, summ, meanyld, meanp &
          , sd, dsd, INTmethodSCS, DEBUGmsgs, maxlen, brd, region, season)
      end if
    end do
    do j = 1,4
      stdvar(j,lacn) = 0.
      stdvar(j,lacn) = covary(j,0,0,305,maxlen,1,covari,covd,sd,lacn,last,precise,1)
    end do
  lplus = (precise - 1)/2.
  do 8 j=1,4
    do 7 k=1,4
!                               Compute 305-d covariances from test day
!                                            Get persistency statistics
      corr305(j,k) = 0.d0
      dV1(j,k) = 0.d0
      dVd(j,k) = 0.d0
      do 7 i=1,305
!                                               Save CPU if precise > 1
        do 7 l=1,305,precise
          vijkl = vary(i,j,1,2,2,1,l,k,1,2,2,1,sd,lacn,last,maxlen)
          vijkl = vijkl*min(precise,306-l)
          corr305(j,k) = corr305(j,k) + vijkl
          dV1(j,k) = dV1(j,k) + i*vijkl
 7            dVd(j,k) = dVd(j,k) + i*(l+lplus)*vijkl
          if ( dim0flag == 1 ) then
            dim0(4*(lacn-1)+j) = dV1(j,j) / stdvar(j,lacn)
            d0(lacn,j) = dim0(4*(lacn-1)+j)
          end if
! 8    d0(j) = dim0(j)
 8 continue
  do 10 j=1,4
    meanp(lacn,j) = meanp(lacn,j) - d0(lacn,j)*meanyld(j,305,lacn)
    varp(lacn,j) = dVd(j,j) - 2*dV1(j,j)*d0(lacn,j) + stdvar(j,lacn)*d0(lacn,j)**2
    do 10 k=1,4
      do 10 l=1,maxlen
 10         covd(j,k,l,305,lacn) = (covd(j,k,l,305,lacn) - &
                 covari(j,k,l,305,lacn)*d0(lacn,j)) / dsqrt(varp(lacn,j))
  do 12 i=1,305
    tdhi = dyield(i,1,lacn) + tdregr
    tdlo = dyield(i,1,lacn) - tdregr
 12     qpq = qpq + (i - d0(lacn,1))**2
  if ( maxprnt > 0 .and. ncall == 0 ) then
    if ( lacn == last ) print *,'=======================' &
        ,'=============================================='
    if ( lacn == last ) print *,'Breed ',brd1,' Lact  Mean, ' &
        ,'St.Dev.   Max.DCR   Lactation Corr.        mid-DIM'
    if ( lacn < last ) plus = ' '
  end if
  do 18 j=1,4
    lacdif = meanyld(j,305,lacn) / meanyld(j,305,last)
    zero = 0.
    sddif = dsqrt(stdvar(j,lacn) / stdvar(j,last))
    if ( maxprnt > 0 .and. ncall == 0 ) then
      if ( UNITSout .eq. 'P' ) then
        print 17,  &
        trt(j+4),lacn,plus,brdyld(brd,j,2)*X100(j) &
          *lacdif*lb,brdsd(brd,j)*X100(j)*sddif*lb,100./rmonth(j,lacn) &
          ,(corr305(j,k)/dsqrt(stdvar(j,lacn)*stdvar(k,lacn)) &
          ,k=1,4),d0(lacn,j)
      else
        print 17,  &
          trt(j+4),lacn,plus,brdyld(brd,j,2)*X100(j) &
          *lacdif,brdsd(brd,j)*X100(j)*sddif,100./rmonth(j,lacn) &
          ,(corr305(j,k)/dsqrt(stdvar(j,lacn)*stdvar(k,lacn)) &
          ,k=1,4),d0(lacn,j)
      end if
    end if
 17   format(a7,i4,a1,2f8.0,f9.1,1x,4f6.2,f8.0)
 18   continue
!     print *,' qpq = ',qpq
      !print *, lacn, ': corr305: ', corr305
      !print *, lacn, ': stdvar: ', stdvar(:,lacn)
 19   continue
  if ( dim0flag == 1 ) then
    print *,'New tipping points (1): ', dim0(1:4)
    print *,'New tipping points (2): ', dim0(5:8)
  end if
  if ( oldbrd .ne. brd ) then
!    if ( DEBUGmsgs > 0 ) print *, '[bestpred]: Breed changed from ', &
!      breed(oldbrd), ' to ', breed(brd), ' in call ', ncall+1, '.'
    oldbrd = brd
  end if
  if ( oldregion .ne. region ) oldregion = region
  if ( oldseason .ne. season ) oldseason = season
!!  if ( maxprnt <= 0 ) go to 20
  if ( maxprnt > 0 .and. ncall == 0 ) then 
    print *,'===========================================' &
              ,'=========================='
  end if
! -------------------------------------------------- DEFINE PLOT LETTERS
  if ( GRAFplot(1) > 0 .or. GRAFplot(2) > 0 .or. GRAFplot(3) > 0 .or. GRAFplot(4) > 0 ) then
    print*,'Compute DCR, predict yield of cow ---, graph contemps ...'
    print*,'Supervised =S    ampm1of2=A  2of3=B  1of3=C' &
           ,' ver=V  LER=L dup=D'
    print*,'OwnerSample=O  OSampm1of2=P  2of3=Q  1of3=R' &
           ,'  ownerLER=U, bad=X'
    print*,'Capital letters if milk was sampled, small if milk-only'
  end if
!
! ------------------------------------ SECTION DONE FOR EVERY LACTATION
20   ncall = ncall + 1
!                                          Assign parities to age group

  lacn = min(parity,last)
  if(lacn <= 0) lacn = last
  read(cowid,'(a1)') brd1
!                                  Check for 0 tests or length > maxlen
  if ( length > maxlen ) length = maxlen
  if ( tstdays <= 0 ) then
    ntd = 1
    dim(1) = maxlen + 1
  else
    ntd = tstdays
  end if
!  print *, "ntd    :", ntd
!  print *, "tstdays:", tstdays
  if ( ntd > maxtd ) ntd = maxtd
!                                  Allow for unsorted test day segments
!                                           Exclude duplicate test days
!                                              Check DIM and MRD ranges
  dimvec = 0
  MRDi = 0
  do 23 i = 1,ntd
!	if ( i >= 20 ) then
!        print *, 'cow ', cowid, ', test ', i, ' milk:', yield(1,1,i), ' fat :', yield(1,2,i) &
!               , ' prot:', yield(1,3,i), ' SCS :', yield(1,4,i), ' MRD: ', MRD(i)
!	end if
    if ( dim(i) <= maxlen ) then
      if ( dim(i) < 1 ) then
        print *,'Bad DIM',dim(i),' cow ',cowid,' test',i
      else
        if ( dim(i) - MRD(i) < 0 ) MRD(i) = dim(i)
        MRDi = MRD(i)
        do 22 k=1,MRDi
          if ( dimvec(dim(i)-k+1) > 0 ) then
            if ( MRDi == 1 ) then
              ndup = ndup + 1
              if ( ndup <= 10 ) print * &
                ,'Warning: cow ',cowid,' has duplicate tests ' &
                ,'lact',parity,' day',dim(i)
            else
              overlaps = overlaps + 1
              if ( overlaps <= 10 ) print * &
                ,'Warning: cow ',cowid,' LER segments overlap' &
                ,' lact',parity,' day',dim(i)-k+1
              MRD(i) = k - 1
              go to 23
            end if
          else
            dimvec(dim(i)-k+1) = i
!!!            print *, 'dimvec(', dim(i)-k+1, '):', dimvec(dim(i)-k+1)
          end if
 22       continue
        end if
      end if
 23   continue
!                                            Adjustments for 3X milking
  read(fresh,'(i4)') yearfr
  !print *, '[bestpred]: dim: ', dim
  !print *, 'cowid: ', cowid, ' (', parity, ')'
  !print *, '[bestpred]: meanyld before 3X adjustment', meanyld
!  print *, '[bestpred]: herdavg before 3X adjustment', herdavg
!  print *, '[bestpred]: ntd before 3X adjustment', ntd
!  print *, '[bestpred]: dim before 3X adjustment', dim
  call adjust3X &
       (ntd,dim,length,yearfr,parity,Xmilk,meanyld, &
        test3X,fact3X,part3X,use3X,maxtd,last,maxlen)
!  print *, '[bestpred]: herdavg after 3X adjustment', herdavg
  do j=1,4
!                        If there is no herd average use the breed mean
    if ( herdavg(1,j) <= 0.d0 ) herdavg(1,j) = brdyld(brd,j,lacn)
!                                               Convert SCS mean to sum
    herd305(1,j) = herdavg(1,j)*X305(j)
!                                       Adjust herd mean for 3X milking
    herd305(1,j) = herd305(1,j)*fact3X(j)
!                                       Deviations are already adjusted
    if ( herd305(1,j) <= 0.d0 ) then
      if(yearfr <  2000) herd305(1,j) = brdyld(brd,j,1)*X305(j)
      if(yearfr >= 2000) herd305(1,j) = brdyld(brd,j,2)*X305(j)
    endif
    hratio(j) = herd305(1,j) / meanyld(j,305,lacn)
    !print *, '[bestpred]: herd305 after 3X adjustment', herd305
!   Create the herd standard curve using a breed mean.
    do i = 1,maxlen
      std(1,j,i) = ymean(j,i,maxlen,1,dyield,lacn,last,hratio,herd305(1,j))
      ! Accumulate lactation totals based on the standard curve.
      if ( i <= 305 ) sum305(j) = sum305(j) + std(1,j,i)
      if ( i <= 365 ) sum365(j) = sum365(j) + std(1,j,i)
      if ( i <= laclen ) sum999(j) = sum999(j) + std(1,j,i)
      if ( i <= length ) sumpart(j) = sumpart(j) + std(1,j,i)
    end do

!   If a herd average was provided by the user adjust the standard curve calculated
!   using the breed mean to a curve based on the herd average by adjusting the curve
!   up or down, as appropriate.
!!! This looks like it never gets invoked because of the code at ca. 797. It's
!!! also not clear to me as of 04/21/2015 exactly what is supposed to happen
!!! here.
    if ( herdavg(1,j) <= 0.d0 ) then
      herd305(1,j) = herd305(1,j) * ( herd305(1,j) / sum305(j) )
    end if
!!!

    herd365(1,j) = herd305(1,j) * ( sum365(j) / sum305(j) )
    herd999(1,j) = herd305(1,j) * ( sum999(j) / sum305(j) )
    herdpart(1,j) = herd305(1,j) * ( sumpart(j) / sum305(j) )

!                                         Actual variances change with:
!                                             lactation, breed, and age
    varfac(j) = 0.
    varfac(j) = brdsd(brd,j)*305./dsqrt(stdvar(j,lacn))
!   minyield and maxyield are 4-by-1 vectors that contain the min and
!   maxyield for each trait for use in plotting.
    maxyield(j) = ymean(j,60,maxlen,1,dyield,lacn,last &
      ,hratio,herd305) / (agefac(j)*part3X(j))
    minyield(j) = ymean(j,maxlen,maxlen,1,dyield,lacn,last &
       ,hratio,herd305) / (agefac(j)*part3X(j))
!   Multiply the components min and max by 10. These values are used only
!   for plotting so we don't have to keep them in their original units.
    if ( j == 2 .or. j == 3 ) then
        minyield(j) = minyield(j)*minyield(1)*10
        maxyield(j) = maxyield(j)*minyield(1)*10
    end if
    if ( j == 4 ) then
        minyield(j) = minyield(j)*10
        maxyield(j) = maxyield(j)*10
    end if
!   Note that we have to swap the high and low values for SCS.
    if ( j == 4 ) then
        k = minyield(j)
        minyield(j) = maxyield(j)
        maxyield(j) = k
    end if
  end do
!  print *, '[bestpred]: herdavg after all adjustments', herdavg
!  print *, '[bestpred]: herd305 after all adjustments', herd305
! maxyld and minyld are the minimum and maximum yields for MILK
! used by PVR's original plotting code.
  maxyld = ymean(1,60,maxlen,1,dyield,lacn,last &
         ,hratio,herd305(1,:)) / (agefac(1)*part3X(1))
  minyld = ymean(1,maxlen,maxlen,1,dyield,lacn,last &
         ,hratio,herd305(1,:)) / (agefac(1)*part3X(1))
  size = 0

!if ( DEBUGmsgs > 1 ) then
!  print *, '[bestpred]: herdavg :', herdavg(1,:)
!  print *, '[bestpred]: herd305 :', herd305(1,:)
!  print *, '[bestpred]: herd365 :', herd365(1,:)
!  print *, '[bestpred]: herd999 :', herd999(1,:)
!  print *, '[bestpred]: herdpart: ', herdpart(1,:)
!  print *, '[bestpred]: minyield: ', minyield
!  print *, '[bestpred]: maxyield: ', maxyield
!end if
! ................................... BEGIN MULTI-TRAIT BEST PREDICTION
! Always do MT if mtrait is 3 (MFP) or 4 (MFPS).
if ( mtrait > 1 ) then
  if ( DEBUGmsgs > 0 ) print *,'Doing MT prediction, mtrait= ', mtrait
!                                                    Edit test day data

!  print *, '[bestpred]: dev before adjusting for 3X, etc.: ', dev

  size = 0
  ntests = 0
  usetd = 0
  do 40 i=1,ntd
    if(dim(i) <   1) go to 40
    if(dim(i) > maxlen) go to 40
    if(dimvec(dim(i)) .ne. i) go to 40
    usetd = usetd + 1
    dimsort(usetd) = dim(i)
    if(sample(i) > Xmilk(i)) sample(i) = Xmilk(i)
    if(MRD(i) > 1) weigh(i) = Xmilk(i)
    !!! 08/28/2007 JBC -- This fixes the problem with out-of-bounds indexing
    !!!                   into dimsort. I think this is okay since dimsort is
    !!!                   used only for milk.
    if ( j == 1 ) dimsort(i) = dim(i)
    do 35 j=1,4
      if(weigh(i) < 1) yield(1,j,i) = zero
      if(sample(i) < 1) then
        yield(1,2,i) = zero
        yield(1,3,i) = zero
!                                            Allow SCS with 0 samples
        if(yield(1,4,i) > zero) sample(i) = 1
      end if
      if(yield(1,j,i) <= zero) go to 35
!                                                     Change % to lbs
      if ( j == 2 ) yield(1,j,i) = yield(1,1,i)*yield(1,j,i)/100.d0
      if ( j == 3 ) yield(1,j,i) = yield(1,1,i)*yield(1,j,i)/100.d0
      if(yield(1,1,i) > maxyld) maxyld = yield(1,1,i)
      if(yield(1,1,i) < minyld) minyld = yield(1,1,i)
      if ( j > mtrait ) go to 35
      if ( size == maxy ) then
        print *,' Too much data (size > maxy) for cow ', cowid
        go to 35
      end if
      size = size + 1
      !print *, "i: ", i, " dim(i): ", dim(i), " j: ", j, " ntests(j): ", ntests(j)
      if(dim(i) <= maxlen) then
        ntests(j) = ntests(j) + 1
        if ( j == 1 ) ntests(1) = ntests(1) + MRD(i) - 1
      end if
      q(1) = dim(i)
      q(2) = j
      q(3) = super(i)
      q(4) = Xmilk(i)
      if(j == 1) q(5) = weigh(i)
      if(j > 1) q(5) = sample(i)
      q(6) = MRD(i)
      qq(:,size) = q(:)

!                                                   Adjust cow for 3X,
!                                                 Contemp mean for age
      dev(size,1) = yield(1,j,i)*test3X(j,i) - &
        ymean(j,dim(i),maxlen,MRD(i),dyield,lacn,last &
        ,hratio,herd305(1,:))/agefac(j)

! Put the unadjusted TD value into tempTD for later use
      if ( yield(1,j,i) /= 0.0 ) tempTD(j,dim(i)) = yield(1,j,i)
!      tempDEV(1,j,dim(i)) = dev(size,1)
!                                                GET TEST DAY VARIANCE
      do  k=1,size
        kk = qq(2,k)
        var(size,k) = vary(q(1),q(2),q(3),q(4),q(5),q(6), &
          qq(1,k),qq(2,k),qq(3,k),qq(4,k),qq(5,k),qq(6,k), &
          sd,lacn,last,maxlen) &
            * varfac(j)*varfac(kk) / (agefac(j)*agefac(kk))
        var(k,size) = var(size,k)
      end do

! ----------------------------------------- Store covariances for graph
      do m = 1,maxlen/dailyfreq
        do l = 1,4
          curves(m,size,l) = vary(q(1),q(2),q(3),q(4),q(5),q(6), &
             m*dailyfreq,l,1,2,2,1,sd,lacn,last,maxlen) &
             * varfac(j)*varfac(l) / (agefac(j)*agefac(l))
        end do
      end do

! Look up covariance with 305
! covary() returns the covariance between TD DIM and lactation yields.
! Note that we're calculating for three endpoints: 305 d, 365 d, and a
! user-defined lactation length (maxlen). If maxlen is equal to either
! 305 or 365 we'll only do the calculation once.
      do 34 k=1,4
!       305-d lactation
        cov305(k,size) = covary(k,j,dim(i),305,maxlen,MRD(i),covari,covd, &
          sd,lacn,last,precise,1) * varfac(j)*varfac(k) / agefac(j)
        covp305(size,k) = cov305(k,size)
!       365-d lactation
        cov365(k,size) = covary(k,j,dim(i),365,maxlen,MRD(i),covari,covd, &
          sd,lacn,last,precise,1) * varfac(j)*varfac(k) / agefac(j)
        covp365(size,k) = cov365(k,size)
!       If laclen is either 305 or 365 don't redo those calculations.

        if ( laclen .ne. 305 .and. laclen .ne. 365 ) then
          cov999(k,size) = covary(k,j,dim(i),laclen,maxlen,MRD(i),covari, &
            covd,sd,lacn,last,precise,1) * varfac(j)*varfac(k) / agefac(j)
          covp999(size,k) = cov999(k,size)
        end if
        dcov(k,size) = covary(k,j,dim(i),305,maxlen,MRD(i),covari,covd, &
           sd,lacn,last,precise,2) * varfac(j)*varfac(k) / agefac(j)
        dcovp(size,k) = dcov(k,size)

!                                       Sum covariances for part record
        covsum(k,size) = 0.d0
        if ( part == 0 ) go to 34
        do 33 m = 1,length
 33       covsum(k,size) = covsum(k,size) + vary(q(1),q(2),q(3), &
            q(4),q(5),q(6),m,k,1,2,2,1,sd,lacn,last,maxlen) &
            * varfac(j)*varfac(k) / agefac(j)
 34     continue

!                                       Save sorted milk tests for graph
!        call binsort(dimsort,size,maxtd)
        call binsort(dimsort,usetd,maxtd)
!        if ( j == 1 ) usetd = size
!        do k=1,size
        do k=1,usetd
          dimmilk(k) = dimsort(k)
        end do
 35     continue
 40   continue

! print *, '[bestpred]: dev after adjusting for 3X, etc. : ', dev

! dev(size,1) = yield(1,j,i)*test3X(j,i) - &
!               ymean(j,dim(i),maxlen,MRD(i),dyield,lacn,last &
!               ,hratio,herd305(1,:))/agefac(j)

 if ( DEBUGmsgs > 1 ) then 
     print *, '[bestpred]: Test day data for milk used to calculate deviations.'
     print *, ''
     print *, 'dev               yield    test3X   ymean    hratio   herd305    agefac'
     print *, '-----------------------------------------------------------------------'
     do i = 1, ntd
             j = 1
             print 1701, dev(i,1), yield(1,j,i), test3X(j,i), ymean(j,dim(i),maxlen,MRD(i),dyield,lacn,last,hratio,herd305(1,:)), &
                 hratio(j), herd305(1,j), agefac(j)
     end do
     1701    format(f48.2, 1x, f7.2, 1x, f7.2, 1x, f7.2, 1x, f7.2, 1x, f7.2, 1x, f7.2 )
 end if

!  print *, '[bestpred.f90]: dimsort', dimsort

!                                                 Check for 0 test days
  if ( size == 0 ) then
    do l=1,4
      vari(l,l) = 0.0
      dvari(l,l) = 0.0
      vari305(l,l) = 0.0
      vari365(l,l) = 0.0
      vari999(l,l) = 0.0
      partrec(l,1) = 0.0
      multi(l,1) = 0.0
      multi305(l,1) = 0.0
      multi365(l,1) = 0.0
      multi999(l,1) = 0.0
      persist(l,1) = 0.0
    end do
    graph1 = 0.0
  else
! -------------------------------- Invert variance -- do only once.
    call invrt2(var,maxy,size,cowid)
! -------------------------------- Multiply covariances, compute multi-trait
!                                  estimates,correlations, expansion factors
!   We need to do this 2 or 3 times, once for each distinct value of laclen
!   305-d lactation
!    print *, '[bestpred]: dev before call mult(...)     : ', dev
!    print *, '[bestpred]: multi305 before call mult(...): ', multi305
    call mult(covar305,cov305,var,4,size,size,4,maxy,maxy)
    call mult(vari305,covar305,covp305,4,size,4,4,maxy,4)
    call mult(multi305,covar305,dev,4,size,1,4,maxy,1)
!    print *, '[bestpred]: dev after call mult(...)     : ', dev
!    print *, '[bestpred]: multi305 after call mult(...): ', multi305
!   365-d lactation
    call mult(covar365,cov365,var,4,size,size,4,maxy,maxy)
    call mult(vari365,covar365,covp365,4,size,4,4,maxy,4)
    call mult(multi365,covar365,dev,4,size,1,4,maxy,1)
!   laclen-d lactation -- note that we avoid 2 matrix multiplications
!   if laclen == 305 or laclen == 365.
    if ( laclen == 305 ) then
        multi999 = multi305
    else if ( laclen == 365 ) then
        multi999 = multi365
    else
        call mult(covar999,cov999,var,4,size,size,4,maxy,maxy)
        call mult(vari999,covar999,covp999,4,size,4,4,maxy,4)
        call mult(multi999,covar999,dev,4,size,1,4,maxy,1)
    end if
!                                  Compute graph points, partial records
!   Do only once because these values do not depend on covariances between
!   TD and lactation yields.
    call mult(curvem,curves,var,maxlen/dailyfreq,size,size,maxlen,maxy,maxy)
    call mult(graph1,curvem,dev,maxlen/dailyfreq,size,1,maxlen,maxy,1)
    do j = 1,mtrait
      do l = 1,maxlen/dailyfreq
        do i = 1,size
          curve(l,i) = curves(l,i,j)
        end do
      end do
      call mult(curvem,curve,var,maxlen/dailyfreq,size,size,maxlen,maxy,maxy)
      call mult(graph1,curvem,dev,maxlen/dailyfreq,size,1,maxlen,maxy,1)
      grafall(1,:,j) = 0.
      do l = 1,maxlen,dailyfreq
        grafall(1,l,j) = graph1(l,1)
      end do
      DAILYbp(1,j,:) = -999.0 ! This preserves the missing value indicator
      do l = 1,maxlen,dailyfreq
        DAILYbp(1,j,l) = graph1(l,1)
      end do
    end do

! Partials need to be calculated only once.
    if ( part == 1 ) then
      call mult(partial,covsum,var,4,size,size,4,maxy,maxy)
      call mult(partrec,partial,dev,4,size,1,4,maxy,1)
    end if

! For the moment persistency will stay defined on a 305 d basis. As a
! result, we will only calculate and return 305 d persistencies.
!                                                   Estimate persistency
    call mult(covar305,dcov,var,4,size,size,4,maxy,maxy)
    call mult(dvari,covar305,dcovp,4,size,4,4,maxy,4)
    call mult(persist,covar305,dev,4,size,1,4,maxy,1)
  end if

! Now we're going to fill up the YLDvec in the following order:
! [305 m,f,p,s|365 m,f,p,s|999 m,f,p,s|part m,f,p,s]
  do l=1,4
    !print *, '[bestpred]: [MT] trait: ', l
    YLDvec(1,l) = multi305(l,1) + herd305(1,l)
    !print *,'[MT]    multi305 : ', multi305(l,1)
    !print *,'[MT]    herd305  : ', herd305(1,l)
    !print *,'[MT]    YLDvec   : ', YLDvec(1,l)
    YLDvec(1,4+l) = multi365(l,1) + herd365(1,l)
    YLDvec(1,8+l) = multi999(l,1) + herd999(1,l)
!                                     Part records only for multi-trait
!    if ( part == 0 .or. mtrait < l ) then
     if ( part == 0 .or. mtrait < l .or. length == 0 ) then ! From J. Boyer
        YLDvec(1,12+l) = 0.d0
        go to 45
    end if
!                                                   herd multiplicative
    YLDvec(1,12+l) = partrec(l,1) + meanyld(l,length,lacn) * hratio(l)

!                                                Multi-trait reliability
! Note that REL is currently calculated for a 305-d lactation. We may need to
! look at this later b/c REL for long lactations could be lower than a 305-d
! lactation based on the same number of TD.
 45 RELyld(l) = vari305(l,l) / (stdvar(l,lacn)*varfac(l)**2)
    regrel(l) = RELyld(l)
!                                      Adjust scale so monthly DCR = 100
    DCRvec(l) = 100.*RELyld(l) / rmonth(l,lacn)
    PERSvec(l) = persist(l,1)
    RELpers(l) = dvari(l,l) / varfac(l)**2
    xmulti(l) = 1.d0 / RELyld(l)
  end do
  n = 0
end if

! ...................................... COMPUTE SINGLE-TRAIT ESTIMATES
! We also want to run ST if mtrait is 3 (SCS ST) or 1 (all ST)
if ( mtrait < 4 ) then
  !if ( DEBUGmsgs > 0 ) print *,'Doing ST prediction, mtrait= ', mtrait
  ! We only want to do the calculations that we need. If we've
  ! already done MT MFP we don't need to do ST MFP.
  if ( mtrait == 1 ) then
    st_start = 1
  else if ( mtrait == 3 ) then
    st_start = 4
  else
    print *, '[ERROR]: Unrecognized value of mtrait (1|3|4): ', mtrait
    !!!Took this out to get a clean build. I don't think it's a problem. JBC.
    !goto 55
  end if
  !print *, 'st_start: ', st_start
  do 55 j=st_start,4
    ntests(j) = 0
    size = 0
!                                        Variance matrix for each trait
    do 50 i=1,ntd
      if(dim(i) <   1) go to 50
      if(dim(i) > maxlen) go to 50
      if(dimvec(dim(i)) .ne. i) go to 50
      if(sample(i) > Xmilk(i)) sample(i) = Xmilk(i)
      if(MRD(i) > 1) weigh(i) = Xmilk(i)
      if(weigh(i) < 1) yield(1,j,i) = zero
      if(sample(i) < 1) then
        yield(1,2,i) = zero
        yield(1,3,i) = zero
!                                            Allow SCS with 0 samples
        if(yield(1,4,i) > zero) sample(i) = 1
      end if
      if(yield(1,j,i) <= zero) go to 50
!                                                     Change % to lbs
      if(j == 2) yield(1,j,i) = yield(1,1,i)*yield(1,j,i)/100.d0
      if(j == 3) yield(1,j,i) = yield(1,1,i)*yield(1,j,i)/100.d0
      if(yield(1,1,i) > maxyld) maxyld = yield(1,1,i)
      if(yield(1,1,i) < minyld) minyld = yield(1,1,i)
    !if( j > mtrait ) go to 50
      if ( size == maxy ) then
        print *,' Too much data (size > maxy) for cow ',cowid
        go to 50
      end if
      size = size + 1
      dimsort(size) = dim(i)
      if(dim(i) <= maxlen) then
        ntests(j) = ntests(j) + 1
        if(j == 1) ntests(1) = ntests(1) + MRD(i) - 1
      end if
      q(1) = dim(i)
      q(2) = j
      q(3) = super(i)
      q(4) = Xmilk(i)
      if(j == 1) q(5) = weigh(i)
      if(j > 1) q(5) = sample(i)
      q(6) = MRD(i)
      qq(:,size) = q(:)
!                                              Adjust test days for 3X
!                                                 Contemp mean for age
      dev(size,1) = yield(1,j,i)*test3X(j,i) -  &
        ymean(j,dim(i),maxlen,MRD(i),dyield,lacn,last &
        ,hratio,herd305(1,:)) / agefac(j)

      do k=1,size
        kk = qq(2,k)
        var(size,k) = vary(q(1),q(2),q(3),q(4),q(5),q(6), &
          qq(1,k),qq(2,k),qq(3,k),qq(4,k),qq(5,k),qq(6,k), &
          sd,lacn,last,maxlen) * &
          varfac(j)*varfac(kk)/(agefac(j)*agefac(kk))
        var(k,size) = var(size,k)
      end do

! Put the unadjusted TD value into tempTD for later use
      if ( yield(1,j,i) /= 0.0 ) tempTD(j,dim(i)) = yield(1,j,i)
!      tempDEV(1,j,dim(i)) = dev(size,1)
! ----------------------------------------- Store covariances for graph
      do m=1,maxlen/dailyfreq
        ! curve contains the covariance between each TD and the
        ! sum of current DIM.
        curve(m,size) = vary(q(1),q(2),q(3),q(4),q(5),q(6), &
          m*dailyfreq,j,1,2,2,1,sd,lacn,last,maxlen) &
          * varfac(j)**2 / agefac(j)**2
      end do
! ----------------------------------------- Look up covariance with yields
!     305-d lactation
      cov305(1,size) = covary(j,j,dim(i),305,maxlen,MRD(i),covari,covd, &
        sd,lacn,last,precise,1) * varfac(j)**2 / agefac(j)
      covp305(size,1) = cov305(1,size)
!     365-d lactation
      cov365(1,size) = covary(j,j,dim(i),365,maxlen,MRD(i),covari,covd, &
        sd,lacn,last,precise,1) * varfac(j)**2 / agefac(j)
      covp365(size,1) = cov365(1,size)
!     If laclen is either 305 or 365 don't redo those calculations.
      if ( laclen .ne. 305 .and. laclen .ne. 365 ) then
        cov999(1,size) = covary(j,j,dim(i),laclen,maxlen,MRD(i),covari, &
          covd,sd,lacn,last,precise,1) * varfac(j)**2 / agefac(j)
        covp999(size,1) = cov999(1,size)
      end if
      dcov(1,size) = covary(j,j,dim(i),305,maxlen,MRD(i),covari,covd, &
        sd,lacn,last,precise,2) * varfac(j)**2 / agefac(j)
      dcovp(size,1) = dcov(1,size)

!                                 Summarize deviations within lactation
      Zdev(i) = dev(size,1)*305/(X305(j) * dsqrt(var(size,size)))

!                                       Sum covariances for part record
! I copied this from the MT section above and made a few changes so that
! it fits into the ST context, such as changing the varfac()s. The array
! indexing is correct -- the references are to covsum(1,size) rather than
! covsum(j,size) because this is the ST run and only the first row of
! covsum is populated in a pass.
      covsum(1,size) = 0.d0
      if ( part == 0 ) go to 50
      do m = 1,length
        covsum(1,size) = covsum(1,size) + vary(q(1),q(2),q(3), &
          q(4),q(5),q(6),m,j,1,2,2,1,sd,lacn,last,maxlen) &
          * varfac(j)**2 / agefac(j)
      end do
    50 continue

    call binsort(dimsort,size,maxtd)
!                                                    Examine smoothness
    bump(j) = 0.d0
    bumpsd(j) = 0.d0
    varb(j) = 0.d0
    if ( size > 1 ) then
      do k=2,size
        Zdif = Zdev(dimvec(dimsort(k))) - Zdev(dimvec(dimsort(k-1)))
        bigbump(j) = max(bump(j),dabs(Zdif))
        bumpsd(j) = bumpsd(j) + Zdif**2
        varb(j) = varb(j) + 2.*(1. - var(k,k-1) / &
          dsqrt(var(k,k)*var(k-1,k-1)))
      end do
!                                            Sum(changes) / Exp(changes)
      bump(j) = bumpsd(j)/varb(j)
      !!! Deal with very smakll or very large bumps
      if ( bump(j) < 0 ) bump(j) = 0.
      if ( bump(j) > 99. ) bump(j) = 99.
    end if
!                                       Save sorted milk tests for graph
    if ( j == 1 ) then
      usetd = size
      do k=1,size
        dimmilk(k) = dimsort(k)
      end do
!                                         Single-trait milk already done
      !if ( mtrait < 2 ) size = 0
      regrel(1) = RELyld(1)
    end if
!                                           Best prediction calculations
    xsingl(j) = 0.d0
    lacwt(j) = 0.d0
    partrec(j,1) = 0.0

!  if(size == 0) go to 55
    if ( size == 0 ) then
      vari(1,1) = 0.0
      dvari(1,1) = 0.0
      vari305(1,1) = 0.0
      vari365(1,1) = 0.0
      vari999(1,1) = 0.0
      partrec(1,1) = 0.0
      single305(1,1) = 0.0
      single365(1,1) = 0.0
      single999(1,1) = 0.0
      pers(1,1) = 0.0
      graph1 = 0.0
    else
! -------------------------------- Invert variance -- do only once.
      call invrt2(var,maxy,size,cowid)
! -------------------------------- Multiply covariances, compute multi-trait
!                                  estimates,correlations, expansion factors
!     We need to do this 2 or 3 times, once for each dustinct value of laclen
!     305-d lactation
      call mult(covar305,cov305,var,1,size,size,4,maxy,maxy)
      call mult(vari305,covar305,covp305,1,size,1,4,maxy,4)
      call mult(single305,covar305,dev,1,size,1,4,maxy,1)
!     365-d lactation
      call mult(covar365,cov365,var,1,size,size,4,maxy,maxy)
      call mult(vari365,covar365,covp365,1,size,1,4,maxy,4)
      call mult(single365,covar365,dev,1,size,1,4,maxy,1)
!     laclen-d lactation -- avoid 2 multiplications if laclen = 305 or 365.
      if ( laclen == 305 ) then
        single999 = single305
      else if ( laclen == 365 ) then
        single999 = single365
      else
        call mult(covar999,cov999,var,1,size,size,4,maxy,maxy)
        call mult(vari999,covar999,covp999,1,size,1,4,maxy,4)
        call mult(single999,covar999,dev,1,size,1,4,maxy,1)
      end if
!                                  Compute graph points, partial records
!     Do only once because these values do not depend on covariances between
!     TD and lactation yields.
      call mult(curvem,curve,var,maxlen/dailyfreq,size,size,maxlen,maxy,maxy)
      call mult(graph1,curvem,dev,maxlen/dailyfreq,size,1,maxlen,maxy,1)

      grafall(1,:,j) = graph1(:,1)
!     Write the daily BP for each trait to the correct row of DAILYbp
       DAILYbp(1,j,:) = -999.0
       where(graph1(:,1) /= 0.0) DAILYbp(1,j,:) = graph1(:,1)
!      DAILYbp(1,j,:) = graph1(:,1)

!                                                 Calculate the partials
!                                                  (actual DIM records)
      if ( part == 1 ) then
        call mult(partial,covsum,var,1,size,size,4,maxy,maxy)
        call mult(partrec,partial,dev,1,size,1,4,maxy,1)
      end if
!                                                 Cov(yield,persistency)
      call mult(qCVC1,covar305,dcovp,1,size,1,4,maxy,4)
!                                               Single-trait persistency
      call mult(covar305,dcov,var,1,size,size,4,maxy,maxy)
      call mult(dvari,covar305,dcovp,1,size,1,4,maxy,4)
      call mult(pers,covar305,dev,1,size,1,4,maxy,1)
    end if

!   Now we're going to fill up the YLDvec in the following order:
!   [305 m,f,p,s|365 m,f,p,s|999 m,f,p,s|part m,f,p,s]
    YLDvec(1,j) = single305(1,1) + herd305(1,j)
    YLDvec(1,4+j) = single365(1,1) + herd365(1,j)
    YLDvec(1,8+j) = single999(1,1) + herd999(1,j)
    if ( part == 0 .or. length == 0 ) then
      YLDvec(1,12+j) = 0.d0
    else
    !print *,'[ST]    YLDvec 305 : ', YLDvec(1,1)
    !print *,'[ST]    YLDvec 305 : ', YLDvec(1,2)
    !print *,'[ST]    YLDvec 305 : ', YLDvec(1,3)
    !print *,'[ST]    YLDvec 305 : ', YLDvec(1,4)

    !                                               herd multiplicative
      YLDvec(1,12+j) = partrec(1,1) + meanyld(j,length,lacn) * hratio(j)
    end if

!                                               Single-trait reliability
!   Note that REL is currently calculated for a 305-d lactation. We may need to
!   look at this later b/c REL for long lactations could be lower than a 305-d
!   lactation based on the same number of TD.

    RELyld(j) = vari305(1,1)/(stdvar(j,lacn)*varfac(j)**2)
    regrel(j) = RELyld(j)

!                                        Compute single-trait estimates,
!                                        correlations, expansion factors
    if(nbump > 0) then
      Fvalue = Fval(nbmp,min(size,10))
      if(bump(j) > Fvalue) then
!                                          Reduce REL of unusual records
        RELyld(j) = RELyld(j) * Fvalue / bump(j)
        DCRvec(j) = DCRvec(j) * Fvalue / bump(j)
        if(bump(j)/Fvalue > maxbump) then
!                                                  Print unusual records
          print *,'------------------------------------------------'
          print *,trt(j),' Regular REL =',regrel(j)
          print *,trt(j),' Reduced REL =',RELyld(j)
          print *,'  Unusual data for cow ',cowid,'   F=',bump(j)
          if ( maxprnt > 0 ) maxprnt = max(maxprnt,ncall)
          maxbump = bump(j)/Fvalue
        end if
      end if
    end if
!                                     Set DCR = 100 for monthly testing
    DCRvec(j) = 100. * RELyld(j) / rmonth(j,lacn)
!                                         Weights and expansion factors
    xsingl(j) = (1.D0 / regrel(j) - 1.)*expamt + 1.
    lacwt(j) = (1. - repty(j)) / (xsingl(j) - repty(j))
    PERSvec(j) = pers(1,1)
    RELpers(j) = dvari(1,1) / varfac(j)**2
!                                                  Yield, pers R matrix
    Rvar(1,1) = (1./regrel(j) - repty(j)) &
      *stdvar(j,lacn)*varfac(j)**2
    Rvar(2,2) = (1./RELpers(j) - reptp(j))*varfac(j)**2
    Rcorr = qCVC1(1,1) / dsqrt(vari(1,1)*dvari(1,1))
    Rvar(2,1) = Rcorr * dsqrt(Rvar(1,1)*Rvar(2,2))
!                                             Force R to be positive def
    if(Rcorr> .99) Rvar(2,1)= .99*dsqrt(Rvar(1,1)*Rvar(2,2))
    if(Rcorr<-.99) Rvar(2,1)=-.99*dsqrt(Rvar(1,1)*Rvar(2,2))
    if(j == 1) Rvecm = 1.d0 / Rvar(1,1)
    if(j == 2) Rvecf = 1.d0 / Rvar(1,1)
    if(j == 3) Rvecp = 1.d0 / Rvar(1,1)
    if(j == 4) Rvecs = 1.d0 / Rvar(1,1)
    Rvar(1,2) = Rvar(2,1)
!                                                             R inverse
    call invrt2(Rvar,2,2,cowid)
    n = 0
    do l = 1,2
      do m = l,2
        n = n + 1
        Rvec(n) = Rvar(l,m)
      end do
    end do
55  continue
else
  print *, '[ERROR]: Unrecognized value of mtrait (1|3|4): ', mtrait
end if
! .................................. OUTPUT DCR AND ACTUAL OR ME YIELDS
!                                       DCRc is mean of fat,protein DCR
!                                         De-adjust for age, parity, 3X
!                                                  Actual 305-d records
  if ( DEBUGmsgs > 1 ) then
    print *,'Correlation among 305-d yields for ', lacn, ' parity'
    print *, corr305
  end if

    do j=1,4
      YLDvec(2,j) = YLDvec(1,j) / (agefac(j)*fact3X(j))
      YLDvec(2,j+4) = YLDvec(1,j+4) / (agefac(j)*fact3X(j))
      YLDvec(2,j+8) = YLDvec(1,j+8) / (agefac(j)*fact3X(j))
      herd305(2,j) = herd305(1,j) / (agefac(j)*fact3X(j))
      herd365(2,j) = herd365(1,j) / (agefac(j)*fact3X(j))
      herd999(2,j) = herd999(1,j) / (agefac(j)*fact3X(j))
!                                                       partial records
      YLDvec(2,j+12) = YLDvec(1,j+12) / (agefac(j)*part3X(j))
    end do
!  end if

  ! All internal calculations are done in kg, so we only need to convert
  ! back to pounds if the user wants pounds.
  if ( UNITSout .eq. 'P' ) then
    do j=1,3
      YLDvec(:,j) = YLDvec(:,j) * lb
      YLDvec(:,j+4) = YLDvec(:,j+4) * lb
      YLDvec(:,j+8) = YLDvec(:,j+8) * lb
      YLDvec(:,j+12) = YLDvec(:,j+12) * lb
      herd305(:,j) = herd305(:,j) * lb
      herd365(:,j) = herd365(:,j) * lb
      herd999(:,j) = herd999(:,j) * lb
    end do
  end if
!                                         Lactation SCS is sum / length
  YLDvec(:,4) = YLDvec(:,4)/305.
  YLDvec(:,8) = YLDvec(:,8)/365.
  YLDvec(:,12) = YLDvec(:,12)/laclen
  if ( length > 0 ) then
    YLDvec(:,16) = YLDvec(:,16)/length
  else
    YLDvec(:,16) = 0.
  end if
  herd305(:,4) = herd305(:,4) / 305.
  herd365(:,4) = herd365(:,4) / 365.
  herd999(:,4) = herd999(:,4) / laclen

  if ( mtrait == 4 ) then
    MTorST = 'MT'
  else if ( mtrait == 1 ) then
    MTorST = 'ST'
  else if ( mtrait == 3 ) then
    MTorST(1:3) = 'MT'
    MTorST(4) = 'ST'
  else
    print *, '[ERROR]: Unrecognized value of mtrait (1|3|4): ', mtrait
  end if

  DCRm = DCRvec(1)
  DCRc = (DCRvec(2) + DCRvec(3))/2.
  DCRs = DCRvec(4)
  if ( ntests(3) == 0 ) DCRc = DCRvec(2)
!                                            Expanded yield and persist
  do j = 1,4
    if ( mtrait > 1 ) then
      do l = 1,mtrait
        Yvec(1,j) = herd305(1,j) + (YLDvec(1,j) - herd305(1,j)) / regrel(j)
        Yvec(2,j) = herd305(2,j) + (YLDvec(2,j) - herd305(2,j)) / regrel(j)
      end do
      if ( mtrait == 3 .and. j == 4 ) then
        Yvec(1,j) = PERSvec(j) / RELpers(j)
        Yvec(2,j) = Yvec(1,j)
      end if
    else
      !!! Added to deal with lactations with no component samples
      if ( RELpers(j) == 0 ) then
        Yvec(1,j) = PERSvec(j)
        Yvec(2,j) = Yvec(1,j)
      else
        Yvec(1,j) = PERSvec(j) / RELpers(j)
        Yvec(2,j) = Yvec(1,j)
      end if
    end if
  end do
  Yvec(1,4) = Yvec(1,4) / 305.
  Yvec(2,4) = Yvec(2,4) / 305.
!
  write(BLUPout(5),60) (Rvec(n),n=1,3)
  write(BLUPout(6),60) (Rvec(n),n=4,6)
  write(BLUPout(7),60) (Rvec(n),n=7,9)
  write(BLUPout(8),60) Rvecm,zero0,zero0
  write(BLUPout(9),60) Rvecf,zero0,zero0
  write(BLUPout(10),60) Rvecp,zero0,zero0
 60   format(3a8)
!                                                length of output record
  do k = 5,10
    BLUPn(k) = 24
  end do
!  print *, '[bestpred]: grafall: ', grafall(1,1:305,1)
!                                                      contemporary mean
  do j = 1, 4
    tempSTD(1,j,:) = 0.
    tempSTD(1,j,:) = grafall(1,:,j) + std(1,j,:)
    !where(std(1,j,:) > -999.0) tempSTD(1,j,:) = grafall(1,:,j) + std(1,j,:)
  end do
!  print *, tempSTD(1,1,1:15)
!  print *, std(1,1,1:15)

!  print *, '[bestpred]: tempSTD: ', tempSTD(1,1,1:305)

!!!
!!! These all need to be checked carefully again. AdV has reported that the actual TD yield is being
!!! written to the ME YD yield column, and plots of the DCRexamples curves made using the Python
!!! program appear to have individual cow BP curves that are too large (shifted up the Y axis).
!!!

  ! De-adjust for age and 3X
  do j = 1,4
    where(std(1,j,:) > -999.0) std(2,j,:) = std(1,j,:) / (agefac(j)*fact3X(j))
    where(grafall(1,:,j) > -999.0) grafall(2,:,j) = grafall(1,:,j) / (agefac(j)*fact3X(j))
    where(tempSTD(1,j,:) > -999.0) tempSTD(2,j,:) = tempSTD(1,j,:) / (agefac(j)*fact3X(j))
    where(DAILYbp(1,j,:) > -999.0) DAILYbp(2,j,:) = DAILYbp(1,j,:) / (agefac(j)*fact3X(j))
    where(yield(1,j,:) > -999.0) yield(2,j,:) = yield(1,j,:) / (agefac(j)*fact3X(j))
  end do

  ! All internal calculations are done in kg, so we only need to convert
  ! back to pounds if the user wants pounds. Note that herd305 is adjusted
  ! upstream and SHOULD NOT also be adjusted here,
  if ( UNITSout .eq. 'P' ) then
    where(std(:,1:3,:) > -999.0) std(:,1:3,:) = std(:,1:3,:) * lb
    where(grafall(:,:,1:3) > -999.0) grafall(:,:,1:3) = grafall(:,:,1:3) * lb
    where(tempSTD(:,1:3,:) > -999.0) tempSTD(:,1:3,:) = tempSTD(:,1:3,:) * lb
    where(DAILYbp(:,1:3,:) > -999.0) DAILYbp(:,1:3,:) = DAILYbp(:,1:3,:) * lb
    where(yield(:,1:3,:) > -999.0) yield(:,1:3,:) = yield(:,1:3,:) * lb
    where(tempTD(1:3,:) > -999.0) tempTD(1:3,:) = tempTD(1:3,:) * lb
!    where(tempDEV(:,1:3,:) > -999.0) tempDEV(:,1:3,:) = tempDEV(:,1:3,:) * lb
    herdavg(:,1:3) = herdavg(:,1:3) * lb
  end if
!  print *, "MEs     :", YLDvec(1,:)
!  print *, "Actuals :", YLDvec(2,:)

  ! Write yield and curve data output.
  if ( WRITEcurve == 1 ) then
    call write_curve_data(mtrait, CURVEfile, SHORTtrait, cowid, STorMTlabel, &
    maxlen, DAILYbp, tempTD, std, parity, GRAFplot, CURVEsmall, CURVEsingle)
  end if
  if ( WRITEdata == 1 ) then
    !print *, ntests
    call write_yield_data(DATAfile,cowid,STorMTlabel,trt, ntests,  &
    X100,YLDvec,herd305,DCRvec,MTorST,PERSvec,RELpers,parity, CURVEsingle)
  end if

!  print *, '[bestpred]: herd305: ', herd305

  !if ( ncall > max(maxprnt,maxshow) ) go to 90
! ------------------------------------------ SECTION FOR GRAPHIC DISPLAY
!if ( mtrait > 1 ) then
    mstart = 1
    if(usetd > 0) then
!                                        Get 3X factors across lactation
      do i=1,usetd
        mstop = dimmilk(i)
        milk3X = test3X(1,dimvec(dimmilk(i)))
        !do m=mstart,mstop,6
        if ( plotfreq > 0 ) then
          do m=mstart,mstop,plotfreq
            freq3X((m+5)/plotfreq) = milk3X
          end do
        end if
        mstart = mstop + 1
      end do
    end if
!                                     Capital letters for standard tests
    do 68 i=1,ntd
      if(dim(i) <   1) go to 68
      if(dim(i) > maxlen) go to 68
      if(MRD(i) > 1) then
        plot(i) = 'L'
        if(super(i) == 2) plot(i) = 'U'
        if(super(i) == 6) plot(i) = 'U'
      else
        plot(i) = 'S'
        if(weigh(i) < Xmilk(i)) plot(i) = 'B'
        if(weigh(i)*2 == Xmilk(i)) plot(i) = 'A'
        if(weigh(i)*2 < Xmilk(i)) plot(i) = 'C'
        if(super(i) == 2 .or. super(i) == 6) then
          plot(i) = 'O'
          if(weigh(i) < Xmilk(i)) plot(i) = 'Q'
          if(weigh(i)*2 == Xmilk(i)) plot(i) = 'P'
          if(weigh(i)*2 < Xmilk(i)) plot(i) = 'R'
        end if
      end if
      if(super(i) > 7) plot(i) = 'V'
      if(super(i) == 0) plot(i) = 'X'
      if(super(i) == 4) plot(i) = 'X'
      if(dimvec(dim(i)) .ne. i) plot(i) = 'D'
!                                    Small letters for milk-only tests
      if(yield(1,2,i) <= zero) then
        k = index(capital,plot(i))
        plot(i) = small(k:k)
      end if
 68     continue
! .......................... Graph lactation curve for cow and contemps
  do i=1,maxlen/24
    fatpct(i) = '    '
    propct(i) = '    '
    scs(i)    = '    '
  end do
  !!! This if...then doesn't work correctly with g95. No error is printed,
  !!! but it evaluates to TRUE even if maxprnt = 0.
  !
  ! GRAFplot(l) is used as an index into several structures in the following loop,
  ! and is used to determine whether mature equivalent values (1) or actual values
  ! (2) should be plotted, if plots are requested at all. I considered adding a
  ! value of 3 that would produce both plots in a single run, but decided that I
  ! can't come up with a good readon for doing that.
  !
!  print *, 'ncall     :', ncall
!  print *, 'maxprnt   : ', maxprnt
!  print *, 'maxshow   : ', maxshow
!  print *, 'max(mp,ms): ', max(maxprnt,maxshow)
  if ( ncall <= maxprnt ) then
    do l = 1,4
      !if ( GRAFplot(l) > 0 .and. plotfreq > 0 ) then
      if ( GRAFplot(l) > 0 .and. maxshow > 0 .and. ncall <= maxshow ) then
        ! Do a little error-trapping here in case bestpred() is called by some
        ! third-party code that doesn't check GRAFplot carefully. If an invalid
        ! value (not 1 or 2) is provided then default to 2 (actual values).
        if ( GRAFplot(l) > 2 ) then
          if ( LOGon == 1 .and. LOGfreq > 0 ) then
            write(int2str, FMT='(I5)') l
            LogMessage = 'GRAFplot ('//int2str//') has an invalid value of '
            write(int2str, FMT='(I5)') GRAFplot(l)
            LogMessage = LogMessage//int2str//'.'//'Using the default of 2.'
            call log_message(LOGfile,LogMessage)
          end if
          GRAFplot(l) = 2
        end if
        ! stdslice and tempSTDslice are used to store the vector of yields for the
        ! current trait because the minval and maxval intrinsics cannot operate
        ! on array slices. Fat, protein, and SCS are multiplied by 10 so that they
        ! plot reasonably well without a major overhaul of the drawing routine.
        stdslice = 0.
        stdslice(1,:) = std(GRAFplot(l),l,:)
        tempSTDslice = 0.
        tempSTDslice(1,:) = tempSTD(GRAFplot(l),l,:)
        if ( l > 1 ) then
          stdslice = stdslice * 10
          tempSTDslice = tempSTDslice * 10
        end if
        ! Hack to deal with very large SCS standard curve values for first 5 d of lactation.
        if ( l == 4 ) then
          stdslice(1,1:5) = 0.
          tempSTDslice(1,1:5) = 0.
        end if
        !minyield(l) = min(minval(stdslice,mask=stdslice>0),minval(tempSTDslice,mask=tempSTDslice>0))
        minyield(l) = 0
        maxyield(l) = max(maxval(stdslice,mask=stdslice>0),maxval(tempSTDslice,mask=tempSTDslice>0))
        do i=1,ntd
          if ( l == 1 ) then
            !if (yield(GRAFplot(l),l,i) < minyield(l) .and. yield(GRAFplot(l),l,i) > 0) &
            !  minyield(l) = yield(GRAFplot(l),l,i)
            if (yield(GRAFplot(l),l,i) > maxyield(l)) maxyield(l) = yield(GRAFplot(l),l,i)
            ! 04/21/2015 -- Something's still screwy hardcode for now.
            if ( maxyield(l) > 150 ) maxyield(l) = 150
          else
            !if ( yield(GRAFplot(l),l,i)*10 < minyield(l) .and. yield(GRAFplot(l),l,i) > 0 ) &
            !  minyield(l) = yield(1,l,i)*10
            if ( yield(GRAFplot(l),l,i)*10 > maxyield(l) ) maxyield(l) = yield(GRAFplot(l),l,i)*10
          end if
        end do ! i
        print *,'-------------------------------------------------------' &
              ,'---------'
        if ( GRAFplot(l) == 1 ) then
          print 70,GRAFname(l),'(ME)',cowid,fresh(1:4),fresh(5:6),fresh(7:8)
        else
          print 70,GRAFname(l),'(Actual)',cowid,fresh(1:4),fresh(5:6),fresh(7:8)
        end if
70       format(a4,' lactation curve ',a,' for cow ',a17,',  fresh ',a4,2(1x,a2))
!                                      Plot predicted ---  and mean ...
        do 73 j = maxyield(l)/YIELDunits(l), minyield(l)/YIELDunits(l), -YIELDsteps(l)
          do 71 i = 1,maxlen/plotfreq
            plt = ' '
            if ( i == maxlen/plotfreq ) plt = ':'
            if ( i == 305/plotfreq ) plt = '.'
            if ( i == length/plotfreq ) plt = '|'
            if ( l > 1 ) then
              testday = std(GRAFplot(l),l,i*plotfreq)*10.0
            else
              testday = std(GRAFplot(l),l,i*plotfreq)
            end if
            if ( testday/YIELDunits(l) == j ) plt = '.'
            if ( l > 1 ) then
              testday = tempSTD(GRAFplot(l),l,i*plotfreq)*10.0
            else
              testday = tempSTD(GRAFplot(l),l,i*plotfreq)
            end if
            if ( testday / YIELDunits(l) == j ) plt = '-'
 71         tests(i) = plt
!                                Plot tests by Super, Owner, Ampm, LER
          do 72 i = 1,ntd
            if (dim(i) <   1) go to 72
            if (yield(GRAFplot(l),2,i) > zero) then
              write(percent,'(f4.1)') yield(GRAFplot(l),2,i)*100 / yield(GRAFplot(l),1,i)
              fatpct(dim(i)/24 + 1) = percent
              write(percent,'(f4.1)') yield(GRAFplot(l),3,i)*100 / yield(GRAFplot(l),1,i)
              propct(dim(i)/24 + 1) = percent
            end if
            if (yield(GRAFplot(l),4,i) > zero) write(percent,'(f4.1)') yield(GRAFplot(l),4,i)
            if (yield(GRAFplot(l),4,i) > zero) scs(dim(i)/24 + 1) = percent
            if ( l > 1 ) then
              testday = yield(GRAFplot(l),l,i)*10
            else
              testday = yield(GRAFplot(l),l,i)
            end if
!                                         Plot at center of days in mean
            middle(i) = dim(i) - (MRD(i)-1)/2
            if ( testday/YIELDunits(l) == j ) then
              !tests(middle(i)/plotfreq) = plot(i)
              if ( middle(i) / plotfreq <= 0 ) then
                tests(middle(i)) = plot(i)
              else if ( middle(i) / plotfreq <= 60 ) then
                tests(middle(i)/plotfreq) = plot(i)
              else
                continue
              end if
            end if
 72       continue
        if ( l == 1 ) then
          print 741,j*YIELDunits(l),tests
 741      format(i4,1x,999a1)
        else
          print 742,float(j*YIELDunits(l))/10.,tests
 742      format(f3.1,1x,999a1)
        end if
 73     continue
        print '(99i5)',(i*30,i=0,maxlen/30)
 75     continue
        print 78,'Fat%',(fatpct(i),i=1,maxlen/24)
        print 78,'Pro%',(propct(i),i=1,maxlen/24)
        print 78,'SCS ',(scs(i)   ,i=1,maxlen/24)
 78     format(1x,a4,99a4)
      end if
    end do

!!!
! Gnarly debugging
!!!
!    print *, 'YLDvec : ', YLDvec
!    print *, 'PERSvec: ', PERSvec
!!!

    print *, ""
    print 790, cowid, fresh(1:4), fresh(5:6), fresh(7:8), parity
 790  format(' Lactation summary for cow ',a17,',  fresh ', a4, 2(1x,a2), ', parity ', i2 )

    print 791,length, laclen
 791  format(/,52x,                 'Cntmp     Adjust       Persistency     Data' &
        ,/,16x,i3,        '-d    305-d    365-d    ',i3,'-d    305-d     Factors    Esti-   Relia-    Coll' &
        ,'ectn',/                                                                           &
        ,' TRAIT Tests    Yield    Yield    Yield    Yield     Mean   Age,     3X  mate    bility    Rating')
    ! Mature equivalents
    do l=1,4
      print 851,trt(l+4),ntests(l), &
        YLDvec(1,l+12)*X100(l),YLDvec(1,l)*X100(l), &
        YLDvec(1,l+4)*X100(l),YLDvec(1,l+8)*X100(l), &
        herd305(1,l)*X100(l),agefac(l),fact3X(l), &
        PERSvec(l),RELpers(l)*100.,DCRvec(l),MTorST(l)
    end do
    ! Actuals
    do l=1,4
      print 851,trt(l),ntests(l), &
        YLDvec(2,l+12)*X100(l),YLDvec(2,l)*X100(l), &
        YLDvec(2,l+4)*X100(l),YLDvec(2,l+8)*X100(l), &
        herd305(2,l)*X100(l),agefac(l),fact3X(l), &
        PERSvec(l),RELpers(l)*100.,DCRvec(l),MTorST(l)
    end do
 851 format(1x,a7,i4,5f9.0,2f7.2,f7.2,f7.0,'% ',f7.0,1x,a2)
    print *,' '
    print *,' '
  end if ! End of IF statement to display plots
  if ( ncall <= maxshow .and. ncall > maxprnt .and. maxshow > 0 ) then

    print *, ""
    print 792, cowid, fresh(1:4), fresh(5:6), fresh(7:8), parity
 792  format(' Lactation summary for cow ',a17,',  fresh ', a4, 2(1x,a2), ', parity ', i2 )


    print 793,length, laclen
 793  format(/,52x,                 'Cntmp     Adjust       Persistency     Data' &
        ,/,16x,i3,        '-d    305-d    365-d    ',i3,'-d    305-d     Factors    Esti-   Relia-    Coll' &
        ,'ectn',/                                                                           &
        ,' Trait Tests    Yield    Yield    Yield    Yield     Mean   Age,     3X  mate    bility    Rating')
    ! Mature equivalents
    do j=1,4
      print 852,trt(j+4),ntests(j), &
        YLDvec(1,j+12)*X100(j),YLDvec(1,j)*X100(j), &
        YLDvec(1,j+4)*X100(j),YLDvec(1,j+8)*X100(j), &
        herd305(1,j)*X100(j),agefac(j),fact3X(j), &
        PERSvec(j),RELpers(j)*100.,DCRvec(j),MTorST(j)
    end do
    ! Actuals
    do j=1,4
      print 852,trt(j),ntests(j), &
        YLDvec(2,j+12)*X100(j),YLDvec(2,j)*X100(j), &
        YLDvec(2,j+4)*X100(j),YLDvec(2,j+8)*X100(j), &
        herd305(2,j)*X100(j),agefac(j),fact3X(j), &
        PERSvec(j),RELpers(j)*100.,DCRvec(j),MTorST(j)
    end do
 852 format(1x,a7,i4,5f9.0,2f7.2,f7.2,f7.0,'% ',f7.0,1x,a2)
    print *,' '
    print *,' '
  end if
 90    do 91 i=1,maxlen
 91      dimvec(i) = 0
!                                      store predicted curves
    do 95 i=1,maxlen/plotfreq
      i6 = i*dailyfreq
      do 92 j=1,ntd
        if(dim(j) <   1) go to 92
        if(dim(j) > maxlen) go to 92
!                                     suppress predictions for test days
 92        continue
! ........... i6= days in milk, grafall= daily yield, std= contemp yield
 95      continue

!  print *, '[bestpred.f90]: YLDvec: ', YLDvec

  if ( maxprnt < 0 ) then
    print *,'======================================================'
    print *,'Outputs from bestpred for cow=',cowid,' fresh=',fresh
    print *,'       tstdays=',tstdays,' parity=',parity, &
                     ' length=',length, ' maxprnt=',maxprnt
    print *, '305d ME yields    : ', YLDvec(1,1:4)
    print *, '305d Actual yields: ', YLDvec(2,1:4)
    print *, 'DCR               : ', DCRvec
    print *, 'Persistencies     : ', PERSvec
    print *, 'Persist REL       : ', RELpers
    print *,'======================================================'
  end if

!  print *, "ME yields    :", YLDvec(1,:)
!  print *, "Actual yields:", YLDvec(2,:)
!  print *, "DCR M        :", DCRm
!  print *, "DCR C        :", DCRc
!  print *, "DCR S        :", DCRs

  if ( DEBUGmsgs > 1 ) print *,'Deallocating memory'
  deallocate(curve)
  deallocate(curvem)
  deallocate(curves)
  deallocate(graph1)
  deallocate(grafall)
  deallocate(freq3X)
  deallocate(std)
  deallocate(stdslice)
  deallocate(tempTD)
  deallocate(tempSTD)
  deallocate(tempSTDslice)
!  deallocate(tempDEV)
  deallocate(fatpct)
  deallocate(propct)
  deallocate(scs)
  deallocate(tests)
  deallocate(dimvec)
  return
end subroutine bestpred
!--------------------------------------------------------------------
  function vary(dim1,trait1,super1,Xmilk1,sample1,MRD1, &
                dim2,trait2,super2,Xmilk2,sample2,MRD2, &
                sd,lacn,last,maxlen)
!
!                       Compute covariance of any two test day yields
!
  integer trait1,trait2,dim1,dim2,super1,super2,Xmilk1,Xmilk2, &
      sample1,sample2,MRD1,MRD2,i,j,ibegin,iend, &
      jbegin,jend,lacn,last,maxlen
  real*8  vary,corr,sd(4,maxlen,last),dimdif,Idiag
!                                    Assign DCR for supervision codes
!                            0   1   2   3   4   5   6   7  8  9
  real, save :: dcr(0:9)=(/  0., 1.,.77,.97, 1., 1.,.77,.97,1.,1./)
  real err1,err2
  real sqc1,sqc2
!***                     milk fat  prot scs   phenotypic correlations
  real, save :: mfpcorr(4,4,2)=reshape((/ &
                      1., .67, .85,-.08,  &
                     .67,  1., .77,-.14,  &
                     .85, .77,  1.,-.10,  &
                    -.08,-.14,-.10,  1.,  &
!                                                    later lactations
                      1., .67, .85,-.08,  &
                     .67,  1., .77,-.14,  &
                     .85, .77,  1.,-.10,  &
                    -.08,-.14,-.10,  1./),(/4,4,2/))
!               former   milk fat  prot scs   phenotypic correlations
! real, save :: mfpcorr(4,4,2)=reshape((/ &
!                     1., .67, .91,-.05, &
!                    .67,  1., .74,-.08, &
!                    .91, .74,  1.,-.02, &
!                   -.05,-.08,-.02,  1., &
!                                                    later lactations
!                     1., .67, .88,-.12, &
!                    .67,  1., .73,-.17, &
!                    .88, .73,  1.,-.11, &
!                   -.12,-.17,-.11,  1./),(/4,4,2/))
  real   Icorr,Vmid,Vlast,Vscs
  vary = 0.d0

!                                          Set MRD range for LER tests
  ibegin = dim1 - MRD1 + 1
  iend = dim1
  jbegin = dim2 - MRD2 + 1
  jend = dim2
!                                    Use center day for fat, prot, SCS
  if ( trait1 > 1 ) then
    ibegin = dim1 - (MRD1 - 1)/2
    iend = ibegin
  end if
  if ( trait2 > 1 ) then
    jbegin = dim2 - (MRD2 - 1)/2
    jend = jbegin
  end if
!                                   Compute test day (co)variance or
!                                       mean (co)variance if MRD > 1
  do 40 i=ibegin,iend
    do 40 j=jbegin,jend
!                                        Average and difference of DIM
      dimdif = iabs(i - j)
!         avgdim = (i + j)/2.d0
!         Icorr= 0.
!                                        Positive definite submatrices
!
      Vmid = (i - i**2/365.)*(j - j**2/365.)/91.25**2
      Vlast= i*j / 365.**2
      Vscs = .995**dimdif
!!!                                         Restrictions for DIM > 365
!      if(i > 365) Vmid = 0.
!      if(j > 365) Vmid = 0.
!      if(i > 365) Vlast = j / 365.
!      if(j > 365) Vlast = i / 365.
      if ( trait1 > 3 .or. trait2 > 3 ) Vscs = 0.
      if ( i == j ) then
        Icorr = mfpcorr(trait1,trait2,lacn)
        Vmid  = 1.
        Vlast = 1.
        if ( trait1 == trait2) Vscs = 1.
        Idiag = 1.
      else
        Idiag = 0.
      end if
!                                           Positive definite function
!      if ( lacn == 1 ) &
!          corr = (.665 * .995**dimdif &
!                + .256 * Vlast        &
!                + .085 * Vscs)        &
!                * mfpcorr(trait1,trait2,lacn)

!!! New 09/05/2007
      if ( lacn == 1 ) then
        if ( trait1 < 4 .and. trait2 < 4 ) then
!          corr = (0.210*Idiag            &
!                 + 0.790*0.998**dimdif ) &
!                 * mfpcorr(trait1,trait2,lacn)
!!! New 12/14/2007
          corr = ( 0.214*Idiag           &
                 + 0.786*0.998**dimdif ) &
                 * mfpcorr(trait1,trait2,lacn)
        else if ( trait1 == 4 .and. trait2 == 4 ) then
          corr = ( 0.199*Idiag           &
                 + 0.801*0.998**dimdif ) &
                 * mfpcorr(trait1,trait2,lacn)
        else
          !!! Calculate correlation between S and other trait
          sqc1 = sqrt(0.214*Idiag + 0.786*0.998**dimdif)
          sqc2 = sqrt(0.199*Idiag + 0.801*0.998**dimdif)
          corr = sqc1 * sqc2 * mfpcorr(trait1,trait2,lacn)
        end if
      end if

!                                                     Later lactations
!       if(lacn > 1) &
!          corr = (.737 * .992**dimdif &
!                + .158 * Vmid         &
!                + .167 * Vlast)       &
!                * mfpcorr(trait1,trait2,lacn)

!!! New 09/05/2007
      if ( lacn > 1 ) then
        if ( trait1 < 4 .and. trait2 < 4 ) then
          !!! Calculate correlations among DIM and/or
          !!! traits for M,F,P
!          corr = (0.557                  &
!                 + 0.098*Idiag           &
!                 + 0.345*0.987**dimdif ) &
!                 * mfpcorr(trait1,trait2,lacn)
          corr = ( 0.132*Idiag           &
                 + 0.868*0.997**dimdif ) &
                 * mfpcorr(trait1,trait2,lacn)
        else if ( trait1 == 4 .and. trait2 == 4 ) then
          !!! Calculate correlations among DIM for SCS
          corr = ( 0.199*Idiag            &
                 + 0.801*0.998**dimdif ) &
                 * mfpcorr(trait1,trait2,lacn)
        else
          !!! Calculate correlation between S and other trait
          sqc1 = sqrt(0.132*Idiag + 0.868*0.997**dimdif)
          sqc2 = sqrt(0.199*Idiag + 0.801*0.998**dimdif)
          corr = sqc1 * sqc2 * mfpcorr(trait1,trait2,lacn)
        end if
      end if

!!!                                         Restrictions for DIM > 365
!       if(lacn == 1 .and. maxlen > 365) corr = .85 * .997**dimdif &
!         *mfpcorr(trait1,trait2,lacn)
!       if(lacn > 1 .and. maxlen > 365) corr = .85 * .995**dimdif &
!         *mfpcorr(trait1,trait2,lacn)
!      if(trait1 == trait2 .and. i == j) &
!        corr = mfpcorr(trait1,trait2,lacn)
      if ( trait1 == trait2 .and. i == j ) corr = 1.0
!
!                                      Old formulas from Canadian data
!                                          (1996 Norman et al abstract)
!         if(i > 20 .and. j > 20) then
!           corr = .807 + .00105*avgdim - .00000383*avgdim*avgdim
!    *           - .00150*dimdif - .00000119*dimdif*dimdif
!    *           - .0000104*avgdim*dimdif
!         else
!                                      Separate function for DIM < 21
!           corr = .384 + .0412*avgdim - .000778*avgdim*avgdim
!    *           - .0212*dimdif - .000140*dimdif*dimdif
!    *           + .000664*avgdim*dimdif
!         end if
!         corr = 1.0 + .3*(1.*Xmilk1/max(sample1,sample2) - 1.)
!         corr = corr * mfpcorr(trait1,trait2,lacn)
!         if(corr > 1.0) corr = 1.0
!
!                                           Increase variance for AMPM
!                                and covariance for traits on same day
      if ( i == j ) then
        corr = corr + .3 * mfpcorr(trait1,trait2,lacn) * &
                         (1.*Xmilk1/max(sample1,sample2) - 1.)
!                                          Large diagonals if UNUSABLE
        if ( trait1 == trait2 ) then
          if ( dcr(super1) == 0.) corr = 1000000.
        end if
      end if
!                              Increase (co)variance for OWNER-SAMPLER
      err1 = 0.
      err2 = 0.
      if(dcr(super1) > 0.) err1 = (1./dcr(super1) - 1.)*.62
      if(dcr(super2) > 0.) err2 = (1./dcr(super2) - 1.)*.62
!
      corr = corr + sqrt(err1*err2) *mfpcorr(trait1,trait2,lacn)
!                                            former method (published)
!         if(super1 == 2 .and. super2 == 2) corr = corr + .18
!
!                          Set correlations of MILK, FAT, PROTEIN, SCS
!                                      Get COV matrix from CORR matrix
!
 40       vary = vary + corr*sd(trait1,i,lacn)*sd(trait2,j,lacn)
!
!                                                Adjust to DAILY basis
  vary = vary/(iend - ibegin + 1)/(jend - jbegin + 1)
  return
end function vary
!--------------------------------------------------------------------
  function covary(trt305,trait,dim,laclen,maxlen,MRD,covari,covd &
                    ,sd,lacn,last,precise,stat)
!                                                Compute covariance of
!                                         test day and lactation yield
  integer trt305,trait,dim,laclen,maxlen,MRD &
        ,lacn,last,precise,stat &
        ,i,j,k,m,mbegin,mend,prec
  real*8  vary,covary,sum,sumd,var,sd(4,maxlen,last) &
        ,covari(4,4,maxlen,maxlen,last),covd(4,4,maxlen,maxlen,last),iplus
  covary = 0.d0
  if(trait == 0) then
!                         Store covariances of any test day observation
!                         with true yield up to day laclen for trait trt305
!                         Store DIM covariances needed for persistancy
    iplus = (precise - 1)/2.
    do 30 j=1,4
      do 30 k=1,maxlen
        sum = 0.d0
        sumd = 0.d0
!                                               Save CPU if precise > 1
        !do 20 i=1,305,precise
        do 20 i=1,maxlen,precise
          var = vary(i,trt305,1,2,2,1,k,j,1,2,2,1,sd,lacn,last,maxlen)
          !prec = min(precise,306-i)
          prec = min(precise,maxlen+1-i)
          sum = sum + var*prec
          sumd = sumd + var*(i+iplus)*prec
          if (i <= 305 .and. k <= 305 .and. trt305 == j) covary = covary + var*prec
          covari(trt305,j,k,i,lacn) = sum
          covd(trt305,j,k,i,lacn) = sumd
 20       continue
!                                    Accumulate variance of 305-d yield
        !if (k < 306 .and. trt305 == j) covary = covary + sum
 30     continue
  else
!                                   Look up covariances stored in table
!                                   or mean of covariances for MRD > 1
    mbegin = dim - MRD + 1
    mend = dim
    if(trait > 1) then
      mbegin = dim - (MRD - 1)/2
      mend = mbegin
    end if
    if ( stat == 1 ) then
      do m = mbegin, mend
        covary = covary + covari(trt305,trait,m,laclen,lacn)
      end do
    else
!                                         Covariances for persistancy
      do m = mbegin, mend
        covary = covary + covd(trt305,trait,m,laclen,lacn)
      end do
    end if
    covary = covary/(mend - mbegin + 1)
  end if
  return
end function covary
!---------------------------------------------------------------------
!                                   Adjust milk, fat, prot, SCS for 3X
subroutine adjust3X &
!                                                         input
          (ntd,dim,length,yearfr,parity,Xmilk,meanyld, &
!                                                          output
           test3X,fact3X,part3X, &
!                                                         control
           use3X,maxtd,last,maxlen)
!                                test3X has factors to adjust each test
!                                fact3X has factors to adjust 305 yield
!                                part3X has factors to adjust part yld
!                                meanyld is cumulative daily yield
  integer :: ntd, dim(maxtd), length, yearfr, parity  &
             , freq(maxlen), use3X, maxtd, last, i, j &
             , Xmilk(maxtd), dimsort(maxtd)           &
             , begin1, begin2, end1, end2, ntd2, yrfr &
             , lacn, lacn3, maxlen
  real   :: year0, year1
  real*8 :: meanyld(4,maxlen,last),test3X(4,maxtd) &
         ,fact3X(4),part3X(4) &
         ,yld3X(4),lac3X(4),a3X(4)
!
!                                           3X adjustment factors from
!                     (milk,fat,prot,scs)   Karaca, 1998 MS Iowa State
  real, save :: new3X(4,3)=reshape((/.12, .09, .10, .00, &
                 .14, .10, .11, .00, &
                 .14, .10, .11, .00/),(/4,3/))
!                                     Old 3X factors from Kendrick,1953
  real, save :: old3X(4,3)=reshape((/.20, .20, .20, .00, &
                   .17, .17, .17, .00, &
                   .15, .15, .15, .00/),(/4,3/))
  real  adj3X(4)
!                                          Parity groups for 3X factors
  lacn = min(parity,last)
  if (lacn <= 0) lacn = last
  lacn3 = min(parity,3)
  if (lacn <= 0) lacn3 = 3
  dimsort = 0
!                                        Calculate amount of adjustment
  do j=1,4
    if (use3X == 0) adj3X(j) = 0.
    if (use3X == 1) adj3X(j) = old3X(j,lacn3)
    if (use3X == 2) adj3X(j) = new3X(j,lacn3)
    if (use3X == 3) then
!                                                      3-year phase-in
      yrfr = min(yearfr,1999)
      yrfr = max(yrfr,1996)
      year0 = float(1999 - yrfr)/(1999 - 1996)
      year1 = 1. - year0
      adj3X(j) = year0*old3X(j,lacn3) + year1*new3X(j,lacn3)
      end if
!                                             Convert to 3X multiplier
    adj3X(j) = 1.0 / (1.0 + adj3X(j))
    fact3X(j) = 1.0
    part3X(j) = 1.0
    lac3X(j) = 0.d0
    yld3X(j) = 0.d0
  end do ! j
  if(ntd == 0) return
  ntd2 = 0
  do i=1,ntd
    if( dim(i) <= maxlen ) then
      dimsort(i) = dim(i)
      freq(dim(i)) = Xmilk(i)
      ntd2 = ntd2 + 1
    end if
!                                             Factors for each test day
    do j=1,4
      test3X(j,i) = 1.0
      if(Xmilk(i) > 2) test3X(j,i) = adj3X(j)
    end do ! j
  end do ! i
!  print *, 'ntd : ', ntd
!  print *, 'ntd2: ', ntd2
  !print *, '[adjust3X]: dimsort before binsort: ', dimsort
  call binsort(dimsort,ntd2,maxtd)
  !print *, '[adjust3X]: dimsort after dimsort: ', dimsort
!                                                  Assume 2X if no data
!  print *, 'dimsort: ', dimsort
  if ( ntd2 == 0 ) return
  if ( dimsort(1) > maxlen ) return
  if ( dimsort(1) <= 0 ) return ! From J. Boyer
  !print *, '[adjust3X]: ntd2', ntd2
  do i=1,ntd2
    end1 = min(dimsort(i),305)
    end2 = min(dimsort(i),length)
    !print *,'[adjust3X]: end1: ', end1
    !print *,'[adjust3X]: end2: ', end2
    !print *,'[adjust3X]: length: ', length
!                                         Expected yields between tests
    do j=1,4
      a3X(j) = 1.0
      !print *, '[adjust3X] i: ', i
      !print *, '[adjust3X]: dimsort: ', dimsort(i)
      !print *, '[adjust3X]: freq: ', freq(dimsort(i))
      if (freq(dimsort(i)) > 2) a3X(j) = adj3X(j)
      if (i>1) lac3X(j) = lac3X(j) - meanyld(j,begin1,lacn)/a3X(j)
      lac3X(j) = lac3X(j) + meanyld(j,end1,lacn)/a3X(j)
      if (i>1) yld3X(j) = yld3X(j) - meanyld(j,begin2,lacn)/a3X(j)
      yld3X(j) = yld3X(j) + meanyld(j,end2,lacn)/a3X(j)
    enddo ! j
    begin1 = end1
    begin2 = end2
  enddo ! i
  do j=1,4
    if (end1 < 305) lac3X(j) = lac3X(j) &
        + (meanyld(j,305,lacn) &
        -  meanyld(j,end1,lacn))/a3X(j)
    if (end1 < length) yld3X(j) = yld3X(j) &
       + (meanyld(j,length,lacn) &
       -  meanyld(j,end2,lacn))/a3X(j)
!                                           Adjustment for 305-d record
      if(lac3X(j) == 0.) then
        fact3X(j) = 1.
      else
        fact3X(j) = meanyld(j,305,lacn)/lac3X(j)
      endif
!                                            Adjustment for part record
      if(yld3X(j) == 0.) then
        part3X(j) = 1.
      else
        part3X(j) = meanyld(j,length,lacn)/yld3X(j)
      endif
  enddo ! j
  do i=1,ntd2
    freq(dim(i)) = 0
  enddo ! i
  return
end subroutine adjust3X
!---------------------------------------------------------------------
!                                            Expected daily yield from
!                                            Standard lactation curves
!
function ymean(trait,dim,maxlen,MRD,dyield,lacn,last,hratio,herd305)
  real*8 ymean,hratio(4),herd305(4),sum
  integer trait,dim,MRD,i,ibegin,iend,lacn,last,maxlen
  real dyield(maxlen,4,last)
  sum = 0.d0

  if ( trait < 1 .or. trait > 4 ) then
    print *, '[ERROR]: An invalid value of trait, ', trait, ', was passed to aipldcr/ymean()!'
    print *, '[ERROR]: Subsequent calculations may be incorrect.'
  else
    !                                   Average over previous days for MRD
    ibegin = dim - MRD + 1
    iend = dim
    if(trait > 1) then
        ibegin = dim - (MRD - 1)/2
        iend = ibegin
        end if
    do i = ibegin,iend
        sum = sum + dyield(i,trait,lacn)
    end do
    !print *, '[ymean]: sum of daily yields for trait ', trait, ': ', sum

    ymean = sum/(iend-ibegin+1)
    !print *, '[ymean]: mean daily yield for trait ', trait, ': ', ymean

    !                                           Adjust for 305-d herd mean
    !                                                       multiplicative
    ymean = ymean*hratio(trait)
    !print *, '[ymean]: adjusted mean daily yield for trait ', trait, ': ', ymean
  end if
  return
end function ymean
!--------------------------------------------------------------------
SUBROUTINE MULT(A,B,C,N1,N2,N3,M1,M2,M3)
  REAL*8 A(M1,M3),B(M1,M2),C(M2,M3),CKJ !,SUM
  INTEGER I,J,K,N1,N2,N3,M1,M2,M3
!                             MATRIX MULTIPLICATION TO GET  A = B * C
  DO 30 J=1,N3
    DO 5 I=1,N1
5         A(I,J)=0.D0
    DO 10 K=1,N2
      CKJ=C(K,J)
      DO 10 I=1,N1
10          A(I,J)=A(I,J) + B(I,K)*CKJ
30        CONTINUE
  RETURN
END SUBROUTINE MULT
!--------------------------------------------------------------------
SUBROUTINE INVRT2(A,IA,N,cowid)
  REAL*8 A(IA,IA),ZERO,SAVE,AIK
  INTEGER I,J,K,N,IA
  character*17 cowid
!             SUBROUTINE TO INVERT A NON-SYMMETRIC MATRIX OF ORDER N.
!             MATRIX MUST BE NON-SINGULAR.  ALGORITHM WORKS BEST FOR
!             POSITIVE DEFINITE MATRICES (NO ROW INTERCHANGES).
!             NOTE :   MATRIX A IS OVERWRITTEN WITH ITS INVERSE
  IF(N.GT.IA) GO TO 201
  ZERO=1.D-12
  DO 80 I=1,N
!                                               CHECK FOR SINGULARITY
     IF(DABS(A(I,I)) < ZERO) THEN
!      A(I,I)=ZERO
       PRINT *,'Matrix not positive definite (diagonal=0) in subr invrt2'
       print *,'CowID = ',cowid
       if(dabs(a(i,i)) <= 0.d0) stop
       END IF
!                                                TRANSFORM THE MATRIX
     SAVE=1.D0/A(I,I)
     DO 50 J=1,N
50         A(J,I)=A(J,I)*SAVE
     A(I,I)=SAVE
     DO 70 K=1,N
       IF(A(I,K) .EQ. 0.D0) GO TO 70
       IF(K.EQ.I)GO TO 70
       AIK=A(I,K)
       DO 60 J=1,N
60           A(J,K)=A(J,K)-A(J,I)*AIK
       A(I,K)=-1.D0*SAVE*AIK
70         CONTINUE
80      CONTINUE
  RETURN
201   PRINT *,'MATRIX DIMENSION IN INVRT2 LARGER THAN DECLARED'
  STOP
END SUBROUTINE INVRT2
! --------------------------------------------------------------------
SUBROUTINE BINSORT(LIST,N,NMAX)
!                                      BINARY SORT OF INTEGER LIST
!                                      USING REPEATED SUBSET MERGES
!                                      ASCENDING SORT ORDER
!                                      USES 2N*LOG2(N) 'IF' STATEMENTS
!                                      USES 2N*LOG2(N) ADDITIONS
  INTEGER LIST(NMAX),N,NMAX,IWK(NMAX),NSRT,I,J,ISTOP2,NSEG
  INTEGER ISTOP1,IP1,IP2,IBIG
!                                              IS LIST ALREADY SORTED?
  IF(N < 2) RETURN
  DO 5 I=2,N
    IF(LIST(I) .LT. LIST(I - 1)) GO TO 6
 5      CONTINUE
  RETURN
 6    NSEG=0
  IBIG=1000000000
  NSRT=1
!                                                      BEGIN SORT LOOP
10    IF(NSRT .GE. N) RETURN
  IP1=1
  IP2=1 + NSRT
  ISTOP1=IP2
  ISTOP2=IP2 + NSRT
  IF(ISTOP2 .GT. N) ISTOP2=N + 1
  DO 40 I=1,N
!                                      COMPARE ELEMENTS OF TWO SUBSETS
    IF(LIST(IP1) .GT. LIST(IP2)) GO TO 25
      IWK(I)=LIST(IP1)
      IP1=IP1 + 1
      IF(IP1 .LT. ISTOP1) GO TO 40
!                                               FIRST SUBSET COMPLETED
        IP1=IP1 - 1
        LIST(IP1)=IBIG
        GO TO 30
25        IWK(I)=LIST(IP2)
      IP2=IP2 + 1
      IF(IP2 .LT. ISTOP2) GO TO 40
!                                              SECOND SUBSET COMPLETED
        IP2=IP2 - 1
        LIST(IP2)=IBIG
30        NSEG=NSEG + 1
      IF(NSEG .LT. 2) GO TO 40
!                                               BOTH SUBSETS COMPLETED
        IP1=IP1 + NSRT + 1
        IP2=IP2 + NSRT + 1
        ISTOP1=IP1 + NSRT
        ISTOP2=IP2 + NSRT
        NSEG=0
        IF(ISTOP2 .GT. N) ISTOP2=N + 1
        IF(IP2 .LE. N) GO TO 40
!                                       FOR NO ELEMENTS IN LAST SUBSET
          DO 35 J=IP1,N
35              IWK(J)=LIST(J)
          GO TO 45
40      CONTINUE
45    DO 50 I=1,N
50      LIST(I)=IWK(I)
  IF(LIST(N) .GE. IBIG) GO TO 100
  NSRT=NSRT*2
  GO TO 10
100   PRINT 101
101   FORMAT('0SOME ELEMENTS TOO BIG, INCREASE IBIG')
  STOP
END SUBROUTINE BINSORT
!--------------------------------------------------------------------
! Valid methods are:
! ALL:  L: Linear interpolation for MFPS
! MFP:  W: Woods curves for MFP
!       C: Season-specific Woods curves for MFP
!       R: Region-specific Woods curves for MFP
!       T: Region- and season-specific Woods curves for MFP
! SCS:  G: Smooth curves (Morant and Gnanasakthy) for SCS
!       D: Season-specific M&G curves for SCS
!       S: Region-specific M&G curves for SCS
!       U: Region- and season-specific M&G curves for SCS
subroutine interpolate(trait, month, dyield, lacn, dyld, sum, meanyld, meanp &
    , sd, dsd, method, DEBUGmsgs, maxlen, breed, region, season)
!                            Interpolate between monthly means and sd
    character :: method
    integer :: month, lacn, i, DEBUGmsgs, maxlen, trait, breed
    integer,parameter :: last = 2
    real :: dyield(maxlen,4,last), dyld(12,4,last), dsd(12,4,last)
    real*8 :: meanyld(4,maxlen,last), sd(4,maxlen,last), meanp(2,4), sum
    real, save :: woods_means(3,3,2,6), woods_sd(3,3,2,6)
    real, save :: mandg_means(4,1,2,6), mandg_sd(4,1,2,6)
    ! You do not have to specify a region -- this is used only for John and Dan's project on regional
    ! adjustments.
    integer :: region, season
    !integer, optional :: region, season
    real, save :: regional_woods_means(3,3,2,7), regional_woods_sd(3,3,2,7)
    real, save :: regional_mandg_means(4,1,2,7), regional_mandg_sd(4,1,2,7)
    real, save :: calving_woods_means(3,3,2,4), calving_woods_sd(3,3,2,4)
    real, save :: calving_mandg_means(4,1,2,4), calving_mandg_sd(4,1,2,4)
    ! 3 parms -> 3 traits -> 2 parities -> four seasons ->seven regions
    real, save :: seasonal_woods_means(3,3,2,4,7), seasonal_woods_sd(3,3,2,4,7)
    ! 4 parms -> 3 traits -> 2 parities -> four seasons ->seven regions
    real, save :: seasonal_mandg_means(4,1,2,4,7), seasonal_mandg_sd(4,1,2,4,7)

    ! All curve parms updated 12/10/2007 by JBC

    data woods_means      /13.7957, 0.1904, 0.00268, & ! Ayrshire
                            0.6496, 0.1217, 0.0016,  &
                            0.513, 0.1244, 0.00152,  &
                           18.195, 0.1997, 0.00421,  &
                            0.9498, 0.1033, 0.00293, &
                            0.7401, 0.1052, 0.00267, &
                           15.0619, 0.1619, 0.00182, & ! Brown Swiss
                            0.6908, 0.1165, 0.00116, &
                            0.4847, 0.1516, 0.00121, &
                           23.3824, 0.13, 0.00274,   &
                            1.1891, 0.0569, 0.00185, &
                            0.8221, 0.0961, 0.00186, &
                           15.0366, 0.146, 0.00214,  & ! Guernsey
                            0.6033, 0.1407, 0.00136, &
                            0.4972, 0.1203, 0.00129, &
                           21.2814, 0.12, 0.00297,   &
                            0.962, 0.0907, 0.00204,  &
                            0.7883, 0.0677, 0.00184, &
                           13.0097, 0.2673, 0.00262, & ! Holstein
                            0.7842, 0.1199, 0.0013,  &
                            0.4625, 0.2033, 0.00161, &
                           22.0087, 0.2155, 0.00357, &
                            1.2874, 0.0731, 0.00213, &
                            0.8539, 0.1317, 0.00232, &
                           11.5338, 0.2022, 0.00222, & ! Jersey
                            0.5027, 0.1954, 0.00161, &
                            0.4096, 0.1797, 0.00148, &
                           17.399, 0.1766, 0.00316,  &
                            0.7902, 0.1634, 0.00249, &
                            0.7146, 0.1231, 0.00219, &
                           13.6246, 0.1877, 0.00231, & ! Milking Shorthorn
                            0.6878, 0.0805, 0.00108, &
                            0.5389, 0.0914, 0.00133, &
                           22.4251, 0.1434, 0.00348, &
                            1.2134, 0.0154, 0.00199, &
                            0.8634, 0.0551, 0.00224  /

!                          Taken from Cole et al. (2007)
    data woods_sd         / 3.8737, 0.1027, 0.000381,  & ! Ayrshire
                            0.2259, 0.034, 9.3e-05,    &
                            0.1229, 0.0875, -4e-05,    &
                            5.8246, 0.0935, 0.00136,   &
                            0.3686, 0.0121, 0.001,     &
                            0.2229, 0.0216, 0.000421,  &
                            4.3308, 0.0905, 0.000127,  & ! Brown Swiss
                            0.2985, 0.0025, -6e-05,    &
                            0.1313, 0.098, -0.00011,   &
                            7.0287, 0.0562, 0.000583,  &
                            0.495, -0.025, 0.000459,   &
                            0.2418, 0.0281, 6e-05,     &
                            4.7574, 0.065, 0.000308,   & ! Guernsey
                            0.3265, -0.0553, -0.00058, &
                            0.1574, 0.035, -0.00023,   &
                            6.5111, 0.0527, 0.000974,  &
                            0.402, -0.0215, 0.000432,  &
                            0.2496, -0.0133, 9.1e-05,  &
                            5.3807, 0.054, -0.0003,    & ! Holstein
                            0.3798, -0.0612, -0.00061, &
                            0.1803, 0.0168, -0.00072,  &
                            8.7545, 0.0282, 0.000439,  &
                            0.5536, -0.0526, 0.000304, &
                            0.327, -0.0363, -0.00022,  &
                            3.8016, 0.074, -0.00012,   & ! Jersey
                            0.2245, 0.045, -0.00023,   &
                            0.1255, 0.0761, -0.0004,   &
                            5.114, 0.0797, 0.000618,   &
                            0.2911, 0.069, 0.000659,   &
                            0.2027, 0.0399, 7.7e-05,   &
                            3.3264, 0.1834, 0.000397,  & ! Milking Shorthorn
                            0.3183, -0.0446, -0.00092, &
                            0.1407, 0.0583, -0.00022,  &
                            7.7328, 0.0241, 3.4e-05,   &
                            0.5617, -0.124, -0.00082,  &
                            0.3494, -0.1302, -0.00139  /

!                          Means and SD for SCS are modelled using curve C$
!                          from Morant and Gnanasankthy (19889).
    data mandg_means         /  1.7911, -0.00233, 2.157e-06, 14.2048, & ! Ayrshire
                                1.7716, -0.00791, -1e-05, 17.8304,    &
                                1.6126, -0.00475, -5.28e-06, 17.8126, & ! Brown Swiss
                                2.2231, -0.00594, -8.51e-06, 13.5807, &
                                2.3295, -0.00316, -4.36e-06, 14.8032, & ! Guernsey
                                2.4706, -0.00426, -5.2e-06, 8.8076,   &
                                1.9798, -0.00344, -3.42e-06, 16.8829, & ! Holstein
                                2.5072, -0.00431, -4.59e-06, 8.9804,  &
                                2.159, -0.00354, -3.68e-06, 19.653,   & ! Jersey
                                2.2129, -0.00572, -8.02e-06, 14.9249, &
                                1.8987, -0.00308, -1.56e-06, 14.1516, & ! Milking Shorthorn
                                1.2378, -0.0115, -3e-05, 20.7716      /

    data mandg_sd            /  1.6442, -0.00055, -3.14e-06, 1.476,   & ! Ayrshire
                                2.1741, 0.00243, 4.193e-06, -1.3917,  &
                                1.9034, 0.00141, 4.422e-06, 0.2612,   & ! Brown Swiss
                                2.2805, 0.00285, 6.958e-06, -2.7067,  &
                                1.647, -0.00064, -1.92e-06, 2.6782,   & ! Guernsey
                                2.243, 0.00131, 1.45e-06, -3.0148,    &
                                1.9551, 0.000248, -1.22e-07, -1.5154, & ! Holstein
                                2.4849, 0.00229, 3.454e-06, -6.3911,  &
                                1.7515, -5e-05, -3.91e-07, 3.1323,    & ! Jersey
                                2.3131, 0.00134, 1.457e-06, -2.5612,  &
                                1.6281, -0.00075, -2.87e-06, 2.2243,  & ! Milking Shorthorn
                                2.4802, 0.00184, 2.61e-06, -1.5926    /

    !!! Structures for John and Dan's regional ajustments project.
    !!! Order of regiona: 1, 2, 3, 4, 5, 6, 7
    !!! Curve parms updated 01/05/2009 by JBC

    data regional_woods_means / 14.3926,  0.2274,  0.00235,  & ! Mideast
                                 0.7313,  0.1248,  0.00135,  &
                                 0.4705,  0.1922,  0.00169,  &
                                22.2123,  0.1879,  0.00320,  &
                                 1.1424,  0.0778,  0.00209,  &
                                 0.8122,  0.1287,  0.00238,  &
                                14.8625,  0.2272,  0.00242,  & ! Midwest
                                 0.9103,  0.0824,  0.00127,  &
                                 0.4949,  0.1838,  0.00167,  &
                                23.4841,  0.1855,  0.00329,  &
                                 1.4331,  0.0378,  0.00206,  &
                                 0.8820,  0.1143,  0.00235,  &
                                12.3776,  0.2689,  0.00253,  & ! Mountain-Prarie
                                 0.7714,  0.1153,  0.00138,  &
                                 0.4534,  0.2028,  0.00173,  &
                                21.9752,  0.2050,  0.00347,  &
                                 1.2546,  0.0663,  0.00222,  &
                                 0.8425,  0.1297,  0.00251,  &
                                15.3781,  0.2134,  0.00236,  & ! Northeast
                                 0.8739,  0.0818,  0.00122,  &
                                 0.5023,  0.1729,  0.00163,  &
                                23.7647,  0.1838,  0.00338,  &
                                 1.3960,  0.0405,  0.00211,  &
                                 0.8706,  0.1159,  0.00244,  &
                                11.5174,  0.3077,  0.00272,  & ! Northwest
                                 0.7262,  0.1452,  0.00155,  &
                                 0.3808,  0.2722,  0.00209,  &
                                23.5710,  0.2103,  0.00338,  &
                                 1.3562,  0.0672,  0.00221,  &
                                 0.9152,  0.1331,  0.00248,  &
                                11.8785,  0.2757,  0.00271,  & ! Southeast
                                 0.6265,  0.1553,  0.00159,  &
                                 0.4371,  0.2141,  0.00190,  &
                                18.9444,  0.2326,  0.00374,  &
                                 0.9578,  0.1113,  0.00236,  &
                                 0.7772,  0.1393,  0.00259,  &
                                13.1788,  0.2591,  0.00245,  & ! Southwest
                                 0.6960,  0.1496,  0.00155,  &
                                 0.4636,  0.2114,  0.00179,  &
                                23.6043,  0.2093,  0.00347,  &
                                 1.2305,  0.0962,  0.00244,  &
                                 0.8819,  0.1450,  0.00264   /

    data regional_woods_sd / 5.9724,  0.0232,  -0.00027,  & ! Mideast
                             0.3922, -0.0707,  -0.00056,  &
                             0.1862,  0.00435, -0.00061,  &
                            10.1147, -0.0304,   0.000135, &
                             0.5498, -0.0690,   0.000234, &
                             0.3541, -0.0777,  -0.00041,  &
                             5.2770,  0.0629,  -0.00002,  & ! Midwest
                             0.4132, -0.0879,  -0.00060,  &
                             0.1744,  0.0239,  -0.00049,  &
                             8.7326,  0.0256,   0.000569, &
                             0.6522, -0.1038,   0.000143, &
                             0.3248, -0.0434,  -0.00013,  &
                             6.3785,  0.0125,  -0.00044,  & ! Mountain-Prarie
                             0.4100, -0.0807,  -0.00052,  &
                             0.2123, -0.0324,  -0.00086,  &
                             9.2347,  0.00953,  0.000518, &
                             0.5346, -0.0597,   0.000425, &
                             0.3292, -0.0444,   3.731E-6, &
                             5.3432,  0.0543,  -0.00007,  & ! Northeast
                             0.4076, -0.1033,  -0.00077,  &
                             0.1647,  0.0321,  -0.00044,  &
                             9.1670,  0.00654,  0.000422, &
                             0.5932, -0.0942,   0.000070, &
                             0.3201, -0.0454,  -0.00016,  &
                             6.3084,  0.00125, -0.00088,  & ! Northwest
                             0.3826, -0.0511,  -0.00054,  &
                             0.2093, -0.0294,  -0.00127,  &
                            10.3325, -0.0281,  -0.00018,  &
                             0.6363, -0.1005,  -0.00016,  &
                             0.4078, -0.1140,  -0.00105,  &
                             4.9953,  0.1081,   0.000272, & ! Southeast
                             0.2721,  0.0707,   0.000514, &
                             0.1328,  0.1433,   0.000401, &
                             8.9977,  0.0576,   0.00101,  &
                             0.4988, -0.00346,  0.000839, &
                             0.3199,  0.0165,   0.000614, &
                             7.3700, -0.0472,  -0.00100,  & ! Southwest
                             0.4305, -0.1142,  -0.00106,  &
                             0.2533, -0.0845,  -0.00138,  &
                             9.6158, -0.0218,  -0.00023,  &
                             0.5436, -0.0701,  -0.00008,  &
                             0.3420, -0.0673,  -0.00069   /

    data regional_mandg_means / 2.4457, -0.00058,  7.297E-6,  14.3740,  & ! Mideast
                                2.4899, -0.00398,  4.112E-7,  10.4446,  &
                                2.1401, -0.00120,  7.27E-6,   14.6920,  & ! Midwest
                                2.2544, -0.00560, -3.97E-6,   12.6155,  &
                                2.0842, -0.00133,  5.295E-6,  13.4341,  & ! Mountain-Prarie
                                2.4834, -0.00334,  4.874E-6,   8.9723,  &
                                2.0111, -0.00059,  9.516E-6,  13.7254,  & ! Northeast
                                2.3053, -0.00451,  1.963E-7,  10.6196,  &
                                1.7327,  0.000917, 0.000021,  15.2489,  & ! Northwest
                                1.9541, -0.00482,  3.818E-6,  13.4681,  &
                                2.1773, -0.00298,  4.274E-7,  15.7907,  & ! Southeast
                                2.3705, -0.00566, -5.19E-6,   12.0497,  &
                                1.9580, -0.00153,  7.828E-6,  15.4836,  & ! Southwest
                                1.7987, -0.00824, -0.00001,   14.6025   /

    data regional_mandg_sd    / 1.8239, -0.00018,  -8.96E-7,   0.5750,  & ! Mideast
                                2.3113,  0.00171,   1.395E-6, -3.2624,  &
                                1.7833, -0.00045,  -3E-6,      1.1579,  & ! Midwest
                                2.3692,  0.00251,   4.185E-6, -3.5947,  &
                                1.8310, -0.00070,  -2.77E-6,   1.0971,  & ! Mountain-Prarie
                                2.3503,  0.00189,   2.984E-6, -2.9599,  &
                                1.6817, -0.00048,  -2.13E-6,   1.6159,  & ! Northeast
                                2.3439,  0.00229,   3.028E-6, -3.8106,  &
                                1.7818, -0.00007,   1.613E-7,  1.1441,  & ! Northwest
                                2.3364,  0.00168,   2.282E-7, -3.1963,  &
                                1.8424, -0.00117,  -5.34E-6,   0.1396,  & ! Southeast
                                2.2432,  0.000769, -1.65E-6,  -2.3455,  &
                                1.7124, -0.00095,  -4.31E-6,   1.1682,  & ! Southwest
                                2.2315,  0.00214,   4.173E-6, -3.9374   /

    ! Season of calving-specific curves
    ! 1: spring (MAM), 2: summer (JJA), 3: fall (SON), 4: winte r (DJF),
    data calving_woods_means  / 15.786,  0.20756, 0.002180946,  & ! Spring
                                 0.8792, 0.08042, 0.001074334,  &
                                 0.4975, 0.17969, 0.001523427,  &
                                26.5014, 0.15818, 0.003100708,  &
                                 1.5492, 0.0106,  0.001743349,  &
                                 0.9469, 0.09402, 0.002090426,  &
                                12.817,  0.26053, 0.002503766,  & ! Summer
                                 0.623,  0.17205, 0.001678702,  &
                                 0.3915, 0.24431, 0.001985798,  &
                                21.6622, 0.20213, 0.003274469,  &
                                 1.0064, 0.12036, 0.002396856,  &
                                 0.7096, 0.17242, 0.002655959,  &
                                12.4642, 0.27854, 0.002774397,  & ! Fall
                                 0.765,  0.13252, 0.001667579,  &
                                 0.4583, 0.21308, 0.001981516,  &
                                20.6226, 0.23288, 0.003704124,  &
                                 1.1985, 0.09865, 0.002613839,  &
                                 0.8157, 0.1535,  0.002853117,  &
                                14.6238, 0.23418, 0.002464894,  & ! Winter
                                 0.9842, 0.05995, 0.00112779,   &
                                 0.5479, 0.16034, 0.001543881,  &
                                23.5231, 0.20455, 0.003615274,  &
                                 1.5954, 0.02206, 0.002092818,  &
                                 0.9895, 0.09734, 0.002365836   /

     data calving_woods_sd  / 5.8618,  0.03394, -0.00025351,   & ! Spring
                              0.413,  -0.08783, -0.000635602,  &
                              0.1764,  0.02674, -0.000490834,  &
                              9.7166, -0.00127,  0.000306529,  &
                              0.6186, -0.09193,  0.000112934,  &
                              0.3253, -0.03589, -0.000133069,  &
                              5.037,   0.06901, -0.000087756,  & ! Summer
                              0.3303, -0.03968, -0.000467502,  &
                              0.1648,  0.03672, -0.000505643,  &
                              8.0277,  0.04631,  0.000611987,  &
                              0.4704, -0.0316,   0.000357456,  &
                              0.2811, -0.00176,  0.000075761,  &
                              5.065,   0.07221,  0.000008603,  & ! Fall
                              0.3952, -0.07742, -0.000567239,  &
                              0.1667,  0.03627, -0.000421237,  &
                              8.3162,  0.04254,  0.000648105,  &
                              0.5382, -0.05281,  0.000400812,  &
                              0.3083, -0.02284,  0.000031604,  &
                              7.1243, -0.02568, -0.000802699,  & ! Winter
                              0.4791, -0.12932, -0.000998077,  &
                              0.2209, -0.03916, -0.001020224,  &
                             10.2597, -0.02345,  0.000041742,  &
                              0.6656, -0.11465, -0.000143019,  &
                              0.3822, -0.08907, -0.000616674   /

    data calving_mandg_means  / 2.4468,  0.00232,  0.000002592, -4.7834,  & ! Spring
                                2.2135, -0.0053,  -0.000003797, 10.424,   &
                                2.1349, -0.00034,  0.000012141, 16.3896,  & ! Summer
                                2.2706, -0.00485,  0.00000048,  14.5693,  &
                                1.8605, -0.00182,  0.000007401, 14.7473,  & ! Fall
                                2.0393, -0.00632, -0.000003872, 13.3113,  &
                                2.0166, -0.00176,  0.000004812, 15.2227,  & ! Winter
                                2.1226, -0.00653, -0.000007269, 11.6598   /

    data calving_mandg_sd     / 1.8669,  0.0003,   0.000000128, -0.486,   & ! Spring
                                2.4679,  0.00326,  0.000007042, -5.2771,  &
                                1.6629, -0.00091, -0.000003794,  2.2219,  & ! Summer
                                2.1873,  0.00186,  0.000002439, -1.684,   &
                                1.6253, -0.00147, -0.000005573,  2.408,   & ! Fall
                                2.21,    0.00148,  0.000001208, -3.0899,  &
                                1.8903, -0.00022, -0.00000259,   0.2213,  & ! Winter
                                2.4468,  0.00232,  0.000002592, -4.7834   /

    ! These are acutally region and season-specific
    data seasonal_woods_means / 17.0437,  0.18043,  0.001958839,  0.821,   0.08534,  0.000904043, 0.5049,  0.16738,  0.001356394,  & ! Mideast spring
                                26.0936,  0.1352,   0.002705453,  1.3223,  0.02258,  0.001492581, 0.8895,  0.08699,  0.001812505,  &
                                13.3792,  0.24137,  0.002305127,  0.5465,  0.19682,  0.00166602,  0.3809,  0.24462,  0.001911111,  & ! summer
                                19.5537,  0.20419,  0.003089742,  0.7718,  0.16541,  0.002440179, 0.5796,  0.20743,  0.002740285,  &
                                12.7817,  0.2633,   0.002697487,  0.6924,  0.14698,  0.001671972, 0.4712,  0.19896,  0.001910598,  & ! fall
                                19.9451,  0.22473,  0.003557596,  1.0545,  0.11032,  0.002452382, 0.7835,  0.15106,  0.002744022,  &
                                14.6635,  0.22746,  0.002540184,  0.9201,  0.07053,  0.001193681, 0.5391,  0.15809,  0.001604175,  & ! winter
                                23.9416,  0.18406,  0.003425372,  1.4994,  0.01883,  0.001923146, 1.0034,  0.08053,  0.002194885,  &
                                17.1634,  0.18562,  0.002092785,  1.0171,  0.04776,  0.000949897, 0.5319,  0.15828,  0.001384069,  & ! Midwest spring
                                27.6632,  0.1379,   0.002925842,  1.7207, -0.01831,  0.001579947, 0.9792,  0.07628,  0.001927055,  &
                                13.9049,  0.23551,  0.00234641,   0.7131,  0.13742,  0.001466017, 0.4094,  0.22656,  0.001830986,  & ! summer
                                21.7122,  0.18918,  0.003111794,  1.053,   0.10086,  0.002225952, 0.7043,  0.16346,  0.002524498,  &
                                13.3241,  0.26235,  0.002709852,  0.8606,  0.10522,  0.001520375, 0.4807,  0.2,      0.001905471,  & ! fall
                                20.5657,  0.22847,  0.003677559,  1.2683,  0.08206,  0.002514697, 0.821,   0.1466,   0.002801255,  &
                                15.3486,  0.22487,  0.002526249,  1.0936,  0.03917,  0.001118781, 0.569,   0.15062,  0.001554824,  & ! winter
                                24.0732,  0.19237,  0.003531736,  1.7434, -0.00178,  0.001974758, 1.0254,  0.08073,  0.002232456,  &
                                13.48,    0.24614,  0.002416709,  0.7788,  0.10706,  0.001289647, 0.4377,  0.20952,  0.001734488,  & ! Mountain-Prarie spring
                                24.7864,  0.17381,  0.003310573,  1.4259,  0.02199,  0.001811255, 0.8909,  0.10848,  0.002269458,  &
                                11.0897,  0.29039,  0.002583227,  0.5916,  0.17715,  0.001622558, 0.378,   0.24583,  0.001919248,  & ! summer
                                21.0423,  0.20017,  0.003200416,  0.9305,  0.12975,  0.002353393, 0.7007,  0.16926,  0.002577369,  &
                                12.3453,  0.27211,  0.002550454,  0.7697,  0.12319,  0.001507992, 0.4881,  0.18713,  0.001675102,  & ! fall
                                20.7794,  0.22079,  0.003558818,  1.2052,  0.08645,  0.002505003, 0.8542,  0.13419,  0.002694594,  &
                                12.3247,  0.27493,  0.002647804,  1.0122,  0.04754,  0.001037794, 0.524,   0.16627,  0.001539437,  & ! winter
                                22.4089,  0.21126,  0.003669395,  1.5537,  0.01994,  0.002081158, 0.9397,  0.10521,  0.002446673,  &
                                17.3334,  0.17923,  0.002109677,  0.9898,  0.04343,  0.000891967, 0.5329,  0.15238,  0.001416768,  & ! Northeast spring
                                27.1048,  0.14382,  0.00307422,   1.6179, -0.01069,  0.001632561, 0.9743,  0.07464,  0.001989632,  &
                                14.2892,  0.22621,  0.002346541,  0.6714,  0.14559,  0.00152212,  0.4146,  0.21817,  0.001830084,  & ! summer
                                22.6572,  0.18332,  0.003218272,  1.0299,  0.10732,  0.002327028, 0.7101,  0.16187,  0.002601094,  &
                                13.8242,  0.24569,  0.00260921,   0.8155,  0.10672,  0.001465873, 0.4798,  0.19138,  0.001852328,  & ! fall
                                21.882,   0.20969,  0.003559749,  1.2871,  0.07441,  0.002481578, 0.819,   0.14302,  0.00279859,   &
                                16.1024,  0.20661,  0.002430498,  1.0514,  0.03661,  0.001049722, 0.5859,  0.13583,  0.001485865,  & ! winter
                                23.3869,  0.20113,  0.003710054,  1.6458,  0.0084,   0.002098594, 0.9678,  0.09562,  0.00245217,   &
                                12.6563,  0.2818,   0.002550693,  0.7421,  0.132,    0.001368242, 0.3856,  0.26772,  0.002041277,  & ! Northwest spring
                                24.4051,  0.2039,   0.003408765,  1.3796,  0.05116,  0.001958237, 0.9211,  0.12784,  0.002382937,  &
                                10.3989,  0.33446,  0.002893651,  0.5465,  0.21928,  0.002003455, 0.308,   0.32923,  0.002481978,  & ! summer
                                25.6466,  0.17892,  0.00306657,   1.1873,  0.0936,   0.002256451, 0.9057,  0.13212,  0.002422058,  &
                                10.7767,  0.3308,   0.002937251,  0.7666,  0.14426,  0.001721803, 0.4011,  0.26494,  0.002150768,  & ! fall
                                21.5479,  0.23335,  0.003523797,  1.3701,  0.07506,  0.002428398, 0.8852,  0.14547,  0.002648453,  &
                                12.5084,  0.28149,  0.002484083,  0.9335,  0.0752,   0.001080822, 0.4548,  0.21892,  0.0016656,    & ! winter
                                22.2512,  0.23444,  0.00363641,   1.4916,  0.05147,  0.002214229, 0.9353,  0.13205,  0.002506354,  &
                                14.5979,  0.2105,   0.002184077,  0.6966,  0.10431,  0.000991857, 0.4816,  0.17614,  0.001460786,  & ! Southeast spring
                                23.706,   0.1583,   0.00309825,   0.9422,  0.08206,  0.001783547, 0.8103,  0.10255,  0.001955669,  &
                                 9.5258,  0.32981,  0.002969289,  0.4057,  0.26797,  0.002160685, 0.3052,  0.3076,   0.002390987,  & ! summer
                                15.4084,  0.26238,  0.003500045,  0.6176,  0.20737,  0.002597081, 0.4896,  0.24657,  0.002973388,  &
                                10.3741,  0.3215,   0.003153717,  0.5908,  0.18893,  0.002129566, 0.4229,  0.23542,  0.002329659,  & ! fall
                                16.9558,  0.26902,  0.00409573,   0.913,   0.1407,   0.002771545, 0.778,   0.154,    0.002935373,  &
                                13.706,   0.24077,  0.002598673,  0.8768,  0.0678,   0.001121388, 0.5523,  0.15148,  0.00153486,   & ! winter
                                21.8619,  0.2129,   0.003949707,  1.344,   0.02959,  0.002058304, 1.0155,  0.07257,  0.002296539,  &
                                14.8289,  0.22363,  0.002142373,  0.7232,  0.1324,   0.001302076, 0.4812,  0.19787,  0.001597059,  & ! Southwest spring
                                27.6048,  0.16465,  0.003095213,  1.4083,  0.05211,  0.001999085, 0.9802,  0.1108,   0.002238183,  &
                                12.266,   0.27712,  0.002590054,  0.559,   0.2047,   0.001909643, 0.4081,  0.24662,  0.002087066,  & ! summer
                                22.0625,  0.21975,  0.003477694,  1.0025,  0.14193,  0.002651575, 0.7549,  0.18211,  0.002838828,  &
                                11.1821,  0.30965,  0.002922253,  0.64,    0.18333,  0.002004642, 0.4167,  0.24631,  0.002203753,  & ! fall
                                20.9626,  0.24612,  0.003817073,  1.1073,  0.13567,  0.002920786, 0.8003,  0.17937,  0.003092096,  &
                                13.9436,  0.24499,  0.002323514,  0.8428,  0.09995,  0.001214844, 0.5212,  0.17937,  0.001524638,  & ! winter
                                24.4486,  0.20658,  0.003521926,  1.448,   0.05851,  0.002244196, 1.0056,  0.1126,   0.002427227   /

    data seasonal_woods_sd    /  6.5683, -0.00245, -0.000428644,  0.4835, -0.13323, -0.00104866,  0.1823,  0.00977, -0.000568653,  & ! Mideast spring
                                12.4,    -0.09267, -0.000359023,  0.6876, -0.13845, -0.000297868, 0.3646, -0.09634, -0.000659871,  &
                                 5.1678,  0.06002, -0.000099431,  0.3139, -0.01998, -0.000365023, 0.1641,  0.02916, -0.000630759,  & ! summer
                                 8.2773,  0.0195,   0.000424932,  0.3455,  0.03734,  0.000652363, 0.2768, -0.01836, -0.00011554,   &
                                 5.0377,  0.06899,  0.00008701,   0.3451, -0.03623, -0.0002192,   0.1721,  0.02509, -0.000381079,  & ! fall
                                 9.755,  -0.02935,  0.000017215,  0.4872, -0.04114,  0.000365591, 0.3279, -0.06003, -0.000277041,  &
                                 8.5379, -0.07793, -0.000902578,  0.4685, -0.11399, -0.000739536, 0.2766, -0.10407, -0.001258177,  & ! winter
                                10.5661, -0.04217,  0.000150384,  0.6388, -0.10467,  0.000122222, 0.4154, -0.12927, -0.00080389,   &
                                 5.9187,  0.03116, -0.000228245,  0.429,  -0.09732, -0.000646677, 0.1846,  0.00897, -0.000591038,  & ! Midwest spring
                                 9.7367, -0.0076,   0.000326253,  0.6968, -0.12279,  0.000001051, 0.3326, -0.05228, -0.000228775,  &
                                 4.9892,  0.0718,  -0.000053291,  0.3478, -0.05849, -0.000595285, 0.1581,  0.04121, -0.000495201,  & ! summer
                                 7.7351,  0.04507,  0.000556914,  0.512,  -0.05923,  0.000226085, 0.2608,  0.00626,  0.000093652,  &
                                 4.8101,  0.08616,  0.000134041,  0.4166, -0.09417, -0.000676219, 0.1659,  0.03308, -0.000431246,  & ! fall
                                 7.6829,  0.0647,   0.000905919,  0.5521, -0.05718,  0.000489766, 0.2924, -0.01287,  0.000158556,  &
                                 5.749,   0.04119, -0.00011823,   0.4389, -0.09559, -0.000546699, 0.1912,  0.00187, -0.00057227,   & ! winter
                                 9.4117,  0.00424,  0.000429862,  0.702,  -0.1234,  -0.000013947, 0.3781, -0.0899,  -0.000514065,  &
                                 6.527,   0.00583, -0.000442552,  0.365,  -0.04717, -0.000268499, 0.2073, -0.02785, -0.00085675,   & ! Mountain-Prarie spring
                                 9.8855, -0.00368,  0.000584482,  0.5887, -0.08866,  0.000175392, 0.3302, -0.04178,  0.000025062,  &
                                 5.8363,  0.03431, -0.000176439,  0.3413, -0.03823, -0.000248429, 0.191,  -0.00409, -0.000552619,  & ! summer
                                 9.1862, -0.00919,  0.000138938,  0.4058, -0.01141,  0.000378593, 0.3297, -0.06117, -0.000352166,  &
                                 5.3554,  0.05192, -0.000226653,  0.361,  -0.04852, -0.00014658,  0.1778,  0.01074, -0.00051082,   & ! fall
                                 8.3202,  0.03716,  0.000676335,  0.5358, -0.0632,   0.000416827, 0.3142, -0.03131,  0.000177215,  &
                                 8.0015, -0.04149, -0.000907973,  0.6202, -0.20042, -0.001637429, 0.2949, -0.11835, -0.001541368,  & ! winter
                                 9.5373,  0.0044,   0.000488444,  0.5279, -0.04715,  0.000596027, 0.3264, -0.03975,  0.00011573,   &
                                 5.4274,  0.05373, -0.000010561,  0.4436, -0.12556, -0.000905878, 0.1624,  0.04137, -0.000311459,  & ! Northeast spring
                                 9.5881, -0.00667,  0.000336985,  0.5816, -0.08918,  0.000130456, 0.3211, -0.04844, -0.000254919,  &
                                 4.7883,  0.07709,  0.000035393,  0.3387, -0.06718, -0.00063237,  0.144,   0.05703, -0.000365465,  & ! summer
                                 8.1997,  0.03128,  0.000590386,  0.4769, -0.04948,  0.000221783, 0.2647,  0.00142,  0.000132645,  &
                                 5.1791,  0.0584,  -0.000091588,  0.3833, -0.09275, -0.000765624, 0.1541,  0.04386, -0.000407036,  & ! fall
                                 9.7802, -0.01968,  0.000071394,  0.5899, -0.09993, -0.000100359, 0.3284, -0.05855, -0.000293425,  &
                                 6.2801,  0.01202, -0.000346654,  0.466,  -0.13229, -0.000877658, 0.2001, -0.01744, -0.000736926,  & ! winter
                                 9.4144,  0.00297,  0.000466678,  0.6186, -0.10188,  0.000076075, 0.3499, -0.07049, -0.000314133,  &
                                 5.8899,  0.0267,  -0.000612112,  0.3281, -0.00605, -0.000075087, 0.1654,  0.04321, -0.000640281,  & ! Northwest spring
                                 9.5348, -0.00166,  0.000145392,  0.549,  -0.05564,  0.00037681,  0.3113, -0.03105, -0.000205712,  &
                                 6.8712, -0.03825, -0.00129269,   0.3353, -0.02129, -0.00049686,  0.2088, -0.04617, -0.001555002,  & ! summer
                                11.3249, -0.05965, -0.00050513,   0.535,  -0.07284, -0.000305407, 0.3685, -0.09432, -0.001064287,  &
                                 5.5608,  0.04464, -0.000457772,  0.4218, -0.07377, -0.000622598, 0.1984, -0.01043, -0.001087772,  & ! fall
                                 9.4062, -0.00023,  0.000072898,  0.6309, -0.08916,  0.000111886, 0.4131, -0.11247, -0.000920949,  &
                                 6.4084,  0.00075, -0.000882701,  0.4391, -0.09473, -0.000861067, 0.2497, -0.06931, -0.001451231,  & ! winter
                                11.8055, -0.07161, -0.000605527,  0.7458, -0.14133, -0.000438468, 0.5175, -0.18209, -0.001616789,  &
                                 5.2337,  0.10148,  0.000250833,  0.2847,  0.04683,  0.000167887, 0.1254,  0.15806,  0.000357397,  & ! Southeast spring
                                11.2702,  0.01867,  0.000859437,  0.4894,  0.00996,  0.0010064,   0.304,   0.05043,  0.000960867,  &
                                 3.4875,  0.20222,  0.00085582,   0.1559,  0.21273,  0.001390713, 0.0891,  0.24145,  0.000908561,  & ! summer
                                 6.6725,  0.1244,   0.001216735,  0.3593,  0.06038,  0.000789503, 0.1623,  0.18264,  0.001449845,  &
                                 5.1401,  0.08855,  0.000189592,  0.3177,  0.02364,  0.000388872, 0.1453,  0.11183,  0.000378497,  & ! fall
                                 8.723,   0.05442,  0.00097472,   0.4842, -0.00419,  0.000833904, 0.3538, -0.02197,  0.000327112,  &
                                 7.1007,  0.01541, -0.000269955,  0.3734, -0.00675,  0.000087373, 0.2008,  0.03922, -0.00010623,   & ! winter
                                 9.6232,  0.03582,  0.000877707,  0.6617, -0.07828,  0.000409743, 0.4406, -0.0812,  -0.000093528,  &
                                 6.7194, -0.01624, -0.00063189,   0.3899, -0.0835,  -0.000713068, 0.2266, -0.04879, -0.001029372,  & ! Southwest spring
                                10.1637, -0.03579, -0.00033737,   0.5485, -0.07269, -0.000121807, 0.3556, -0.07466, -0.000767174,  &
                                 5.5576,  0.03506, -0.000261567,  0.3329, -0.04566, -0.000531772, 0.2122, -0.0282,  -0.000808764,  & ! summer
                                 6.8617,  0.07502,  0.000612173,  0.4045,  0.00957,  0.000581181, 0.2528,  0.01853,  0.000040613,  &
                                 5.7162,  0.03559, -0.000170228,  0.3776, -0.06596, -0.000504688, 0.2014, -0.01246, -0.000636067,  & ! fall
                                 7.7535,  0.05152,  0.000581445,  0.4688, -0.01772,  0.00052925,  0.2883, -0.00879,  0.000013632,  &
                                10.2449, -0.15036, -0.001938473,  0.5489, -0.1894,  -0.00176384,  0.3184, -0.15659, -0.002068627,  & ! winter
                                12.9285, -0.11996, -0.001144001,  0.6702, -0.1384,  -0.000764201, 0.4316, -0.1454,  -0.001458735   /

    data seasonal_mandg_means /  2.6658,  0.00064,  0.000009754, 11.0239,  & ! Mideast spring
                                 2.5144, -0.0038,  -0.00000174,   6.0578,  &
                                 2.7328,  0.00164,  0.000016318, 12.8311,  & ! Mideast summer
                                 2.7724, -0.00262,  0.000004644, 10.8673,  &
                                 2.1655, -0.00229,  0.000002976, 15.8377,  & ! Mideast fall
                                 2.1199, -0.00605, -0.000005217, 16.158,   &
                                 2.2645, -0.00192,  0.00000201,  17.7088,  & ! Mideast winter
                                 2.5043, -0.0041,   0.000000446,  7.9796,  &
                                 2.4953,  0.00115,  0.000014321, 10.9857,  & ! Midwest spring
                                 2.4155, -0.00455, -0.000002164,  8.4987,  &
                                 2.1957, -0.00059,  0.000010277, 17.3705,  & ! Midwest summer
                                 2.4622, -0.0038,   0.000003404, 14.8041,  &
                                 1.824,  -0.00281,  0.000003382, 15.9011,  & ! Midwest fall
                                 2.0599, -0.00646, -0.000004627, 13.3347,  &
                                 2.141,  -0.00194,  0.000003381, 13.8371,  & ! Midwest winter
                                 2.1351, -0.00705, -0.000010201, 12.9145,  &
                                 2.303,  -0.00108,  0.000002302, 10.5964,  & ! Mountain-Prarie spring
                                 2.7033, -0.00247,  0.000004624,  5.5784,  &
                                 2.0921, -0.00117,  0.000007635, 18.108,   & ! Mountain-Prarie summer
                                 2.5527, -0.00182,  0.000012721, 12.6629,  &
                                 1.6707, -0.0035,  -0.000002243, 14.413,   & ! Mountain-Prarie fall
                                 2.1674, -0.00518, -0.00000127,  10.6217,  &
                                 2.385,   0.00137,  0.000017873, 10.674,   & ! Mountain-Prarie winter
                                 2.5755, -0.00364,  0.000003547,  6.0502,  &
                                 2.3028,  0.00127,  0.000015224, 10.3572,  & ! Northeast spring
                                 2.3953, -0.00434, -0.000001902,  8.5863,  &
                                 2.1384,  0.00073,  0.000015133, 15.6513,  & ! Northeast summer
                                 2.5191, -0.00212,  0.000011082, 11.0002,  &
                                 1.6689, -0.00268,  0.000002541, 15.1274,  & ! Northeast fall
                                 2.1254, -0.00534, -0.000001012, 13.0412,  &
                                 1.9272, -0.00191,  0.00000384,  14.1052,  & ! Northeast winter
                                 2.2042, -0.0059,  -0.000005551,  9.8707,  &
                                 1.6659, -0.00049,  0.000013293, 15.2826,  & ! Northwest spring
                                 1.6891, -0.00728, -0.000009414, 15.6219,  &
                                 1.7417,  0.00106,  0.000023385, 16.9648,  & ! Northwest summer
                                 2.263,  -0.00275,  0.00001232,  11.7284,  &
                                 1.8314,  0.00326,  0.000031697, 10.5912,  & ! Northwest fall
                                 1.8939, -0.00527,  0.000003062, 14.2684,  &
                                 1.7464,  0.00046,  0.000016128, 17.1519,  & ! Northwest winter
                                 1.7709, -0.00613, -0.000001217, 14.4485,  &
                                 2.5734, -0.00103,  0.000006634, 11.4676,  & ! Southeast spring
                                 2.4998, -0.00598, -0.000010064,  7.448,   &
                                 2.5833,  0.00039,  0.000013884, 13.2926,  & ! Southeast summer
                                 2.781,  -0.0026,   0.000005393, 14.1729,  &
                                 2.184,  -0.00126,  0.000010979, 13.0096,  & ! Southeast fall
                                 2.367,  -0.00423,  0.000004275, 12.1613,  &
                                 1.601,  -0.0079,  -0.000019697, 23.7518,  & ! Southeast winter
                                 2.1465, -0.00791, -0.000013538,  9.9888,  &
                                 2.0074, -0.0017,   0.000005392, 15.579,   & ! Southwest spring
                                 1.6635, -0.00823, -0.000011205, 16.063,   &
                                 1.9784, -0.0016,   0.000008436, 16.0774,  & ! Southwest summer
                                 1.6515, -0.00975, -0.000016698, 17.7184,  &
                                 1.9485, -0.0007,   0.000013126, 13.5226,  & ! Southwest fall
                                 1.9055, -0.00748, -0.00000679,  12.6665,  &
                                 1.9546, -0.00146,  0.000007293, 16.1875,  & ! Southwest winter
                                 1.9712, -0.00715, -0.000007514, 11.8503   /

    data seasonal_mandg_sd    /  2.0313,  0.00139,  0.000004835, -2.1678,  & ! Mideast spring
                                 2.4006,  0.00194,  0.000000776, -4.7794,  &
                                 1.6143, -0.00151, -0.000006503,  2.8182,  & ! Mideast summer
                                 2.16,    0.00081, -0.000002399, -0.4493,  &
                                 1.704,  -0.0009,  -0.000003001,  1.3246,  & ! Mideast fall
                                 2.2038,  0.00154,  0.000002282, -2.9371,  &
                                 1.9548,  0.0002,   0.000000058,  0.3223,  & ! Mideast winter
                                 2.5335,  0.00271,  0.000004396, -6.0361,  &
                                 1.9038,  0.00062,  0.000001112, -0.5132,  & ! Midwest spring
                                 2.4965,  0.00334,  0.000006957, -5.0526,  &
                                 1.598,  -0.00131, -0.000005838,  2.8493,  & ! Midwest summer
                                 2.2584,  0.0024,   0.000004198, -2.023,   &
                                 1.6511, -0.00128, -0.000005202,  2.2125,  & ! Midwest fall
                                 2.2587,  0.00187,  0.000002587, -3.4419,  &
                                 1.9981,  0.00031, -0.000001577, -0.4307,  & ! Midwest winter
                                 2.4673,  0.00241,  0.000002577, -4.4914,  &
                                 2.1205,  0.00118,  0.000002887, -2.8497,  & ! Mountain-Prarie spring
                                 2.6515,  0.00361,  0.000006148, -6.7892,  &
                                 1.721,  -0.001,   -0.000003397,  2.6333,  & ! Mountain-Prarie summer
                                 2.2439,  0.00239,  0.000006825, -1.2853,  &
                                 1.4706, -0.0034,  -0.000012277,  4.4711,  & ! Mountain-Prarie fall
                                 2.1913,  0.00138,  0.000003727, -2.0401,  &
                                 1.9725,  0.00008,  0.000000313,  0.0847,  & ! Mountain-Prarie winter
                                 2.4792,  0.00179,  0.000001882, -4.2953,  &
                                 1.7181, -0.00032, -0.000002433,  0.8495,  & ! Northeast spring
                                 2.4219,  0.00279,  0.00000514,  -4.7571,  &
                                 1.653,   0.00002,  0.000001196,  2.1161,  & ! Northeast summer
                                 2.162,   0.00157,  0.000000182, -1.5636,  &
                                 1.5257, -0.00143, -0.000004878,  2.6719,  & ! Northeast fall
                                 2.34,    0.00262,  0.000005206, -4.0166,  &
                                 1.8218, -0.00015, -0.000002027,  0.5417,  & ! Northeast winter
                                 2.454,   0.00231,  0.000002006, -5.1689,  &
                                 1.7977, -0.00005, -0.00000074,  -0.6732,  & ! Northwest spring
                                 2.3711,  0.0015,  -0.000001974, -3.1833,  &
                                 1.6548, -0.00072, -0.000002121,  3.042,   & ! Northwest summer
                                 2.4643,  0.0033,   0.000006883, -4.7796,  &
                                 1.6864, -0.0004,   0.000000087,  2.231,   & ! Northwest fall
                                 2.1458,  0.00024, -0.000004839, -1.5646,  &
                                 1.922,  -0.00008, -0.000001537,  0.3732,  & ! Northwest winter
                                 2.2199, -0.00002, -0.000007649, -1.5974,  &
                                 1.9993,  0.00001, -0.00000105,  -1.211,   & ! Southeast spring
                                 2.3744,  0.00161,  0.000001344, -3.8538,  &
                                 1.8732, -0.0008,  -0.000003594,  0.0965,  & ! Southeast summer
                                 2.3392,  0.00258,  0.000007033, -2.4519,  &
                                 1.518,  -0.00392, -0.000016616,  3.3022,  & ! Southeast fall
                                 1.9548, -0.00144, -0.000010889,  0.4838,  &
                                 1.8866, -0.00124, -0.000007125, -0.6012,  & ! Southeast winter
                                 2.3742,  0.0003,  -0.000006127, -5.2165,  &
                                 1.7784, -0.00032, -0.000002015, -0.073,   & ! Southwest spring
                                 2.4043,  0.00392,  0.000011894, -6.0768,  &
                                 1.6845, -0.00101, -0.00000439,   1.3717,  & ! Southwest summer
                                 1.9745,  0.00082, -0.000000201, -1.0723,  &
                                 1.6148, -0.00167, -0.000006673,  2.5655,  & ! Southwest fall
                                 2.123,   0.00114,  0.000000098, -3.3872,  &
                                 1.7716, -0.00078, -0.000004093,  0.6427,  & ! Southwest winter
                                 2.4701,  0.00328,  0.000007631, -5.9746   /

    dyield(:,trait,lacn) = 0.
    sd(trait,:,lacn) = 0.
    meanyld(trait,:,lacn) = 0.

    !
    ! Check to see if we've been given valid values for methods when we're
    ! interpolating. If we haven't, revert to linear interpolation.
    !
    if ( trait < 4 ) then
      if ( method .ne. 'L' .and. method .ne. 'W' .and. method .ne. 'R' .and. method .ne. 'T' .and. method .ne. 'C' ) then
        print *, "[WARNING]: Invalid interpolation method, ", method, &
          " provided for trait ", trait, ". Using linear interpolation (L)."
        method = 'L'
      end if
    else
      if ( method .ne. 'L' .and. method .ne. 'G' .and. method .ne. 'S' .and. method .ne. 'U' .and. method .ne. 'D' ) then
        print *, "[WARNING]: Invalid interpolation method, ", method, &
          " provided for SCS. Using linear interpolation (L)."
        method = 'L'
      end if
    end if

    ! If an invalid breed is provided default to Holstein.
    if ( breed < 1 .or. breed > 6  ) then
      print *, "[WARNING]: Invalid breed, ", breed, &
        " provided for trait ", trait, ". Using Holstein (4)."
      breed = 4
    end if

    select case(method)
        !
        ! L = linear interpolation
        !
        case ('L')
          !            Interpolate between monthly means and sd
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i=1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = ((i-month*30+15)*dyld(month+1,trait,lacn) &
              + ((month+1)*30-15-i)*dyld(month,trait,lacn))/30.
            !            Constant yield assumed after day 365
            if(i > 365) dyield(i,trait,lacn) = dyld(12,trait,lacn)
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if(i <= 305) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = ((i-month*30+15)*dsd(month+1,trait,lacn) + &
            ((month+1)*30-15-i)*dsd(month,trait,lacn))/30.
            !                   Constant SD assumed after day 365
            if(i > 365) sd(trait,i,lacn) = dsd(12,trait,lacn)
          end do
        !
        ! W = Woods curves
        !
        case('W')
          ! Interpolate for MFP mean and SD using the Woods curves
          ! recommended by Mahinda Dematawewa et al. (2007, in press).
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = woods_means(1,trait,lacn,breed) &
              * float(i)**woods_means(2,trait,lacn,breed) &
              * exp(-float(i)*woods_means(3,trait,lacn,breed))
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = woods_sd(1,trait,lacn,breed) &
              * float(i)**woods_sd(2,trait,lacn,breed) &
              * exp(-float(i)*woods_sd(3,trait,lacn,breed))
          end do
        !
        ! R = Region-specific Woods curves
        !
        case('R')
          ! Interpolate for MFP mean and SD using the region-specific Woods curves
          ! calculated by Dan Null using the dataset of Dematawewa et al. (2008).
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = regional_woods_means(1,trait,lacn,region) &
              * float(i)**regional_woods_means(2,trait,lacn,region) &
              * exp(-float(i)*regional_woods_means(3,trait,lacn,region))
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = regional_woods_sd(1,trait,lacn,region) &
              * float(i)**regional_woods_sd(2,trait,lacn,region) &
              * exp(-float(i)*regional_woods_sd(3,trait,lacn,region))
          end do
          if ( DEBUGmsgs > 0 ) then
            print *, '[R][', lacn, ']: mean curve:', regional_woods_means(:,trait,lacn,region)
            print *, '[R][', lacn, ']: SD curve  :', regional_woods_sd(:,trait,lacn,region)
          end if
        !
        ! C = Season of calving-specific Woods curves
        !
        case('C')
          ! Interpolate for MFP mean and SD using the region-specific Woods curves
          ! calculated by Dan Null using the dataset of Dematawewa et al. (2008).
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = calving_woods_means(1,trait,lacn,season) &
              * float(i)**calving_woods_means(2,trait,lacn,season) &
              * exp(-float(i)*calving_woods_means(3,trait,lacn,season))
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = calving_woods_sd(1,trait,lacn,season) &
              * float(i)**calving_woods_sd(2,trait,lacn,season) &
              * exp(-float(i)*calving_woods_sd(3,trait,lacn,season))
          end do
        !
        ! T = Region- and season of calving-specific Woods curves
        !
        case('T')
          ! Interpolate for MFP mean and SD using the region- and season-specific Woods
          ! curves calculated by Dan Null using the dataset of Dematawewa et al. (2008).
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = seasonal_woods_means(1,trait,lacn,season,region) &
              * float(i)**seasonal_woods_means(2,trait,lacn,season,region) &
              * exp(-float(i)*seasonal_woods_means(3,trait,lacn,season,region))
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = seasonal_woods_sd(1,trait,lacn,season,region) &
              * float(i)**seasonal_woods_sd(2,trait,lacn,season,region) &
              * exp(-float(i)*seasonal_woods_sd(3,trait,lacn,season,region))
          end do
          if ( DEBUGmsgs > 0 ) then
            print *, '[T][', trait, ', ', lacn, ', ', season, ',', region, ']: mean curve:', &
              seasonal_woods_means(:,trait,lacn,season,region)
            print *, '[T][', trait, ', ', lacn, ', ', season, ',', region, ']: SD curve  :', &
              seasonal_woods_sd(:,trait,lacn,season,region)
          end if
        !
        ! G = Morant and Gnanasankthy (1989) curves for SCS
        !
        case('G')
          ! Interpolate using Morant and Gnanasankthy's curve C4 for
          ! means and SD of SCS (y = a - (b*dim) + (c*dim**2)/2 + d/dim).
          ! Computed values are extremelay large for small values of DIM
          ! (<10-d), so DIM are shifted by 10-d to avoid this boundary condition.
          if ( trait < 4 .or. trait > 4 ) then
              print *, "[ERROR]: Invalid interpolation method, ", method, &
                       "provided for trait ", trait,                      &
                       ". Morant and Gnanasankthy curves are provided only for SCS."
              return
          end if
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = mandg_means(1,1,lacn,breed) &
              - mandg_means(2,1,lacn,breed) * float(i+10)      &
              + mandg_means(3,1,lacn,breed) * float(i+10)**2   &
              + mandg_means(4,1,lacn,breed) / float(i+10)
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = mandg_sd(1,1,lacn,breed)             &
              - mandg_sd(2,1,lacn,breed) * float(i+10)              &
              + ( mandg_sd(3,1,lacn,breed) * float(i+10)**2 ) / 2   &
              + mandg_sd(4,1,lacn,breed) / float(i+10)
          end do
        !
        ! S = Region-specific Morant and Gnanasankthy (1989) curves for SCS
        !
        case('S')
          ! Interpolate using Morant and Gnanasankthy's curve C4 for
          ! means and SD of SCS (y = a - (b*dim) + (c*dim**2)/2 + d/dim).
          ! Computed values are extremelay large for small values of DIM
          ! (<10-d), so DIM are shifted by 10-d to avoid this boundary condition.
          if ( trait < 4 .or. trait > 4 ) then
              print *, "[ERROR]: Invalid interpolation method, ", method, &
                       "provided for trait ", trait,                      &
                       ". Region-specific Morant and Gnanasankthy curves are provided only for SCS."
              return
          end if
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = regional_mandg_means(1,1,lacn,region) &
              - regional_mandg_means(2,1,lacn,region) * float(i+10)      &
              + regional_mandg_means(3,1,lacn,region) * float(i+10)**2   &
              + regional_mandg_means(4,1,lacn,region) / float(i+10)
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = regional_mandg_sd(1,1,lacn,region)             &
              - regional_mandg_sd(2,1,lacn,region) * float(i+10)              &
              + ( regional_mandg_sd(3,1,lacn,region) * float(i+10)**2 ) / 2   &
              + regional_mandg_sd(4,1,lacn,region) / float(i+10)
          end do
          if ( DEBUGmsgs > 0 ) then
            print *, '[S][', lacn, ']: mean curve:', regional_mandg_means(:,1,lacn,region)
            print *, '[S][', lacn, ']: SD curve  :', regional_mandg_sd(:,1,lacn,region)
          end if
        !
        ! S = Region-specific Morant and Gnanasankthy (1989) curves for SCS
        !
        case('D')
          ! Season of calving specific Morant and Gnanasankthy's curves for SCS
          if ( trait < 4 .or. trait > 4 ) then
              print *, "[ERROR]: Invalid interpolation method, ", method, &
                       "provided for trait ", trait,                      &
                       ". Region-specific Morant and Gnanasankthy curves are provided only for SCS."
              return
          end if
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = calving_mandg_means(1,1,lacn,season) &
              - calving_mandg_means(2,1,lacn,season) * float(i+10)      &
              + calving_mandg_means(3,1,lacn,season) * float(i+10)**2   &
              + calving_mandg_means(4,1,lacn,season) / float(i+10)
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = calving_mandg_sd(1,1,lacn,season)             &
              - calving_mandg_sd(2,1,lacn,season) * float(i+10)              &
              + ( calving_mandg_sd(3,1,lacn,season) * float(i+10)**2 ) / 2   &
              + calving_mandg_sd(4,1,lacn,season) / float(i+10)
          end do
        !
        ! U = Region and season of calving specific Morant and Gnanasankthy (1989) curves for SCS
        !
        case('U')
          ! Region- and season-specific Morant and Gnanasankthy's curves for SCS
          if ( trait < 4 .or. trait > 4 ) then
              print *, "[ERROR]: Invalid interpolation method, ", method, &
                       "provided for trait ", trait,                      &
                       ". Region-specific Morant and Gnanasankthy curves are provided only for SCS."
              return
          end if
          sum = 0.d0
          meanp(lacn,trait) = 0.d0
          do i = 1,maxlen
            month = min(max(1,(i+15)/30),11)
            dyield(i,trait,lacn) = seasonal_mandg_means(1,1,lacn,season,region) &
              - seasonal_mandg_means(2,1,lacn,season,region) * float(i+10)      &
              + seasonal_mandg_means(3,1,lacn,season,region) * float(i+10)**2   &
              + seasonal_mandg_means(4,1,lacn,season,region) / float(i+10)
            sum = sum + dyield(i,trait,lacn)
            meanyld(trait,i,lacn) = sum
            if ( i <= maxlen ) meanp(lacn,trait) = meanp(lacn,trait) + dyield(i,trait,lacn)*i
            sd(trait,i,lacn) = seasonal_mandg_sd(1,1,lacn,season,region)             &
              - seasonal_mandg_sd(2,1,lacn,season,region) * float(i+10)              &
              + ( seasonal_mandg_sd(3,1,lacn,season,region) * float(i+10)**2 ) / 2   &
              + seasonal_mandg_sd(4,1,lacn,season,region) / float(i+10)
          end do
          if ( DEBUGmsgs > 0 ) then
            print *, '[U][', lacn, ', ', season, ',', region, ']: mean curve:', &
              seasonal_mandg_means(:,1,lacn,season,region)
            print *, '[U][', lacn, ', ', season, ',', region, ']: SD curve  :', &
              seasonal_mandg_sd(:,1,lacn,season,region)
          end if
        !
        ! Otherwise, the user gave us an invalid method. D'oh! Since we're
        ! not (yet) doing any explicit error handling here I think that the
        ! subroutine will end up returning an array of NaN, but I'm not
        ! sure of that.
        !
        case default
            print *, "[ERROR]: Invalid interpolation method, ", method, &
              ", for trait ", trait, "!"
    end select

    return
end subroutine interpolate
!--------------------------------------------------------------------
subroutine write_curve_data(mtrait, CURVEfile, SHORTtrait, cowid, STorMTlabel, &
  maxlen, DAILYbp, tempTD, std, parity, GRAFplot, CURVEsmall, CURVEsingle)

  character*64,intent(in) :: CURVEfile
  character(1),intent(in) :: SHORTtrait(4)
  character,intent(in)    :: cowid*17
  character(2),intent(in) :: STorMTlabel(2)
  integer,intent(in)      :: mtrait, maxlen, parity
  real*8,intent(in)       :: DAILYbp(2,4,maxlen), tempTD(4,maxlen), std(2,4,maxlen)!, tempDEV(2,4,maxlen)
  integer                 :: j, l
  integer, intent(in)     :: GRAFplot(4), CURVEsmall, CURVEsingle
  character*128 :: LONGcurvefile
  character*8  :: curr_date

  ! First, form the complete file name. If CURBEsingle is 1 then we're going to write a separate
  ! file for each animal.
  call DATE_AND_TIME(curr_date)
  if ( CURVEsingle == 0 ) then
    LONGcurvefile = trim(CURVEfile)//'.'//STorMTlabel(1)//'.'//curr_date
  else
    LONGcurvefile = trim(CURVEfile)//'.'//trim(cowid)//'.'//STorMTlabel(1)//'.'//curr_date
  end if
!  print *, '[bestpred]: CURVEsingle = ', CURVEsingle
!  print *, '[bestpred]: Cow ID = ', cowid
!  print *, '[bestpred]: Cow lactation file name = ', LONGcurvefile

  ! Then open LONGcurvefile for writing.
  OPEN (FILE=LONGcurvefile, UNIT=99, POSITION='APPEND', ACTION='WRITE')

  ! Write the curves to the file.
  do j = 1,4
    if ( GRAFplot(j) /= 0 ) then
      do l = 1,maxlen
        ! FIXED: cowid, lacn, trait, DIM, TD yield
        ! Actual: daily BP, standard curve
        ! ME: daily BP, standard curve
        if ( CURVEsmall == 1 ) then
          write (99,*) cowid, ' ', parity, ' ', j, ' ', l, ' ', tempTD(j,l), ' ', &
            DAILYbp(2,j,l), ' ', std(2,j,l)
        else
          write (99,*) cowid, ' ', parity, ' ', j, ' ', l, ' ', tempTD(j,l), ' ', &
            DAILYbp(2,j,l), ' ', std(2,j,l), ' ', DAILYbp(1,j,l), ' ', std(1,j,l), ' '
        end if
      end do
    end if
  end do
  CLOSE(99)
end subroutine write_curve_data
!--------------------------------------------------------------------
subroutine write_yield_data(DATAfile, cowid, STorMTlabel, trt, ntests,  &
  X100, YLDvec, herd305, DCRvec, MTorST, PERSvec, RELpers, parity,      &
  CURVEsingle)

  character(len=*), intent(in)  :: DATAfile
  character, intent(in)     :: cowid*17
  character(2), intent(in)  :: STorMTlabel(2),MTorST(4)
  integer, intent(in)       :: ntests(4),parity
  character(7), intent(in)  :: trt(-3:8)
  real, intent(in)          :: X100(4)
  real*8, intent(in)        :: YLDvec(2,16),herd305(2,4),DCRvec(4),PERSvec(4),RELpers(4)

  character*128             :: LONGdatafile
  character*8               :: curr_date
  integer                   :: j, CURVEsingle

  ! First, form the complete file name.
  call DATE_AND_TIME(curr_date)
  if ( CURVEsingle == 0 ) then
    LONGdatafile = trim(DATAfile)//'.'//STorMTlabel(1)//'.'//curr_date
  else
    LONGdatafile = trim(DATAfile)//'.'//trim(cowid)//'.'//STorMTlabel(1)//'.'//curr_date
  end if
!  print *, '[bestpred]: Cow data file name = ', LONGdatafile

  ! Then open LONGdatafile for writing.
  OPEN (FILE=LONGdatafile, UNIT=99, POSITION='APPEND', ACTION='WRITE')

  ! Write lactation yields,reliabilities, etc. to a text file.
  ! CowID,Lac#,Trait,Tests,ME: 305-d 365-d laclen-d LTD Contemp,Actual: 305-d
  ! 365-d laclen-d LTD Contemp,DCR,Method,Pers,Rel
  do j = 1,4
    write (UNIT=99,FMT=999) cowid , parity, trt(j), ntests(j), &
      YLDvec(1,j)*X100(j),YLDvec(1,j+4)*X100(j), YLDvec(1,j+8)*X100(j), &
      YLDvec(1,j+12)*X100(j),herd305(1,j)*X100(j),YLDvec(2,j)*X100(j), &
      YLDvec(2,j+4)*X100(j),YLDvec(2,j+8)*X100(j),YLDvec(2,j+12)*X100(j), &
      herd305(2,j)*X100(j),DCRvec(j), MTorST(j),PERSvec(j),RELpers(j)*100
    999 format(a17,1x,i2,1x,a7,i4,10f7.0,f5.0,a3,1x,f5.2,f4.0)
  end do
  CLOSE(99)
end subroutine write_yield_data
!--------------------------------------------------------------------
! Assign herds to regions of the country.  Regions are the same as
! those used in:
! E. Hare, H. D. Norman and J. R. Wright. 2004. Duration of herd
! participation in Dairy Herd Improvement milk recording in the United
! States. J. Dairy Sci. 87:2743-2747.
! Input:  2-digit state (from the herd code on the Format 4) as a string
! Output: Region code as an integer.
subroutine state_to_region(herdstate, region, DEBUGmsgs)
  character*2, intent(in)  :: herdstate
  integer, intent(inout) :: region
  integer, intent(in) :: DEBUGmsgs
  !
  ! 1 Mideast      Delaware, Kentucky, Maryland, North Carolina, Tennessee,
  !                Virginia, and West Virginia
  ! 2 Midwest      Illinois, Indiana, Iowa, Michigan, Minnesota, Missouri,
  !                Ohio, and Wisconsin
  ! 3 Mountain-Prairie     Colorado, Kansas, Montana, Nebraska, North
  !                        Dakota, South Dakota, Utah, and Wyoming
  ! 4 Northeast    Connecticut, Maine, Massachusetts, New Hampshire, New
  !                Jersey, New York, Pennsylvania, Rhode Island, and
  !                Vermont
  ! 5 Northwest    Alaska, Idaho, Oregon, and Washington
  ! 6 Southeast    Alabama, Arkansas, Florida, Georgia, Louisiana,
  !                Mississippi, Oklahoma, Puerto Rico, South Carolina,
  !                and Texas
  ! 7 Southwest    Arizona, California, Hawaii, Nevada, and New Mexico
  SELECT CASE(herdstate)
    CASE('50', '51', '52', '54', '55', '61', '63')
        region = 1  ! Mideast
    CASE('31', '32', '33', '34', '35', '41', '42', '43')
        region = 2  ! Midwest
    CASE('45', '46', '47', '48', '81', '83', '84', '87')
        region = 3  ! Mountain-Prarie
    CASE('11', '12', '13', '14', '15', '16', '21', '22', '23')
        region = 4  ! Northeast
    CASE('82', '91', '92', '96')
        region = 5  ! Northwest
    CASE('56', '57', '58', '64', '65', '71', '72', '73', '74', '94')
        region = 6  ! Southeast
    CASE('85', '86', '88', '93', '95')
        region = 7  ! Southwest
    CASE DEFAULT
        region = 2  ! Midwest
        WRITE(*,*) 'Could not assign the state  ', herdstate, ' to a region. Defaulting to ', region, '.'
  END SELECT
  if ( DEBUGmsgs > 0 ) print *, '[state_to_region]: State ', herdstate, ' is in region ', region , '.'
  return
end subroutine state_to_region
!--------------------------------------------------------------------
subroutine date_to_season(freshdate, season, DEBUGmsgs)
  character*8, intent(in)  :: freshdate     !YYYYMMDD
  integer, intent(inout) :: season
  integer, intent(in) :: DEBUGmsgs
  !
  ! 1: spring (MAM), 2: summer (JJA), 3: fall (SON), 4: winter (DJF),
  SELECT CASE(freshdate(5:6))
    CASE('03', '04', '05')
        season = 1  ! Spring
    CASE('06', '07', '08')
        season = 2  ! Summer
    CASE('09', '10', '11')
        season = 3  ! Fall
    CASE('12', '01', '02')
        season = 4  ! Winter
    CASE DEFAULT
        season = 1  ! Spring
        WRITE(*,*) 'Could not assign the month  ', freshdate(5:6), ' to a season. Defaulting to ', season, '.'
  END SELECT
  if ( DEBUGmsgs > 0 ) print *, '[date_to_season]: Month ', freshdate(5:6), ' is in season ', season , '.'
  return
end subroutine date_to_season
!--------------------------------------------------------------------
