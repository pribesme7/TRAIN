c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program train
      implicit none
c-----------------------------------------------------------------------
c program for bunch train interaction simulation
c-----------------------------------------------------------------------
c
c Conventions:
c beam circulates clockwise in ring-1, anti-clockwise in ring-2
c the start of both machines is either IP3 or IP4
c IP5 must always be a BB interaction point
c slots (half-buckets) count clockwise from 0 at IP5
c buckets (filling scheme) count anti-clockwise from 0 at IP5 for ring-1
c i.e.  slotnumber = mdslt - 2*bucketnumber   mod mdslt (at start)
c                                     clockwise from 0 at IP5 for ring-2
c i.e.  slotnumber = 2*bucketnumber (at start)
c bunches count anti-clockwise from 1 at IP5 for ring-1
c                    clockwise from 1 at IP5 for ring-2 
c BB points count clockwise from 1 at leftmost parasitic at IP5
c
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2), para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur

      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     + eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      
      double precision survey_x, survey_y, survey_z,
     +   survey_x_ip, survey_y_ip, survey_z_ip,
     +   survey_sep_x, survey_sep_y 
      common /surveyf/
     +    survey_x(mcol,2), survey_y(mcol,2), survey_z(mcol,2),
     +    survey_x_ip(9,2), survey_y_ip(9,2), survey_z_ip(9,2),
     +    survey_sep_x(mcol), survey_sep_y(mcol)
      save /surveyf/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr

c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer nint
      real time1, time2
c-----------------------------------------------------------------------

      debug = 0
      print *,"next: hajimeru, A2"
      call hajimeru
 
c ask user for global data
      print *,"next: dialog2  "
      call dialog2
 
c read collision schedule
      print *, ' '
      print *, 'collision schedule:',mucoll
      call collsch1(mucoll)
      print *, ' '
 
c read optics table
      if (debug .gt. 0)  print *, 'Reading optics files . . .'
      call rdoptc(1, nint)
      if (nint .gt. mcol)  then
        print *, 'fatal: no. of interactions points ', nint,
     +  ' > ', mcol
        stop
      endif

      ninter = nint
      call rdoptc(2, nint)
      if (nint .ne. ninter)  then
        print *,
     +  'number of interaction points differ in beams 1 + 2: ',
     +  ninter, nint
        stop
      endif

c find pit names and optical data
      if (debug .gt. 0)  print *, 'call mkpits2'
      call mkpits2
 
c assemble groups of collision points
      if (debug .gt. 0)  print *, 'call assoc'
      call assoc
 
c read survey files
      if (debug .gt. 0)  print *, 'Reading survey files . . .'
      call rdsurvey(1, nint)
      if (nint .ne. ninter)  then
        print *,
     +  'number of interaction points differ in SURVEY B1: ',
     +  ninter, nint
        stop
      endif

      call rdsurvey(2, nint)
      if (nint .ne. ninter)  then
        print *,
     +  'number of interaction points differ in SURVEY B2: ',
     +  ninter, nint
        stop
      endif
      
      if (debug .gt. 0)  print *, 'Calculating survey separation . . .'
      CALL calcsurvey

c initialize tables for transfer maps
      print *, 'Reading map files 2 . . .'
      call rdmaps2
c set up equivalence classes + collision scheme
      if (debug .gt. 0)  print *, 'call prcoll'
      call prcoll
c prepare bunch numbers
      if (debug .gt. 0)  print *, 'call set_cuem'
      call set_cuem
     
 
c initial guess for displacements in collision points
      if (debug .gt. 0)  print *, 'call orbit02'
      call orbit02
      
c output orbits
      if (all_write)  then
        call orbw
        close(orbout)
      endif
c calculate eigenvectors at observation point
      if (nturns .gt. 0)  then
        if (debug .gt. 0)  print *, 'call eicoll'
        call eicoll
      endif

      if (c_orbit)  then
c analyse tunes
        if (debug .gt. 0)  print *, 'call mktunemad2'
        call mktunemad2(.false.)
        
        if (all_write)  call print2(.false., .true.)
        
c attempt to find the closed orbits
        if (debug .gt. 0)  print *, 'call orbitb2'
        call orbitb2
        if (debug .gt. 0)  print *, 'call print2'
        if (all_write)  call print2(.true., .false.)

        if(c_tunes .gt. 0)  then
c analyse tunes
          if (debug .gt. 0)  print *, 'call mktune'
          call mktunemad2(.true.)

          if (all_write)  call print2(.true., .true.)
        endif

c--- multi-turn tracking
        if (nturns .gt. 0)  then
c--- initialize track orbits for all bunches
          if (debug .gt. 0)  print *, 'call initrack'
          call initrack
c track bunches
          print *, 'start strong-strong tracking and timing ..'
          if (debug .gt. 0)  print *, 'call bottrack'
!          call timex(time1)
          do c_turn = 1, nturns
            open(13, file='kick.out', status = 'UNKNOWN')
            call bottrack
            close(13)
          enddo
 !         call timex(time2)
          print *, 'time for ', nturns, ' turns: ', time2-time1,
     +    ' sec'
c print result
          if (debug .gt. 0)  print *, 'call pbunch'
          call pbunch
        endif
      endif
      stop
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine assoc
      implicit none
      integer i,j,ip,m,n,kpar
c-----------------------------------------------------------------------
c associates the collision points with pits and sets global counters.
c
c The algorithm was copied from E. Keil's "orbit8" program.
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol,twompar
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
c both sides parasitic interactions
      parameter (twompar = 2*mpar)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter,lnpar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
     
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/

      common /selecip/ onepit,ipstart,intip
      save /selecip/ 
      logical onepit
      integer intip
      character*4 ipstart
c-----------------------------------------------------------------------
      integer headon, parasit,nintip
      integer inpar
      integer lnpar(8)
c-----------------------------------------------------------------------
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/ 
c associate collision point with an IP by requiring that it is
c less than 500 m from the IP
      do 10 n = 1, npit
         ncoll(n) = 0
   10 continue
      
      do n = 1,8
        lnpar(n) = 0
      enddo
 
      do 30 i = 1, ninter
         kpar = parasit(name(i,1),1)
         if (headon(name(i,1)) .gt. 0 .or. kpar .gt. 0) then
           if (kpar .gt. 0) then
             ip = 0
             do 20 m = 1, npit
               if (abs(s(i,1)-si(m)) .lt. 500.0) ip = m
               
   20        continue
             
c check association
             if (ip .eq. 0) then
               print *, 'Error - ',
     +              'Collision point ', i, ' not associated with pit'
               stop
             else
               ncoll(ip) = ncoll(ip) + 1

               read(name(ipit(ip),1)(5:5),'(i1)') inpar
               lnpar(inpar) = ncoll(ip)
             endif
 
             if (ncoll(ip) .gt. mlocal) then
               print *, 'Error - ',
     +         'at pit', ip, 'too many collision points:', ncoll(ip)
               stop
             endif
           endif
 
         endif
   30 continue

      do i = 1,8
         lnpar(i) = lnpar(i)/2
         if (lnpar(i) .ne. 0) then
           print *,'total no. of one-side parasitic', 
     + ' in IP',i,': ', lnpar(i)
         endif
      enddo
  
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine bbtrac2(bnum1, bnum2, intp, fact)
c-----------------------------------------------------------------------
c--- tracks one bunch pair over int. point - orig. orbit
c  input
c  bnum1    bunch number ring-1
c  bnum2    bunch number ring-2
c  intp     interaction point
c  fact     factor for long-range
c-----------------------------------------------------------------------
      implicit none
      integer bnum1, bnum2, intp
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      double precision fact,xs,ys,ccp,sx,sy,phix0, phiy0, phix, phiy
      double precision sigx, sigy

      double precision kick
      common/ kck/ kick(2,mbunch,mcol,2)
      save /kck/
 
      xs = z1a(1,bnum1,intp) - z2a(1,bnum2,intp)
      ys = z1a(3,bnum1,intp) - z2a(3,bnum2,intp)
      ccp = fact * xisign * sqrt(bcurr1(bnum1) * bcurr2(bnum2))
     +      / (ech * frev)
      sx = sqrt(sigx(bnum1,intp,1)**2 + sigx(bnum2,intp,2)**2)
      sy = sqrt(sigy(bnum1,intp,1)**2 + sigy(bnum2,intp,2)**2)
      call bbtr2(ccp, sx, sy, xs, ys, phix0, phiy0)
      xs = xs + ztr(1,bnum1,1) - ztr(1,bnum2,2)
      ys = ys + ztr(3,bnum1,1) - ztr(3,bnum2,2)
      call bbtr2(ccp, sx, sy, xs, ys, phix, phiy)
      ztr(2,bnum1,1) = ztr(2,bnum1,1) + phix - phix0
      ztr(4,bnum1,1) = ztr(4,bnum1,1) + phiy - phiy0
      ztr(2,bnum2,2) = ztr(2,bnum2,2) - phix + phix0
      ztr(4,bnum2,2) = ztr(4,bnum2,2) - phiy + phiy0
      
        write(13,*) bnum1, ztr(2,bnum1,1),ztr(4,bnum1,1),bnum2,
     +ztr(2,bnum2,2),ztr(4,bnum2,2)
     
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine bbtr2(ccp, sx, sy, xs, ys, phix, phiy)
      implicit none
      double precision cbx,cby,ccp,crx,cry,exk,exkc,explim,fk,
     +phix,phiy,r,r2,rho2,rk,sx,sx2,sy,sy2,tk,xb,xr,xs,yb,yr,ys
c-----------------------------------------------------------------------
c transport map for beam-beam element
c-----------------------------------------------------------------------
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
c     if x > explim, exp(-x) is outside machine limits.
      parameter         (explim = 150.0d0)
      double precision newtol
      parameter (newtol = 1.0d-12)
c-----------------------------------------------------------------------
 
c factor for beam-beam kick
      phix = 0
      phiy = 0
      fk = two * tradius * ccp / gamma

      if (fk .eq. 0.d0)  return
      sx2 = sx * sx
      sy2 = sy * sy
 
c limit formulas for sigma(x) = sigma(y)
c      if (sx2 .eq. sy2) then
      if (abs(sx2-sy2) .le. newtol) then 
         rho2 = xs * xs + ys * ys
 
c do nothing for xs = ys = 0
c general case
         if (rho2 .ne. 0.0) then
            tk = rho2 / (two * sx2)
            if (tk .gt. explim) then
               exk  = 0.0
               exkc = one
               phix = xs * fk / rho2
               phiy = ys * fk / rho2
            else
               exk  = exp(-tk)
               exkc = one - exk
               phix = xs * fk / rho2 * exkc
               phiy = ys * fk / rho2 * exkc
            endif
         endif
c case sigma(x) > sigma(y)
      else
         r2 = two * (sx2 - sy2)
         if (sx2 .gt. sy2) then
            r  = sqrt(r2)
            rk = fk * sqrt(pi) / r
            xr = abs(xs) / r
            yr = abs(ys) / r
            call errf(xr, yr, crx, cry)
            tk = (xs * xs / sx2 + ys * ys / sy2) / two
            if (tk .gt. explim) then
               exk = 0.0
               cbx = 0.0
               cby = 0.0
            else
               exk = exp(-tk)
               xb  = (sy / sx) * xr
               yb  = (sx / sy) * yr
               call errf(xb, yb, cbx, cby)
            endif
 
c case sigma(x) < sigma(y)
         else
            r  = sqrt(-r2)
            rk = fk * sqrt(pi) / r
            xr = abs(xs) / r
            yr = abs(ys) / r
            call errf(yr, xr, cry, crx)
            tk = (xs * xs / sx2 + ys * ys / sy2) / two
            if (tk .gt. explim) then
               exk = 0.0
               cbx = 0.0
               cby = 0.0
            else
               exk = exp(-tk)
               xb  = (sy / sx) * xr
               yb  = (sx / sy) * yr
               call errf(yb, xb, cby, cbx)
            endif
         endif
         phix = rk * (cry - exk * cby) * sign(one, xs)
         phiy = rk * (crx - exk * cbx) * sign(one, ys)
 
      endif
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision function bc_fun(x)
      implicit none
      double precision x
      real rn32
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
c      bc_fun = 1 + 2 * (rn32() - 0.5) * x
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine bottrack
c-----------------------------------------------------------------------
c--- tracks all bunches one full turn from initial position
c-----------------------------------------------------------------------
      implicit none
      integer i,j,j2,n,np,nb
      integer outbunch
      double precision fact
      logical pit_number
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      
      double precision temp(2,mbunch,2)
      common /tempc/ temp
      save /tempc/
      do i = 1,nbunch
        do j = 1,2
          do n = 1,2
            temp(j,i,n) = 0 
          enddo
        enddo
      enddo 
c--- loop over buckets
      do i = 1, mdslt
c--- loop over bunches, ring-1
        do j = 1, nbunch
          if (present(j,1) .ne. 0)  then
c--- increase all bunch positions, track if new position is BB
            np = tcount(j,1)
            tcount(j,1) = np + 1
            if (tcount(j,1) .eq. mdslt) tcount(j,1) = 0
            n = tcount(j,1)
            if (maskmi(n) .ne. 0)  then
              call gaptrack(j, maskmp(np), maskmp(n), 1)
            endif
          endif
        enddo
c--- loop over bunches, ring-2
        do j = 1, nbunch
          if (present(j,2) .ne. 0)  then
c--- decrease all bunch positions, track if new position is BB
            np = tcount(j,2)
            tcount(j,2) = np - 1
            if (tcount(j,2) .lt. 0) tcount(j,2) = mdslt-1
            n = tcount(j,2)
            if (maskmi(n) .ne. 0)  then
              call gaptrack(j, maskmn(np), maskmn(n), -1)
            endif
          endif
        enddo
c--- track bunch pairs, change orbit
        do j = 1, nbunch
          n = tcount(j,1)
          if (present(j,1) .ne. 0 .and. maskmi(n) .ne. 0)  then
            j2 = ibnch1(j, maskmi(n))
            if (j2 .ne. 0)  then
              if (present(j2,2) .ne. 0) then
                if (pit_number(maskmi(n)))  then
                  fact = hofact
                else
                  fact = xifact
                endif
                if (maskmi(n) .eq. outpos) then
                  ztr(1,j,1) = ztr(1,j,1) + orb_amp(1,j,1)
                  ztr(3,j,1) = ztr(3,j,1) + orb_amp(2,j,1)
                  ztr(1,j2,2) = ztr(1,j2,2) + orb_amp(1,j2,2)
                  ztr(3,j2,2) = ztr(3,j2,2) + orb_amp(2,j2,2)
                  orb_amp(1,j,1) = zero
                  orb_amp(2,j,1) = zero
                  orb_amp(1,j2,2) = zero
                  orb_amp(2,j2,2) = zero
                endif
                call bbtrac2(j, j2, maskmi(n), fact)
                temp(1,j,1) = temp(1,j,1) + ztr(2,j,1)
                temp(2,j,1) = temp(2,j,1) + ztr(4,j,1)
                temp(1,j2,2) = temp(1,j2,1) + ztr(2,j2,2)
                temp(2,j2,2) = temp(2,j2,2) + ztr(4,j2,2)
                if (maskmi(n) .eq. outpos) call trstat(j, j2)
              endif
            endif
            nb = outbunch(j)
            if (nb .gt. 0 .and. maskmi(n) .eq. outpos)
     +      call wtrack(j, j2, nb)
          endif
        enddo
      enddo   

c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine check(flag, a)
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      integer flag,j
      double precision a(6,6), v(6), loc(6,6)
      save loc
      if (flag .eq. 0)  then
        call mxone(loc, 6, 6)
      elseif (flag .eq. 1)  then
        call mxmpy(a, loc, loc, 6, 6)
      elseif (flag .eq. 2)  then
        print '(a,2x,1p,6e12.4)','initial: ', (a(j,1), j = 1, 6)
        call mxbyv(loc, a, v)
        print '(a,2x,1p,6e12.4)','final:   ', (v(j), j = 1, 6)
      else
        print *, 'check: wrong flag: ', flag
        stop
      endif
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine collide(mbuck, mdslt, mcol, b2off, npar, a, mask,
     +maski, colcnt_f, list_f, part_f, colcnt_b, list_b, part_b)
      implicit none
      integer mbuck, mdslt, mcol, npar(8), b2off, a(2,0:mbuck-1),
     +mask(0:mdslt-1), maski(0:mdslt-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1)
      integer i, j, k, fpart, bpart, ppos,mpar,counttemp
      counttemp = 0
      do fpart = 0, mbuck-1

        if (a(1,fpart) .ne. 0)  then
c--- ppos is the position in half buckets with respect to IP5
c    of bunch fpart of beam_1 when it hits bunch 0 of beam 2;
c    b2off is the off_set of beam_2 (in half-buckets) with respect to beam_1
c    at IP5 (negative towards IP4, positive towards IP6).
          ppos = b2off - 2 * fpart
          do i = 0, mdslt-1
            k = mod(ppos, mdslt)
            if (k .lt. 0) k = k + mdslt
            if (mask(k) .ne. 0)  then

              call selipnum(npar, mcol,mdslt,k,mpar) ! One side number of parasitic collisions per IP
              j = k + i -2*mpar

              if (j .lt. 0)     j = j + mdslt
              if (j .ge. mdslt) j = j - mdslt
              bpart = j / 2

              if (a(2,bpart) .ne. 0) then
                colcnt_f(fpart) = colcnt_f(fpart) + 1
                colcnt_b(bpart) = colcnt_b(bpart) + 1
                list_f(colcnt_f(fpart),fpart) = maski(k)
                list_b(colcnt_b(bpart),bpart) = maski(k)
                part_f(colcnt_f(fpart),fpart) = bpart
                part_b(colcnt_b(bpart),bpart) = fpart
              endif
            endif
            ppos = ppos + 1
          enddo
        endif

      enddo

c--- sort for collision point

      do fpart = 0, mbuck-1
        if (colcnt_f(fpart) .gt. 0)  then
          call mysort(colcnt_f(fpart), list_f(1,fpart), part_f(1,fpart))
        endif
      enddo
      do bpart = 0, mbuck-1
        if (colcnt_b(bpart) .gt. 0)  then
          call mysort(colcnt_b(bpart), list_b(1,bpart), part_b(1,bpart))
        endif
      enddo

      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine collsch1(unit)
      implicit none
      integer groups, lastnb
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
      double precision bc1,bc2,epsxb1,epsxb2,epsyb1,epsyb2
c   for bunch currents
      common /bcc/    bc1(mbunch), bc2(mbunch), epsxb1(mbunch), 
     +     epsxb2(mbunch), epsyb1(mbunch), epsyb2(mbunch)
c   ztr trajectory
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist      
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      SAVE /flagsi/
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer unit
      character * 100 comment
      integer i, j, k, tmp(2,100)
      real             xmp(2,100)
      integer ii, iii
      integer n
 
      i = 0
      n = 1
   10 continue
      if (emitt_f .eq. 8)  then
         read (unit,*,end = 20) ii,collsk(1,i),collsk(2,i),xmp(1,2),
     +       xmp(1,3),xmp(1,4),xmp(1,5),xmp(1,6),xmp(1,7)
      else
         read (unit,*,end = 20) ii,collsk(1,i),collsk(2,i),xmp(1,2),
     +       xmp(1,3)
      endif
          iii = i + 1
          write(77,*) iii,collsk(1,i),collsk(2,i)
          if ((collsk(1,i).ne.0) .or. (collsk(2,i).ne.0)) then
             bc1(n) = xmp(1,2)
             bc2(n) = xmp(1,3)
             if (collsk(1,i).ne.0) then
                epsxb1(n) = 1.0
                epsyb1(n) = 1.0
             else
                epsxb1(n) = 0.0
		epsyb1(n) = 0.0
             endif
	     
             if (collsk(2,i).ne.0) then
                epsxb2(n) = 1.0
                epsyb2(n) = 1.0
             else
                epsxb2(n) = 0.0
		epsyb2(n) = 0.0
             endif

	     if (emitt_f .eq. 8)  then
	       epsxb1(n) = xmp(1,4)
	       epsxb2(n) = xmp(1,5)
	       epsyb1(n) = xmp(1,6)
	       epsyb2(n) = xmp(1,7)        
               write(79,*) iii,collsk(1,i),collsk(2,i),n,bc1(n),bc2(n),
     +         epsxb1(n),epsxb2(n),epsyb1(n),epsyb2(n)
             else
	       write(79,*) iii,collsk(1,i),collsk(2,i),n,bc1(n),bc2(n)
	     endif
             n = n + 1
          endif
          i = i + 1
      goto 10
   20 continue
      if (i .ne. mbuck)  then
        print *, 'fatal: collision count, slots: ', i, mbuck,n
       stop
      endif
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ddcopy(d_in, d_out,d_count)
c--- copies double precision to double precision arrays
      double precision d_in(*), d_out(*)
      integer i, d_count
      do i = 1, d_count
        d_out(i) = d_in(i)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine drcopy(d_in, d_out, d_count)
c--- copies double to real precision arrays
      double precision d_in(*)
      real             d_out(*)
      integer i, d_count
      do i = 1, d_count
        d_out(i) = d_in(i)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rdcopy(d_in, d_out,d_count)
c--- copies real to double precision arrays
      real             d_in(*)
      double precision d_out(*)
      integer i, d_count
      do i = 1, d_count
        d_out(i) = d_in(i)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dzero(d_in, d_count)
c--- zeros double precision arrays
      double precision d_in(*)
      integer i, d_count
      do i = 1, d_count
        d_in(i) = 0.d0
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine decode(token, jtok, ntok, number, rval, sval)
      implicit none
      integer idig,ilen,ise,ive,jtok,length,npl,ntok,num
      double precision rval,sig,val
c-----------------------------------------------------------------------
c decode one table field according to format code and length
c-----------------------------------------------------------------------
 
      logical number
      character*1 token(*)
      character*8 sval
      logical dig, pnt
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
c-----------------------------------------------------------------------
c skip leading blanks
   10 if (jtok .le. ntok  .and.  token(jtok) .eq. ' ') then
        jtok = jtok + 1
        go to 10
      endif
 
      rval = zero
      sval = ' '
 
      if (jtok .le. ntok) then
         if (number) then
c any numeric character ?
            if (index('0123456789+-.', token(jtok)) .ne. 0) then
               val = zero
               sig = one
               ive = 0
               ise = 1
               npl = 0
               dig = .false.
               pnt = .false.
 
c sign
               if (token(jtok) .eq. '+') then
                  jtok = jtok + 1
               else if (token(jtok) .eq. '-') then
                  sig = - one
                  jtok = jtok + 1
               endif
 
c digit or decimal point?
 20            idig = index('0123456789', token(jtok))
 
               if (idig .ne. 0) then
                  val = ten * val + float(idig-1)
                  dig = .true.
                  if (pnt) npl = npl + 1
                  jtok = jtok + 1
                  go to 20
               else if (token(jtok) .eq. '.') then
                  pnt = .true.
                  jtok = jtok + 1
                  go to 20
               endif
 
c exponent ?
               if (jtok .le. ntok  .and.
     +             index('DEde', token(jtok)) .ne. 0) then
                  jtok = jtok + 1
                  dig = .false.
 
                  if (token(jtok) .eq. '+') then
                     jtok = jtok + 1
                  else if (token(jtok) .eq. '-') then
                     ise = -1
                     jtok = jtok + 1
                  endif
 
 30               idig = index('0123456789', token(jtok))
 
                  if (idig .ne. 0) then
                     ive = 10 * ive + idig - 1
                     dig = .true.
                     jtok = jtok + 1
                     go to 30
                  endif
               endif
 
c return value
               rval = sig * val * ten ** (ise * ive - npl)
               if (token(jtok) .eq. ',') jtok = jtok + 1
            endif
c string
         else
            ilen = 0
 
            if (token(jtok) .eq. '"') then
               jtok = jtok + 1
               length = len(sval)
               num = 0
 
 40            if (jtok .le. ntok) then
                  if (token(jtok) .eq. '"') then
                     jtok = jtok + 1
                     return
                  else if (num .lt. length) then
                     ilen = ilen + 1
                     sval(ilen:ilen) = token(jtok)
                  endif
 
                  jtok = jtok + 1
                  go to 40
               endif
 
c otherwise, string is up to next blank or comma
            else
 50            if (jtok .lt. ntok) then
                  if (ilen .lt. len(sval)) then
                     ilen = ilen + 1
                     sval(ilen:ilen) = token(jtok)
                  endif
 
                  jtok = jtok + 1
 
                  if (token(jtok) .eq. ',') then
                     jtok = jtok + 1
                     return
                  else if (token(jtok) .eq. ' ') then
                     return
                  endif
 
                  go to 50
               endif
            endif
         endif
      endif
 
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dialog2
      implicit none
      double precision seed
      integer i
      character * 120 nxline
c-----------------------------------------------------------------------
c ask user for optional data
c-----------------------------------------------------------------------
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol,twompar
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
c both sides parasitic interactions
      parameter (twompar = 2*mpar)

      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      character*120 filename, chst,filenm
      common /filesch/ filenm,chst
      save / filesch/
      integer idum(12)
      integer idumtemp
      
      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/     
c-----------------------------------------------------------------------
c open collision schedule input file
      filename = nxline()
      filenm = filename
      open(mucoll, file = filename, status = 'UNKNOWN')
      chst = nxline()
      read (chst, *) idum
      write(*,*) chst
      all_write= idum(1) .ne. 0
      w_coll   = all_write .and. idum(2) .ne. 0
      w_frequ  = all_write .and. idum(3) .ne. 0
      w_equ    = all_write .and. idum(4) .ne. 0
      w_set    = all_write .and. idum(5) .ne. 0
      w_alt    = all_write .and. idum(6) .ne. 0
      c_orbit  = idum(7) .ne. 0
      f_second = idum(8) .ne. 0
      w_detail = idum(9) .ne. 0
      beamc_f  = idum(10)
      emitt_f  = idum(11)
      chst = nxline()
      read (chst, *) idum
      f_coll   = idum(1) .ne. 0
      nturns   = idum(2)
      debug    = idum(3)
      outbcnt  = min(idum(4),max_list)
      outpos   = idum(5)
      outnorm  = idum(6)
      xifact   = idum(7)
      hofact   = idum(8)
      amp_bunch = idum(9)
      amp_fac  = idum(10)
      lumi_hist = idum(11) .ne. 0
      b2_off = idum(12)
      chst = nxline()
      read (chst, *) idum
c--- list of output bunches
      do i = 1, outbcnt
        outblist(i) = idum(i)
      enddo
c open luminosity histogram file (from tracking)
      if (lumi_hist)
     +open(lumilist, file = 'hist.list', status = 'UNKNOWN')
c open orbit output file
      filename = nxline()
      if (all_write)
     +open(orbout, file = filename, status = 'UNKNOWN')
c open list output file
      filename = nxline()
      if (all_write)
     +open(mulist, file = filename, status = 'UNKNOWN')
      chst = nxline()
      read (chst, *) c_tunes
      print *, 'coherent TUNE requested for ', c_tunes, ' bunches'
      chst = nxline()
      read(chst, *) ampx, ampy
      print *, 'start amplitude sigmas x(1,2) + y(1,2): ', ampx, ampy
      chst = nxline()
      read(chst, *) sigb, sigem
      print *, 'beam current + emittance sigma: ', sigb, sigem
      bcfile = .false.
      chst = nxline()
      read(chst, *) seed
      if (seed .lt. 1.e5  .or.  seed .gt. 1.e10) then
        seed = 12345678.
        print *, 'Seed should be between 1.0e5 and 1.0e10 - set to: ',
     +  seed
      else
        print *, 'random seed: ', seed
      endif
      iseed = seed
c Extra elements on the maps, the default ones are the crab cavities and the triplets
      chst = nxline()
      read (chst,*) idumtemp
      extraelem = idumtemp .ne. 0
      print *, 'Use extra elements? ', extraelem     
      end
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dialog
      implicit none
      double precision seed
      integer i
      character * 120 nxline
c-----------------------------------------------------------------------
c ask user for optional data
c-----------------------------------------------------------------------
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol,twompar
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
c both sides parasitic interactions
      parameter (twompar = 2*mpar)

      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      common /selecip/ onepit, ipstart,intip
      save / selecip/
      
      integer intip
      character*120 filename, chst,filenm
      common /filesch/ filenm,chst
      save / filesch/
      integer idum(12)
      logical onepit
      character*4 ipstart
c-----------------------------------------------------------------------
c open collision schedule input file
      filename = nxline()
      filenm = filename
      open(mucoll, file = filename, status = 'UNKNOWN')
      chst = nxline()
      read (chst, *) idum
      write(*,*) chst
      all_write= idum(1) .ne. 0
      w_coll   = all_write .and. idum(2) .ne. 0
      w_frequ  = all_write .and. idum(3) .ne. 0
      w_equ    = all_write .and. idum(4) .ne. 0
      w_set    = all_write .and. idum(5) .ne. 0
      w_alt    = all_write .and. idum(6) .ne. 0
      c_orbit  = idum(7) .ne. 0
      f_second = idum(8) .ne. 0
      w_detail = idum(9) .ne. 0
      beamc_f  = idum(10)
      emitt_f  = idum(11)
      chst = nxline()
      read (chst, *) idum
      f_coll   = idum(1) .ne. 0
      nturns   = idum(2)
      debug    = idum(3)
      outbcnt  = min(idum(4),max_list)
      outpos   = idum(5)
      outnorm  = idum(6)
      xifact   = idum(7)
      hofact   = idum(8)
      amp_bunch = idum(9)
      amp_fac  = idum(10)
      lumi_hist = idum(11) .ne. 0
      b2_off = idum(12)
      chst = nxline()
      read (chst, *) idum
c--- list of output bunches
      do i = 1, outbcnt
        outblist(i) = idum(i)
      enddo
c open luminosity histogram file (from tracking)
      if (lumi_hist)
     +open(lumilist, file = 'hist.list', status = 'UNKNOWN')
c open orbit output file
      filename = nxline()
      if (all_write)
     +open(orbout, file = filename, status = 'UNKNOWN')
c open list output file
      filename = nxline()
      if (all_write)
     +open(mulist, file = filename, status = 'UNKNOWN')
      chst = nxline()
      read (chst, *) c_tunes
      print *, 'coherent TUNE requested for ', c_tunes, ' bunches'
      chst = nxline()
      read(chst, *) ampx, ampy
      print *, 'start amplitude sigmas x(1,2) + y(1,2): ', ampx, ampy
      chst = nxline()
      read(chst, *) sigb, sigem
      print *, 'beam current + emittance sigma: ', sigb, sigem
      bcfile = .false.
      chst = nxline()
      read(chst, *) seed
      if (seed .lt. 1.e5  .or.  seed .gt. 1.e10) then
        seed = 12345678.
        print *, 'Seed should be between 1.0e5 and 1.0e10 - set to: ',
     +  seed
      else
        print *, 'random seed: ', seed
      endif
      iseed = seed
      print *,'before seed'
c     call rndm(iseed)
      print *,'after  seed'
     
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine disp(r, aux, d)
      implicit none
      integer i,irank,j
      real d(6)
      double precision a(4,5),aux(6),r(6,6)
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
c   initial values for dispersion or its first derivative by delta
c   only the first four components of the vectors are set
c-----------------------------------------------------------------------
 
      do 20 i = 1, 4
         do 10 j = 1, 4
            a(i,j) = r(i,j)
 10      continue
 
         a(i,i) = a(i,i) - 1.0
         a(i,5) = - aux(i)
 20   continue
 
!      if (debug .gt. 0)  print *, 'call solver in disp'
      call solver(a, 4, 1, irank)
 
      if (irank .ge. 4) then
         do 30 i = 1, 4
            d(i) = a(i,5)
 30      continue
      else
         print *, 'Unable to compute dispersion - set to zero'
 
         do 40 i = 1, 4
            d(i) = 0.0
 40      continue
      endif
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine eicoll
      implicit none
      integer i, iobs
      logical fsec
      double precision q1, q2, rt(6,6),orbit(6)
c-----------------------------------------------------------------------
c find eigenvectors of one-turn matrix at observation point
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      fsec = .true.
      do iobs = 1, outbcnt
c initialise map for forward bunch
        call mxone(tr1(1,1,1), 6, 6)
        call dzero(tt1(1,1,1,1), 216)
        call rdcopy(z1(1,outblist(iobs),outpos), orbit, 6)
c track through sectors outpos to end and catenate map
        do i = outpos, ninter
           call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i), orbit, rt)
           call mapcat(fsec, rt, st1(1,1,1,i),
     +        tr1(1,1,1), tt1(1,1,1,1), tr1(1,1,1), tt1(1,1,1,1))
        enddo
 
c track through sectors 0 to outpos-1 and catenate map
        do i = 0, outpos-1
           call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i), orbit, rt)
           call mapcat(fsec, rt, st1(1,1,1,i),
     +        tr1(1,1,1), tt1(1,1,1,1), tr1(1,1,1), tt1(1,1,1,1))
        enddo
c        call macheck(tr1(1,1,1))
        call ddcopy(tr1(1,1,1), eiv1(1,1,iobs), 36)
        call eigen(eiv1(1,1,iobs), 6, q1, q2)
c initialise map for backward bunch
        call mxone(tr2(1,1,1), 6, 6)
        call dzero(tt2(1,1,1,1), 216)
        call rdcopy(z2(1,outblist(iobs),outpos), orbit, 6)
        call rdcopy(z2(1,outblist(iobs),outpos), orb0_2(1,iobs), 6)
 
c track through sectors outpos-1 to 0 and catenate map
        do i = outpos-1, 0, -1
           call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i), orbit, rt)
           call mapcat(fsec, rt, st2(1,1,1,i),
     +        tr2(1,1,1), tt2(1,1,1,1), tr2(1,1,1), tt2(1,1,1,1))
        enddo
 
c track through sectors ninter to outpos and catenate map
        do i = ninter, outpos, -1
           call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i), orbit, rt)
           call mapcat(fsec, rt, st2(1,1,1,i),
     +        tr2(1,1,1), tt2(1,1,1,1), tr2(1,1,1), tt2(1,1,1,1))
        enddo
c        call macheck(tr2(1,1,1))
        call ddcopy(tr2(1,1,1), eiv2(1,1,iobs), 36)
        call eigen(eiv2(1,1,iobs), 6, q1, q2)
      enddo
c-----------------------------------------------------------------------
      end
      subroutine eigen(a, n, q1, q2)
c-----------------------------------------------------------------------
c return eigenvalues and eigenvectors of a n*n matrix
c-----------------------------------------------------------------------
      implicit none
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer n, ilo, ihi, ierr
      double precision a(n,n)
      double precision q1, q2
      integer mvary
      parameter (mvary = 6)
c-----------------------------------------------------------------------
      double precision reval(mvary), aival(mvary), v(mvary,mvary)
      double precision d(mvary,2)
c-----------------------------------------------------------------------
c compute eigenvalues and eigenvectors
      ilo = 1
      ihi = n
      call orthes(mvary, n, ilo, ihi, a, d)
      call ortran(mvary, n, ilo, ihi, a, d, v)
      call hqr2(mvary, n, ilo, ihi, a, reval, aival, v, ierr)
      
      if (ierr .ne. 0) then
        print *, 'Error in "train" program - ',
     +        'Unable to find eigenvalues for matrix'
c       stop
      endif
      q1 = atan2(aival(2), reval(2)) / (2.0 * pi)
      q2 = atan2(aival(4), reval(4)) / (2.0 * pi)

c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine gaptrack(bnum, inpt1, inpt2, dir)
c-----------------------------------------------------------------------
c--- tracks one bunch over interaction gap
c-----------------------------------------------------------------------
c  input
c  bnum     bunch number
c  inpt1    start of tracking
c  inpt2    end of tracking
c  dir      direction: 1 ring-1, -1 ring-2
c-----------------------------------------------------------------------
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      integer i, bnum, inpt1, inpt2, dir
      double precision rdummy(6,6)
      if (dir .gt. 0)  then
c--- forward beam
        if (inpt2 .gt. inpt1)  then
c--- track through sectors inpt1 -> inpt2-1
          do i = inpt1, inpt2-1
            call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i),
     +      ztr(1,bnum,1), rdummy)
          enddo
        else
c--- track through inpt1 -> last sector -> sector 0 -> inpt2-1
          do i = inpt1, ninter
            call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i),
     +      ztr(1,bnum,1), rdummy)
          enddo
          do i = 0, inpt2-1
            call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i),
     +      ztr(1,bnum,1), rdummy)
          enddo
        endif
      else
c--- backward beam
        if (inpt1 .gt. inpt2)  then
c--- track through sector inpt1-1 -> inpt2
          do i = inpt1-1, inpt2, -1
            call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i),
     +      ztr(1,bnum,2), rdummy)
          enddo
        else
c--- track through sector inpt1-1 -> sector 0 -> last sector -> inpt2
          do i = inpt1-1, 0, -1
            call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i),
     +      ztr(1,bnum,2), rdummy)
          enddo
          do i = ninter, inpt2, -1
            call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i),
     +      ztr(1,bnum,2), rdummy)
          enddo
        endif
      endif
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine hajimeru
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      integer i
      root2 = sqrt(two)
      lumicnt = 0
      lumiav = 0
      n_parasit = 0
      do i = 1, 8
        iact(i) = 0
        ippos(i) = -1000.
      enddo
      do i = 0, mbunch
        write(code(i), '(i4)') i
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function headon(name)
      implicit none
      integer lastnb
      character *(*) name
      if (name(1:4) .eq. 'MKIP' .and. lastnb(name) .eq. 5)  then
        read (name(5:5), '(i1)') headon
      else
        headon = 0
      endif
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine initrack
c-----------------------------------------------------------------------
c--- initializes strong-strong tracking
c-----------------------------------------------------------------------
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      integer i, j, low, up, step
      real rg32cut
      double precision vx, vy
      double precision sigx, sigy
      character * 1000 fn
 
c--- open bunch output units
      fn = 'bunch.'
      do i = 1, outbcnt
        write(fn(7:10), '(i4.4)') outblist(i)
        open(mtrack+i,file=fn, status = 'UNKNOWN')
      enddo
c--- set factor for luminosity
      lumifact = 1.d0 / (4 * sqrt(sigx(1,outpos,1) * sigy(1,outpos,1)
     + * sigx(1,outpos,2) * sigy(1,outpos,2)))
c--- loop over two beams, reset map kicks
      do i = 0, ninter-1
        call dzero(sk1(1,i), 6)
      enddo
      do i = ninter, 1, -1
        call dzero(sk2(1,i), 6)
      enddo
c--- reset c.o.
      call dzero(ztr, 12*mbunch)
c--- add amplitudes to closed orbit
      if (amp_bunch .eq. 0)  then
        low = 1
        up = nbunch
        step = 1
      elseif (amp_bunch .lt. 0)  then
        low = -amp_bunch
        up  = nbunch
        step = -amp_bunch
      else
        low = amp_bunch
        up = amp_bunch
        step = 1
      endif
      do j = low, up, step
        do i = 1, 2
          vx = sigx(j,outpos,i) * ampx(i)
          vy = sigy(j,outpos,i) * ampy(i)
          if (amp_fac .eq. 0)  then
            vx = vx * rg32cut(3.)
            vy = vy * rg32cut(3.)
          endif
          orb_amp(1,j,i) = vx
          orb_amp(2,j,i) = vy
        enddo
      enddo
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine macheck(a)
      implicit none
      double precision q1, q2, a(6,6),s(6,6),b(6,6),c(6,6),d(6,6)
      integer i,j
      if (.false.)  then
        do i = 1,6
          do j = 1,6
            s(j,i) = 0.d0
            c(j,i) = a(i,j)
          enddo
        enddo
        do i = 1,3
          s(2*i-1,2*i) = -1.d0
          s(2*i,2*i-1) =  1.d0
        enddo
        print *, 'symplecticity check:'
        call mxmpy(s,a,b,6,6)
        call mxmpy(c,b,d,6,6)
        do i = 1, 6
          print '(1p,6e12.4)', (d(i,j),j=1,6)
        enddo
      endif
      do i = 1, 6
        print '(1p,6e12.4)', (a(i,j),j=1,6)
      enddo
      call eigen(a, 6, q1, q2)
      do i = 1, 6
        print '(1p,6e12.4)', (a(i,j),j=1,6)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mapcat(fsec, rb, tb, ra, ta, rd, td)
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      integer i1,i2,i3
      double precision ra,rb,rd,rw,ta,tb,td,ts,tw
c-----------------------------------------------------------------------
c   concatenate two transport maps
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      logical           fsec
      dimension         rb(6,6), tb(36,6), ra(6,6), ta(6,6,6)
      dimension         rd(6,6), td(6,6,6)
 
      dimension         rw(6,6), tw(6,6,6), ts(36,6)
c-----------------------------------------------------------------------

!     write(*,*) 'ra in mapcat: ',ra
!     write(*,*) 'rb in mapcat: ',rb
c transfer matrix
      do  i2 = 1, 6
        do  i1 = 1, 6
        rw(i1,i2) = rb(i1,1)*ra(1,i2) + rb(i1,2)*ra(2,i2)
     +            + rb(i1,3)*ra(3,i2) + rb(i1,4)*ra(4,i2)
     +            + rb(i1,5)*ra(5,i2) + rb(i1,6)*ra(6,i2)
        enddo
      enddo
!     write(*,*) 'rw (1) in mapcat: ',rw
!     write(38,*) 'rw (1) in mapcat: ',rw
 
c second order terms
      if (fsec) then
        do  i3 = 1, 6
          do  i1 = 1, 36
          ts(i1,i3) = tb(i1,1)*ra(1,i3) + tb(i1,2)*ra(2,i3)
     +              + tb(i1,3)*ra(3,i3) + tb(i1,4)*ra(4,i3)
     +              + tb(i1,5)*ra(5,i3) + tb(i1,6)*ra(6,i3)
          enddo
        enddo
 
        do i2 = 1, 6
         do i3 = i2, 6
          do i1 = 1, 6
          tw(i1,i2,i3) =
     +        rb(i1,1)*ta(1,i2,i3) + rb(i1,2)*ta(2,i2,i3)
     +      + rb(i1,3)*ta(3,i2,i3) + rb(i1,4)*ta(4,i2,i3)
     +      + rb(i1,5)*ta(5,i2,i3) + rb(i1,6)*ta(6,i2,i3)
     +      + ts(i1,   i2)*ra(1,i3) + ts(i1+ 6,i2)*ra(2,i3)
     +      + ts(i1+12,i2)*ra(3,i3) + ts(i1+18,i2)*ra(4,i3)
     +      + ts(i1+24,i2)*ra(5,i3) + ts(i1+30,i2)*ra(6,i3)
          tw(i1,i3,i2) = tw(i1,i2,i3)
          enddo
         enddo
        enddo
      endif
 
c copy result
      call ddcopy(rw, rd, 36)
      if (fsec) call ddcopy(tw, td, 216)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mkpits2
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol,twompar
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
c both sides parasitic interactions
      parameter (twompar = 2*mpar)

      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam,ipstart
      common /pitf/ si(mpit)
      save /pitf/
      
      common /selecip/ onepit,ipstart,intip
      save /selecip/ 
      integer intip
      integer i,m
      logical onepit
      integer isinpit,pitnum, nampos
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
c-----------------------------------------------------------------------
c Fills a table with parameters at the pits with suffix i and
c 1 dimension, the pit number.
c The square roots of beta-functions are stored.
c
c The algorithm was copied from E. Keil's "orbit8" program.
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      integer headon
c generate tables at IP's of beta functions, etc.
      frev   = clight / circum
      gamma  = gev / tmass
 
      npit = 0
 
      do 10 i = 1, ninter
         m = headon(name(i,1))
         if (m .gt. 0) then
            if (npit .ge. mpit) then
               print *, 'Error in "train" program - ',
     +              'too many pits in optics table'
               stop
            endif
 
            npit = npit + 1
            ipit(npit) = i
            pitnam(npit) = name(i,1)(3:)
            si(npit) = s(i,1)
            print *, 'active pit: ', m
            iact(m) = 1
            actlist(npit) = m
         endif
   10 continue

      if (npit .eq. 0) then
         print *, 'Error in "train" program - ',
     +        'No IP found in optics table'
         stop
      endif

      print *, 'no. of active pits: ', npit
      

c-----------------------------------------------------------------------
      end
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mkpits
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol,twompar
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
c both sides parasitic interactions
      parameter (twompar = 2*mpar)

      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam,ipstart
      common /pitf/ si(mpit)
      save /pitf/
      
      common /selecip/ onepit,ipstart,intip
      save /selecip/ 
      integer intip
      integer i,m
      logical onepit
      integer isinpit,pitnum, nampos
      
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
c-----------------------------------------------------------------------
c Fills a table with parameters at the pits with suffix i and
c 1 dimension, the pit number.
c The square roots of beta-functions are stored.
c
c The algorithm was copied from E. Keil's "orbit8" program.
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      integer headon
c generate tables at IP's of beta functions, etc.
      frev   = clight / circum
      gamma  = gev / tmass
 
      npit = 0
 
      do 10 i = 1, ninter
         m = headon(name(i,1))
         if (m .gt. 0) then
            if (npit .ge. mpit) then
               print *, 'Error in "train" program - ',
     +              'too many pits in optics table'
               stop
            endif
 
            npit = npit + 1
            ipit(npit) = i
            pitnam(npit) = name(i,1)(3:)
            si(npit) = s(i,1)
            print *, 'active pit: ', m
            iact(m) = 1
         endif
   10 continue

      if (npit .eq. 0) then
         print *, 'Error in "train" program - ',
     +        'No IP found in optics table'
         stop
      endif

      print *, 'no. of active pits: ', npit
      

c-----------------------------------------------------------------------
      end
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mktune(flag)
      implicit none
      integer i,j,k,l,last,m
      double precision aux1,aux2,cosmux,cosmuy,rp1,rp2,sinmu2,
     +sinmux,sinmuy,temp1,temp2,q1f,q2f,q1b,q2b
c-----------------------------------------------------------------------
c analyse tunes
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
cc maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
cc velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
cc pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
cc global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
cc dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
cc   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
c 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
cc   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
cc maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical flag
      dimension aux1(6), aux2(6)
      dimension rp1(6,6), rp2(6,6)
c-----------------------------------------------------------------------
      double precision ma(4,4)

      print *, 'Building second-order maps . . .'
 
c track second-order maps
      if (flag) then
         call trackb(.true., .false., -1)
         last = min(c_tunes, nbunch)
      else
         call track0(.true., .false.)
         last = 1
      endif
 
c loop over bunches
      m = ninter + 1
      do 90 l = 1, last
         if (l .eq. 1 .or. l .eq. last .or. mod(l,100) .eq. 0)
     +   print *, 'Analysing bunch pair ', code(l), ' . . .'
c compute inital first-order dispersion
         do i = 1, 6
           do j = 1, 6
             rp1(j,i) = tr1(j,i,l)
             rp2(j,i) = tr2(j,i,l)
           enddo
         enddo
         
         call eigen(rp1, 6, q1f, q2f)
         call eigen(rp2, 6, q1b, q2b)

         call disp(tr1(1,1,l), tr1(1,6,l), d1(1,l,0))
         call disp(tr2(1,1,l), tr2(1,6,l), d2(1,l,m))
         d1(5,l,0) = 0.0
         d1(6,l,0) = 1.0
         d2(5,l,m) = 0.0
         d2(6,l,m) = 1.0
 
c derivative of transfer matrix w.r.t. delta(p)/p
         do 30 i = 1, 6
            aux1(i) = 0.0
            aux2(i) = 0.0
 
            do 20 k = 1, 6
               temp1 = 0.0
               temp2 = 0.0
 
               do 10 j = 1, 6
                  temp1 = temp1 + tt1(i,j,k,l) * d1(j,l,0)
                  temp2 = temp2 + tt2(i,j,k,l) * d2(j,l,m)
 10            continue
 
               aux1(i) = aux1(i) + temp1 * d1(k,l,0)
               rp1(i,k) = 2.0 * temp1
               aux2(i) = aux2(i) + temp2 * d2(k,l,m)
               rp2(i,k) = 2.0 * temp2
 20         continue
 30      continue
 
c compute inital second-order dispersion
         call disp(tr1(1,1,l), aux1, dd1(1,l,0))
         call disp(tr2(1,1,l), aux2, dd2(1,l,m))
         dd1(5,l,0) = 0.0
         dd1(6,l,0) = 0.0
         dd2(5,l,m) = 0.0
         dd2(6,l,m) = 0.0
 
c horizontal motion

         cosmux = (tr1(1,1,l) + tr1(2,2,l)) / 2.0
         qx1(l) = acos(cosmux) / (2.0 * pi)

         sinmu2 = - tr1(1,2,l) * tr1(2,1,l) -
     +        (tr1(1,1,l) - tr1(2,2,l))**2 / 4.0
         if (sinmu2 .lt. 0.0) sinmu2 = toler
         sinmux = sign(sqrt(sinmu2), tr1(1,2,l))
         qxp1(l) = - (rp1(1,1) + rp1(2,2)) / (4.0 * pi * sinmux)
c         q11(l) = q1f
 
         cosmux = (tr2(1,1,l) + tr2(2,2,l)) / 2.0
         qx2(l) = acos(cosmux) / (2.0 * pi)
         sinmu2 = - tr2(1,2,l) * tr2(2,1,l) -
     +        (tr2(1,1,l) - tr2(2,2,l))**2 / 4.0
         if (sinmu2 .lt. 0.0) sinmu2 = toler
         sinmux = sign(sqrt(sinmu2), tr2(1,2,l))
         qxp2(l) = - (rp2(1,1) + rp2(2,2)) / (4.0 * pi * sinmux)
c         q12(l) = q1b
 
c vertical motion
         cosmuy = (tr1(3,3,l) + tr1(4,4,l)) / 2.0
         qy1(l) = acos(cosmuy) / (2.0 * pi)
         sinmu2 = - tr1(3,4,l) * tr1(4,3,l) -
     +        (tr1(3,3,l) - tr1(4,4,l))**2 / 4.0
         if (sinmu2 .lt. 0.) sinmu2 = toler
         sinmuy = sign(sqrt(sinmu2), tr1(3,4,l))
         qyp1(l) = - (rp1(3,3) + rp1(4,4)) / (4.0 * pi * sinmuy)
c         q21(l) = q2f
 
         cosmuy = (tr2(3,3,l) + tr2(4,4,l)) / 2.0
         qy2(l) = acos(cosmuy) / (2.0 * pi)
         sinmu2 = - tr2(3,4,l) * tr2(4,3,l) -
     +        (tr2(3,3,l) - tr2(4,4,l))**2 / 4.0
         if (sinmu2 .lt. 0.) sinmu2 = toler
         sinmuy = sign(sqrt(sinmu2), tr2(3,4,l))
         qyp2(l) = - (rp2(3,3) + rp2(4,4)) / (4.0 * pi * sinmuy)
c         q22(l) = q2b
 
c track dispersion
         if (flag) then
            call trackb(.true., .true., l)
         else
            call track0(.true., .true.)
         endif
 90   continue
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine maprint(n, t, m)
      implicit none
      integer n, i, j
      character *(*) t
      double precision m(6,6)
      print *, 'number: ', n, ' ', t
      do i = 1, 6
        print '(1p,6d12.4)', (m(i,j), j = 1, 6)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mxbyv(amat, avec, target)
      implicit none
      integer i,j
      double precision amat,avec,target,temp
c-----------------------------------------------------------------------
c multiply matrix times vector
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      dimension         amat(6,6), avec(6), target(6)
      dimension         temp(6)
c-----------------------------------------------------------------------
      do 20 i = 1, 6
        temp(i) = 0.0
        do 10 j = 1, 6
          temp(i) = temp(i) + amat(i,j) * avec(j)
   10   continue
   20 continue
 
      call ddcopy(temp, target, 6)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mxmpy(fact1, fact2, target, m, n)
      implicit none
      integer i,j,k,m,n
      double precision fact1,fact2,target,temp
c-----------------------------------------------------------------------
c multiply two matrices
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
c-----------------------------------------------------------------------
      dimension fact1(m,m), fact2(m,m), target(m,m)
      dimension temp(mvary)
c-----------------------------------------------------------------------
      do 90 k = 1, n
         do 10 i = 1, n
            temp(i) = 0.
 10      continue
 
         do 30 j = 1, n
            do 20 i = 1, n
               temp(i) = temp(i) + fact1(i,j) * fact2(j,k)
 20         continue
 30      continue
 
         do 40 i = 1, n
            target(i,k) = temp(i)
 40      continue
 90   continue
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mxone(target, m, n)
      implicit none
      integer i,j,m,n
      double precision target
c-----------------------------------------------------------------------
c   set matrix to unity
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      dimension target(m,m)
c-----------------------------------------------------------------------
      do 20 j = 1, n
         do 10 i = 1, n
            target(i,j) = 0.0
 10      continue
 
         target(j,j) = 1.0
 20   continue
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mxtrns(r, a, l, n, m)
      implicit none
      integer j,k,l,m,n
      double precision a,r,vec
c-----------------------------------------------------------------------
c tranform block at l by r
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
c-----------------------------------------------------------------------
      dimension r(6,6), a(m,m), vec(mdim)
c-----------------------------------------------------------------------
      do 150 j = 1, n
         do 130 k = 1, mdim
            vec(k) = r(k,1) * a(l+1,j) + r(k,2) * a(l+2,j) +
     +               r(k,3) * a(l+3,j) + r(k,4) * a(l+4,j)
 130     continue
 
         do 140 k = 1, mdim
            a(l+k,j) = vec(k)
 140     continue
 150  continue
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function mylist(kelem,nlist,klist)
c-----------------------------------------------------------------------
c
c---Purpose:    finds number in sorted list (ascending)
c               with binary search.
c
c---Input
c   KELEM           number to be looked up
c   NLIST           length of table
c   KLIST           table
c
c---Output function value:
c                   = 0: number not in table
c                   > 0: position in table
c
c---Author :    HG      date: 17.5.79     last revision: 20.6.84
c
c-----------------------------------------------------------------------
      implicit none
      integer kelem, nlist, klist(*)
      integer ipos, kpos, last, m, n
      ipos=0
      last=0
      n=nlist
      if(n.gt.0)  then
         kpos=0
   10    m=(n+1)/2
         last=kpos+m
         if (kelem.lt.klist(last))  then
            n=m
            last=last-1
            if (n.gt.1) goto 10
         elseif (kelem.gt.klist(last))  then
            kpos=last
            n=n-m
            if (n.gt.0) goto 10
         else
            ipos=last
         endif
      endif
      mylist = ipos
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character * 120 function nxline()
      implicit none
      integer lastnb
      character * 120 text
   10 continue
      read (5, '(a)', end = 100) text
      print *, text(:lastnb(text))
      if (text(1:1) .eq. '#')  goto 10
      nxline = text
      return
  100 print *, 'EOF on input before all data read'
      stop
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine orbw
      implicit none
      integer i,j
      double precision vdist, hdist
      double precision sigx, sigy
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      do j = 0, ninter
        if (j .eq. 0) then
          i = ninter
        else
          i = j
        endif
        hdist = abs(z2(1,1,j) - z1(1,1,j)) / sigx(1,i,2)
        vdist = abs(z2(3,1,j) - z1(3,1,j)) / sigy(1,i,2)
        write(orbout, '(i3, 6f12.6)') j, hdist, vdist,
     +  z1(1,1,j)/sigx(1,i,1), z1(3,1,j)/sigy(1,i,1),
     +  z2(1,1,j)/sigx(1,i,2), z2(3,1,j)/sigy(1,i,2)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine orbit02
      implicit none
      integer i,irank,itmax,itra,j,m
      double precision a,err
c-----------------------------------------------------------------------
c find closed orbits for unperturbed bunches,
c fill positions at collision points
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      parameter (itmax = 5)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      dimension a(4,5)
c-----------------------------------------------------------------------
      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps),eqvoptmap
      logical extraelem
      character*12 mname(mmaps)
      save /extel/
      integer ipit,ncoll
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
c-----------------------------------------------------------------------

      m = nmaps + 1

      do 10 i = 1, 5
         z1(i,1,0) = zero
         z2(i,1,m) = zero
 10   continue
      z1(6,1,0) = delta
      z2(6,1,m) = delta
c initial guess
      print *, 'Solving c.o. without beam-beam interactions . . .'
 
c iteration for closed orbits
      do 200 itra = 1, itmax  

c track orbits and transfer matrices
         call track02(.false., .false.)
         err = 0.0     
 
c solve for forward beam
         do 120 i = 1, 4
            do 110 j = 1, 4
               a(i,j) = tr1(i,j,1)
            write(51,*) i,j,tr1(i,j,1)
 110        continue

            a(i,i) = a(i,i) - 1.0
            a(i,5) = z1(i,1,m) - z1(i,1,0)
            err = max(abs(a(i,5)), err)
 120     continue
         
      if (debug .gt. 0)  print *, 'call solver in orbit0 forward beam'
         call solver(a, 4, 1, irank)
         if (irank .lt. 4) stop
 
         do 140 i = 1, 4
            z1(i,1,0) = z1(i,1,0) - a(i,5)
 140     continue
 
c solve for backward beam
         do 170 i = 1, 4
            do 160 j = 1, 4
               a(i,j) = tr2(i,j,1)
            write(52,*) i,j,tr2(i,j,1)
 160        continue
 
            a(i,i) = a(i,i) - 1.0
            a(i,5) = z2(i,1,0) - z2(i,1,m)
            err = max(abs(a(i,5)), err)
 170     continue
 
      if (debug .gt. 0)  print *, 'call solver in orbit0 backward beam'
         call solver(a, 4, 1, irank)
         if (irank .lt. 4) stop
 
         do 190 i = 1, 4
            z2(i,1,m) = z2(i,1,m) - a(i,5)
 190     continue
 
c message
         print *, 'Iteration ', itra, ', error = ', err
 
c convergence test
         if (err .lt. toler) return
 200  continue
     
c---- no convergence.
      print *, 'Closed orbits did not converge in orbit02', itmax,
     +     ' iterations'
      stop
c-----------------------------------------------------------------------
 910  format(' '/' Looking for unperturbed closed orbits:')
 920  format(16x,'forward beam',22x,'backward beam'/' iter.',
     +     '  x       px      y       py    ',2x,
     +     '  x       px      y       py      error'/
     +     6x,'  [mm]    [mrad]  [mm]    [mrad]',
     +     2x,'  [mm]    [mrad]  [mm]    [mrad]')
 930  format(1x,i4,1x,3p,4f8.5,2x,4f8.5,1p,e14.6)
c-----------------------------------------------------------------------
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine orbit0
      implicit none
      integer i,irank,itmax,itra,j,m
      double precision a,err
c-----------------------------------------------------------------------
c find closed orbits for unperturbed bunches,
c fill positions at collision points
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      parameter (itmax = 5)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      dimension a(4,5)
c-----------------------------------------------------------------------
c initialize start points for one bunch per direction
      m = ninter + 1
      do 10 i = 1, 5
         z1(i,1,0) = zero
         z2(i,1,m) = zero
 10   continue
      z1(6,1,0) = delta
      z2(6,1,m) = delta
 
c initial guess
      print *, 'Solving c.o. without beam-beam interactions . . .'
 
c iteration for closed orbits
      do 200 itra = 1, itmax
 
c track orbits and transfer matrices
         call track0(.false., .false.)
         err = 0.0
 
c solve for forward beam
         do 120 i = 1, 4
            do 110 j = 1, 4
               a(i,j) = tr1(i,j,1)
            write(51,*) i,j,tr1(i,j,1)
 110        continue

            a(i,i) = a(i,i) - 1.0
            a(i,5) = z1(i,1,m) - z1(i,1,0)
            err = max(abs(a(i,5)), err)
 120     continue
 
      if (debug .gt. 0)  print *, 'call solver in orbit0 forward beam'
         call solver(a, 4, 1, irank)
         if (irank .lt. 4) stop
 
         do 140 i = 1, 4
            z1(i,1,0) = z1(i,1,0) - a(i,5)
 140     continue
 
c solve for backward beam
         do 170 i = 1, 4
            do 160 j = 1, 4
               a(i,j) = tr2(i,j,1)
            write(52,*) i,j,tr2(i,j,1)
 160        continue
 
            a(i,i) = a(i,i) - 1.0
            a(i,5) = z2(i,1,0) - z2(i,1,m)
            err = max(abs(a(i,5)), err)
 170     continue
 
      if (debug .gt. 0)  print *, 'call solver in orbit0 backward beam'
         call solver(a, 4, 1, irank)
         if (irank .lt. 4) stop
 
         do 190 i = 1, 4
            z2(i,1,m) = z2(i,1,m) - a(i,5)
 190     continue
 
c message
         print *, 'Iteration ', itra, ', error = ', err
 
c convergence test
         if (err .lt. toler) return
                    
        
 200  continue
 
c---- no convergence.
      print *, 'Closed orbits did not converge in ', itmax,
     +     ' iterations'
      stop
c-----------------------------------------------------------------------
 910  format(' '/' Looking for unperturbed closed orbits:')
 920  format(16x,'forward beam',22x,'backward beam'/' iter.',
     +     '  x       px      y       py    ',2x,
     +     '  x       px      y       py      error'/
     +     6x,'  [mm]    [mrad]  [mm]    [mrad]',
     +     2x,'  [mm]    [mrad]  [mm]    [mrad]')
 930  format(1x,i4,1x,3p,4f8.5,2x,4f8.5,1p,e14.6)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine orbitb2
      implicit none
      integer i,iout,irank,itmax,itra,j,k,m
      double precision a,dz1,dz2,err
c-----------------------------------------------------------------------
c find closed orbits for all bunches
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.50d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      parameter (itmax = 500)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
c-----------------------------------------------------------------------
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      integer hd
      dimension a(4,5)
c-----------------------------------------------------------------------
      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/

      integer ipit,ncoll
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/

c initialize start points for all bunches
      m = nmaps + 1
      print *, 'Solving c.o. with beam-beam interactions . . .'
 
      do 20 j = 2, nbunch
         do 10 i = 1, 6
            z1(i,j,0) = z1(i,1,0)
            z2(i,j,m) = z2(i,1,m)
 10      continue
  20  continue

c iteration for closed orbits
      do 400 iout = 1, itmax
         do 300 itra = 1, itmax
c track orbits and transfer matrices
            call trackb2(f_second, .false., -1)
            err = 0.0
c solve for forward bunch
            do 200 k = 1, nbunch
               do 120 i = 1, 4
                  do 110 j = 1, 4
                     a(i,j) = tr1(i,j,k)
 110              continue
 
                  a(i,i) = a(i,i) - 1.0
                  a(i,5) = z1(i,k,m) - z1(i,k,0)
                  err = max(abs(a(i,5)), err)
 120           continue
      if (debug .gt. 1)  print *, 'call solver in orbitb'
               call solver(a, 4, 1, irank)
               if (irank .lt. 4) stop
 
               do 130 i = 1, 4
                  z1(i,k,0) = z1(i,k,0) - a(i,5)
 130           continue
 
c solve for backward bunch
               do 170 i = 1, 4
                  do 160 j = 1, 4
                     a(i,j) = tr2(i,j,k)
 160              continue
 
                  a(i,i) = a(i,i) - 1.0
                  a(i,5) = z2(i,k,0) - z2(i,k,m)
                  err = max(abs(a(i,5)), err)
 170           continue
 
      if (debug .gt. 1)  print *, 'call solver in orbitb'
               call solver(a, 4, 1, irank)
               if (irank .lt. 4) stop
 
               do 180 i = 1, 4
                  z2(i,k,m) = z2(i,k,m) - a(i,5)
 180           continue
 200        continue

            print *, '   Inner iteration ', itra, ', error = ', err
c convergence test
            if (err .lt. toler) goto 310
 300     continue
 
c---- no convergence.
         print *, 'Closed orbits did not converge in ', itmax,
     +        ' iterations'
         stop
 
 310     err = 0.0
         
         do 340 i = 1, ninter
            do 330 j = 1, nbunch
               
               do 320 k = 1, 6
                  dz1 = half*(z1b(k,j,i)+z1a(k,j,i))-z1(k,j,i)
                  z1(k,j,i) = z1(k,j,i) + dz1
                  dz2 = half*(z2b(k,j,i)+z2a(k,j,i))-z2(k,j,i)
                  z2(k,j,i) = z2(k,j,i) + dz2
                  err = max(err, abs(dz1), abs(dz2))
                  
 320           continue

 330        continue
 340     continue
 
c in case of success, compute second-order map
         print *, 'Outer iteration ', iout, ', error = ', err

         if (err .lt. toler) return
            
         
 400  continue
c-----------------------------------------------------------------------
 920  format(' '/
     +     ' Initial guess for orbits with beam-beam interactions:')
 930  format(17x,'forward beam',24x,'backward beam'/
     +     ' bunch ',2('  x       px      y       py        ')/
     +     7x,2('  [mm]    [mrad]  [mm]    [mrad]    '))
 940  format(' ',a5,' ',3p,4f8.6,4x,4f8.6)
 950  format(' '/' Iteration = ',i2,'/',i2,', error = ',1p,e16.6)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine orbitb
      implicit none
      integer i,iout,irank,itmax,itra,j,k,m
      double precision a,dz1(6),dz2(6),err
c-----------------------------------------------------------------------
c find closed orbits for all bunches
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      parameter (itmax = 500)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      
      integer ipit,ncoll,l
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      
      character*120 filename,chst
      common /filesch/ filename,chst
      save / filesch/

c-----------------------------------------------------------------------
      integer hd,nstart
      dimension a(4,5) 
c-----------------------------------------------------------------------
c initialize and print start points for all bunches
      m = ninter + 1
      print *, 'Solving c.o. with beam-beam interactions . . .'
 
      do 20 j = 2, nbunch
         do 10 i = 1, 6
            z1(i,j,0) = z1(i,1,0)
            z2(i,j,m) = z2(i,1,m)
 10      continue
  20   continue

      do l = 1, npit
        open(73 + l,file='orb.'// name(ipit(l),1)(3:), 
     +  status = 'UNKNOWN')
      enddo 
c iteration for closed orbits
      do 400 iout = 1, itmax
         do 300 itra = 1, itmax
c track orbits and transfer matrices
            call trackb(f_second, .false., -1)
            err = 0.0
 
c solve for forward bunch
            do 200 k = 1, nbunch
               do 120 i = 1, 4
                  do 110 j = 1, 4
                     a(i,j) = tr1(i,j,k)
 110              continue
 
                  a(i,i) = a(i,i) - 1.0
                  a(i,5) = z1(i,k,m) - z1(i,k,0)
                  err = max(abs(a(i,5)), err)
 120           continue
      if (debug .gt. 1)  print *, 'call solver in orbitb'
               call solver(a, 4, 1, irank)
               if (irank .lt. 4) stop
 
               do 130 i = 1, 4
                  z1(i,k,0) = z1(i,k,0) - a(i,5)
 130           continue
 
c solve for backward bunch
               do 170 i = 1, 4
                  do 160 j = 1, 4
                     a(i,j) = tr2(i,j,k)
 160              continue
 
                  a(i,i) = a(i,i) - 1.0
                  a(i,5) = z2(i,k,0) - z2(i,k,m)
                  err = max(abs(a(i,5)), err)
 170           continue
 
      if (debug .gt. 1)  print *, 'call solver in orbitb'
               call solver(a, 4, 1, irank)
               if (irank .lt. 4) stop
 
               do 180 i = 1, 4
                  z2(i,k,m) = z2(i,k,m) - a(i,5)
 180           continue
 200        continue
 

            hd = nstart()
            do l = 1, npit
              write(73 + l,*) iout,itra,z1(1,hd,ipit(l)),
     + z2(1,hd,ipit(l)),z1(3,hd,ipit(l)),z2(3,hd,ipit(l)),
     + z1(2,hd,ipit(l)),z2(2,hd,ipit(l)),z1(4,hd,ipit(l)),
     + z2(4,hd,ipit(l))
                         
            enddo 
            
c convergence test
            if (err .lt. toler) goto 310
 300     continue
 
c---- no convergence.
         print *, 'Closed orbits did not converge in ', itmax,
     +        ' iterations'
         stop
 
 310     err = 0.0
         
         do 340 i = 1, ninter
            do 330 j = 1, nbunch
               
               do 320 k = 1, 6
                  dz1(k) = half*(z1b(k,j,i)+z1a(k,j,i))-z1(k,j,i)
                  z1(k,j,i) = z1(k,j,i) + dz1(k)
                  dz2(k) = half*(z2b(k,j,i)+z2a(k,j,i))-z2(k,j,i)
                  z2(k,j,i) = z2(k,j,i) + dz2(k)
                  err = max(err, abs(dz1(k)), abs(dz2(k)))
                  
 320           continue
 330        continue
 340     continue
 
c in case of success, compute second-order map
         print *, 'Outer iteration ', iout, ', error = ', err
         if (err .lt. toler) then
           do l = 1, npit
             close(73+l)
           enddo 


           return
         endif 
 400  continue
c-----------------------------------------------------------------------
 920  format(' '/
     +     ' Initial guess for orbits with beam-beam interactions:')
 930  format(17x,'forward beam',24x,'backward beam'/
     +     ' bunch ',2('  x       px      y       py        ')/
     +     7x,2('  [mm]    [mrad]  [mm]    [mrad]    '))
 940  format(' ',a5,' ',3p,4f8.6,4x,4f8.6)
 950  format(' '/' Iteration = ',i2,'/',i2,', error = ',1p,e16.6)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function outbunch(bunch)
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      integer bunch, i
      outbunch = 0
      do i = 1, outbcnt
        if (bunch .eq. outblist(i)) then
          outbunch = i
          return
        endif
      enddo
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function parasit(name, nbeam)
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      integer nbeam
      character *(*) name
      integer i, lastnb, nampos
      if (name(1:4) .eq. 'MKIP' .and. lastnb(name) .gt. 5) then
        if (nbeam .eq. 1)  then
          i = nampos(name, para_names, n_parasit)
          if (i .gt. 0)  then
            parasit = i
          else
            n_parasit = n_parasit + 1
            para_names(n_parasit) = name
            parasit = n_parasit
          endif
        else
          parasit = 1
        endif
      else
        parasit = 0
      endif
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pbunch
c-----------------------------------------------------------------------
c--- prints global tracking results
c-----------------------------------------------------------------------
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      print *, 'average luminosity: ', lumiav / lumicnt
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      logical function pit_number(n)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer n, i
      pit_number = .false.
      do i = 1, npit
        if (n .eq. ipit(i)) pit_number = .true.
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prcdmp
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      integer i, j, cnt, regp, reg_list(10)
      if (w_coll)  then
        open(11,file='coll.count', status = 'UNKNOWN')
        write(11,'(a)')
     +  '#  slot, f + b # of bunches, bunch numbers, bunch ref.'
        do i = 0, mbuck-1
          if ((collsk(1,i).ne.0).or.(collsk(2,i).ne.0))  then
            write(11, '(3i10)')
     +      i, colcnt_f(i), colcnt_b(i)
            write(11, '(20i5)') (list_f(j,i), j = 1, colcnt_f(i))
            write(11, '(20i5)') (part_f(j,i), j = 1, colcnt_f(i))
          endif
        enddo
        close(11)
      endif
      if (w_frequ)  then
        open(12, file='freq_f.count', status = 'UNKNOWN')
        write(12,'(a)')
     +  '#  bforward: class, # coll., # bunches, (bunch numbers)'
        do i = 1, ordl_f
          write(12, '(3i10)') i, cordl_f(i), nordl_f(i)
          if (w_detail)
     +    write(12, '(20i5)') (lordl_f(j,i), j = 1, nordl_f(i))
        enddo
        close(12)
        open(12, file='freq_b.count', status = 'UNKNOWN')
        write(12,'(a)')
     +  '#  backward: class, # coll., # bunches, (bunch numbers)'
        do i = 1, ordl_b
          write(12, '(3i10)') i, cordl_b(i), nordl_b(i)
          if (w_detail)
     +    write(12, '(20i5)') (lordl_b(j,i), j = 1, nordl_b(i))
        enddo
        close(12)
      endif
      if (w_equ)  then
        open(13, file='equ_f.count', status = 'UNKNOWN')
        write(13,'(a)')
     +  '#  forward: equ. class, # coll., # bunches, bunch numbers'
        do i = 1, equl_f
          write(13, '(3i10)') i, cequl_f(i), nequl_f(i)
          write(13, '(20i5)') (lequl_f(j,i), j = 1, nequl_f(i))
        enddo
        close(13)
        open(13, file='equ_b.count', status = 'UNKNOWN')
        write(13,'(a)')
     +  '#  backward: equ. class, # coll., # bunches, bunch numbers'
        do i = 1, equl_b
          write(13, '(3i10)') i, cequl_b(i), nequl_b(i)
          write(13, '(20i5)') (lequl_b(j,i), j = 1, nequl_b(i))
        enddo
        close(13)
      endif
      if (w_set)  then
        open(13, file='set_f.list', status = 'UNKNOWN')
        write(13,'(a)')
     +  '#  bucket numbers in forward set'
        write(13, '(20i5)') (set_f(j), j = 1, nset_f)
        close(13)
        open(13, file='set_b.list', status = 'UNKNOWN')
        write(13,'(a)')
     +  '#  bucket numbers in backward set'
        write(13, '(20i5)') (set_b(j), j = 1, nset_b)
        close(13)
      endif
      if (w_alt)  then
        open(13, file='alt.list', status = 'UNKNOWN')
        write(13,'(a)')
     +  '#  coll. point i, bunch j, ibnch1(j,i), ibnch2(j,i)'
        do i = 1, mcol
          do j = 1, nbunch
            if (ibnch1(j,i) .ne. 0 .or. ibnch2(j,i) .ne. 0)  then
              write(13, '(4i10)') i, j, ibnch1(j,i), ibnch2(j,i)
            endif
          enddo
        enddo
        close(13)
        open(13, file='reg.list', status = 'UNKNOWN')
        write(13,'(a)')
     +  '#  list of regular bunches'
        regp = 0
        do j = 1, nbunch
          cnt = 0
          do i = 1, mcol
            if (ibnch1(j,i) .ne. 0) cnt = cnt + 1
          enddo
          if (cnt .eq. ninter) then
            regp = regp + 1
            reg_list(regp) = j
            if (regp .eq. 10)  then
              write(13, '(10i6)') reg_list
              regp = 0
            endif
          endif
        enddo
        if (regp .ne. 0) write(13, '(10i6)') (reg_list(j), j = 1, regp)
        close(13)
      endif
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prcoll
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
c   a           injection scheme (supposed symmetric to IP5)
c   mask        half-slot mask for a
c   maski       coll. point number at half-slot
c   rb          half-slot bunch scheme

 
      integer a(2,0:mbuck-1), mask(0:mdslt-1), maski(0:mdslt-1)
      integer cp, i, j, k, n, k1, k2
      integer low_f, high_f, low_b, high_b
      integer mylist
      integer temppar
c      equivalence (a(1), collsk(1))
c      equivalence (a(2),collsk(2))
      k = 0
      do i = 0, mbuck-1
        if (collsk(1,i) .ne. 0)  k = k + 1
      enddo
      print *, 'total bunches forward beam: ', k

      k = 0
      do i = 0, mbuck-1
        if (collsk(2,i) .ne. 0)  k = k + 1
      enddo
      print *, 'total bunches backward beam: ', k


      call getmask2(mdslt, npar, circum, iact, ippos, mask,
     +maski,colpnt)
      print *, 'total number of collision points: ', colpnt
      call collide(mbuck, mdslt, mcol, b2_off, npar, collsk, mask, 
     +maski,colcnt_f, list_f, part_f, colcnt_b, list_b, part_b)
      k = 0
      do i = 0, mbuck-1
        if (colcnt_f(i) .ne. colcnt_b(i)) then
          k = k + 1
        endif
      enddo
      print *, 'no. of bunches with different counts f/b: ', k
      call equ_class(mbuck, mcol, collsk, colcnt_f, list_f, low_f, 
     +high_f,equl_f, cequl_f, nequl_f, lequl_f,
     +ordl_f, cordl_f, nordl_f, lordl_f)
      print *, 'low_f, high_f: ', low_f, high_f
      call equ_class(mbuck, mcol, collsk, colcnt_b, list_b, low_b, 
     +high_b,equl_b, cequl_b, nequl_b, lequl_b,
     +ordl_b, cordl_b, nordl_b, lordl_b)
      print *, 'low_b, high_b: ', low_b, high_b
      ntotal_f = 0
      do i = 1, ordl_f
        ntotal_f = ntotal_f + nordl_f(i)
      enddo
      ctotal_f = 0
      do i = 1, equl_f
        ctotal_f = ctotal_f + nequl_f(i)
      enddo
      ntotal_b = 0
      do i = 1, ordl_b
        ntotal_b = ntotal_b + nordl_b(i)
      enddo
      ctotal_b = 0
      do i = 1, equl_b
        ctotal_b = ctotal_b + nequl_b(i)
      enddo
      print *, 'forward  total, classes: ', ntotal_f, ctotal_f
      print *, 'backward total, classes: ', ntotal_b, ctotal_b
      if (f_coll)  then
        print *, '*===* full collision flag: all bunches taken *===*'
        do k = 0, mbuck-1
          hitlist_f(k) = collsk(1,k)
          hitlist_b(k) = collsk(2,k)
        enddo
      else
        call get_hits(mbuck, mcol, part_f, part_b,
     +  colcnt_f, colcnt_b,
     +  hitlist_f, hitlist_b, equl_f, nequl_f, lequl_f)
      endif
      nset_f = 0
      nset_b = 0
      do k = 0, mbuck-1
        if (hitlist_f(k) .ne. 0) then
          nset_f = nset_f + 1
          set_f(nset_f) = k
        endif
        if (hitlist_b(k) .ne. 0) then
          nset_b = nset_b + 1
          set_b(nset_b) = k
        endif
      enddo
      print *, 'bunch subset counts 1 + 2: ', nset_f, nset_b
      if (nset_f .gt. mbunch)  then
        print *, 'warning: no. of bunches ', nset_f,
     +  ' cut at ', mbunch
      endif
c--- make schedule: ibnch1(i,j) is the bunch in beam 2 colliding
c    with bunch i from beam 1 at collision point j
c    ibnch2 is the same with 1 and 2 inverted
      do i = 1, mbunch
        do j = 1, mcol
          ibnch1(i,j) = 0
          ibnch2(i,j) = 0
        enddo
      enddo
      nbunch = min(nset_f, mbunch)
      if (nbunch .ne. min(nset_b, mbunch))  then
        print *, '+++++++ warning: bunch numbers differ: ',
     +  nbunch, '(1),', min(nset_b, mbunch), '(2)'
      endif
      do i = 1, colpnt
        cp = i - 1
        do j = 1, nbunch
          k = set_f(j)
          n = mylist(cp, colcnt_f(k), list_f(1,k))
          if (n .ne. 0)  then
            ibnch1(j,i) = mylist(part_f(n,k), nset_b, set_b)
          endif
        enddo
        do j = 1, min(nset_b, mbunch)
          k = set_b(j)
          n = mylist(cp, colcnt_b(k), list_b(1,k))
          if (n .ne. 0)  then
            ibnch2(j,i) = mylist(part_b(n,k), nset_f, set_f)
          endif
        enddo
      enddo
      do i = 0, mdslt-1
        
        k = i - npar(5)
        
        if (k .lt. 0) k = k + mdslt
        if (i .eq. 0 .or. maski(i) .gt. 0) then
          maskmi(k) = maski(i) + 1
        endif
      enddo
      do i = 0, mdslt-1
        if (maskmi(i) .gt. 0) maskm(maskmi(i)) = i
      enddo
      maskm(ninter+1) = maskm(1)
      do i = 1, ninter
        k1 = maskm(i)
        k2 = maskm(i+1)
        if (k2 .lt. k1)  k2 = k2 + mdslt
        do j = k1, k2
          k = j
          if (j .gt. mdslt) k = j - mdslt
          if (maskmi(k) .gt. 0)  then
            maskmp(k) = maskmi(k)
            maskmn(k) = maskmi(k)
          elseif (i .eq. ninter)  then
            maskmp(k) = i
            maskmn(k) = 1
          else
            maskmp(k) = i
            maskmn(k) = i+1
          endif
        enddo
      enddo
      call prsup
      print *,'after prsup'
      call prcdmp
      print *,'after prcdmp'
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsup
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      integer add, lmm, suppac, nominal,
     +fpip(mpit), fmult(mpit), flist(mbuck), collcnt(mbuck)
      parameter (lmm = 10)
      integer i, j, k, cd, cc, lcc(lmm)
 
c--- forward
      suppac = 0
      do i = 1, mbuck
        flist(i) = 0
        collcnt(i) = 0
      enddo
      nominal = 0
      do i = 1, nbunch
        do j = 1, mcol
          if (ibnch1(i,j) .ne. 0) collcnt(i) = collcnt(i) + 1
        enddo
      if (collcnt(i) .eq. ninter) nominal = nominal + 1
      enddo
      do i = 1, npit
      add = 2 ** (i-1)
      fmult(i) = 0
      fpip(i) = 0
      k = ipit(i)
      if (w_detail) print *, 'forward super-pacman bunches at pit: ',
     +pitnam(i)
        cc = 0
        
        do j = 1, nbunch
          
          if (ibnch1(j,k) .eq. 0) then
            flist(j) = flist(j) + add
            if (cc .eq. lmm) then
              cc = 0
              if (w_detail) print '(10i6)', lcc
            endif
            cc = cc + 1
            lcc(cc) = set_f(j)
          endif
        enddo
        if (w_detail .and. cc .gt. 0) print '(10i6)', (lcc(j),j=1,cc)
      enddo
      do j = 1, nbunch
        cd = flist(j)
        if (cd .gt. 0) suppac = suppac + 1
        cc = 0
        do i = 1, npit
          if (mod(cd,2) .ne. 0) then
            fpip(i) = fpip(i) + 1
            cc = cc + 1
          endif
          cd = cd / 2
        enddo
        if (cc .gt. 0)  fmult(cc) = fmult(cc) + 1
      enddo
      print *, 'forward nominal bunches:              ', nominal
      print *, 'forward super-pacman total:           ', suppac
      print *, 'super-pacman counts per active pit:   ',
     +(fpip(i), i = 1, npit)
      print *, 'super-pacman multiplicity counts (1):     ',
     +(fmult(i), i = 1, npit)
c--- backward
      suppac = 0
      do i = 1, mbuck
        flist(i) = 0
        collcnt(i) = 0
      enddo
      nominal = 0
      
      do i = 1, nbunch
        do j = 1, mcol
          if (ibnch2(i,j) .ne. 0) collcnt(i) = collcnt(i) + 1
        enddo
      if (collcnt(i) .eq. ninter) nominal = nominal + 1
      enddo
      do i = 1, npit
      add = 2 ** (i-1)
      fmult(i) = 0
      fpip(i) = 0
      k = ipit(i)
      if (w_detail) print *, 'backward super-pacman bunches at pit:',
     +pitnam(i)
        cc = 0
        do j = 1, nbunch
          if (ibnch2(j,k) .eq. 0) then
            flist(j) = flist(j) + add
            if (cc .eq. lmm) then
              cc = 0
              if (w_detail) print '(10i6)', lcc
            endif
            cc = cc + 1
            lcc(cc) = set_b(j)
          endif
        enddo
        if (w_detail .and. cc .gt. 0) print '(10i6)', (lcc(j),j=1,cc)
      enddo
      do j = 1, nbunch
        cd = flist(j)
        if (cd .gt. 0) suppac = suppac + 1
        cc = 0
        do i = 1, npit
          if (mod(cd,2) .ne. 0) then
            fpip(i) = fpip(i) + 1
            cc = cc + 1
          endif
          cd = cd / 2
        enddo
        if (cc .gt. 0)  fmult(cc) = fmult(cc) + 1
      enddo
      print *, 'backward nominal bunches:             ', nominal
      print *, 'backward super-pacman total:          ', suppac
      print *, 'super-pacman counts per active pit:   ',
     +(fpip(i), i = 1, npit)
      print *, 'super-pacman multiplicity counts (2):     ',
     +(fmult(i), i = 1, npit)
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine print2(flag, tunes)
      implicit none
      integer m
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical flag, tunes,exist
      integer i


      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/

      print *, 'Printing . . .'
 
      if (flag) then
         write (mulist, 910)
      else
         write (mulist, 920)
      endif
 
c table 1. tunes and chromaticities
      if (tunes) call prtune(flag)

c perturbed beams
      if (flag) then
        if (tunes)  then
c table 2. offset, slope, and dispersion per pit
         open(108,file = 'fort.orbitshift')
         call prsclb2('vdisp','Vertical dispersion in mm',
     +        1, d1, 3, ten3p)
         call prsclb2('vdisp','Vertical dispersion in mm',
     +        2, d2, 3, ten3p)
 
         call prsclb2('vsecdisp',
     +   'Vertical second-order dispersion in m', 1, dd1, 3, one)
         call prsclb2('vsecdisp',
     +   'Vertical second-order dispersion in m', 2, dd2, 3, one)
 
         call prsclb2('hdisp', 'Horizontal dispersion in mm',
     +        1, d1, 1, ten3p)
         call prsclb2('hdisp', 'Horizontal dispersion in mm',
     +        2, d2, 1, ten3p)
 
         call prsclb2('hsecdisp',
     +   'Horizontal second-order dispersion in m', 1, dd1, 1, one)
         call prsclb2('hsecdisp',
     +   'Horizontal second-order dispersion in m', 2, dd2, 1, one)
         
         call pravdisp('fort.av_disp',ten3p) 
c         call prshft
c         call prsprd
          call prlumi
          call pravlumi
        else
c table 2. offset, slope, and dispersion per pit
         inquire(file = 'fort.orbitshift',exist=exist)
         if (exist) then
            open(108,file = 'fort.orbitshift',status="old", 
     +      position="append", action="write")
         else
            open(108,file = 'fort.orbitshift',status="new", 
     +      action="write")
         endif

         call prsclb2('voff', 'Vertical offset in microns',
     +        1, z1, 3, ten6p)
         call prsclb2('voff', 'Vertical offset in microns',
     +        2, z2, 3, ten6p)
 
         call prsclb2('vslope', 'Vertical slope in mrad',
     +        1, z1, 4, ten3p)
         call prsclb2('vslope', 'Vertical slope in mrad',
     +        2, z2, 4, ten3p)
 
         call prsclb2('hoff', 'Horizontal offset in microns',
     +        1, z1, 1, ten6p)
         call prsclb2('hoff', 'Horizontal offset in microns',
     +        2, z2, 1, ten6p)
 
         call prsclb2('hslope', 'Horizontal slope in mrad',
     +        1, z1, 2, ten3p)
         call prsclb2('hslope', 'Horizontal slope in mrad',
     +        2, z2, 2, ten3p)
         close(108)

c table 2. offset, slope, and dispersion per extra element in the maps
         if (extraelem) then
            call prextra('voff', 'Vertical offset in microns',
     +      1,z1,3,ten6p)
            call prextra('voff', 'Vertical offset in microns',
     +      2,z2,3,ten6p)

            call prextra('vslope', 'Vertical slope in mrad',
     +      1,z1,4,ten3p)
            call prextra('vslope', 'Vertical slope in mrad',
     +      2,z2,4,ten3p)
              
            call prextra('hoff', 'Horizontal offset in microns',
     +        1, z1, 1, ten6p)
            call prextra('hoff', 'Horizontal offset in microns',
     +        2, z2, 1, ten6p)
 
            call prextra('hslope', 'Horizontal slope in mrad',
     +        1, z1, 2, ten3p)
            call prextra('hslope', 'Horizontal slope in mrad',
     +        2, z2, 2, ten3p)
         endif 

c table 2. offset, slope, and dispersion per observation point
         call prsobs2('voff', 'Vertical offset in microns',
     +        1, z1, 3, ten6p)
         call prsobs2('voff', 'Vertical offset in microns',
     +        2, z2, 3, ten6p)
 
         call prsobs2('vslope', 'Vertical slope in mrad',
     +        1, z1, 4, ten3p)
         call prsobs2('vslope', 'Vertical slope in mrad',
     +        2, z2, 4, ten3p)
 
         call prsobs2('hoff', 'Horizontal offset in microns',
     +        1, z1, 1, ten6p)
         call prsobs2('hoff', 'Horizontal offset in microns',
     +        2, z2, 1, ten6p)
 
         call prsobs2('hslope', 'Horizontal slope in mrad',
     +        1, z1, 2, ten3p)
         call prsobs2('hslope', 'Horizontal slope in mrad',
     +        2, z2, 2, ten3p)

         
c table 3. separations, crossing angles, etc. per crossing point
         do i = 1, npit
           call prsep2(msep, i, 3)
           call prsep2(msep, i, 1)
         enddo
	endif
c unperturbed beams
      else
 
c table 2. offsets, slopes and dispersions
         call prscl02('Vertical offset in microns',
     +        z1, z2, 3, ten6p)
         call prscl02('Vertical slope in mrad',
     +        z1, z2, 4, ten3p)
         call prscl02('Vertical dispersion in mm',
     +        d1, d2, 3, ten3p)
         call prscl02('Vertical second-order dispersion in m',
     +        dd1, dd2, 3, one)
         call prscl02('Horizontal offset in microns',
     +        z1, z2, 1, ten6p)
         call prscl02('Horizontal slope in mrad',
     +        z1, z2, 2, ten3p)
         call prscl02('Horizontal dispersion in mm',
     +        d1, d2, 1, ten3p)
         call prscl02('Horizontal second-order dispersion in m',
     +        dd1, dd2, 1, one)
      endif
c-----------------------------------------------------------------------
 910  format(' '/' '/'Results including beam-beam collisions:')
 920  format(' '/' '/'Results excluding beam-beam collisions:')
      end
c-----------------------------------------------------------------------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine print(flag, tunes)
      implicit none
      integer m
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical flag, tunes,exist
      integer i
c-----------------------------------------------------------------------

      print *, 'Printing . . .'
      m = ninter + 1
 
      if (flag) then
         write (mulist, 910)
      else
         write (mulist, 920)
      endif
 
c table 1. tunes and chromaticities
      
      if (tunes) call prtune(flag)
      
c perturbed beams
      if (flag) then
        if (tunes)  then
c table 2. offset, slope, and dispersion per pit
         open(108,file = 'fort.orbitshift')
         call prsclb('vdisp','Vertical dispersion in mm',
     +        1, d1, 3, ten3p)
         call prsclb('vdisp','Vertical dispersion in mm',
     +        2, d2, 3, ten3p)
 
         call prsclb('vsecdisp',
     +   'Vertical second-order dispersion in m', 1, dd1, 3, one)
         call prsclb('vsecdisp',
     +   'Vertical second-order dispersion in m', 2, dd2, 3, one)
 
         call prsclb('hdisp', 'Horizontal dispersion in mm',
     +        1, d1, 1, ten3p)
         call prsclb('hdisp', 'Horizontal dispersion in mm',
     +        2, d2, 1, ten3p)
 
         call prsclb('hsecdisp',
     +   'Horizontal second-order dispersion in m', 1, dd1, 1, one)
         call prsclb('hsecdisp',
     +   'Horizontal second-order dispersion in m', 2, dd2, 1, one)
 
c         call prshft
c         call prsprd
          call prlumi
          call pravlumi
        else
c table 2. offset, slope, and dispersion per pit
         inquire(file = 'fort.orbitshift',exist=exist)
         if (exist) then
            open(108,file = 'fort.orbitshift',status="old", 
     +      position="append", action="write")
         else
            open(108,file = 'fort.orbitshift',status="new", 
     +      action="write")
         endif

         call prsclb('voff', 'Vertical offset in microns',
     +        1, z1, 3, ten6p)
         call prsclb('voff', 'Vertical offset in microns',
     +        2, z2, 3, ten6p)
 
         call prsclb('vslope', 'Vertical slope in mrad',
     +        1, z1, 4, ten3p)
         call prsclb('vslope', 'Vertical slope in mrad',
     +        2, z2, 4, ten3p)
 
         call prsclb('hoff', 'Horizontal offset in microns',
     +        1, z1, 1, ten6p)
         call prsclb('hoff', 'Horizontal offset in microns',
     +        2, z2, 1, ten6p)
 
         call prsclb('hslope', 'Horizontal slope in mrad',
     +        1, z1, 2, ten3p)
         call prsclb('hslope', 'Horizontal slope in mrad',
     +        2, z2, 2, ten3p)
         close(108)
c table 2. offset, slope, and dispersion per observation point
         call prsobs('voff', 'Vertical offset in microns',
     +        1, z1, 3, ten6p)
         call prsobs('voff', 'Vertical offset in microns',
     +        2, z2, 3, ten6p)
 
         call prsobs('vslope', 'Vertical slope in mrad',
     +        1, z1, 4, ten3p)
         call prsobs('vslope', 'Vertical slope in mrad',
     +        2, z2, 4, ten3p)
 
         call prsobs('hoff', 'Horizontal offset in microns',
     +        1, z1, 1, ten6p)
         call prsobs('hoff', 'Horizontal offset in microns',
     +        2, z2, 1, ten6p)
 
         call prsobs('hslope', 'Horizontal slope in mrad',
     +        1, z1, 2, ten3p)
         call prsobs('hslope', 'Horizontal slope in mrad',
     +        2, z2, 2, ten3p)
         
c table 3. separations, crossing angles, etc. per crossing point
         do i = 1, npit
           call prsep(msep, i, 3)
           call prsep(msep, i, 1)
         enddo
	endif
c unperturbed beams
      else

 
c table 2. offsets, slopes and dispersions
         call prscl0('Vertical offset in microns',
     +        z1, z2, 3, ten6p)
         call prscl0('Vertical slope in mrad',
     +        z1, z2, 4, ten3p)
         call prscl0('Vertical dispersion in mm',
     +        d1, d2, 3, ten3p)
         call prscl0('Vertical second-order dispersion in m',
     +        dd1, dd2, 3, one)
         call prscl0('Horizontal offset in microns',
     +        z1, z2, 1, ten6p)
         call prscl0('Horizontal slope in mrad',
     +        z1, z2, 2, ten3p)
         call prscl0('Horizontal dispersion in mm',
     +        d1, d2, 1, ten3p)
         call prscl0('Horizontal second-order dispersion in m',
     +        dd1, dd2, 1, one)
      endif
c-----------------------------------------------------------------------
 910  format(' '/' '/'Results including beam-beam collisions:')
 920  format(' '/' '/'Results excluding beam-beam collisions:')
      end
c-----------------------------------------------------------------------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prlumi
      implicit none
      integer i,j,k,l,m
      double precision zz1,zz2
c-----------------------------------------------------------------------
c print relative luminosities
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      double precision sigx, sigy
      double precision fact, beamav(6,mpit,2),zline(mcol)
      double precision zmax(mcol), zmin(mcol), zsum(mcol), zdif(mcol)
      character * 120 temp1
c-----------------------------------------------------------------------
 
      do i = 1, npit
        temp1 = 'lumi.' // name(ipit(i),1)(3:)
        open(ustart+i, file = temp1)
      enddo
 
      write (mulist, 910)
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do l = 1,2
        do i = 1, mpit
          do j = 1, 6
            beamav(j,i,l) = 0
          enddo
        enddo
      enddo
      do i = 1, mpit
        do j = 1, 6
          do l = 1, nbunch
            beamav(j,i,1) = beamav(j,i,1) + z1(j,l,ipit(i))
          enddo
          beamav(j,i,1) = beamav(j,i,1) / nbunch
        enddo
      enddo
      do i = 1, mpit
        do j = 1, 6
          do l = 1, nbunch
            beamav(j,i,2) = beamav(j,i,2) + z2(j,l,ipit(i))
          enddo
          beamav(j,i,2) = beamav(j,i,2) / nbunch
        enddo
      enddo
      do l = 1, nbunch
         do i = 1, npit
            m = ipit(i)
            k = ibnch1(l,m)
            if (k .eq. 0)  then
              zline(i) = 0
            else
              fact = bcurr1(l) * bcurr2(k) / bcurr**2
              zz1 = z1(1,l,m)-z2(1,k,m)
              zz2 = z1(3,l,m)-z2(3,k,m)
              zline(i) = fact *
     +        exp(-(zz1**2 + zz2**2)/(4 * sigx(l,m,1) * sigy(l,m,1)))
            endif
            if (l .eq. 1)  then
              zmax(i) = zline(i)
              zmin(i) = zline(i)
              zsum(i) = zline(i)
            else
              zmax(i) = max(zmax(i), zline(i))
              zmin(i) = min(zmin(i), zline(i))
              zsum(i) = zsum(i) + zline(i)
            endif
         enddo
         if (l .eq. 1 .or. l .eq. nbunch)
     +   write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +        code(l), set_f(l), (' ', zline(i), i = 1, npit), '  '
         do i = 1, npit
           write(ustart+i, '(i5, 4x, f12.6)') set_f(l), zline(i)
         enddo
      enddo
      do i = 1, npit
         close(ustart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
      enddo
 
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/'Relative luminosities:')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pravlumi
      implicit none
      integer i,j,k,l,m
      double precision zz1,zz2, sigx, sigy
c-----------------------------------------------------------------------
c print relative luminosities after elimination of average distance
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      double precision fact, beamav(6,mpit,2),zline(mcol)
      double precision zmax(mcol), zmin(mcol), zsum(mcol), zdif(mcol)
      character * 120 temp1
c-----------------------------------------------------------------------
      do i = 1, npit
        temp1 = 'av_lumi.' // name(ipit(i),1)(3:)
        open(ustart+i, file = temp1)
      enddo
 
      write (mulist, 910)
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do l = 1,2
        do i = 1, mpit
          do j = 1, 6
            beamav(j,i,l) = 0
          enddo
        enddo
      enddo
      do i = 1, mpit
        do j = 1, 6
          do l = 1, nbunch
            beamav(j,i,1) = beamav(j,i,1) + z1(j,l,ipit(i))
          enddo
          beamav(j,i,1) = beamav(j,i,1) / nbunch
        enddo
      enddo
      do i = 1, mpit
        do j = 1, 6
          do l = 1, nbunch
            beamav(j,i,2) = beamav(j,i,2) + z2(j,l,ipit(i))
          enddo
          beamav(j,i,2) = beamav(j,i,2) / nbunch
        enddo
      enddo
      do l = 1, nbunch
         do i = 1, npit
            m = ipit(i)
            k = ibnch1(l,m)
            if (k .eq. 0)  then
              zline(i) = 0
            else
              fact = bcurr1(l) * bcurr2(k) / bcurr**2
              zz1 = (z1(1,l,m)-beamav(1,i,1))-(z2(1,k,m)-beamav(1,i,2))
              zz2 = (z1(3,l,m)-beamav(3,i,1))-(z2(3,k,m)-beamav(3,i,2))
              zline(i) = fact *
     +        exp(-(zz1**2 + zz2**2)/(4 * sigx(l,m,1) * sigy(l,m,1)))
            endif
            if (l .eq. 1)  then
              zmax(i) = zline(i)
              zmin(i) = zline(i)
              zsum(i) = zline(i)
            else
              zmax(i) = max(zmax(i), zline(i))
              zmin(i) = min(zmin(i), zline(i))
              zsum(i) = zsum(i) + zline(i)
            endif
         enddo
         if (l .eq. 1 .or. l .eq. nbunch)
     +   write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +        code(l), set_f(l), (' ', zline(i), i = 1, npit), '  '
         do i = 1, npit
           write(ustart+i, '(i5, 4x, f12.6)') set_f(l), zline(i)
         enddo
      enddo
      do i = 1, npit
         close(ustart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
      enddo
 
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/'Relative luminosities (after shift):')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prscl02(header, z1, z2, index, scale)
      implicit none
      integer i,index
      real z1,z2
      double precision scale,zline
c-----------------------------------------------------------------------
c print value table per interaction point for both unperturbed beams
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      character*(*) header
      dimension z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1)
      dimension zline(mbunch)

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps) 
      save /extel/
      integer itemp,eqvoptmap
c-----------------------------------------------------------------------
      write (mulist, 910) header
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do 10 i = 1, npit
         itemp = eqvoptmap(ipit(i))
         zline(i) = scale * z1(index,1,itemp)
 10   continue
 
      write (mulist, 940)
     +     'F   ', (' ', zline(i), i = 1, npit), '  '
 
      do 20 i = 1, npit
         itemp = eqvoptmap(ipit(i))
         zline(i) = scale * z2(index,1,itemp)
 20   continue
 
      write (mulist, 940)
     +     'B   ', (' ', zline(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/a,':')
 920  format(' ')
 930  format(5x,16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prscl0(header, z1, z2, index, scale)
      implicit none
      integer i,index
      real z1,z2
      double precision scale,zline
c-----------------------------------------------------------------------
c print value table per interaction point for both unperturbed beams
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      character*(*) header
      dimension z1(6,mbunch,0:mcol+1), z2(6,mbunch,0:mcol+1)
      dimension zline(mbunch)
c-----------------------------------------------------------------------
      write (mulist, 910) header
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do 10 i = 1, npit
         zline(i) = scale * z1(index,1,ipit(i))
 10   continue
 
      write (mulist, 940)
     +     'F   ', (' ', zline(i), i = 1, npit), '  '
 
      do 20 i = 1, npit
         zline(i) = scale * z2(index,1,ipit(i))
 20   continue
 
      write (mulist, 940)
     +     'B   ', (' ', zline(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/a,':')
 920  format(' ')
 930  format(5x,16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsclb2(fhead, header, bd, z, index, scale)
      implicit none
      integer i,index,k,l,bd,lastnb
      real z
      double precision scale,zdif,zline,zmax,zmin,zsum,muline
c-----------------------------------------------------------------------
c print value table per bunch and all pits for one beam
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/ 
c------------------------------------------------------------------------
      double precision qx,qy,q
      common /globtun/ qx,qy
      save/ globtun/

      double precision kick
      common/ kck/ kick(2,mbunch,mcol,2)
      save /kck/
      
      double precision orbdev
      double precision betaline, beta
      double precision kicktemp
      double precision templist(nbunch)
      double precision diflist(mcol)
      double precision temp(2,mbunch,2)
      common /tempc/ temp
      save /tempc/
c-----------------------------------------------------------------------
      character*(*) fhead, header
      character * 8  dir
      character * 1000 temp1, temp2, temp3
      dimension z(6,mbunch,0:mmaps+1)
      dimension zline(mbunch), zmax(mbunch), zmin(mbunch), zsum(mbunch)
      dimension zdif(mbunch),muline(mbunch),betaline(mbunch),
     +beta(mbunch)

      character*120 filename, chst,filenm
      common /filesch/ filenm,chst
      save / filesch/ 
      integer nstart,filestart,j
      character*4 kstr
      double precision difmatr(1,nbunch,nbunch),
     +tunematr(2,nbunch,nbunch)
      double precision chromatr(2,nbunch,nbunch)
      double precision enj,enl,betanew,sigmaj,sigmal,sigmapj,sigmapl,xj,
     +xl,xpj,xpl,ql,qj,qpl,qpj,r,phi
      character*3 part
      double precision dmatmax,rms
      double precision difarray(nbunch*(nbunch-1)/2)
      integer diflen,exptemp 
      real(kind = 8) :: cum,cumn,cumtemp,maxtol
      logical cumflag
      

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/
      integer itemp,eqvoptmap
c-----------------------------------------------------------------------

      if (bd .eq. 1)  then
        dir = 'forward'
        temp1 = fhead(:lastnb(fhead)) // '_f'

      else
        dir = 'backward'
        temp1 = fhead(:lastnb(fhead)) // '_b'
      endif
      write (mulist, 910) header, dir
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      
      do 10 i = 1, npit
         itemp = eqvoptmap(ipit(i))
         zmax(i) = scale * z(index,1,itemp)
         zmin(i) = zmax(i)
         zsum(i) = 0.0
 10   continue
      do i = 1, npit
        temp2 = temp1(:lastnb(temp1)) // '.' // name(ipit(i),1)(3:)
        open(ustart+i, file = temp2)
      enddo

      do 30 l = 1, nbunch
         do 20 i = 1, npit
            itemp = eqvoptmap(ipit(i))
            zline(i) = scale * z(index,l,itemp)
            zmax(i) = max(zmax(i), zline(i))
            zmin(i) = min(zmin(i), zline(i))
            zsum(i) = zsum(i) + zline(i)

 20      continue
         if (bd .eq. 1)  then
           k = set_f(l)
         else
           k = set_b(l)
         endif
         
         if (l .eq. 1 .or. l .eq. nbunch)  then
           write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +     code(l), k, (' ', zline(i), i = 1, npit), '  '
         endif
         do i = 1, npit

           write(ustart+i, '(i5, 4x, f12.6,4x,f12.6)') k, zline(i)

         enddo
 30   continue

      do 40 i = 1, npit
         close(ustart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
 40   continue    

c Select the maximum variation in closed orbit of all interaction points
c mcol, corresponding to the collision points where a long range interaction 
c takes place
      
      dmatmax = 0
      if (temp2(2:4) .eq. 'off') then

       call matrzero(difmatr,1,nbunch,nbunch)
       call dzero(difarray,nbunch*(nbunch-1)/2)

       i = 1
       
       temp2 = 'fort.ps.'//temp1(:lastnb(temp1))//name(ipit(i),1)(3:)
       open(764 + i,file = temp2)

       if (bd .eq. 1) then
         part = '_f.'
       else
         part = '_b.'
       endif


       temp3 = 'fort.ip'//part//temp1(1:1)//'.'//name(ipit(i),1)(3:)

       open(1785 + i,file = temp3)
       do j = nstart(),nbunch
         do l = nstart(),nbunch

            if (temp1(1:1) .eq. 'h') then
               enj = epsx(j,bd)*gev/tmass
               enl = epsx(l,bd)*gev/tmass
               betanew = betx(ipit(i),bd)

            else
               enj = epsy(j,bd)*gev/tmass
               enl = epsy(l,bd)*gev/tmass
               betanew = bety(ipit(i),bd)
                 

            endif
              
            sigmaj = sqrt((enj*betanew)/gamma)
            sigmal = sqrt((enl*betanew)/gamma)
            sigmapj = sqrt(enj/(betanew*gamma))
            sigmapl = sqrt(enl/(betanew*gamma))
            itemp = eqvoptmap(ipit(i))  
            xj = z(index,j,itemp)/sigmaj
            xl = z(index,l,itemp)/sigmal
            xpj = z(index+1,j,itemp)/sigmapj
            xpl = z(index+1,l,itemp)/sigmapl

            difmatr(i,j,l) = sqrt((xj-xl)**2 + (xpj-xpl)**2)

   
            if (difmatr(i,j,l) .gt. dmatmax) dmatmax = 
     +      difmatr(i,j,l)

            if (difmatr(i,j,l) .eq. 0) then
              if (bd .eq. 1) then
                 write(1785 + i,*) ipit(i),set_f(j),sigmaj,
     +           sigmapj,sqrt(xj**2+xpj**2),atan2(xpj,xj)
              else
                 write(1785 + i,*)ipit(i),set_b(j),sigmaj,
     +           sigmapj,sqrt(xj**2+xpj**2),atan2(xpj,xj)
              endif      
            endif
          enddo
          r = sqrt(xj**2+xpj**2)
          phi = atan2(xpj,xj)
          if (bd .eq. 1) then
             write(764 + i,*) set_f(j),xj,xpj,r,phi
          else
             write(764 + i,*) set_b(j),xj,xpj,r,phi
          endif
       enddo
       close(764 + i)
       close(1785 + i)

       diflen = 0 
       do i = nstart(),nbunch
         do j = i + 1, nbunch
           diflen = diflen + 1
           difarray(diflen) = difmatr(1,i,j)
         enddo
       enddo
       cumtemp = 0 
       cum = 0 
       
       maxtol =  10.0d0**10         
       do i = 1, diflen
         cum = cum +  difarray(i)**2
         if (cum .ge. maxtol) then
            cumtemp = cum/diflen + cumtemp
            cum = 0
            cumflag = .True.
         endif
         
       enddo

       if (cumflag) then 
          cumn = cum/diflen + cumtemp
       else
          cumn = cum/diflen
       endif

       rms = sqrt(cumn)
       
       open(108,file='fort.'//temp2)
       
       write(108,*) 'MAXVAL orbit difference offset PS ',temp1(1:7), 
     + dmatmax
       write(108,*) 'RMS orbit difference offset ',temp1(1:7),
     + rms
       close(108)
       write(*,*) 'MAXVAL orbit difference offset PS ',temp1(1:7), 
     + dmatmax
       write(*,*) 'RMS orbit difference offset ',temp1(1:7),
     + rms
      endif  
      
            
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)

c-----------------------------------------------------------------------
 910  format(' '/a,' for ',a,' beam:')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end  

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsclb(fhead, header, bd, z, index, scale)
      implicit none
      integer i,index,k,l,bd,lastnb
      real z
      double precision scale,zdif,zline,zmax,zmin,zsum,muline
c-----------------------------------------------------------------------
c print value table per bunch and all pits for one beam
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/ 
c------------------------------------------------------------------------
      double precision qx,qy,q
      common /globtun/ qx,qy
      save/ globtun/

      double precision kick
      common/ kck/ kick(2,mbunch,mcol,2)
      save /kck/
      
      double precision orbdev
      double precision betaline, beta
      double precision kicktemp
      double precision templist(nbunch)
      double precision diflist(mcol)
      double precision temp(2,mbunch,2)
      common /tempc/ temp
      save /tempc/
c-----------------------------------------------------------------------
      character*(*) fhead, header
      character * 8  dir
      character * 1000 temp1, temp2, temp3
      dimension z(6,mbunch,0:mcol+1)
      dimension zline(mbunch), zmax(mbunch), zmin(mbunch), zsum(mbunch)
      dimension zdif(mbunch),muline(mbunch),betaline(mbunch),
     +beta(mbunch)

      character*120 filename, chst,filenm
      common /filesch/ filenm,chst
      save / filesch/ 
      integer nstart,filestart,j
      character*4 kstr
      double precision difmatr(1,nbunch,nbunch),
     +tunematr(2,nbunch,nbunch)
      double precision chromatr(2,nbunch,nbunch)
      double precision enj,enl,betanew,sigmaj,sigmal,sigmapj,sigmapl,xj,
     +xl,xpj,xpl,ql,qj,qpl,qpj,r,phi
      character*3 part
      double precision dmatmax,rms
      double precision difarray(nbunch*(nbunch-1)/2)
      integer diflen,exptemp 
      real(kind = 8) :: cum,cumn,cumtemp,maxtol
      logical cumflag
      
       
c-----------------------------------------------------------------------

      if (bd .eq. 1)  then
        dir = 'forward'
        temp1 = fhead(:lastnb(fhead)) // '_f'

      else
        dir = 'backward'
        temp1 = fhead(:lastnb(fhead)) // '_b'
      endif
      write (mulist, 910) header, dir
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do 10 i = 1, npit
         zmax(i) = scale * z(index,1,ipit(i))
         zmin(i) = zmax(i)
         zsum(i) = 0.0
 10   continue

      do i = 1, npit
        temp2 = temp1(:lastnb(temp1)) // '.' // name(ipit(i),1)(3:)
        open(ustart+i, file = temp2)
      enddo
      do 30 l = 1, nbunch
         do 20 i = 1, npit
            zline(i) = scale * z(index,l,ipit(i))
            zmax(i) = max(zmax(i), zline(i))
            zmin(i) = min(zmin(i), zline(i))
            zsum(i) = zsum(i) + zline(i)

 20      continue
         if (bd .eq. 1)  then
           k = set_f(l)
         else
           k = set_b(l)
         endif
         if (l .eq. 1 .or. l .eq. nbunch)  then
           write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +     code(l), k, (' ', zline(i), i = 1, npit), '  '
         endif
         do i = 1, npit

           write(ustart+i, '(i5, 4x, f12.6,4x,f12.6)') k, zline(i)

         enddo
 30   continue
 
      do 40 i = 1, npit
         close(ustart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
 40   continue
     
c Select the maximum variation in closed orbit of all interaction points
c mcol, corresponding to the collision points where a long range interaction 
c takes place

      dmatmax = 0
      if (temp2(2:4) .eq. 'off') then

        
       call matrzero(difmatr,1,nbunch,nbunch)
       call dzero(difarray,nbunch*(nbunch-1)/2)

       i = 1
       
       temp2 = 'fort.ps.'//temp1(:lastnb(temp1))//name(ipit(i),1)(3:)
       open(764 + i,file = temp2)

       if (bd .eq. 1) then
         part = '_f.'
       else
         part = '_b.'
       endif


       temp3 = 'fort.ip'//part//temp1(1:1)//'.'//name(ipit(i),1)(3:)

       open(1785 + i,file = temp3)
       do j = nstart(),nbunch
         do l = nstart(),nbunch

            if (temp1(1:1) .eq. 'h') then
               enj = epsx(j,bd)*gev/tmass
               enl = epsx(l,bd)*gev/tmass
               betanew = betx(ipit(i),bd)

            else
               enj = epsy(j,bd)*gev/tmass
               enl = epsy(l,bd)*gev/tmass
               betanew = bety(ipit(i),bd)
                 

            endif
              
            sigmaj = sqrt((enj*betanew)/gamma)
            sigmal = sqrt((enl*betanew)/gamma)
            sigmapj = sqrt(enj/(betanew*gamma))
            sigmapl = sqrt(enl/(betanew*gamma))
              
            xj = z(index,j,ipit(i))/sigmaj
            xl = z(index,l,ipit(i))/sigmal
            xpj = z(index+1,j,ipit(i))/sigmapj
            xpl = z(index+1,l,ipit(i))/sigmapl

            difmatr(i,j,l) = sqrt((xj-xl)**2 + (xpj-xpl)**2)

   
            if (difmatr(i,j,l) .gt. dmatmax) dmatmax = 
     +      difmatr(i,j,l)


            if (difmatr(i,j,l) .eq. 0) then
              if (bd .eq. 1) then
                 write(1785 + i,*) ipit(i),set_f(j),sigmaj,
     +           sigmapj,sqrt(xj**2+xpj**2),atan2(xpj,xj)
              else
                 write(1785 + i,*)ipit(i),set_b(j),sigmaj,
     +           sigmapj,sqrt(xj**2+xpj**2),atan2(xpj,xj)
              endif      
            endif
          enddo
          r = sqrt(xj**2+xpj**2)
          phi = atan2(xpj,xj)
          if (bd .eq. 1) then
             write(764 + i,*) set_f(j),xj,xpj,r,phi
          else
             write(764 + i,*) set_b(j),xj,xpj,r,phi
          endif
       enddo
        write(*,*) 'MAXVAL orbit difference offset PS ',temp1(1:7), 
     + dmatmax

       close(764 + i)
       close(1785 + i)

       diflen = 0 
       do i = nstart(),nbunch
         do j = i + 1, nbunch
           diflen = diflen + 1
           difarray(diflen) = difmatr(1,i,j)
         enddo
       enddo
       cumtemp = 0 
       cum = 0 
       
       maxtol =  10.0d0**10         
       do i = 1, diflen
         cum = cum +  difarray(i)**2
         if (cum .ge. maxtol) then
            cumtemp = cum/diflen + cumtemp
            cum = 0
            cumflag = .True.
         endif
         
       enddo

       if (cumflag) then 
          cumn = cum/diflen + cumtemp
       else
          cumn = cum/diflen
       endif

       rms = sqrt(cumn)
       
       open(108,file='fort.'//temp2)
       
       write(108,*) 'MAXVAL orbit difference offset PS ',temp1(1:7), 
     + dmatmax
       write(108,*) 'RMS orbit difference offset',temp1(1:7),
     + rms
       close(108)
       write(*,*) 'MAXVAL orbit difference offset PS ',temp1(1:7), 
     + dmatmax
       write(*,*) 'RMS orbit difference offset',temp1(1:7),
     + rms
      endif  
      
    

      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)

c-----------------------------------------------------------------------
 910  format(' '/a,' for ',a,' beam:')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsobs2(fhead, header, bd, z, index, scale)
      implicit none
      integer i,index,k,l,bd,lastnb
      real z
      double precision scale,zdif,zline,zmax,zmin,zsum
c-----------------------------------------------------------------------
c print value table per bunch at observation points for one beam
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      character*(*) fhead, header
      character * 8  dir
      character * 1000 temp1, temp2, temp3
      dimension z(6,mbunch,0:mmaps+1)
      dimension zline(mbunch), zmax(mbunch), zmin(mbunch), zsum(mbunch)
      dimension zdif(mbunch)

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/
      integer itemp,eqvoptmap
c-----------------------------------------------------------------------

      if (bd .eq. 1)  then
        dir = 'forward'
        temp1 = fhead(:lastnb(fhead)) // '_f'
      else
        dir = 'backward'
        temp1 = fhead(:lastnb(fhead)) // '_b'
      endif
      temp3 = '.obs.'
      write (mulist, 910) header, dir
      write (mulist, 920)
      write (mulist, 930) (' ', outblist(i), i = 1, outbcnt), '  '
      write (mulist, 920)
 
      do 10 i = 1, outbcnt
         itemp = eqvoptmap(outblist(i))
         zmax(i) = scale * z(index,1,itemp)
         zmin(i) = zmax(i)
         zsum(i) = 0.0
 10   continue
      do i = 1, outbcnt
        write(temp3(5:8), '(i4.4)') outblist(i)
        temp2 = temp1(:lastnb(temp1)) // temp3(:lastnb(temp3))
        open(ustart+i, file = temp2)
      enddo
      do 30 l = 1, nbunch
         do 20 i = 1, outbcnt
            itemp = eqvoptmap(outblist(i))
            zline(i) = scale * z(index,l,itemp)
            zmax(i) = max(zmax(i), zline(i))
            zmin(i) = min(zmin(i), zline(i))
            zsum(i) = zsum(i) + zline(i)
 20      continue
         if (bd .eq. 1)  then
           k = set_f(l)
         else
           k = set_b(l)
         endif
         if (l .eq. 1 .or. l .eq. nbunch)  then
           write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +     code(l), k, (' ', zline(i), i = 1, outbcnt), '  '
         endif
         do i = 1, outbcnt
           write(ustart+i, '(i5, 4x, f12.6)') k, zline(i)
         enddo
 30   continue
 
      do 40 i = 1, outbcnt
         close(ustart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
 40   continue
 
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, outbcnt), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, outbcnt), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, outbcnt), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, outbcnt), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/a,' for ',a,' beam:')
 920  format(' ')
 930  format('bunch bucket',16('   ',a,4x,i4,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end



c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsobs(fhead, header, bd, z, index, scale)
      implicit none
      integer i,index,k,l,bd,lastnb
      real z
      double precision scale,zdif,zline,zmax,zmin,zsum
c-----------------------------------------------------------------------
c print value table per bunch at observation points for one beam
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      character*(*) fhead, header
      character * 8  dir
      character * 1000 temp1, temp2, temp3
      dimension z(6,mbunch,0:mcol+1)
      dimension zline(mbunch), zmax(mbunch), zmin(mbunch), zsum(mbunch)
      dimension zdif(mbunch)
c-----------------------------------------------------------------------

      if (bd .eq. 1)  then
        dir = 'forward'
        temp1 = fhead(:lastnb(fhead)) // '_f'
      else
        dir = 'backward'
        temp1 = fhead(:lastnb(fhead)) // '_b'
      endif
      temp3 = '.obs.'
      write (mulist, 910) header, dir
      write (mulist, 920)
      write (mulist, 930) (' ', outblist(i), i = 1, outbcnt), '  '
      write (mulist, 920)
 
      do 10 i = 1, outbcnt
         zmax(i) = scale * z(index,1,outblist(i))
         zmin(i) = zmax(i)
         zsum(i) = 0.0
 10   continue
      do i = 1, outbcnt
        write(temp3(5:8), '(i4.4)') outblist(i)
        temp2 = temp1(:lastnb(temp1)) // temp3(:lastnb(temp3))
        open(ustart+i, file = temp2)
      enddo
      do 30 l = 1, nbunch
         do 20 i = 1, outbcnt
            zline(i) = scale * z(index,l,outblist(i))
            zmax(i) = max(zmax(i), zline(i))
            zmin(i) = min(zmin(i), zline(i))
            zsum(i) = zsum(i) + zline(i)
 20      continue
         if (bd .eq. 1)  then
           k = set_f(l)
         else
           k = set_b(l)
         endif
         if (l .eq. 1 .or. l .eq. nbunch)  then
           write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +     code(l), k, (' ', zline(i), i = 1, outbcnt), '  '
         endif
         do i = 1, outbcnt
           write(ustart+i, '(i5, 4x, f12.6)') k, zline(i)
         enddo
 30   continue
 
      do 40 i = 1, outbcnt
         close(ustart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
 40   continue
 
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, outbcnt), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, outbcnt), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, outbcnt), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, outbcnt), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/a,' for ',a,' beam:')
 920  format(' ')
 930  format('bunch bucket',16('   ',a,4x,i4,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsep2(unit, pitnum, index)
      implicit none
      integer i, index, k, unit, pitnum
      double precision zz1,zz2,norm
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      double precision sigx, sigy
      character * 120 filename,filename2

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/
      integer itemp,eqvoptmap
      double precision avsep
      double precision sig1(mbunch),sig2(mbunch)

      if (index .eq. 1)  then
        filename = 'hsep_sig.' // pitnam(pitnum)
        norm = sigx(1,ipit(pitnum),1)
        do i = 1,nbunch
           sig1(i) = sigx(i,ipit(pitnum),1)
           sig2(i) = sigx(i,ipit(pitnum),2)
        enddo
      else
        filename = 'vsep_sig.' // pitnam(pitnum)
        norm = sigy(1,ipit(pitnum),1)
        do i = 1,nbunch
           sig1(i) = sigy(i,ipit(pitnum),1)
           sig2(i) = sigy(i,ipit(pitnum),2)
        enddo
      endif
      open(unit, file = filename, status = 'UNKNOWN')
     
      avsep = 0
      do i = 1, nbunch
        itemp = eqvoptmap(ipit(pitnum))
        k = ibnch1(i,ipit(pitnum))

        if (k .ne. 0)  then
          zz1 = z1(index,i,itemp)
          zz2 = z2(index,k,itemp)
          avsep = avsep + (zz2-zz1)
          write(unit, '(i5, 2x, f14.8)') set_f(i), (zz2-zz1)/norm
        endif
      enddo
      close(unit)

      avsep = avsep/nbunch 

      if (index .eq. 1)  then
        filename = 'hsep_mu.' // pitnam(pitnum)
        filename2 = 'hsep_mu_norm.' // pitnam(pitnum) 
      else
        filename = 'vsep_mu.' // pitnam(pitnum)
        filename2 = 'vsep_mu_norm.' // pitnam(pitnum) 
      endif
      open(unit, file = filename, status = 'UNKNOWN')
      open(unit+2, file = filename2, status = 'UNKNOWN')

      do i = 1, nbunch
        itemp = eqvoptmap(ipit(pitnum))
        k = ibnch1(i,ipit(pitnum))
        if (k .ne. 0)  then
          zz1 = z1(index,i,itemp)
          zz2 = z2(index,k,itemp)
          write(unit, '(i5, 2x, f14.8)') set_f(i), 1.d6*(zz2-zz1)
          write(unit+2,*) set_f(i),abs(zz2-zz1-avsep),
     +sig1(i),sig2(k)! Absolute value of the separation of both beams in micrometers once the average separation has been corrected.

        endif
      enddo
      close(unit)
      close(unit+2)
     
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsep(unit, pitnum, index)
      implicit none
      integer i, index, k, unit, pitnum
      double precision zz1,zz2,norm
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      double precision sigx, sigy
      character * 120 filename, filename2
 
      if (index .eq. 1)  then
        filename = 'hsep_sig.' // pitnam(pitnum)
        norm = sigx(1,ipit(pitnum),1)
      else
        filename = 'vsep_sig.' // pitnam(pitnum)
        norm = sigy(1,ipit(pitnum),1)
      endif
      open(unit, file = filename, status = 'UNKNOWN')
 
      do i = 1, nbunch
        k = ibnch1(i,ipit(pitnum))
        if (k .ne. 0)  then
          zz1 = z1(index,i,ipit(pitnum))
          zz2 = z2(index,k,ipit(pitnum))
          write(unit, '(i5, 2x, f14.8)') set_f(i), (zz2-zz1)/norm
        endif
      enddo
      close(unit)
      if (index .eq. 1)  then
        filename = 'hsep_mu.' // pitnam(pitnum)
      else
        filename = 'vsep_mu.' // pitnam(pitnum)
      endif
      open(unit, file = filename, status = 'UNKNOWN')
      do i = 1, nbunch
        k = ibnch1(i,ipit(pitnum))
        if (k .ne. 0)  then
          zz1 = z1(index,i,ipit(pitnum))
          zz2 = z2(index,k,ipit(pitnum))
          write(unit, '(i5, 2x, f14.8)') set_f(i), 1.d6*(zz2-zz1)
        endif
      enddo
      close(unit)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prshft
      implicit none
      integer i,k,l,m
      double precision denom,denum,zdif,zline,zmax,
     +zmin,zsum
c-----------------------------------------------------------------------
c print CM energy shift
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      dimension zline(mcol)
      dimension zmax(mcol), zmin(mcol), zsum(mcol), zdif(mcol)
c-----------------------------------------------------------------------
      write (mulist, 910)
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do l = 1, nbunch
         do i = 1, npit
            m = ipit(i)
            k = ibnch1(l,m)
            if (k .eq. 0)  then
              zline(i) = 0
            else
              denum = gev * deltap**2 * (d1(3,l,m) - d2(3,k,m))
              denom = 2.0 * epsy0 * bety(m,1) +
     +           deltap**2 * (d1(3,l,m)**2 + d2(3,k,m)**2)
              zline(i) = - ten3p * (z1(3,l,m) - z2(3,k,m)) *
     +           denum / denom
            endif
            if (l .eq. 1)  then
              zmax(i) = zline(i)
              zmin(i) = zline(i)
              zsum(i) = zline(i)
            else
              zmax(i) = max(zmax(i), zline(i))
              zmin(i) = min(zmin(i), zline(i))
              zsum(i) = zsum(i) + zline(i)
            endif
         enddo
         write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +        code(l), set_f(l), (' ', zline(i), i = 1, npit), '  '
      enddo
 
      do i = 1, npit
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
      enddo
 
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/'CM energy shift in MeV:')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prsprd
      implicit none
      integer i,k,l,m
      double precision denom,denum,zdif,zline,zmax,
     +zmin,zsum
c-----------------------------------------------------------------------
c print CM energy spread
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      dimension zline(mcol)
      dimension zmax(mcol), zmin(mcol), zsum(mcol), zdif(mcol)
c-----------------------------------------------------------------------
      write (mulist, 910)
      write (mulist, 920)
      write (mulist, 930) (' ', name(ipit(i),1)(3:), i = 1, npit), '  '
      write (mulist, 920)
 
      do l = 1, nbunch
         do i = 1, npit
            m = ipit(i)
            k = ibnch1(l,m)
            if (k .eq. 0)  then
              zline(i) = 0
            else
              denum = deltap**2 * (d1(3,l,m) + d2(3,l,m))**2 +
     +           4.0 * epsy0 * bety(m,1)
              denom = 2.0 * epsy0 * bety(m,1) +
     +           deltap**2 * (d1(3,l,m)**2 + d2(3,l,m)**2)
              zline(i) = ten3p * gev * deltap *
     +           sqrt(denum / denom)
            endif
            if (l .eq. 1)  then
              zmax(i) = zline(i)
              zmin(i) = zline(i)
              zsum(i) = zline(i)
            else
              zmax(i) = max(zmax(i), zline(i))
              zmin(i) = min(zmin(i), zline(i))
              zsum(i) = zsum(i) + zline(i)
            endif
         enddo
         write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +        code(l), set_f(l), (' ', zline(i), i = 1, npit), '  '
      enddo
 
      do i = 1, npit
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
      enddo
 
      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/'CM energy spread in MeV:')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prtune(flag)
      implicit none
      integer i,l
      double precision zave,zdif,zline,zmax,zmin
c-----------------------------------------------------------------------
c print tunes and chromaticities
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
c-----------------------------------------------------------------------
      logical flag
      dimension zline(5), zmax(5), zmin(5), zave(5), zdif(5)
      character*120 filename, chst,filenm
      common /filesch/ filenm,chst
      save / filesch/
      double precision tunematrx(2,nbunch,nbunch),
     +tunematry(2,nbunch,nbunch),chromatrx(2,nbunch,nbunch),
     +chromatry(2,nbunch,nbunch)
      integer nstart
c-----------------------------------------------------------------------
c perturbed beams
      if (flag) then
 
c forward beam

         open(ustart+1,file='tune_f.list')
         write (mulist, 910) 'forward'
         write (mulist, 920)
         zline(1) = bcurr1(1)
         zline(2) = qx1(1)
         zline(3) = qy1(1)
c-- chromaticity
         zline(4) = qxp1(1)
         zline(5) = qyp1(1)
c-- eigen-tune
!         zline(4) = q11(1)
!         zline(5) = q21(1)
 
         do 10 i = 1, 5
            zmax(i) = zline(i)
            zmin(i) = zline(i)
            zave(i) = 0.0
 10      continue
 
         do 30 l = 1, nbunch
            zline(1) = bcurr1(l)
            zline(2) = qx1(l)
            zline(3) = qy1(l)
c-- chromaticity
            zline(4) = qxp1(l)
            zline(5) = qyp1(l)
c-- eigen-tune
!            zline(4) = q11(l)
!            zline(5) = q21(l)
 
            do 20 i = 1, 5
               zmax(i) = max(zmax(i), zline(i))
               zmin(i) = min(zmin(i), zline(i))
               zave(i) = zave(i) + zline(i)
 20         continue
            if (l .eq. 1 .or. l .eq. nbunch)  then
              write (mulist, '(a5,i6,5('' '',f11.6),''   '')')
     +        code(l), set_f(l), (zline(i), i = 1, 5)
            endif
            write(ustart+1, *) set_f(l),
     +      (zline(i), i = 1, 5) 
 30      continue
         close(ustart+1)
         do 40 i = 1, 5
            zave(i) = zave(i) / nbunch
            zdif(i) = zmax(i) - zmin(i)
 40      continue

         if (filenm(16:24) .eq. '25ns_2760'(1:9) ) then
              nstart = 25
          else if (filenm(16:24) .eq. '25ns_2748'(1:9)) then
              nstart = 13
          else
              nstart = 13
          endif

         call matrzero(tunematrx,2,nbunch,nbunch)
         call matrzero(chromatrx,2,nbunch,nbunch)
         do i = nstart, nbunch
            do l = nstart,nbunch
               tunematrx(1,i,l) = abs(qx1(l) - qx1(i))
               chromatrx(1,i,l) = abs(qxp1(l) - qxp1(i))
               tunematry(1,i,l) = abs(qy1(l) - qy1(i))
               chromatry(1,i,l) = abs(qyp1(l) - qyp1(i))
               tunematrx(2,i,l) = abs(qx2(l) - qx2(i))
               chromatrx(2,i,l) = abs(qxp2(l) - qxp2(i))
               tunematry(2,i,l) = abs(qy2(l) - qy2(i))
               chromatry(2,i,l) = abs(qyp2(l) - qyp2(i))
            
            enddo
         enddo
         open(ustart+3,file='fort.tuneschromashift')
         write(ustart+3,*) 'Maximum horizontal tune shift, b1',
     +   maxval(tunematrx(1,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum vertical tune shift, b1',
     +   maxval(tunematry(1,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum horizontal chromaticity shift, b1',
     +   maxval(chromatrx(1,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum vertical chromaticitye shift, b1',
     +   maxval(chromatry(1,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum horizontal tune shift, b2',
     +   maxval(tunematrx(2,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum vertical tune shift, b2',
     +   maxval(tunematry(2,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum horizontal chromaticity shift, b2',
     +   maxval(chromatrx(2,nstart:nbunch,nstart:nbunch))
         write(ustart+3,*) 'Maximum vertical chromaticitye shift, b2',
     +   maxval(chromatry(2,nstart:nbunch,nstart:nbunch))
         close(ustart+3)
         
         write(*,*) 'Maximum horizontal tune shift, b1',
     +   maxval(tunematrx(1,nstart:nbunch,nstart:nbunch))
         write(*,*) 'Maximum vertical tune shift, b1',
     +   maxval(tunematry(1,nstart:nbunch,nstart:nbunch))
         write(*,*) 'Maximum horizontal chromaticity shift, b1',
     +   maxval(chromatrx(1,nstart:nbunch,nstart:nbunch))
         write(*,*) 'Maximum vertical chromaticitye shift, b1',
     +   maxval(chromatry(1,nstart:nbunch,nstart:nbunch))
         
         write (mulist, 920)
         write (mulist, 930) 'mx  ', (zmax(i), i = 1, 5)
         write (mulist, 930) 'mn  ', (zmin(i), i = 1, 5)
         write (mulist, 930) 'av  ', (zave(i), i = 1, 5)
         write (mulist, 930) 'df  ', (zdif(i), i = 1, 5)
         write (mulist, 920)
 
c backward beam
         open(ustart+1,file='tune_b.list')
         write (mulist, 910) 'backward'
         write (mulist, 920)
         zline(1) = bcurr2(1)
         zline(2) = qx2(1)
         zline(3) = qy2(1)
c-- chromaticity
         zline(4) = qxp2(1)
         zline(5) = qyp2(1)
c-- eigen-tune
!         zline(4) = q12(1)
!         zline(5) = q22(1)
 
         do 60 i = 1, 5
            zmax(i) = zline(i)
            zmin(i) = zline(i)
            zave(i) = 0.0
 60      continue
 
         do 80 l = 1, nbunch
            zline(1) = bcurr2(l)
            zline(2) = qx2(l)
            zline(3) = qy2(l)
c-- chromaticity
            zline(4) = qxp2(l)
            zline(5) = qyp2(l)
c-- eigen-tune
!            zline(4) = q12(l)
!            zline(5) = q22(l)
 
            do 70 i = 1, 5
               zmax(i) = max(zmax(i), zline(i))
               zmin(i) = min(zmin(i), zline(i))
               zave(i) = zave(i) + zline(i)
 70         continue
            if (l .eq. 1 .or. l .eq. nbunch)  then
              write (mulist, '(a5,i6,5('' '',f12.8),''   '')')
     +        code(l), set_b(l), (zline(i), i = 1, 5)
            endif
            write(ustart+1,*) set_b(l), (zline(i), i = 1, 5)
 80      continue
         close(ustart+1)
 
         do 90 i = 1, 5
            zave(i) = zave(i) / nbunch
            zdif(i) = zmax(i) - zmin(i)
 90      continue
         write (mulist, 920)
         write (mulist, 930) 'mx  ', (zmax(i), i = 1, 5)
         write (mulist, 930) 'mn  ', (zmin(i), i = 1, 5)
         write (mulist, 930) 'av  ', (zave(i), i = 1, 5)
         write (mulist, 930) 'df  ', (zdif(i), i = 1, 5)
      else
 
c unperturbed beams
         write (mulist, 940)
         write (mulist, 920)
         write (mulist, 930)
     +        'F   ', bcurr1(1), qx1(1), qy1(1), qxp1(1), qyp1(1)
         write (mulist, 930)
     +        'B   ', bcurr2(1), qx2(1), qy2(1), qxp2(1), qyp2(1)
      endif
 
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/'Tune and chromaticity for ',a,' beam:'/
     +     'bunch bucket     bcurr       qx         qy     ',
     +     '   qx''      qy'' ')
c     +     '      q1          q2  ')
 920  format(' ')
 930  format(a5,5('  ',f12.8),'   ')
 940  format(' '/'Tune and chromaticity:'/
     +     '            bcurr         qx            qy     ',
     +     '       qx''           qy''       ')
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rdmaps2
      implicit none
c-----------------------------------------------------------------------
c fill tables for transfer maps
c-----------------------------------------------------------------------
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      double precision sk1,sk2,sr1,sr2,st1,st2
     
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/      

c-----------------------------------------------------------------------
      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/   
c-----------------------------------------------------------------------
      character*200 templine
      integer i, j, i1, i2, i3, ibeam,iostat
      integer lastnb
      character * 48 t_name
      double precision ss1, ss2
      character *(mcnam) dummy
      integer cont1,cont2
      character*(mcnam)    name2(mcol)
c-----------------------------------------------------------------------
      
      if (extraelem) then ! Read total number of maps. nmaps >= ninter
         t_name = 'train.manf'
         nmaps = 0
         
         open(iunit, file=t_name, status='OLD', iostat=iostat)
         if (iostat .ne. 0)  then
            print *, '<<< TRAIN >>> fatal: MAP file ',
     +      t_name(:lastnb(t_name)), ' not found'
            stop
         endif
112      continue
         read (iunit, *, end = 1120) templine
         if (templine(1:3) .eq. 'LHC') then
            goto 1120 
         else
            nmaps = nmaps + 1
            goto 112
         endif 
1120  continue
      close(iunit)   
               
      else
         nmaps = ninter
         
      endif 
      
      do ibeam = 1, 2
        if (ibeam .eq. 1)  then
          i1 = 0
          i2 = nmaps 
          i3 = 1
          t_name = 'train.manf'
        else
          i1 = nmaps 
          i2 = 0
          i3 = -1
          t_name = 'train.manb'
        endif

        open(iunit, file=t_name,status='OLD', iostat=iostat)
            if (iostat .ne. 0)  then          
              print *, '<<< TRAIN >>> fatal: MAP file ',
     +        t_name(:lastnb(t_name)), ' not found'
              stop
            endif

        do i = i1, i2, i3
           if (ibeam .eq. 1)  then
              read (iunit, *, end=1000) dummy,ss1,(sk1(j,i), j=1,6),
     +	      (sr1(j,1,i), j=1,36),(st1(j,1,1,i), j=1,216) 
              read(dummy(1:lastnb(dummy)),*)  mname(i)
              
           else
              read (iunit, *, end=1000) dummy,ss2,(sk2(j,i), j=1,6),
     +        (sr2(j,1,i), j=1,36),(st2(j,1,1,i), j=1,216)
           endif
           
        enddo
        close(iunit)
      enddo

      
      

      do j = 0,nmaps - 1
         do i = 1,ninter
            if (mname(j)(1:lastnb(mname(j))).eq.
     +name(i,1)(1:lastnb(name(i,1)))) mapmask(j + 1) = 1
         enddo
      enddo

      return
1000  continue
        
      print *, '<<< TRAIN >>> fatal: MAP file',
     +t_name(:lastnb(t_name)), ' empty or too short'
      stop

      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rdmaps
      implicit none
c-----------------------------------------------------------------------
c fill tables for transfer maps
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      double precision sk1,sk2,sr1,sr2,st1,st2
      double precision gk1,gk2,gr1,gr2,gt1,gt2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      common /sectorn/ gk1(6,0:mcol), gr1(6,6,0:mcol),
     +                gt1(6,6,6,0:mcol),
     +                gk2(6,0:mcol), gr2(6,6,0:mcol),
     +                gt2(6,6,6,0:mcol)
      save /sectorn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/

c----------------------------------------------------------------------
      common /selecip/ onepit,ipstart,intip
      save /selecip/ 

      logical onepit
      integer intip
      character*4 ipstart
      character*200 templine
c-----------------------------------------------------------------------
      integer i, j, i1, i2, i3, ibeam, nd, nch, iostat
      integer lastnb
      character * 48 t_name
      character * 48 g_name

      integer junit
      double precision ss1, ss2

      character *(mcnam) dummy

      junit = 12
      

c-----------------------------------------------------------------------
c--- loop over two beams
      do ibeam = 1, 2
        if (ibeam .eq. 1)  then
          i1 = 0
          i2 = ninter
          i3 = 1
          t_name = 'train.manf'
        else
          i1 = ninter
          i2 = 0
          i3 = -1
          t_name = 'train.manb'
        endif
        open(iunit, file=t_name,status='OLD', iostat=iostat)
            if (iostat .ne. 0)  then
              print *, '<<< TRAIN >>> fatal: MAP file ',
     +        t_name(:lastnb(t_name)), ' not found'
              stop
            endif

        if (onepit) then
          i = i1
102       continue
          read(iunit, *, end = 1000)  templine
          if (templine(3:5) .eq. ipstart .or. 
     +    templine(1:2) .ne. 'MK') then
            read (iunit, *, end=1000) dummy,ss1,(sk1(j,i), j=1,6),
     +		(sr1(j,1,i), j=1,36),(st1(j,1,1,i), j=1,216)
            i = i + i3
          if (i .ne. i2) goto 102 
          endif
         

        else
          do i = i1, i2, i3
c            nd = msect
c            nch = mcnam
!            if (ibeam .eq. 1) g_name = 'train.mapf'
!            if (ibeam .eq. 2) g_name = 'train.mapb'
            
!          open(junit, file=g_name,status='OLD', iostat=iostat)
!          if (iostat .ne. 0)  then
!            print *, '<<< TRAIN >>> fatal: MAP file ',
!     +      g_name(:lastnb(g_name)), ' not found'
!            stop
!          endif
!         read(iunit, *, end = 1000) dummy

            if (ibeam .eq. 1)  then
!            read(junit, *, end = 100) dummy,ss1,(gk1(j,i), j=1,6),
!     +          (gr1(j,1,i), j = 1,36),(gt1(j,1,1,i), j=1,216)
!            write(91,*) dummy,ss1,(gk1(j,i), j=1,6),
!     +          (gr1(j,1,i), j = 1,36),(gt1(j,1,1,i), j=1,216)
           
             read (iunit, *, end=1000) dummy,ss1,(sk1(j,i), j=1,6),
     +		(sr1(j,1,i), j=1,36),(st1(j,1,1,i), j=1,216)
C            write(21,*) dummy,(sk1(j,i), j=1,6),
!     +		(sr1(j,1,i), j=1,36),(st1(j,1,1,i), j=1,216)
           

            else
!            read(junit, *, end = 100) dummy,ss2,(gk2(j,i), j=1,6),
!     +          (gr2(j,1,i), j = 1,36),(gt2(j,1,1,i), j=1,216)
!            write(92,*) dummy,ss2,(gk2(j,i), j=1,6),
!     +          (gr2(j,1,i), j = 1,36),(gt2(j,1,1,i), j=1,216)

             read (iunit, *, end=1000) dummy,ss2,(sk2(j,i), j=1,6),
     +       (sr2(j,1,i), j=1,36),(st2(j,1,1,i), j=1,216)

            endif

          enddo
        endif  
!       read(junit, *, end = 100) dummy
!       print *, '<<< TRAIN >>> fatal: MAP file ',
!    +  g_name(:lastnb(g_name)), ' too long'
!       stop
100     close(iunit)
101     close(junit)
      enddo
      return
1000  continue
      print *, '<<< TRAIN >>> fatal: MAP file',
     +t_name(:lastnb(t_name)), ' empty or too short'
      stop
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine read_opth(unit, nbeam, particle, seq_name, gev_l,
     +bcurr_l, xix, xiy, qx, qy, circum, delta, epsx0, epsy0,
     +deltap, ncolumn, c_names,mass,q)
      implicit none
      character *(*) particle, seq_name
      integer unit, ncolumn, i1, i2, get_tokens, nbeam, lastnb
      double precision gev_l, bcurr_l, xiy, xix, qx, qy, circum,
     +delta, epsx0, epsy0, deltap,mass,q
      character * 10000 line
      character * 1000 c_names (100)
      double precision pmass
      parameter         (pmass  = 0.938 272 31 d+00)
      
      mass = pmass
      q = 1.d0
  100 continue
      read(unit, '(a)', end=999) line
      if (line(1:1) .eq. '$')  then
         goto 999
      else if(line(1:1) .eq. '*')  then
         line(1:1) = ' '
         ncolumn = get_tokens(lastnb(line), line, ' *', c_names)

      else if (line(3:10) .eq. 'PARTICLE') then
        i1 = index(line, '"')
        i2 = i1 + index(line(i1+1:), '"')
        particle = line(i1+1:i2-1)
        if (particle .eq. 'ION') then
           read(unit, '(a)', end=999) line
	   if(line(3:6) .eq. 'MASS') then
              read(line(23:), *) mass
           endif
           read(unit, '(a)', end=999) line    
           if (line(3:8) .eq. 'CHARGE') then
              read(line(23:), *) q
           endif
        endif
      elseif (line(3:10) .eq. 'SEQUENCE') then
        i1 = index(line, '"')
        i2 = i1 + index(line(i1+1:), '"')
        seq_name = line(i1+1:i2-1)
      elseif (line(3:8) .eq. 'ENERGY') then
        if (nbeam .eq. 1) read(line(23:), *) gev_l
      elseif (line(3:10) .eq. 'BCURRENT') then
        if (nbeam .eq. 1) read(line(23:), *) bcurr_l
      elseif (line(3:4) .eq. 'Q1') then
        if (nbeam .eq. 1) read(line(23:), *) qx
      elseif (line(3:4) .eq. 'Q2') then
        if (nbeam .eq. 1) read(line(23:), *) qy
      elseif (line(3:5) .eq. 'DQ1') then
        if (nbeam .eq. 1) read(line(23:), *) xix
      elseif (line(3:5) .eq. 'DQ2') then
        if (nbeam .eq. 1) read(line(23:), *) xiy
      elseif (line(3:8) .eq. 'LENGTH') then
        if (nbeam .eq. 1) read(line(23:), *) circum
      elseif (line(3:8) .eq. 'DELTAP') then
        if (nbeam .eq. 1) read(line(23:), *) delta
      elseif (line(3:6) .eq. 'SIGE') then
        if (nbeam .eq. 1) read(line(23:), *) deltap
       elseif (line(3:4) .eq. 'EX') then
        if (nbeam .eq. 1) read(line(23:), *) epsx0
       elseif (line(3:4) .eq. 'EY') then
        if (nbeam .eq. 1) read(line(23:), *) epsy0
      endif
      goto 100
c-----------------------------------------------------------------------
  999 end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rdoptc(nbeam, nint)
      implicit none
      integer nbeam, nint, j, eflag, nampos, irow, it
      double precision betxt,betyt,dxt,dyt,qx,qy,st,
     +xix,xiy,alfxt,xt,alfyt,yt,muxt,muyt
c-----------------------------------------------------------------------
c read an optics table and fills the following arrays:
c
c  - name          name of element
c  - s             distance from start of sequence
c  - x             horizontal orbit offset in mm
c  - betx          horizontal beta-function
c  - alfx          horizontal alpha
c  - dx            horizontal dispersion
c  - y             vertical orbit offset in mm
c  - bety          vertical beta-function
c  - alfy          vertical alpha
c  - dy            vertical dispersion
c
c This table is fabricated in a MAD run with data which have extra
c markers BB in the neighbourhood of the pits at the places where
c separated beam-beam collisions occur for bunch trains.
c The program assumes that LHC is cycled such that all collision
c points close to a pit are also close to it in the optics file.
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy,mux,muy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/


      common /globtun/ qx,qy
      save/ globtun/
      
      common /selecip/ onepit,ipstart,intip
      save /selecip/ 



      logical onepit
      integer intip
      character*4 ipstart

c-----------------------------------------------------------------------
      character*(mcnam) namet
      integer mr, mc
c-- mc below limits the no. of columns to 20 - the OPTICS table has
c   49 possible columns from which the user can select a subset
      parameter (mr = 10000, mc = 20)
      integer i, k, ncolumn, lchar, 
     +n, lastnb, nc
      character *48 t_name
      character particle(2)*(mcnam)
      logical t1e, t2e, s1, s2
      integer headon, parasit, iostat, get_tokens, c_pos(100)
      double precision zz
      character * 10000 tmpline
      character * 1000 tokens(100), c_names (100)
      double precision mass(2),q(2)
c-----------------------------------------------------------------------
      lchar = mcnam
      n = 100
      zz = zero
      eflag = 0
      if (nbeam .eq. 1)  t_name = 'train.optf'
      if (nbeam .eq. 2)  t_name = 'train.optb'
      open(iunit, file=t_name, status='OLD', iostat=iostat)
      if (iostat .ne. 0)  then
        print *, '<<< TRAIN >>> fatal: OPTICS file ',
     +  t_name(:lastnb(t_name)), ' not found'
        stop
      endif
      call read_opth(iunit, nbeam, particle(nbeam), seq_name(nbeam),
     +gev, bcurr, xix, xiy, qx, qy, circum, delta, epsx0, epsy0,
     +deltap, ncolumn, c_names,mass(nbeam),q(nbeam))
      call get_cpos(ncolumn, c_names, c_pos)
      if (nbeam .eq. 1)  then
        if (all_write)  then
          write (mulist, 910) epsx0*1.d9, epsy0*1.d9, deltap*1.d3
          write (mulist, 950) xix, xiy, qx, qy, circum, delta
        endif
      endif
c--- get names, positions, betas etc.
      k = 0
      nint = 0
100   continue
      read (iunit, '(a)', end = 1000) tmpline
      if (tmpline(1:1) .ne. ' ')  goto 100
      nc = get_tokens(lastnb(tmpline), tmpline, ' "', tokens)
      if (nc .ne. ncolumn)  then
        print *, '<<< TRAIN >>> optics file column mismatch: ',
     +  ncolumn, nc
        stop
      endif
      irow = irow + 1
      namet = tokens(c_pos(1))
      it = 1
      read(tokens(c_pos(2)), *) st
      read(tokens(c_pos(3)), *) xt
      read(tokens(c_pos(4)), *) betxt
      read(tokens(c_pos(5)), *) alfxt
      read(tokens(c_pos(6)), *) dxt
      read(tokens(c_pos(7)), *) yt
      read(tokens(c_pos(8)), *) betyt
      read(tokens(c_pos(9)), *) alfyt
      read(tokens(c_pos(10)), *) dyt
      


      if (namet(1:2) .eq. 'IP')  then
        read (namet(3:3), '(i1)')  j
        namet(4:) = ' '
        if (nbeam .eq. 1)  then
          if (j .le. 8 .and. ippos(j) .lt. 0.d0) then
            if (k .eq. 8)  then
              print *, 'rdoptc fatal: more than 8 IPs found'
              stop
            endif
            k = k + 1
            ippos(j) = st
          endif
        endif
      endif
c keep only collisiont points
      if ((headon(namet) .gt. 0 .or. parasit(namet,nbeam) .gt. 0) 
     + .and. (.not. onepit .or. (onepit .and. namet(3:5) .eq.
     + ipstart))) then
        if (nint .ge. mcol) then
          print *, 'Error in "train" program - ', nint,
     +    ' collision points in optics table >= mcol = ', mcol
          stop
        endif
        
        if (nbeam .eq. 1)  then
          nint         = nint + 1
          name(nint,1) = namet
          s(nint,1)    = st
          x(nint,1)    = xt * ten3m
          alfx(nint,1) = alfxt
          betx(nint,1) = betxt
          dx(nint,1)   = dxt
          y(nint,1)    = yt * ten3m
          alfy(nint,1) = alfyt
          bety(nint,1) = betyt
          dy(nint,1)   = dyt
          occur(nint,1)= it
          
        else
          i = nampos(namet, name(1,1), ninter)
          if (i .eq. 0)  then
            print *, 'element ', namet, ' only in backward beam'
            eflag = eflag + 1
          else
            nint      = nint + 1
            name(i,2) = namet
            s(i,2)    = st
            x(i,2)    = xt * ten3m
            betx(i,2) = betxt
            alfx(i,2) = alfxt
            dx(i,2)   = dxt
            y(i,2)    = yt * ten3m
            bety(i,2) = betyt
            alfy(i,2) = alfyt
            dy(i,2)   = dyt
            occur(i,2)= it
            

          endif
        endif
      endif
      goto 100
1000  continue

c store particle types
      if (nbeam .eq. 2)  then
        t1e = particle(1) .eq. 'ELECTRON'
     +        .or. particle(1) .eq. 'POSITRON'
        t2e = particle(2) .eq. 'ELECTRON'
     +        .or. particle(2) .eq. 'POSITRON'
        s1 = particle(1) .eq. 'POSITRON'
     +        .or. particle(1) .eq. 'PROTON'
        s2 = particle(1) .eq. 'POSITRON'
     +        .or. particle(1) .eq. 'PROTON'
        if (t1e .and. .not. t2e .or. t2e .and. .not. t1e) then
          print *, 'Sorry - no different particle types yet'
          stop
        endif
        if (t1e)  then
          tmass = emass
          tradius = erad
	elseif (.not. t1e .and. .not. t2e .and. .not. s1 .and.
     + .not. s2) then
          if ((mass(1) .ne. mass(2)) .or. (abs(q(1)) .ne. 
     +    abs(q(2)))) then
             print *, 'Sorry - no different particle types yet'
             stop
          else
             tmass = mass(1)
             tradius = erad*mass(1)/emass*q(1)**2
          endif   
        else
          tmass = pmass
          tradius = prad
        endif
        if (s1 .and. s2 .or. .not. s1 .and. .not. s2) then
          if (.not. t1e .and. .not. t2e .and. .not. s1 .and.
     +    .not. s2) then
             if (q(1) .eq. q(2)) then
                xisign = 1.d0
             else
                xisign = -1.d0
             endif
          else
             xisign = 1.d0
          endif
        else
          if (.not. t1e .and. .not. t2e .and. .not. s1 .and.
     +    .not. s2) then
             if (q(1) .eq. q(2)) then
                xisign = 1.d0
             else
                xisign = -1.d0
             endif
          else
             xisign = -1.d0
          endif 
        endif
      endif
      if (all_write)
     +write (mulist, '(''particle '', i2, '': '',a)') nbeam,
     +particle(nbeam)
      if (all_write .and. nbeam .eq. 2) 
     +write(*,*) 'Particle 1:',particle(1),', m:',mass(1),
     +'GeV, q:',q(1),' e'
      if (all_write .and. nbeam .eq. 2)
     +write(*,*) 'Particle 2:',particle(2),', m:',mass(2),
     +'GeV, q:',q(2),' e'	 
     
      if (all_write .and. nbeam .eq. 1)
     +write (mulist, 920) gev, bcurr, sigb, iseed
      if (nbeam .eq. 2)  then
        do i = 1, ninter
          j = nampos(name(i,1), name(1,2), ninter)
          if (j .eq. 0)  then
            print *, 'element ', name(i,1), ' only in forward beam'
            eflag = eflag + 1
          endif
        enddo
      endif
      if (eflag .ne. 0)  then
        print *, 'exiting after above ', eflag, ' errors..'
        stop
      endif
      close(iunit)
      print *, 'ninter :', nint
      
c-----------------------------------------------------------------------
 910  format(' '/
     +     ' horizontal emittance            = ',f12.6,' nm'/
     +     ' vertical emittance              = ',f12.6,' nm'/
     +     ' relative momentum spread        = ',f12.6,' 1/1000')
 920  format(
     +     ' energy per particle             = ',f12.6,' GeV'/
     +     ' current per bunch               = ',f12.6,' A'/
     +     ' sigma of bunch current          = ',f12.6,' 1/1000'/
     +     ' seed for random generator       = ',i12)
 950  format(' '/
     +     ' horizontal chromaticity         = ',0p,f12.6/
     +     ' vertical chromaticity           = ',f12.6/
     +     ' horizontal tune                 = ',f12.6/
     +     ' vertical tune                   = ',f12.6/
     +     ' machine circumference           = ',f12.6,' m'/
     +     ' relative momentum error         = ',f12.6,' 1')
c-----------------------------------------------------------------------
999   end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rdsurvey(nbeam, nint)
      implicit none
      integer nbeam, nint, i, j, nampos, irow, nc
      double precision sx, sy, sz
c-----------------------------------------------------------------------
c read an survey file and fill surveyf - survey_x,y,z.
c selected elements must exactly match optics / twiss file
c  -- hence assumes b2 survey is backwards (as optics and maps)
c it also assumes that bb markers are named MKIP# and MKIP#P...
c-----------------------------------------------------------------------
 
      integer iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mpit, mpar, mdslt, mcol, mlocal, mcnam,melm
      integer mmaps
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c input/output units
      parameter (iunit=11)
c various array sizes etc.
      parameter (mcnam = 16)
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     + amp_bunch, amp_fac, b2_off, n_parasit
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     + iseed, amp_bunch, amp_fac, n_parasit, iact(8)
      save /globi/
      

      double precision survey_x, survey_y, survey_z,
     +   survey_x_ip, survey_y_ip, survey_z_ip,
     +   survey_sep_x, survey_sep_y 
      common /surveyf/
     +    survey_x(mcol,2), survey_y(mcol,2), survey_z(mcol,2),
     +    survey_x_ip(9,2), survey_y_ip(9,2), survey_z_ip(9,2),
     +    survey_sep_x(mcol), survey_sep_y(mcol)
      save /surveyf/

c-----------------------------------------------------------------------
      character*(mcnam) namet
      integer headon, parasit, iostat, get_tokens, c_pos(100), lastnb
      character *48 t_name
      character * 10000 tmpline
      character * 1000 tokens(100), c_names (100)
c   --- columns in the file: NAME, X, Y, Z
      integer ncolumn
      parameter (ncolumn=4) 
c-----------------------------------------------------------------------
      common /selecip/ onepit,ipstart,intip
      save /selecip/ 
      logical onepit
      integer intip
      character*4 ipstart
c-----------------------------------------------------------------------

      if (nbeam .eq. 1)  t_name = 'train.surf'
      if (nbeam .eq. 2)  t_name = 'train.surb'
      open(iunit, file=t_name, status='OLD', iostat=iostat)
      if (iostat .ne. 0)  then
        print *, '<<< TRAIN >>> fatal: SURVEY file ',
     +  t_name(:lastnb(t_name)), ' not found'
        stop
      endif

      nint = 0
      irow = 0
100   continue
      read (iunit, '(a)', end = 1000) tmpline
      if (tmpline(1:1) .ne. ' ')  goto 100
      nc = get_tokens(lastnb(tmpline), tmpline, ' "', tokens)
      if (nc .ne. 4)  then
        print *, '<<< TRAIN >>> survey file column mismatch: ',
     +       nc, '>>> NEEDS exactly NAME, X, Y, Z'
        stop
      endif
      irow = irow + 1
      namet = tokens(1)
      read(tokens(2), *) sx
      read(tokens(3), *) sy
      read(tokens(4), *) sz

      if (namet(1:2) .eq. 'IP')  then
        read (namet(3:3), '(i1)')  j
        namet(4:) = ' '
        if (j .le. 8) then
          survey_x_ip(j, nbeam) = sx
          survey_y_ip(j, nbeam) = sy
          survey_z_ip(j, nbeam) = sz
        endif
      endif
      if (headon(namet) .gt. 0 .or. parasit(namet,nbeam) .gt. 0) then
        if (nint .ge. mcol) then
          print *, 'Error in "train" program - ', nint,
     +    ' collision points in survey file >= mcol = ', mcol
          stop
        endif
        i = nampos(namet, name(1,1), ninter)
        if (i .eq. 0)  then
          print *, 'Error in "train" program: element ', namet,
     +             ' from SURVEY not in TWISS file'
          STOP
        else
          nint      = nint + 1
          survey_x(i,nbeam) = sx
          survey_y(i,nbeam) = sy
          survey_z(i,nbeam) = sz
        endif
      endif
      GO TO 100
1000  continue
      close(iunit)
      print *, 'ninter in SURVEY :', nint
999   end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine calcsurvey()
      implicit none
      integer i,j,nb,xsign
      double precision sx(2),sy(2),sz(2), angle
c-----------------------------------------------------------------------
      integer iunit, mulist
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mpit, mpar, mdslt, mcol, mlocal, mcnam,
     +        max_list,melm,mmaps
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c input/output units
      parameter (mulist=23)
c various array sizes etc.
      parameter (mcnam = 16)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c-----------------------------------------------------------------------
      double precision zero, one, two, pi
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter         (pi     = 3.1415926535898d0)
c-----------------------------------------------------------------------
c global counters
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact

      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     + iseed, amp_bunch, amp_fac, n_parasit, iact(8)
      save /globi/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     + amp_bunch, amp_fac, b2_off, n_parasit

c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name

      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     + eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/

      double precision survey_x, survey_y, survey_z,
     +   survey_x_ip, survey_y_ip, survey_z_ip,
     +   survey_sep_x, survey_sep_y 
      common /surveyf/
     +    survey_x(mcol,2), survey_y(mcol,2), survey_z(mcol,2),
     +    survey_x_ip(9,2), survey_y_ip(9,2), survey_z_ip(9,2),
     +    survey_sep_x(mcol), survey_sep_y(mcol)
      save /surveyf/
c-----------------------------------------------------------------------
      print *, 'BBLR with additional design orbit separation:'
      DO i = 1, ninter
         read (name(i,1)(5:5), '(i1)')  j
         do nb = 1,2
           sx(nb) = survey_x(i,nb) - survey_x_ip(j,nb)
           sy(nb) = survey_y(i,nb) - survey_y_ip(j,nb)
           sz(nb) = survey_z(i,nb) - survey_z_ip(j,nb)
         enddo
         xsign = 1
         angle = atan2(-sz(2)-sz(1), sx(2)-sx(1)) - 
     +             two * pi * s(i,1)/circum
100      if (angle > pi) angle = angle - two * pi
         if (angle < -pi) angle = angle + two * pi
         if (abs(angle) > pi) go to 100
         if (abs(angle) > pi/two) xsign = -1
         
         survey_sep_x(i) = sqrt((sx(2)-sx(1)) ** 2 +
     +          (-sz(2)-sz(1)) ** 2) * xsign
         survey_sep_y(i) = (sy(2)-sy(1))
         if (abs(survey_sep_x(i)) > 1e-7 .or. 
     +       abs(survey_sep_y(i)) > 1e-7) then
            print *, '   ', NAME(i,1), ' x=',
     +       survey_sep_x(i), 'y=', survey_sep_y(i)
         end if
      ENDDO
999   end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine set_cuem
      implicit none
      integer n
      double precision bc_fun
c-----------------------------------------------------------------------
c set beam current for all particles, emittances, set bunches active
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
      double precision bc1,bc2,epsxb1,epsxb2,epsyb1,epsyb2
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bcc/    bc1(mbunch), bc2(mbunch), epsxb1(mbunch), 
     +     epsxb2(mbunch), epsyb1(mbunch), epsyb2(mbunch)
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
c-----------------------------------------------------------------------
      real rg32cut
c-----------------------------------------------------------------------
      do n = 1, nbunch
        present(n,1) = 1
        present(n,2) = 1
      enddo
      do n = nbunch+1, mbunch
        present(n,1) = 0
        present(n,2) = 0
      enddo
        write(91,*) 'beamc_f: ',beamc_f
      print *,'Base beam current set to bcurr=',bcurr
      if (beamc_f .eq. 3)  then
	print *,'Using user defined beam currents!'
      endif
      do n = 1, nbunch
        if (beamc_f .eq. 1)  then
c two alternate beam currents
          if (mod(n,2) .eq. 0)  then
            bcurr1(n) = sigb * bcurr
            bcurr2(n) = sigb * bcurr
          else
            bcurr1(n) = bcurr
            bcurr2(n) = bcurr
          endif
        elseif (beamc_f .eq. 2)  then
c use function for beam current
          bcurr1(n) = bcurr * bc_fun(sigb)
          bcurr2(n) = bcurr * bc_fun(sigb)
        elseif (beamc_f .eq. 3)  then

c use function for beam current
!         write(*,*) 'try to read beam currents for ',nbunch, ' bunches'
          bcurr1(n) = bcurr*bc1(n)                                           
          bcurr2(n) = bcurr*bc2(n)                                         
        else
c set random numbers
!          print *,'set currents'
!          bcurr1(n) = max(0.d0, bcurr * (1.0 + sigb * rg32cut(3.)))
!          bcurr2(n) = max(0.d0, bcurr * (1.0 + sigb * rg32cut(3.)))
!          print *,'currents are set '
        endif
	if (debug .gt. 2) print *,'set to id',n,' current of ',bcurr1(n)
      enddo
      print *,'Base emittances set to epsx0=',epsx0,' to epsy0=', epsy0
      if (emitt_f .eq. 8)  then
	print *,'Using user defined emittances!'
      endif
      do n = 1, nbunch
        if (emitt_f .eq. 1)  then
c two alternate beam currents
          if (mod(n,2) .eq. 0)  then
            epsx(n,1) = sigem * epsx0
            epsx(n,2) = sigem * epsx0
            epsy(n,1) = sigem * epsy0
            epsy(n,2) = sigem * epsy0
          else
            epsx(n,1) = epsx0
            epsx(n,2) = epsx0
            epsy(n,1) = epsy0
            epsy(n,2) = epsy0
          endif
	
        elseif (emitt_f .eq. 2)  then
c use function for beam current
          epsx(n,1) = epsx0 * bc_fun(sigem)
          epsx(n,2) = epsx0 * bc_fun(sigem)
          epsy(n,1) = epsy0 * bc_fun(sigem)
          epsy(n,2) = epsy0 * bc_fun(sigem)
        elseif (emitt_f .eq. 3)  then
c first half gets factor
          if (n .le. nbunch/2)  then
            epsx(n,1) = sigem * epsx0
            epsx(n,2) = sigem * epsx0
            epsy(n,1) = sigem * epsy0
            epsy(n,2) = sigem * epsy0
          else
            epsx(n,1) = epsx0
            epsx(n,2) = epsx0
            epsy(n,1) = epsy0
            epsy(n,2) = epsy0
          endif
        elseif (emitt_f .eq. 8)  then
c use function from file
          epsx(n,1) = epsx0 * epsxb1(n)
          epsx(n,2) = epsx0 * epsxb2(n)
          epsy(n,1) = epsy0 * epsyb1(n)
          epsy(n,2) = epsy0 * epsyb2(n)
        else
c set random numbers
          epsx(n,1) = max(0.d0, epsx0 * (1.0 + sigem * rg32cut(3.)))
          epsx(n,2) = max(0.d0, epsx0 * (1.0 + sigem * rg32cut(3.)))
          epsy(n,1) = max(0.d0, epsy0 * (1.0 + sigem * rg32cut(3.)))
          epsy(n,2) = max(0.d0, epsy0 * (1.0 + sigem * rg32cut(3.)))
        endif
	if (debug .gt. 2) print *,'set emitt_x fo id',n,' to ', epsx(n,1)
      enddo
      do  n = 1,nbunch
          write(93,*) n,bcurr1(n),bcurr2(n),
     &                  epsx(n,1),epsx(n,2),epsy(n,1),epsy(n,2)
      enddo
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine track02(fsec, fdis)
      implicit none
      integer i,j
      double precision rt,orbit
c-----------------------------------------------------------------------
c fill initial guess for collision points
c-----------------------------------------------------------------------

      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/

      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical fsec, fdis
      dimension orbit(6), rt(6,6)
c-----------------------------------------------------------------------
      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/
      
c-----------------------------------------------------------------------
c initialise map for forward bunch
      call mxone(tr1(1,1,1), 6, 6)
      if (fsec) call dzero(tt1(1,1,1,1), 216)
      call rdcopy(z1(1,1,0), orbit, 6)
 
c track through sector 0 and catenate map
      call trmap(sk1(1,0), sr1(1,1,0), st1(1,1,1,0), orbit, rt)
      call mapcat(fsec, rt, st1(1,1,1,0),
     +     tr1(1,1,1), tt1(1,1,1,1), tr1(1,1,1), tt1(1,1,1,1))
 
      if (fdis) then
         call trdisp(rt, st1(1,1,1,0), d1(1,1,0), dd1(1,1,0),
     +        d1(1,1,1), dd1(1,1,1))
      endif

c for sectors 1 to nmaps do...
      do 40 i = 1, nmaps
         
c save collision position
         do 30 j = 1, nbunch
            call drcopy(orbit, z1(1,j,i), 6)
 30      continue
c track bunch through sector and catenate map
         call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i), orbit, rt)
         call mapcat(fsec, rt, st1(1,1,1,i),
     +        tr1(1,1,1), tt1(1,1,1,1), tr1(1,1,1), tt1(1,1,1,1))
 
         if (fdis) then
            call trdisp(rt, st1(1,1,1,i), d1(1,1,i), dd1(1,1,i),
     +           d1(1,1,i+1), dd1(1,1,i+1))
         endif
 40   continue
c save end point of orbit
      call drcopy(orbit, z1(1,1,nmaps+1), 6)  

c initialise map for backward bunch
      call mxone(tr2(1,1,1), 6, 6)
      if (fsec) call dzero(tt2(1,1,1,1), 216)
      call rdcopy(z2(1,1,nmaps+1), orbit, 6)
 
c for sectors nmaps to 1 do...  
      do 90 i = nmaps, 1, -1
 
c track through sector and catenate map
         call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i), orbit, rt)
         call mapcat(fsec, rt, st2(1,1,1,i),
     +        tr2(1,1,1), tt2(1,1,1,1), tr2(1,1,1), tt2(1,1,1,1))
 
         if (fdis) then
            call trdisp(rt, st2(1,1,1,i), d2(1,1,i+1), dd2(1,1,i+1),
     +           d2(1,1,i), dd2(1,1,i))
         endif
 
c save collision positions
         do 80 j = 1, nbunch
            call drcopy(orbit, z2(1,j,i), 6)
 80      continue
 90   continue
 
c track backward bunch through sector 0 and catenate map
      call trmap(sk2(1,0), sr2(1,1,0), st2(1,1,1,0), orbit, rt)
      call mapcat(fsec, rt, st2(1,1,1,0),
     +     tr2(1,1,1), tt2(1,1,1,1), tr2(1,1,1), tt2(1,1,1,1))
 
      if (fdis) then
         call trdisp(rt, st2(1,1,1,0), d2(1,1,1), dd2(1,1,1),
     +        d2(1,1,0), dd2(1,1,0))
      endif
c save end point of orbit
      call drcopy(orbit, z2(1,1,0), 6)
 
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine track0(fsec, fdis)
      implicit none
      integer i,j
      double precision rt,orbit
c-----------------------------------------------------------------------
c fill initial guess for collision points
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical fsec, fdis
      dimension orbit(6), rt(6,6)
c-----------------------------------------------------------------------
c initialise map for forward bunch
      call mxone(tr1(1,1,1), 6, 6)
      if (fsec) call dzero(tt1(1,1,1,1), 216)
      call rdcopy(z1(1,1,0), orbit, 6)
 
c track through sector 0 and catenate map
      call trmap(sk1(1,0), sr1(1,1,0), st1(1,1,1,0), orbit, rt)
      call mapcat(fsec, rt, st1(1,1,1,0),
     +     tr1(1,1,1), tt1(1,1,1,1), tr1(1,1,1), tt1(1,1,1,1))
 
      if (fdis) then
         call trdisp(rt, st1(1,1,1,0), d1(1,1,0), dd1(1,1,0),
     +        d1(1,1,1), dd1(1,1,1))
      endif
 
c for sectors 1 to ninter do...
      do 40 i = 1, ninter
 
c save collision position
         do 30 j = 1, nbunch
            call drcopy(orbit, z1(1,j,i), 6)
 30      continue
c track bunch through sector and catenate map
         call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i), orbit, rt)
         call mapcat(fsec, rt, st1(1,1,1,i),
     +        tr1(1,1,1), tt1(1,1,1,1), tr1(1,1,1), tt1(1,1,1,1))
 
         if (fdis) then
            call trdisp(rt, st1(1,1,1,i), d1(1,1,i), dd1(1,1,i),
     +           d1(1,1,i+1), dd1(1,1,i+1))
         endif
 40   continue
 
c save end point of orbit
      call drcopy(orbit, z1(1,1,ninter+1), 6)
 
c initialise map for backward bunch
      call mxone(tr2(1,1,1), 6, 6)
      if (fsec) call dzero(tt2(1,1,1,1), 216)
      call rdcopy(z2(1,1,ninter+1), orbit, 6)
 
c for sectors ninter to 1 do...
      do 90 i = ninter, 1, -1
 
c track through sector and catenate map
         call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i), orbit, rt)
         call mapcat(fsec, rt, st2(1,1,1,i),
     +        tr2(1,1,1), tt2(1,1,1,1), tr2(1,1,1), tt2(1,1,1,1))
 
         if (fdis) then
            call trdisp(rt, st2(1,1,1,i), d2(1,1,i+1), dd2(1,1,i+1),
     +           d2(1,1,i), dd2(1,1,i))
         endif
 
c save collision positions
         do 80 j = 1, nbunch
            call drcopy(orbit, z2(1,j,i), 6)
 80      continue
 90   continue
 
c track backward bunch through sector 0 and catenate map
      call trmap(sk2(1,0), sr2(1,1,0), st2(1,1,1,0), orbit, rt)
      call mapcat(fsec, rt, st2(1,1,1,0),
     +     tr2(1,1,1), tt2(1,1,1,1), tr2(1,1,1), tt2(1,1,1,1))
 
      if (fdis) then
         call trdisp(rt, st2(1,1,1,0), d2(1,1,1), dd2(1,1,1),
     +        d2(1,1,0), dd2(1,1,0))
      endif
 
c save end point of orbit
      call drcopy(orbit, z2(1,1,0), 6)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trackb2(fsec, fdis, bno)
      implicit none
      integer i,j,k, bno, low, up,l
      double precision ccp,rt,sx,sy,tt,xm,ym,orbit
      double precision sigx, sigy
c-----------------------------------------------------------------------
c track orbits for all bunches, or a selected bunch pair
c-----------------------------------------------------------------------

      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      
      double precision survey_x, survey_y, survey_z,
     +   survey_x_ip, survey_y_ip, survey_z_ip,
     +   survey_sep_x, survey_sep_y 
      common /surveyf/
     +    survey_x(mcol,2), survey_y(mcol,2), survey_z(mcol,2),
     +    survey_x_ip(9,2), survey_y_ip(9,2), survey_z_ip(9,2),
     +    survey_sep_x(mcol), survey_sep_y(mcol)
      save /surveyf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical fsec, fdis, pit_number
      dimension orbit(6), rt(6,6), tt(6,6,6)

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/ 
c-----------------------------------------------------------------------
c for all bunches, or bunch pair bno do ..
      if (bno .lt. 0)  then
        low = 1
        up = nbunch
      else
        low = bno
        up = bno
      endif
      do 90 j = low, up 
c initialise map for forward bunch j
         call mxone(tr1(1,1,j), 6, 6)
         if (fsec) call dzero(tt1(1,1,1,j), 216)
         call rdcopy(z1(1,j,0), orbit, 6)
 
c track through sector 0 and catenate map
         call trmap(sk1(1,0), sr1(1,1,0), st1(1,1,1,0), orbit, rt)
         call mapcat(fsec, rt, st1(1,1,1,0),
     +        tr1(1,1,j), tt1(1,1,1,j), tr1(1,1,j), tt1(1,1,1,j))
         if (fdis) then
            call trdisp(rt, st1(1,1,1,0), d1(1,j,0), dd1(1,j,0),
     +           d1(1,j,1), dd1(1,j,1))
         endif
 
c for all sectors 1 to nmaps do ...
         l = 0
         do 10 i = 1, nmaps
c save position before collision
            call drcopy(orbit, z1b(1,j,i), 6)
            if (mapmask(i) .ne. 0) then
               l = l + 1
               if(l .gt. ninter) then 
                  write(*,*) 'Error- Wrong correspondance between ', 
     +            'optic functions and maps.'
               endif

               k = ibnch1(j,l)
               
               if (k .ne. 0) then
                  ccp = xisign * bcurr2(k) / (ech * frev)
              
                  if (pit_number(l))  then
                     ccp = ccp * hofact
                  else
                     ccp = ccp * xifact
                  endif
                  sx = sqrt(sigx(j,l,1)**2 + sigx(k,l,2)**2)
                  sy = sqrt(sigy(j,l,1)**2 + sigy(k,l,2)**2)
                  xm = z2(1,k,i) + survey_sep_x(l)
                  ym = z2(3,k,i) + survey_sep_y(l)
                  call trbb(fsec, ccp, sx, sy, xm, ym, orbit, rt, tt)
                  call mapcat(fsec, rt, tt,
     +            tr1(1,1,j), tt1(1,1,1,j), tr1(1,1,j), tt1(1,1,1,j))
                  if (fdis) then
                     call trdisp(rt, tt, d1(1,j,i), dd1(1,j,i),
     +                 d1(1,j,i), dd1(1,j,i))
                  endif
                  
               endif  
               
               
            endif
c save position after collision
            call drcopy(orbit, z1a(1,j,i), 6) 
                             
c track bunch through sector and catenate map
            call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i), orbit, rt)
            call mapcat(fsec, rt, st1(1,1,1,i),
     +           tr1(1,1,j), tt1(1,1,1,j), tr1(1,1,j), tt1(1,1,1,j))
            if (fdis) then
               call trdisp(rt, st1(1,1,1,i), d1(1,j,i), dd1(1,j,i),
     +              d1(1,j,i+1), dd1(1,j,i+1))
            endif
            
 10      continue
c save end position of orbit
         call drcopy(orbit, z1(1,j,nmaps+1), 6)
 
c initialise map for backward bunch j
         call mxone(tr2(1,1,j), 6, 6)
         if (fsec) call dzero(tt2(1,1,1,j), 216)
         call rdcopy(z2(1,j,nmaps+1), orbit, 6)

c for all sectors ninter to 1 do...
         l = ninter + 1
         do 80 i = nmaps, 1, -1
 
c track through sector and catenate map
            call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i), orbit, rt)
            call mapcat(fsec, rt, st2(1,1,1,i),
     +           tr2(1,1,j), tt2(1,1,1,j), tr2(1,1,j), tt2(1,1,1,j))
 
            if (fdis) then
               call trdisp(rt, st2(1,1,1,i), d2(1,j,i+1), dd2(1,j,i+1),
     +              d2(1,j,i), dd2(1,j,i))
            endif
 
c save position before collision
            call drcopy(orbit, z2b(1,j,i), 6) 

c interaction for backward beam
            
            if (mapmask(i) .ne. 0) then
               l = l - 1 
               
               if(l .lt. 0) then 
                  write(*,*) 'Error- Wrong correspondance between ', 
     +            'optic functions and maps.'
               endif

c interaction for backward beam
               k = ibnch2(j,l)            
               if (k .ne. 0) then
                  ccp = xisign * bcurr1(k) / (ech * frev)
                  if (pit_number(l))  then
                    ccp = ccp * hofact
                  else
                    ccp = ccp * xifact
                  endif
                  sx = sqrt(sigx(j,l,1)**2 + sigx(k,l,2)**2)
                  sy = sqrt(sigy(j,l,1)**2 + sigy(k,l,2)**2)
                  xm = z1(1,k,i) - survey_sep_x(l)
                  ym = z1(3,k,i) - survey_sep_y(l)
                  call trbb(fsec, ccp, sx, sy, xm, ym, orbit, rt, tt)
                  call mapcat(fsec, rt, tt,
     +              tr2(1,1,j), tt2(1,1,1,j), tr2(1,1,j), tt2(1,1,1,j))
 
                  if (fdis) then
                     call trdisp(rt, tt, d2(1,j,i), dd2(1,j,i),
     +                 d2(1,j,i), dd2(1,j,i))
                  endif
               endif
               
            endif
c save position after collision
            call drcopy(orbit, z2a(1,j,i), 6)
 80      continue
c track backward bunch through sector 0 and catenate map
         call trmap(sk2(1,0), sr2(1,1,0), st2(1,1,1,0), orbit, rt)
         call mapcat(fsec, rt, st2(1,1,1,0),
     +        tr2(1,1,j), tt2(1,1,1,j), tr2(1,1,j), tt2(1,1,1,j))
 
         if (fdis) then
            call trdisp(rt, st2(1,1,1,0), d2(1,j,1), dd2(1,j,1),
     +           d2(1,j,0), dd2(1,j,0))
         endif
 
c save end position of orbit
         call drcopy(orbit, z2(1,j,0), 6)
         
 90   continue
     
c-----------------------------------------------------------------------
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trackb(fsec, fdis, bno)
      implicit none
      integer i,j,k, bno, low, up
      double precision ccp,rt,sx,sy,tt,xm,ym,orbit
      double precision sigx, sigy
c-----------------------------------------------------------------------
c track orbits for all bunches, or a selected bunch pair
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      
      double precision survey_x, survey_y, survey_z,
     +   survey_x_ip, survey_y_ip, survey_z_ip,
     +   survey_sep_x, survey_sep_y 
      common /surveyf/
     +    survey_x(mcol,2), survey_y(mcol,2), survey_z(mcol,2),
     +    survey_x_ip(9,2), survey_y_ip(9,2), survey_z_ip(9,2),
     +    survey_sep_x(mcol), survey_sep_y(mcol)
      save /surveyf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      integer ibcnt1,ibcnt2,ibnch1,ibnch2,maskm,maskmi,maskmp,maskmn,
     +        present
c-----------------------------------------------------------------------
c collision schedules
      common /sched/ ibcnt1, ibcnt2, maskm(mcol+1), present(mbunch,2),
     +maskmi(0:mdslt), maskmp(0:mdslt), maskmn(0:mdslt),
     +ibnch1(mbunch,mcol), ibnch2(mbunch,mcol)
      save /sched/
c--- maskm   for collision point i, maskm(i) is the slot number
c--- maskmi  for slot i, maskmi(i) is 0 or the number of the coll. point
c--- maskmp  for slot i, maskmp(i) = number of prev. or current
c            coll. point
c--- maskmn  for slot i, maskmn(i) = number of next or current
c            coll. point
c    present present(i,j) with i = bunch number, j = 1,2 (ring):
c            1 if (still) present, 0 if not
c    ibnch1  for bunch i of ring_1, ibnch1(i,j) is the ring_2 bunch
c            it collides with at collision point j
c    ibnch2  for bunch i of ring_2, ibnch2(i,j) is the ring_1 bunch
c            it collides with at collision point j
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical fsec, fdis, pit_number
      dimension orbit(6), rt(6,6), tt(6,6,6)
c-----------------------------------------------------------------------
c for all bunches, or bunch pair bno do ..
      if (bno .lt. 0)  then
        low = 1
        up = nbunch
      else
        low = bno
        up = bno
      endif
      do 90 j = low, up
 
c initialise map for forward bunch j
         call mxone(tr1(1,1,j), 6, 6)
         if (fsec) call dzero(tt1(1,1,1,j), 216)
         call rdcopy(z1(1,j,0), orbit, 6)
 
c track through sector 0 and catenate map
         call trmap(sk1(1,0), sr1(1,1,0), st1(1,1,1,0), orbit, rt)
         call mapcat(fsec, rt, st1(1,1,1,0),
     +        tr1(1,1,j), tt1(1,1,1,j), tr1(1,1,j), tt1(1,1,1,j))
 
         if (fdis) then
            call trdisp(rt, st1(1,1,1,0), d1(1,j,0), dd1(1,j,0),
     +           d1(1,j,1), dd1(1,j,1))
         endif
 
c for all sectors 1 to ninter do ...
         do 10 i = 1, ninter
 
c save position before collision
            call drcopy(orbit, z1b(1,j,i), 6)
 
c interaction
            k = ibnch1(j,i)
            if (k .ne. 0) then
               ccp = xisign * bcurr2(k) / (ech * frev)
              
               if (pit_number(i))  then
                  ccp = ccp * hofact
               else
                  ccp = ccp * xifact
               endif
               sx = sqrt(sigx(j,i,1)**2 + sigx(k,i,2)**2)
               sy = sqrt(sigy(j,i,1)**2 + sigy(k,i,2)**2)
               xm = z2(1,k,i) + survey_sep_x(i)
               ym = z2(3,k,i) + survey_sep_y(i)
               call trbb(fsec, ccp, sx, sy, xm, ym, orbit, rt, tt)
               call mapcat(fsec, rt, tt,
     +              tr1(1,1,j), tt1(1,1,1,j), tr1(1,1,j), tt1(1,1,1,j))
               if (fdis) then
                  call trdisp(rt, tt, d1(1,j,i), dd1(1,j,i),
     +                 d1(1,j,i), dd1(1,j,i))
               endif
            endif
 
c save position after collision
            call drcopy(orbit, z1a(1,j,i), 6)
 
c track bunch through sector and catenate map
            call trmap(sk1(1,i), sr1(1,1,i), st1(1,1,1,i), orbit, rt)
            call mapcat(fsec, rt, st1(1,1,1,i),
     +           tr1(1,1,j), tt1(1,1,1,j), tr1(1,1,j), tt1(1,1,1,j))
            if (fdis) then
               call trdisp(rt, st1(1,1,1,i), d1(1,j,i), dd1(1,j,i),
     +              d1(1,j,i+1), dd1(1,j,i+1))
            endif
 10      continue
 
c save end position of orbit
         call drcopy(orbit, z1(1,j,ninter+1), 6)
 
c initialise map for backward bunch j
         call mxone(tr2(1,1,j), 6, 6)
         if (fsec) call dzero(tt2(1,1,1,j), 216)
         call rdcopy(z2(1,j,ninter+1), orbit, 6)
 
c for all sectors ninter to 1 do...
         do 80 i = ninter, 1, -1
 
c track through sector and catenate map
            call trmap(sk2(1,i), sr2(1,1,i), st2(1,1,1,i), orbit, rt)
            call mapcat(fsec, rt, st2(1,1,1,i),
     +           tr2(1,1,j), tt2(1,1,1,j), tr2(1,1,j), tt2(1,1,1,j))
 
            if (fdis) then
               call trdisp(rt, st2(1,1,1,i), d2(1,j,i+1), dd2(1,j,i+1),
     +              d2(1,j,i), dd2(1,j,i))
            endif
 
c save position before collision
            call drcopy(orbit, z2b(1,j,i), 6)
 
c interaction for backward beam
            k = ibnch2(j,i)
            if (k .ne. 0) then
               ccp = xisign * bcurr1(k) / (ech * frev)
               if (pit_number(i))  then
                 ccp = ccp * hofact
               else
                 ccp = ccp * xifact
               endif
               sx = sqrt(sigx(j,i,1)**2 + sigx(k,i,2)**2)
               sy = sqrt(sigy(j,i,1)**2 + sigy(k,i,2)**2)
               xm = z1(1,k,i) - survey_sep_x(i)
               ym = z1(3,k,i) - survey_sep_y(i)
               call trbb(fsec, ccp, sx, sy, xm, ym, orbit, rt, tt)
               call mapcat(fsec, rt, tt,
     +              tr2(1,1,j), tt2(1,1,1,j), tr2(1,1,j), tt2(1,1,1,j))
 
               if (fdis) then
                  call trdisp(rt, tt, d2(1,j,i), dd2(1,j,i),
     +                 d2(1,j,i), dd2(1,j,i))
               endif
            endif
 
c save position after collision
            call drcopy(orbit, z2a(1,j,i), 6)
 80      continue
 
c track backward bunch through sector 0 and catenate map
         call trmap(sk2(1,0), sr2(1,1,0), st2(1,1,1,0), orbit, rt)
         call mapcat(fsec, rt, st2(1,1,1,0),
     +        tr2(1,1,j), tt2(1,1,1,j), tr2(1,1,j), tt2(1,1,1,j))
 
         if (fdis) then
            call trdisp(rt, st2(1,1,1,0), d2(1,j,1), dd2(1,j,1),
     +           d2(1,j,0), dd2(1,j,0))
         endif
 
c save end position of orbit
         call drcopy(orbit, z2(1,j,0), 6)
 90   continue
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trbb(fsec, ccp, sx, sy, xm, ym, orbit, re, te)
      implicit none
      double precision cbx,cby,ccp,crx,cry,exk,exkc,explim,fk,orbit,
     +phix,phixx,phixxx,phixxy,phixy,phixyy,phiy,phiyy,phiyyy,r,r2,re,
     +rho2,rho4,rho6,rk,sx,sx2,sy,sy2,te,tk,xb,xm,xr,xs,yb,ym,
     +yr,ys
c-----------------------------------------------------------------------
c transport map for beam-beam element
c-----------------------------------------------------------------------
 
      logical fsec
      dimension orbit(6), re(6,6), te(6,6,6)
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
c     if x > explim, exp(-x) is outside machine limits.
      parameter         (explim = 150.0d0)
      double precision newtol
      parameter (newtol = 1.0d-12)
c-----------------------------------------------------------------------
c initialize
      call mxone(re, 6, 6)
      call dzero(te, 216)
 
c factor for beam-beam kick
      fk = two * tradius * ccp / gamma
      if (fk .eq. 0.d0)  return
      sx2 = sx * sx
      sy2 = sy * sy
      xs  = orbit(1) - xm
      ys  = orbit(3) - ym
 
c limit formulas for sigma(x) = sigma(y)
c      if (sx2 .eq. sy2) then
      if (abs(sx2-sy2) .le. newtol) then
         rho2 = xs * xs + ys * ys
 
c limit case for xs = ys = 0
         if (rho2 .eq. 0.0) then
            re(2,1) = fk / (two * sx2)
            re(4,3) = fk / (two * sx2)
 
c general case
         else
            tk = rho2 / (two * sx2)
            if (tk .gt. explim) then
               exk  = 0.0
               exkc = one
               phix = xs * fk / rho2
               phiy = ys * fk / rho2
            else
               exk  = exp(-tk)
               exkc = one - exk
               phix = xs * fk / rho2 * exkc
               phiy = ys * fk / rho2 * exkc
            endif
 
c orbit kick
            
            orbit(2) = orbit(2) + phix
            orbit(4) = orbit(4) + phiy
 
c first-order effects
            rho4 = rho2 * rho2
            phixx = fk * (- exkc * (xs*xs - ys*ys) / rho4
     +           + exk * xs*xs / (rho2 * sx2))
            phixy = fk * (- exkc * two * xs * ys / rho4
     +           + exk * xs*ys / (rho2 * sx2))
            phiyy = fk * (+ exkc * (xs*xs - ys*ys) / rho4
     +           + exk * ys*ys / (rho2 * sx2))
            re(2,1) = phixx
            re(2,3) = phixy
            re(4,1) = phixy
            re(4,3) = phiyy
 
c second-order effects
            if (fsec) then
               rho6 = rho4 * rho2
               phixxx = fk*xs * (+ exkc * (xs*xs - three*ys*ys) / rho6
     +              - exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2)
     +              - exk * xs*xs / (two * rho2 * sx2**2))
               phixxy = fk*ys * (+ exkc * (three*xs*xs - ys*ys) / rho6
     +              - exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2)
     +              - exk * xs*xs / (two * rho2 * sx2**2))
               phixyy = fk*xs * (- exkc * (xs*xs - three*ys*ys) / rho6
     +              + exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2)
     +              - exk * ys*ys / (two * rho2 * sx2**2))
               phiyyy = fk*ys * (- exkc * (three*xs*xs - ys*ys) / rho6
     +              + exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2)
     +              - exk * ys*ys / (two * rho2 * sx2**2))
               te(2,1,1) = phixxx
               te(2,1,3) = phixxy
               te(2,3,1) = phixxy
               te(4,1,1) = phixxy
               te(2,3,3) = phixyy
               te(4,1,3) = phixyy
               te(4,3,1) = phixyy
               te(4,3,3) = phiyyy
            endif
         endif
 
c case sigma(x) > sigma(y)
      else
         r2 = two * (sx2 - sy2)
         if (sx2 .gt. sy2) then
            r  = sqrt(r2)
            rk = fk * sqrt(pi) / r
            xr = abs(xs) / r
            yr = abs(ys) / r
            call errf(xr, yr, crx, cry)
            tk = (xs * xs / sx2 + ys * ys / sy2) / two
            if (tk .gt. explim) then
               exk = 0.0
               cbx = 0.0
               cby = 0.0
            else
               exk = exp(-tk)
               xb  = (sy / sx) * xr
               yb  = (sx / sy) * yr
               call errf(xb, yb, cbx, cby)
            endif
 
c case sigma(x) < sigma(y)
         else
            r  = sqrt(-r2)
            rk = fk * sqrt(pi) / r
            xr = abs(xs) / r
            yr = abs(ys) / r
            call errf(yr, xr, cry, crx)
            tk = (xs * xs / sx2 + ys * ys / sy2) / two
            if (tk .gt. explim) then
               exk = 0.0
               cbx = 0.0
               cby = 0.0
            else
               exk = exp(-tk)
               xb  = (sy / sx) * xr
               yb  = (sx / sy) * yr
               call errf(yb, xb, cby, cbx)
            endif
         endif
         
c orbit kick
         phix = rk * (cry - exk * cby) * sign(one, xs)
         phiy = rk * (crx - exk * cbx) * sign(one, ys)
         orbit(2) = orbit(2) + phix
         orbit(4) = orbit(4) + phiy
 
c first-order effects
         phixx = (two / r2) * (- (xs * phix + ys * phiy)
     +        + fk * (one - (sy / sx) * exk))
         phixy = (two / r2) * (- (xs * phiy - ys * phix))
         phiyy = (two / r2) * (+ (xs * phix + ys * phiy)
     +        - fk * (one - (sx / sy) * exk))
         re(2,1) = phixx
         re(2,3) = phixy
         re(4,1) = phixy
         re(4,3) = phiyy
         
c second-order effects
         if (fsec) then
            phixxx = (- phix - (xs * phixx + ys * phixy)
     +           + fk * xs * sy * exk / sx**3) / r2
            phixxy = (- phiy - (xs * phixy - ys * phixx)) / r2
            phixyy = (+ phix - (xs * phiyy - ys * phixy)) / r2
            phiyyy = (+ phiy + (xs * phixy + ys * phiyy)
     +           - fk * ys * sx * exk / sy**3) / r2
            te(2,1,1) = phixxx
            te(2,1,3) = phixxy
            te(2,3,1) = phixxy
            te(4,1,1) = phixxy
            te(2,3,3) = phixyy
            te(4,1,3) = phixyy
            te(4,3,1) = phixyy
            te(4,3,3) = phiyyy
         endif
           
      endif
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trdisp(re, te, da, dda, db, ddb)
      implicit none
      integer i,j,k
      real da,dda,db,ddb
      double precision dc,ddc,re,te,temp
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
c track dispersion for one bunch through a sector
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      dimension re(6,6), te(6,6,6), da(6), dda(6), db(6), ddb(6)
      dimension dc(6), ddc(6)
c-----------------------------------------------------------------------
      do 30 i = 1, 6
         dc(i) = 0.0
         ddc(i) = 0.0
 
         do 20 k = 1, 6
            temp = 0.0
 
            do 10 j = 1, 6
               temp = temp + te(i,j,k)*da(j)
 10         continue
 
            dc(i)  = dc(i)  + re(i,k)*da(k)
            ddc(i) = ddc(i) + temp*da(k) + re(i,k)*dda(k)
 20      continue
 30   continue
      call drcopy(dc, db, 6)
      call drcopy(ddc, ddb, 6)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trmap(ek, re, te, orbit, rt)
      implicit none
      integer i,k,l
      double precision ek,orbit,re,rt,sum1,te,temp
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
c track orbit and change reference for re matrix into rt
c-----------------------------------------------------------------------
 
      dimension ek(6), re(6,6), te(6,6,6), orbit(6), rt(6,6)
      dimension temp(6)
c-----------------------------------------------------------------------
      do 30 i = 1, 6
         temp(i) = ek(i)
 
         do 20 k = 1, 6
            sum1 = 0.0
            do 10 l = 1, 6
               sum1 = sum1 + te(i,k,l) * orbit(l)
 10         continue
 
            temp(i) = temp(i) + (re(i,k) + sum1) * orbit(k)
            rt(i,k) = re(i,k) + sum1 + sum1
 20      continue
 30   continue
 
      do 40 i = 1, 6
         orbit(i) = temp(i)
 40   continue
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trstat(bunch1, bunch2)
      implicit none
c--- statistics information at observation point
      integer bunch1, bunch2
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      double precision zz1, zz2, tmp
c-----------------------------------------------------------------------
      lumicnt = lumicnt + 1
      zz1 = ztr(1,bunch1,1) - ztr(1,bunch2,2)
      zz2 = ztr(3,bunch1,1) - ztr(3,bunch2,2)
      tmp = exp(-(zz1**2 + zz2**2)*lumifact)
      lumiav = lumiav + tmp
      if (lumi_hist) write(lumilist, '(1p,e12.4)') tmp
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine wtrack(bunch1, bunch2, lp)
      implicit none
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer bunch1, bunch2, lp, jj, kp, kq
      double precision znt_1(4), znt_2(4), wcs(4)
      double precision sigx, sigy
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (outnorm .eq. 0)  then
        write(mtrack+lp,'(i5,1p,8e14.6)') c_turn,
     +  (ztr(jj,bunch1,1), jj = 1,4), (ztr(jj,bunch2,2), jj = 1,4)
      else
        do kq = 1, 3, 2
          kp = kq + 1
          znt_1(kq) =
     +      eiv1(2,kp,lp) * ztr(1,bunch1,1) -
     +      eiv1(1,kp,lp) * ztr(2,bunch1,1) +
     +      eiv1(4,kp,lp) * ztr(3,bunch1,1) -
     +      eiv1(3,kp,lp) * ztr(4,bunch1,1) +
     +      eiv1(6,kp,lp) * ztr(5,bunch1,1) -
     +      eiv1(5,kp,lp) * ztr(6,bunch1,1)
          znt_1(kp) =			
     +      eiv1(1,kq,lp) * ztr(2,bunch1,1) -
     +      eiv1(2,kq,lp) * ztr(1,bunch1,1) +
     +      eiv1(3,kq,lp) * ztr(4,bunch1,1) -
     +      eiv1(4,kq,lp) * ztr(3,bunch1,1) +
     +      eiv1(5,kq,lp) * ztr(6,bunch1,1) -
     +      eiv1(6,kq,lp) * ztr(5,bunch1,1)
          znt_2(kq) =			
     +      eiv2(2,kp,lp) * ztr(1,bunch2,2) -
     +      eiv2(1,kp,lp) * ztr(2,bunch2,2) +
     +      eiv2(4,kp,lp) * ztr(3,bunch2,2) -
     +      eiv2(3,kp,lp) * ztr(4,bunch2,2) +
     +      eiv2(6,kp,lp) * ztr(5,bunch2,2) -
     +      eiv2(5,kp,lp) * ztr(6,bunch2,2)
          znt_2(kp) =			
     +      eiv2(1,kq,lp) * ztr(2,bunch2,2) -
     +      eiv2(2,kq,lp) * ztr(1,bunch2,2) +
     +      eiv2(3,kq,lp) * ztr(4,bunch2,2) -
     +      eiv2(4,kq,lp) * ztr(3,bunch2,2) +
     +      eiv2(5,kq,lp) * ztr(6,bunch2,2) -
     +      eiv2(6,kq,lp) * ztr(5,bunch2,2)
        enddo
        znt_1(2) = (znt_1(2) * betx(outpos,1)
     +  + znt_1(1) * alfx(outpos,1))/sigx(bunch1,outpos,1)
        znt_1(4) = (znt_1(4) * bety(outpos,1)
     +  + znt_1(3) * alfy(outpos,1))/sigy(bunch1,outpos,1)
        znt_2(2) = (znt_2(2) * betx(outpos,2)
     +  + znt_2(1) * alfx(outpos,2))/sigx(bunch2,outpos,2)
        znt_2(4) = (znt_2(4) * bety(outpos,2)
     +  + znt_2(3) * alfy(outpos,2))/sigy(bunch2,outpos,2)
        znt_1(1) = znt_1(1) / sigx(bunch1,outpos,1)
        znt_1(3) = znt_1(3) / sigy(bunch1,outpos,1)
        znt_2(1) = znt_2(1) / sigx(bunch2,outpos,2)
        znt_2(3) = znt_2(3) / sigy(bunch2,outpos,2)
        wcs(1) = znt_1(1)**2 + znt_1(2)**2
        wcs(2) = znt_1(3)**2 + znt_1(4)**2
        wcs(3) = znt_2(1)**2 + znt_2(2)**2
        wcs(4) = znt_2(3)**2 + znt_2(4)**2
        write(mtrack+lp,'(i5,1p,12e14.6)') c_turn, znt_1, znt_2, wcs
      endif
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ddump(text, array, n)
      implicit none
      character *(*) text
      double precision array(*)
      integer i, n
      print '(a)', text
      print '(1p,6e14.6)', (array(i), i = 1, n)
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine errf(xx, yy, wx, wy)
      implicit none
      integer n,nc,nu
      double precision cc,h,one,q,rx,ry,saux,sx,sy,tn,two,tx,ty,wx,wy,x,
     +xh,xl,xlim,xx,y,yh,ylim,yy
c-----------------------------------------------------------------------
c   modification of wwerf, double precision complex error function,
c   written at CERN by K. Koelbig.
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      parameter (cc     = 1.12837 91670 9551d0)
      parameter (one    = 1.0d0)
      parameter (two    = 2.0d0)
      parameter (xlim   = 5.33d0)
      parameter (ylim   = 4.29d0)
      dimension rx(33), ry(33)
c-----------------------------------------------------------------------
      x = abs(xx)
      y = abs(yy)
 
      if (y .lt. ylim  .and.  x .lt. xlim) then
        q  = (one - y / ylim) * sqrt(one - (x/xlim)**2)
        h  = one / (3.2d0 * q)
        nc = 7 + int(23.0*q)
        xl = h**(1 - nc)
        xh = y + 0.5d0/h
        yh = x
        nu = 10 + int(21.0*q)
        rx(nu+1) = 0.
        ry(nu+1) = 0.
 
        do 10 n = nu, 1, -1
          tx = xh + n * rx(n+1)
          ty = yh - n * ry(n+1)
          tn = tx*tx + ty*ty
          rx(n) = 0.5d0 * tx / tn
          ry(n) = 0.5d0 * ty / tn
   10   continue
 
        sx = 0.
        sy = 0.
 
        do 20 n = nc, 1, -1
          saux = sx + xl
          sx = rx(n) * saux - ry(n) * sy
          sy = rx(n) * sy + ry(n) * saux
          xl = h * xl
   20   continue
 
        wx = cc * sx
        wy = cc * sy
      else
        xh = y
        yh = x
        rx(1) = 0.
        ry(1) = 0.
 
        do 30 n = 9, 1, -1
          tx = xh + n * rx(1)
          ty = yh - n * ry(1)
          tn = tx*tx + ty*ty
          rx(1) = 0.5d0 * tx / tn
          ry(1) = 0.5d0 * ty / tn
   30   continue
 
        wx = cc * rx(1)
        wy = cc * ry(1)
      endif
 
      if (y .eq. 0.) wx = exp(-x**2)
      if (yy .lt. 0.) then
        wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx
        wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy
        if (xx .gt. 0.) wy = -wy
      else
        if (xx .lt. 0.) wy = -wy
      endif
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine solver(augmat, ndim, mdim, irank)
      implicit none
      integer ic,ip,ir,irank,it,mdim,nc,ndim,nr
      double precision augmat,h,pivot
c-----------------------------------------------------------------------
c solve the linear equation  a * x = b
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      dimension augmat(ndim,ndim+mdim)
c-----------------------------------------------------------------------
      nr = ndim
      nc = ndim + mdim
 
      do 100 it = 1, nr
        pivot = 0.
        ip = 0
        do 10 ir = it, nr
          if (abs(augmat(ir,it)) .ge. abs(pivot)) then
            pivot = augmat(ir,it)
            ip = ir
          endif
   10   continue
 
!          write(*,*) augmat
!          write(*,*) ndim,mdim,irank
        if (pivot .eq. 0.0) then
c           write(*,*) augmat
c           write(*,*) ndim,mdim,irank
           print 910,  irank, ndim
           return
        endif
 
        irank = it
 
        do 30 ic = 1, nc
          augmat(ip,ic) = augmat(ip,ic) / pivot
   30   continue
 
        if (ip .ne. it) then
          do 50 ic = 1, nc
            h = augmat(ip,ic)
            augmat(ip,ic) = augmat(it,ic)
            augmat(it,ic) = h
   50     continue
        endif
 
        do 70 ir = 1, nr
          if (ir .ne. it) then
            h = augmat(ir,it)
            do 60 ic = 1, nc
              augmat(ir,ic) = augmat(ir,ic) - h * augmat(it,ic)
   60       continue
          endif
   70   continue
  100 continue
 
      irank = ndim
c-----------------------------------------------------------------------
 910  format(' '/' error: iteration matrix has rank ',i20,
     +     ' and dimension ',i5)
c-----------------------------------------------------------------------
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function nampos(sname,slist,nlist)
c-----------------------------------------------------------------------
c
c   gives position in name list or 0
c
c   input
c   SNAME                   name to be looked up
c   SLIST                   name list
c   NLIST                   no. of names in SLIST
c-----------------------------------------------------------------------
      implicit none
      character *(*) sname,slist(*)
      integer nlist, i
      nampos = 0
      do i = 1, nlist
        if (sname .eq. slist(i)) then
          nampos = i
          return
        endif
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function lastnb(string)
      implicit none
c-----------------------------------------------------------------------
c
c--- find alst non-blank in string
c--- input
c    STRING
c--- output
c    function value = last non-blank (at least 1)
c-----------------------------------------------------------------------
      integer i,lastnb
      character *(*) string
      do i = len(string), 1, -1
        if (string(i:i) .ne. ' ')  goto 20
      enddo
      i = 1
   20 lastnb = i
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine upper(sl)
c***********************************************************************
c
c   Purpose: Converts all characters in SL into upper case.
c
c--- Input/Output
c    SL           string to be modified
c
c   Author: H. Grote / CERN                        date: Dec. 19, 1991
c                                              last mod: Dec. 19, 1991
c
c***********************************************************************
      implicit none
      integer i, k
      character * (*) sl
      character alfbet(2) * 26
      save alfbet
      data alfbet(1) / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data alfbet(2) / 'abcdefghijklmnopqrstuvwxyz'/
 
      do i=1,len(sl)
        k = index(alfbet(2),sl(i:i))
        if (k .ne. 0)  sl(i:i) = alfbet(1)(k:k)
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_hits(mbuck, mcol, part_f,part_b,colcnt_f, colcnt_b,
     +hitlist_f, hitlist_b, equl, nequl, lequl)
      implicit none
      integer mbuck, mcol, hitlist_f(0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +colcnt_f(0:mbuck-1), colcnt_b(0:mbuck-1),
     +equl, nequl(*), lequl(mbuck,*)
      integer i, j, k
c--- fill hitlist with equ. class bunches from beam 1
      do i = 1, equl
        do j = 1, nequl(equl)
          hitlist_f(lequl(j,i)) = 1
        enddo
      enddo
c--- loop until all bunches in 1 + 2 assembled
 10   k = 0
      do i = 0, mbuck-1
        if (hitlist_f(i) .ne. 0)  then
          do j = 1, colcnt_f(i)
            if (hitlist_b(part_f(j,i)) .eq. 0) then
              k = k + 1
              hitlist_b(part_f(j,i)) = 1
            endif
          enddo
        endif
      enddo
      do i = 0, mbuck-1
        if (hitlist_b(i) .ne. 0)  then
          do j = 1, colcnt_b(i)
            if (hitlist_f(part_b(j,i)) .eq. 0) then
              k = k + 1
              hitlist_f(part_b(j,i)) = 1
            endif
          enddo
        endif
      enddo
      if (k .ne. 0)  goto 10
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_cpos(n, names, pos)
      implicit none
      integer mnames
      parameter ( mnames = 10)
      integer n, i, j, pos(*)
      character * 1000 names(*), cname(mnames)
      data cname / 'NAME', 'S', 'X', 'BETX',
     +            'ALFX', 'DX', 'Y', 'BETY', 'ALFY',
     +            'DY' /
      if (n .lt. mnames)  then
        print *, '<<< TRAIN >>> fatal: # columns in optics file = ',
     +  n, ' lower than minimum = ', mnames
        stop
      endif
      do 10 j = 1, n
        do i = 1, mnames
          if (cname(j) .eq. names(i))  then
            pos(j) = i
            goto 10
          endif
        enddo
        print *, '<<< TRAIN >>> fatal: column ', cname(j),
     +  ' not in optics file'
        stop
   10 continue
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function get_tokens(length, string, skip, tokens)
      implicit none
      character *(*) string, skip
      character * 1000 tokens(*)
      integer length, n, i, k
      n = 0
      k = 0
      do i = 1, length
        if (index(skip, string(i:i)) .eq. 0)  then
          if (k .eq. 0)  then
             n = n + 1
             tokens(n) = ' '
          endif
          k = k + 1
          tokens(n)(k:k) = string(i:i)
        else
          k = 0
        endif
      enddo
      get_tokens = n
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine equ_class(mbuck, mcol, a, colcnt, list, low, high,
     +equl, cequl, nequl, lequl, ordl, cordl, nordl, lordl)
      implicit none
      integer mbuck, mcol, colcnt(0:mbuck-1), list(mcol,0:mbuck-1),
     +a(2,0:mbuck-1),
     +equl, nequl(*), cequl(*), lequl(mbuck,*),
     +ordl, nordl(*), cordl(*), lordl(mbuck,*)
      integer i, j, k, low, high
      logical in_list
      equl = 0
c--- get lowest + highest class
      low = 100000
      high = 0
      do i = 0, mbuck-1
        if ((a(1,i) .ne. 0) .or. (a(2,i) .ne. 0))  then
          if (colcnt(i) .gt. high)  high = colcnt(i)
          if (colcnt(i) .lt. low)   low  = colcnt(i)
        endif
      enddo
c--- list classes according to coll. counts
      ordl = 0
      do i = low, high
        k = 0
        do j = 0, mbuck-1
          if (colcnt(j) .eq. i)  then
            k = k + 1
            lordl(k,ordl+1) = j
          endif
        enddo
        if (k .gt. 0)  then
          ordl = ordl + 1
          nordl(ordl) = k
          cordl(ordl) = i
        endif
      enddo
c--- find equivalence classes
      equl = ordl
      do i = 1, ordl
        nequl(i) = 0
        cequl(i) = cordl(i)
        do j = 1, nordl(i)
          k = lordl(j,i)
          if (colcnt(k) .ne. cequl(i))  then
            print *, 'fatal: colcnt(k), cequl(i) = ',
     +      colcnt(k), cequl(i)
            stop
          endif
          if (.not.in_list(mcol, k, colcnt(k), list,
     +    nequl(i), lequl(1,i))) then
            nequl(i) = nequl(i) + 1
            lequl(nequl(i),i) = k
          endif
        enddo
      enddo
      end
      logical function in_list(mcol, n, cnt, list, neq, leq)
      implicit none
      integer mcol, n, cnt, list(mcol,0:mcol), neq, leq(*)
      integer i, j, k
      do 1  i = 1, neq
        k = leq(i)
        do j = 1, cnt
          if (list(j,n) .ne. list(j,k))  goto 1
        enddo
        in_list = .true.
        return
 1    continue
      in_list = .false.
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mysort(n, l, p)
      implicit none
      integer n, j, k, f, l(*), p(*)
 10   continue
      f = 0
      do j = 1, n-1
        if (l(j) .gt. l(j+1))  then
          k = l(j)
          l(j) = l(j+1)
          l(j+1) = k
          k = p(j)
          p(j) = p(j+1)
          p(j+1) = k
          f = 1
        endif
      enddo
      if (f .ne. 0)  goto 10
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine getmask2(mdslt, npar, circ, iact, ip, mask, maski,
     +count)
      implicit none
      integer mdslt, npar(8), mask(0:mdslt-1), maski(0:mdslt-1),
     + ipos(8)
      
      integer i, j, k, count, keep, iact(8),rotate
      double precision circ, hstep, ip(8)
      common /selnpar/ ipos
      save /selnpar/

      hstep = circ / mdslt

      do i = 1, 8
        ipos(i) = ip(i) / hstep + 0.1
      enddo

      keep = rotate(ipos)

      do i = 1, 8
        ipos(i) = ipos(i) - keep
        if (ipos(i) .lt. 0) ipos(i) = ipos(i) + mdslt
      enddo

      do i = 0, mdslt-1
        mask(i) = 0
        maski(i) = 0
      enddo
      do i = 1, 8

        if (iact(i) .ne. 0)  then
          mask(ipos(i)) = 2
          if (npar(i) .ne. 0) then
            do j = 1, npar(i)
              k = ipos(i) - j
              if (k .lt. 0)  k = k + mdslt
              mask(k) = 1
              k = ipos(i) + j
              if (k .ge. mdslt)  k = k - mdslt
              mask(k) = 1
            enddo
          endif  
        endif
      enddo
      count = 0
      do i = 0, mdslt-1
        if (mask(i) .gt. 0)  then
          maski(i) = count
          count = count + 1
        endif
      enddo

      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine getmask(mdslt, npar, circ, iact, ip, mask, maski,
     +count)
      implicit none
      integer mdslt, npar(8), mask(0:mdslt-1), maski(0:mdslt-1),
     + ipos(8)
      
      integer i, j, k, count, keep, iact(8)
      double precision circ, hstep, ip(8)
c-----------------------------------------------------------------------
      common /selecip/ onepit,ipstart,intip
      save /selecip/ 
      common /selnpar/ ipos
      save /selnpar/
      logical onepit
      integer intip
      character*4 ipstart
      integer mbuck
      mbuck =3564
c-----------------------------------------------------------------------

      hstep = circ / mdslt
      do i = 1, 8
        ipos(i) = ip(i) / hstep + 0.1
      enddo

      keep = ipos(5) - npar(5)
      do i = 1, 8
        ipos(i) = ipos(i) - keep
        if (ipos(i) .lt. 0) ipos(i) = ipos(i) + mdslt
      enddo
      do i = 0, mdslt-1
        mask(i) = 0
        maski(i) = 0
      enddo
      do i = 1, 8
        if (iact(i) .ne. 0)  then
          mask(ipos(i)) = 2
          if (npar(i) .ne. 0) then
            do j = 1, npar(i)
              k = ipos(i) - j
              if (k .lt. 0)  k = k + mdslt
              mask(k) = 1
              k = ipos(i) + j
              if (k .ge. mdslt)  k = k - mdslt
              mask(k) = 1
            enddo
          endif  
        endif
      enddo
      count = 0
      do i = 0, mdslt-1
        if (mask(i) .gt. 0)  then
          maski(i) = count
          count = count + 1
        endif
      enddo
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hqr2(ndim, n, ilow, iupp, h, wr, wi, vecs, ierr)
      implicit none
      integer i,ien,ierr,ilow,its,iupp,j,k,l,m,n,na,ndim
      double precision den,epsmch,h,hnorm,p,q,r,ra,s,sa,t,temp,tempi,
     +tempr,vecs,vi,vr,w,wi,wr,x,y,z
c----------------------------------------------------------------------*
c purpose:                                                             *
c   finds eigenvalues and eigenvectors of an unsymmetric real matrix,  *
c   a which has been reduced to upper hessenberg form, h, by the       *
c   subroutine orthes. the orthogonal transformations must be placed   *
c   in the array vecs by subroutine ortran.                            *
c                                                                      *
c   translation of the algol procedure hqr2 in:                        *
c   handbook series linear algebra,                                    *
c   num. math. 16, 181 - 204 (1970) by g. peters and j. h. wilkinson.  *
c input:                                                               *
c   n         (integer) order of the hessenberg matrix h.              *
c   ilow,iupp (integer)                                                *
c   h(ndim,n) (real)    the hessenberg matrix produced by orthes.      *
c   vecs(ndim,n) (real) a square matrix of order n containing the      *
c                       similarity transformation from a to h          *
c output:                                                              *
c   h(ndim,n) (real)    modified.                                      *
c   wr(n)     (real)    real parts of eigenvalues of h (or a).         *
c   wi(n)     (real)    imaginary parts of eigenvalues of h (or a).    *
c   vecs(ndim,n) (real) the unnormalized eigenvectors of a.            *
c                       complex vectors are stored as pairs of reals.  *
c----------------------------------------------------------------------*
 
      dimension h(ndim,n), wr(n), wi(n), vecs(ndim,n)
      parameter (epsmch = 1.0e-12)
 
      ierr = 0
 
c---- store isolated roots.
      do 10 i = 1, n
        if (i .lt. ilow  .or.  i .gt. iupp) then
          wr(i) = h(i,i)
          wi(i) = 0.0
        endif
   10 continue
 
      ien = iupp
      t = 0.0
 
c---- next eigenvalue.
   60 if (ien .ge. ilow) then
        its = 0
        na = ien - 1
 
c---- next iteration; look for single small sub-diagonal element.
   70   continue
          do 80 l = ien, ilow + 1, -1
            if (abs(h(l,l-1)) .le.
     +          epsmch * (abs(h(l-1,l-1)) + abs(h(l,l)))) go to 100
   80     continue
          l = ilow
  100     continue
          x = h(ien,ien)
          if (l .eq. ien) go to 270
          y = h(na,na)
          w = h(ien,na) * h(na,ien)
          if (l .eq. na) go to 280
          if (its .eq. 30) then
            ierr = ien
            go to 9999
          endif
 
c---- form exceptional shift.
          if (its .eq. 10  .or.  its .eq. 20) then
            t = t + x
            do 120 i = ilow, ien
              h(i,i) = h(i,i) - x
  120       continue
            s = abs(h(ien,na)) + abs(h(na,ien-2))
            x = 0.75 * s
            y = x
            w = - 0.4375 * s * s
          endif
          its = its + 1
 
c---- look for two consecutive small sub-diagonal elements.
          do 140 m = ien - 2, l, - 1
            z = h(m,m)
            r = x - z
            s = y - z
            p = (r * s - w) / h(m+1,m) + h(m,m+1)
            q = h(m+1,m+1) - z - r - s
            r = h(m+2,m+1)
            s = abs(p) + abs(q) + abs(r)
            p = p / s
            q = q / s
            r = r / s
            if (m .eq. l) go to 150
            if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. epsmch * abs(p)
     x       * (abs(h(m-1,m-1)) + abs(z) + abs(h(m+1,m+1)))) go to 150
  140     continue
 
  150     continue
          h(m+2,m) = 0.0
          do 160 i = m + 3, ien
            h(i,i-2) = 0.0
            h(i,i-3) = 0.0
  160     continue
 
c---- double qr step involving rows l to ien and columns m to ien.
          do 260 k = m, na
            if (k .ne. m) then
              p = h(k,k-1)
              q = h(k+1,k-1)
              if (k .ne. na) then
                r = h(k+2,k-1)
              else
                r = 0.0
              endif
              x = abs(p) + abs(q) + abs(r)
              if (x .eq. 0.0) go to 260
              p = p / x
              q = q / x
              r = r / x
            endif
            s = sign(sqrt(p**2+q**2+r**2),p)
            if (k .ne. m) then
              h(k,k-1) = - s * x
            else if (l .ne. m) then
              h(k,k-1) = - h(k,k-1)
            endif
            p = p + s
            x = p / s
            y = q / s
            z = r / s
            q = q / p
            r = r / p
 
c---- row modification.
            do 210 j = k, n
              p = h(k,j) + q * h(k+1,j)
              if (k .ne. na) then
                p = p + r * h(k+2,j)
                h(k+2,j) = h(k+2,j) - p * z
              endif
              h(k+1,j) = h(k+1,j) - p * y
              h(k,j) = h(k,j) - p * x
  210       continue
 
c---- column modification.
            j = min(ien,k+3)
            do 230 i = 1, j
              p = x * h(i,k) + y * h(i,k+1)
              if (k .ne. na) then
                p = p + z * h(i,k+2)
                h(i,k+2) = h(i,k+2) - p * r
              endif
              h(i,k+1) = h(i,k+1) - p * q
              h(i,k) = h(i,k) - p
  230       continue
 
c---- accumulate transformations.
            do 250 i = ilow, iupp
              p = x * vecs(i,k) + y * vecs(i,k+1)
              if (k .ne. na) then
                p = p + z * vecs(i,k+2)
                vecs(i,k+2) = vecs(i,k+2) - p * r
              endif
              vecs(i,k+1) = vecs(i,k+1) - p * q
              vecs(i,k) = vecs(i,k) - p
  250       continue
  260     continue
 
c---- go to next iteration.
        go to 70
 
c==== one real root found.
  270   h(ien,ien) = x + t
        wr(ien) = h(ien,ien)
        wi(ien) = 0.0
        ien = na
        go to 60
 
c==== two roots (real pair or complex conjugate) found.
  280   p = (y - x) / 2.0
        q = p**2 + w
        z = sqrt(abs(q))
        x = x + t
        h(ien,ien) = x
        h(na,na) = y + t
 
c---- real pair.
        if (q .gt. 0.0) then
          z = p + sign(z,p)
          wr(na) = x + z
          wr(ien) = x - w / z
          wi(na) = 0.0
          wi(ien) = 0.0
          x = h(ien,na)
          r = sqrt(x**2+z**2)
          p = x / r
          q = z / r
 
c---- row modification.
          do 290 j = na, n
            z = h(na,j)
            h(na,j) = q * z + p * h(ien,j)
            h(ien,j) = q * h(ien,j) - p * z
  290     continue
 
c---- column modification.
          do 300 i = 1, ien
            z = h(i,na)
            h(i,na) = q * z + p * h(i,ien)
            h(i,ien) = q * h(i,ien) - p * z
  300     continue
 
c---- accumulate transformations.
          do 310 i = ilow, iupp
            z = vecs(i,na)
            vecs(i,na) = q * z + p * vecs(i,ien)
            vecs(i,ien) = q * vecs(i,ien) - p * z
  310     continue
 
c---- complex pair.
        else
          wr(na) = x + p
          wr(ien) = x + p
          wi(na) = z
          wi(ien) = -z
        endif
 
c----- go to next root.
        ien = ien - 2
        go to 60
      endif
 
c==== compute matrix norm.
      hnorm = 0.0
      k = 1
      do 520 i = 1, n
        do 510 j = k, n
          hnorm = hnorm + abs(h(i,j))
  510   continue
        k = i
  520 continue
 
c==== back substitution.
      do 690 ien = n, 1, -1
        p = wr(ien)
        q = wi(ien)
        na = ien - 1
 
c---- real vector.
        if (q .eq. 0.0) then
          m = ien
          h(ien,ien) = 1.0
          do 640 i = na, 1, -1
            w = h(i,i) - p
            r = h(i,ien)
            do 610 j = m, na
              r = r + h(i,j) * h(j,ien)
  610       continue
            if (wi(i) .lt. 0.0) then
              z = w
              s = r
            else
              m = i
              if (wi(i) .eq. 0.0) then
                temp = w
                if (w .eq. 0.0) temp = epsmch * hnorm
                h(i,ien) = - r / temp
              else
                x = h(i,i+1)
                y = h(i+1,i)
                q = (wr(i) - p)**2 + wi(i)**2
                t = (x * s - z * r) / q
                h(i,ien) = t
                if (abs(x) .gt. abs(z)) then
                  h(i+1,ien) = - (r + w * t) / x
                else
                  h(i+1,ien) = - (s + y * t) / z
                endif
              endif
            endif
  640     continue
 
c---- complex vector associated with lamda = p - i * q.
        else if (q .lt. 0.0) then
          m = na
          if (abs(h(ien,na)) .gt. abs(h(na,ien))) then
            h(na,na) = - (h(ien,ien) - p) / h(ien,na)
            h(na,ien) = - q / h(ien,na)
          else
            den = (h(na,na) - p)**2 + q**2
            h(na,na) = - h(na,ien) * (h(na,na) - p) / den
            h(na,ien) = h(na,ien) * q / den
          endif
          h(ien,na) = 1.0
          h(ien,ien) = 0.0
          do 680 i = ien - 2, 1, - 1
            w = h(i,i) - p
            ra = h(i,ien)
            sa = 0.0
            do 660 j = m, na
              ra = ra + h(i,j) * h(j,na)
              sa = sa + h(i,j) * h(j,ien)
  660       continue
            if (wi(i) .lt. 0.0) then
              z = w
              r = ra
              s = sa
            else
              m = i
              if (wi(i) .eq. 0.0) then
                den = w**2 + q**2
                h(i,na) = - (ra * w + sa * q) / den
                h(i,ien) = (ra * q - sa * w) / den
              else
                x = h(i,i+1)
                y = h(i+1,i)
                vr = (wr(i) - p)**2 + wi(i)**2 - q**2
                vi = 2.0 * (wr(i) - p) * q
                if (vr .eq. 0.0  .and.  vi .eq. 0.0) then
                  vr = epsmch * hnorm
     +               * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z))
                endif
                tempr = x * r - z * ra + q * sa
                tempi = x * s - z * sa - q * ra
                den = vr**2 + vi**2
                h(i,na) = (tempr * vr + tempi * vi) / den
                h(i,ien) = (tempi * vr - tempr * vi) / den
                if (abs(x) .gt. abs(z) + abs(q)) then
                  h(i+1,na) = (- ra - w * h(i,na) + q * h(i,ien)) / x
                  h(i+1,ien) = (- sa - w * h(i,ien) - q * h(i,na)) / x
                else
                  tempr = - r - y * h(i,na)
                  tempi = - s - y * h(i,ien)
                  den = z**2 + q**2
                  h(i+1,na) = (tempr * z + tempi * q) / den
                  h(i+1,ien) = (tempi * z - tempr * q) / den
                endif
              endif
            endif
  680     continue
        endif
  690 continue
 
c==== vectors of isolated roots.
      do 720 i = 1, n
        if (i .lt. ilow  .or.  i .gt. iupp) then
          do 710 j = i, n
            vecs(i,j) = h(i,j)
  710     continue
        endif
  720 continue
 
c==== multiply by transformation matrix to give eigenvectors of the
c     original full matrix.
      do 790 j = n, ilow, - 1
        m = min(j,iupp)
        if (wi(j) .lt. 0.0) then
          l = j - 1
          do 740 i = ilow, iupp
            y = 0.0
            z = 0.0
            do 730 k = ilow, m
              y = y + vecs(i,k) * h(k,l)
              z = z + vecs(i,k) * h(k,j)
  730       continue
            vecs(i,l) = y
            vecs(i,j) = z
  740     continue
        else if (wi(j) .eq. 0.0) then
          do 760 i = ilow, iupp
            z = 0.0
            do 750 k = ilow, m
              z = z + vecs(i,k) * h(k,j)
  750       continue
            vecs(i,j) = z
  760     continue
        endif
  790 continue
 
 9999 end
      subroutine orthes(ndim, n, ilow, iupp, a, d)
      implicit none
      integer i,ilow,iupp,j,m,n,ndim
      double precision a,d,f,g,h,scale
c----------------------------------------------------------------------*
c purpose:                                                             *
c   converts an unsymmetric real matrix, a, to upper hessenberg form   *
c   applying successive orthogonal transformations.                    *
c                                                                      *
c   translation of the algol procedure orthes in:                      *
c   handbook series linear algebra,                                    *
c   num. math. 12, 349-368 (1968) by r. s. martin and j. h. wilkinson. *
c input:                                                               *
c   n         (integer) order of the matrix a.                         *
c   ilow,iupp (integer) determine a submatrix, set by balanc.          *
c                       may be set to 1 and n respectively.            *
c   a(ndim,n) (real)    input matrix.                                  *
c output:                                                              *
c   a(ndim,n) (real)    the matrix a, converted to upper hessenberg.   *
c                       the lower triangle contains information        *
c                       about the orthogonal transformations.          *
c   d(n)      (real)    further information.                           *
c----------------------------------------------------------------------*
 
      dimension         a(ndim,n), d(n)
c-----------------------------------------------------------------------
      do 90 m = ilow + 1, iupp - 1
        h = 0.0
        d(m) = 0.0
 
c---- find scale factor.
        scale = 0.0
        do 10 i = m, iupp
          scale = scale + abs(a(i,m-1))
   10   continue
        if (scale .ne. 0.0) then
          do 20 i = iupp, m, - 1
            d(i) = a(i,m-1) / scale
            h = h + d(i) * d(i)
   20     continue
 
          g = sign(sqrt(h),d(m))
          h = h + d(m) * g
          d(m) = d(m) + g
 
c---- form (i - (u*ut) / h) * a.
          do 50 j = m, n
            f = 0.0
            do 30 i = iupp, m, - 1
              f = f + d(i) * a(i,j)
   30       continue
            f = f / h
            do 40 i = m, iupp
              a(i,j) = a(i,j) - f * d(i)
   40       continue
 
   50     continue
 
c---- form (i - (u*ut) / h) * a * (i - (u*ut) / h).
          do 80 i = 1, iupp
            f = 0.0
            do 60 j = iupp, m, - 1
              f = f + d(j) * a(i,j)
   60       continue
            f = f / h
            do 70 j = m, iupp
              a(i,j) = a(i,j) - f * d(j)
   70       continue
   80     continue
 
          d(m) = scale * d(m)
          a(m,m-1) = - scale * g
        endif
   90 continue
c-----------------------------------------------------------------------
      end
      subroutine ortran(ndim, n, ilow, iupp, h, d, v)
      implicit none
      integer i,ilow,iupp,j,k,m,n,ndim
      double precision d,h,v,x,y
c----------------------------------------------------------------------*
c purpose:                                                             *
c   accumulate the orthogonal similarity transformation used by        *
c   orthes to reduce a general real matrix a to upper hessenberg form. *
c                                                                      *
c   translation of the algol procedure ortrans in:                     *
c   handbook series linear algebra,                                    *
c   num. math. 16, 181-204 (1970) by g. peters and j. h. wilkinson.    *
c input:                                                               *
c   n         (integer) order of the matrices a and v.                 *
c   ilow,iupp (integer) determine a sub-matrix set by balanc.          *
c                       may be set to 1 and n respectively.            *
c   h(ndim,n) (real)    the matrix resulting from running orthes.      *
c   d(n)      (real)    further information about the transformation.  *
c output:                                                              *
c   v(ndim,n) (real)    the accumulated transformation.                *
c   d(n)      (real)    destroyed.                                     *
c----------------------------------------------------------------------*
 
      dimension         h(ndim,n), d(n), v(ndim,n)
c-----------------------------------------------------------------------
c---- initialize v to identity matrix.
      do 20 i = 1, n
        do 10 j = 1, n
          v(i,j) = 0.0
   10   continue
        v(i,i) = 1.0
   20 continue
 
c---- accumulate transformations.
      do 90 k = iupp - 2, ilow, - 1
        m = k + 1
        y = h(m,k)
        if (y .ne. 0.0) then
          y = y * d(m)
 
          do 30 i = k + 2, iupp
            d(i) = h(i,k)
   30     continue
c
          do 60 j = m, iupp
            x = 0.0
            do 40 i = m, iupp
              x = x + d(i) * v(i,j)
   40       continue
            x = x / y
            do 50 i = m, iupp
              v(i,j) = v(i,j) + x * d(i)
   50       continue
   60     continue
        endif
   90 continue
c-----------------------------------------------------------------------
      end
      real function rg32cut(cut)
      implicit none
      real cut, rg32
      real rvec(10)
      integer klen
      klen = 5
   10 continue
!      call rnorml(rvec,klen)
      rg32cut = rvec(1)
      if (abs(rg32cut) .gt. cut .and. cut .gt. 0.)  goto 10
      end
      double precision function sigx(bunchnum, ipnt, beamnum)
      implicit none
      integer bunchnum, ipnt, beamnum
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      sigx = sqrt(epsx(bunchnum, beamnum) * betx(ipnt,beamnum)
     +           + (deltap * dx(ipnt,beamnum))**2)
      end



      double precision function sigy(bunchnum, ipnt, beamnum)
      implicit none
      integer bunchnum, ipnt, beamnum
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      sigy = sqrt(epsy(bunchnum, beamnum) * bety(ipnt,beamnum)
     +           + (deltap * dy(ipnt,beamnum))**2)
      end



      double precision function orbdev(kick,beta,betast,q,mu)

      double precision pi
c-----------------------------------------------------------------------
c pi:
      parameter         (pi     = 3.1415926535898d0)
c-----------------------------------------------------------------------
      double precision kick, betast,beta,q,mu

      orbdev = kick*sqrt(betast*beta)/(2*sin(pi*q))*cos(mu-pi*q)
      end


      subroutine matrzero(a_in,a_count1,a_count2,a_count3)
      implicit none
      double precision a_in(a_count1,a_count2,a_count3)
      integer i,j,k, a_count1,a_count2,a_count3
      do i = 1, a_count1
         do j = 1, a_count2
            do k = 1,a_count3
              a_in(i,j,k) = 0.d0
            enddo
         enddo
      enddo
      end
      
      subroutine selipnum( npar, mcol,mdslt,k,output)
      implicit none
      integer npar(8),ipos(8)
      integer mcol, mdslt, k, output
      integer right, left, i

      common /selnpar/ ipos
      save /selnpar/
      
c Return the number of parasitic interactions of each ip depending on 
c which pit the interaction point is closer to      
      do i = 1,8
         left = ipos(i) - mcol
         right = ipos(i) + mcol
        
         if ((k .ge. left) .and. (k .le. right)) then
            output = npar(i)
            return
         endif
      enddo
      write(*,*) 'error: Interaction point ',k,
     +' not associated with any pit' 
      stop  
      end




c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mktunemad2(flag)
      implicit none
      integer i,j,k,l,last,m
      double precision aux1,aux2,cosmux,cosmuy,rp1,rp2,sinmu2,
     +sinmux,sinmuy,temp1,temp2,q1f,q2f,q1b,q2b
c-----------------------------------------------------------------------
c analyse tunes
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical flag
      double precision rp11,rp12
      dimension aux1(6), aux2(6),rp11(6,6),rp12(6,6)
      dimension rp1(4,4), rp2(4,4)
      double precision a1(2,2),b1(2,2),a2(2,2),b2(2,2)
      double precision qx1aux,qy1aux,qx2aux,qy2aux 
      double precision dqmin1,dqmin2,deth1,deth2,rb(2,2)
      character*24 filename1, filename2
      double precision auxf,auxf2,auxb,auxb2,rf(2,2)
      double precision betxt(mcol,2),betyt(mcol,2)

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/ 
c-----------------------------------------------------------------------
      print *, 'Building second-order maps . . .'

c Opening file where linear coupling will be stored      

      filename1 = 'fort.closest_tune_app_f'
      open(13, file=filename1, status = 'UNKNOWN')
      filename2 = 'fort.closest_tune_app_b'
      open(14, file=filename2, status = 'UNKNOWN')

c track second-order maps
      if (flag) then
         call trackb2 (.true., .false., -1)
         last = min(c_tunes, nbunch)
         
      else
         call track02(.true., .false.)
         last = 1
      endif
      
      
c loop over bunches
      m = nmaps + 1
      do 90 l = 1, last
         if (l .eq. 1 .or. l .eq. last .or. mod(l,100) .eq. 0)
     +   print *, 'Analysing bunch pair ', code(l), ' . . .'
c copy x,y components of the one turn map into auxiliar variables
         do i = 1, 4
           do j = 1, 4
             rp1(j,i) = tr1(j,i,l)
             rp2(j,i) = tr2(j,i,l)
           enddo
         enddo
         
c Find diagonal boxes of one turn map for beam 1
         call eigenmad(rp1,4,a1,b1,deth1,rf,auxf2)

c Find diagonal boxes of one turn map for beam 2
         call eigenmad(rp2,4,a2,b2,deth2,rb,auxb2)

c Compute the tunes and optic functions         
         call twisspar(a1,qx1(l),betx(1,1),alfx(1,1))
         call twisspar(b1,qy1(l),bety(1,1),alfy(1,1))
         
         call twisspar(a2,qx2(l),betx(m,2),alfx(m,2))
         call twisspar(b2,qy2(l),bety(m,2),alfy(m,2))


c Compute linear coupling
         call lnrcoup(deth1,qx1(l),qy1(l),dqmin1)
         call lnrcoup(deth2,qx2(l),qy2(l),dqmin2)
c aux, dqminf,betxf,betyf,gammaf
         write(13,*) set_f(l),dqmin1,betxt(1,1),betyt(1,1),gamma,rf
         write(14,*) set_b(l),dqmin2,betxt(m,2),betyt(m,2),gamma,rb
c compute initial first-order dispersion
         
         do i = 1, 6
           do j = 1, 6
             rp11(j,i) = tr1(j,i,l)
             rp12(j,i) = tr2(j,i,l)
           enddo
         enddo
         
         call eigen(rp11, 6, q1f, q2f)
         call eigen(rp12, 6, q1b, q2b)

         call disp(tr1(1,1,l), tr1(1,6,l), d1(1,l,0))
         call disp(tr2(1,1,l), tr2(1,6,l), d2(1,l,m))
         d1(5,l,0) = 0.0
         d1(6,l,0) = 1.0
         d2(5,l,m) = 0.0
         d2(6,l,m) = 1.0
c derivative of transfer matrix w.r.t. delta(p)/p
         do 30 i = 1, 6
            aux1(i) = 0.0
            aux2(i) = 0.0
 
            do 20 k = 1, 6
               temp1 = 0.0
               temp2 = 0.0
 
               do 10 j = 1, 6
                  temp1 = temp1 + tt1(i,j,k,l) * d1(j,l,0)
                  temp2 = temp2 + tt2(i,j,k,l) * d2(j,l,m)
 10            continue
 
               aux1(i) = aux1(i) + temp1 * d1(k,l,0)
               rp11(i,k) = 2.0 * temp1
               aux2(i) = aux2(i) + temp2 * d2(k,l,m)
               rp12(i,k) = 2.0 * temp2
 20         continue
 30      continue
 
c compute inital second-order dispersion
         call disp(tr1(1,1,l), aux1, dd1(1,l,0))
         call disp(tr2(1,1,l), aux2, dd2(1,l,m))
         dd1(5,l,0) = 0.0
         dd1(6,l,0) = 0.0
         dd2(5,l,m) = 0.0
         dd2(6,l,m) = 0.0
c horizontal motion
         sinmu2 = - tr1(1,2,l) * tr1(2,1,l) -
     +        (tr1(1,1,l) - tr1(2,2,l))**2 / 4.0
         if (sinmu2 .lt. 0.0) sinmu2 = toler
         sinmux = sign(sqrt(sinmu2), tr1(1,2,l))
         qxp1(l) = - (rp11(1,1) + rp11(2,2))/(4.0 * pi * sinmux)

         sinmu2 = - tr2(1,2,l) * tr2(2,1,l) -
     +        (tr2(1,1,l) - tr2(2,2,l))**2 / 4.0
         if (sinmu2 .lt. 0.0) sinmu2 = toler
         sinmux = sign(sqrt(sinmu2), tr2(1,2,l))
         qxp2(l) = - (rp12(1,1) + rp12(2,2))/(4.0 * pi * sinmux)

c vertical motion
         sinmu2 = - tr1(3,4,l) * tr1(4,3,l) -
     +        (tr1(3,3,l) - tr1(4,4,l))**2 / 4.0
         if (sinmu2 .lt. 0.) sinmu2 = toler
         sinmuy = sign(sqrt(sinmu2), tr1(3,4,l))
         qyp1(l) = - (rp11(3,3) + rp11(4,4))/(4.0 * pi * sinmuy)

         sinmu2 = - tr2(3,4,l) * tr2(4,3,l) -
     +        (tr2(3,3,l) - tr2(4,4,l))**2 / 4.0
         if (sinmu2 .lt. 0.) sinmu2 = toler
         sinmuy = sign(sqrt(sinmu2), tr2(3,4,l))
         qyp2(l) = - (rp12(3,3) + rp12(4,4))/(4.0 * pi * sinmuy)

c track dispersion
         if (flag) then
            call trackb2(.true., .true., l)
         else
            call track02(.true., .true.)
         endif
 90   continue
      close(13)
      close(14)
c-----------------------------------------------------------------------
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine mktunemad(flag)
      implicit none
      integer i,j,k,l,last,m
      double precision aux1,aux2,cosmux,cosmuy,rp1,rp2,sinmu2,
     +sinmux,sinmuy,temp1,temp2,q1f,q2f,q1b,q2b
c-----------------------------------------------------------------------
c analyse tunes
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
c--- equivalenced with collision class arrays to save space
      integer mstart, madd
      parameter (mstart = 15*mbuck+1, madd = mbuck*mcol)
      real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1)
c   ntotal_f           sum of all ordered list bunches = total number
c   ctotal_f           sum of all equ. list bunches = no. of equ. part.
c   nset_f             no. of bunches (from hitlist) in set_f
c   hitlist_f(i)       bunch (slot) mask for all equ. class bunches
c   set_f(i)           bunch (slot) number (i.e. all equ. bunches)
c   colcnt_f(i)        no. of collision points of bunch i
c   list_f(j,i)        collision point numbers (j) of bunch i
c   part_f(j,i)        colliding bunch in backward beam at coll. point
c   ordl        number of ordered collision lists (i.e. number of one's
c               in a)
c   cordl(i)    collision count for lordl(i)
c   nordl(i)    number of bunches (slots) in lordl(i)
c   lordl(j,i)  bunch number j in list lordl(i)
c   equl, cequ, nequl, lequl as ordl etc. above where each class is
c   a different ordered list of collision points (equ. classes)
      integer hitlist_f(0:mbuck-1),
     +colcnt_f(0:mbuck-1), list_f(mcol,0:mbuck-1),
     +part_f(mcol,0:mbuck-1),
     +hitlist_b(0:mbuck-1),
     +colcnt_b(0:mbuck-1), list_b(mcol,0:mbuck-1),
     +part_b(mcol,0:mbuck-1),
     +cequl_f(mcol), nequl_f(mcol), lequl_f(mbuck,mcol),
     +cordl_f(mcol), nordl_f(mcol), lordl_f(mbuck,mcol),
     +cequl_b(mcol), nequl_b(mcol), lequl_b(mbuck,mcol),
     +cordl_b(mcol), nordl_b(mcol), lordl_b(mbuck,mcol)
 
      integer total(60*mbunch*mcol)
      equivalence (total, z1)
      equivalence (hitlist_f(0), total(mbuck+1)),
     +(colcnt_f(0), total(3*mbuck+1)),
     +(hitlist_b(0), total(4*mbuck+1)),
     +(colcnt_b(0), total(6*mbuck+1)),
     +(cequl_f(1),total(7*mbuck+1)), (nequl_f(1),total(8*mbuck+1)),
     +(cordl_f(1),total(9*mbuck+1)), (nordl_f(1),total(10*mbuck+1)),
     +(cequl_b(1),total(11*mbuck+1)), (nequl_b(1),total(12*mbuck+1)),
     +(cordl_b(1),total(13*mbuck+1)), (nordl_b(1),total(14*mbuck+1)),
     +(list_f(1,0), total(mstart)),
     +(part_f(1,0), total(mstart+madd)),
     +(list_b(1,0), total(mstart+2*madd)),
     +(part_b(1,0), total(mstart+3*madd)),
     +(lequl_f(1,1), total(mstart+4*madd)),
     +(lordl_f(1,1), total(mstart+5*madd)),
     +(lequl_b(1,1), total(mstart+6*madd)),
     +(lordl_b(1,1), total(mstart+7*madd))
      save /corbit/
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision sk1,sk2,sr1,sr2,st1,st2
c-----------------------------------------------------------------------
c maps from one interaction point to the next
      common /sector/ sk1(6,0:mmaps), sr1(6,6,0:mmaps),
     +                st1(6,6,6,0:mmaps),
     +                sk2(6,0:mmaps), sr2(6,6,0:mmaps),
     +                st2(6,6,6,0:mmaps)
      save /sector/
      double precision tr1,tr2,tt1,tt2
c-----------------------------------------------------------------------
c maps per bunch for one turn
      common /turn/ tr1(6,6,mbunch), tt1(6,6,6,mbunch),
     +              tr2(6,6,mbunch), tt2(6,6,6,mbunch)
      save /turn/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/
c-----------------------------------------------------------------------
      logical flag
      double precision rp11,rp12
      dimension aux1(6), aux2(6),rp11(6,6),rp12(6,6)
      dimension rp1(4,4), rp2(4,4)
      double precision a1(2,2),b1(2,2),a2(2,2),b2(2,2)
      double precision qx1aux,qy1aux,qx2aux,qy2aux 
      double precision dqmin1,dqmin2,deth1,deth2,rb(2,2)
      character*24 filename1, filename2
      double precision auxf,auxf2,auxb,auxb2,rf(2,2)
      double precision betxt(mcol,2),betyt(mcol,2)
c-----------------------------------------------------------------------
      print *, 'Building second-order maps . . .'

c Opening file where linear coupling will be stored      

      filename1 = 'fort.closest_tune_app_f'
      open(13, file=filename1, status = 'UNKNOWN')
      filename2 = 'fort.closest_tune_app_b'
      open(14, file=filename2, status = 'UNKNOWN')

c track second-order maps
      write(*,*) 'ctunes,nbunch',c_tunes,nbunch
      if (flag) then
         call trackb(.true., .false., -1)
         last = min(c_tunes, nbunch)
      else
         call track0(.true., .false.)
         last = 1
      endif
      
      
c loop over bunches
      m = ninter + 1
      do 90 l = 1, last
         if (l .eq. 1 .or. l .eq. last .or. mod(l,100) .eq. 0)
     +   print *, 'Analysing bunch pair ', code(l), ' . . .'
c copy x,y components of the one turn map into auxiliar variables
         do i = 1, 4
           do j = 1, 4
             rp1(j,i) = tr1(j,i,l)
             rp2(j,i) = tr2(j,i,l)
           enddo
         enddo
         
c Find diagonal boxes of one turn map for beam 1
         call eigenmad(rp1,4,a1,b1,deth1,rf,auxf2)

c Find diagonal boxes of one turn map for beam 2
         call eigenmad(rp2,4,a2,b2,deth2,rb,auxb2)

c Compute the tunes and optic functions         
         call twisspar(a1,qx1(l),betx(1,1),alfx(1,1))
         call twisspar(b1,qy1(l),bety(1,1),alfy(1,1))
         
         call twisspar(a2,qx2(l),betx(m,2),alfx(m,2))
         call twisspar(b2,qy2(l),bety(m,2),alfy(m,2))

c Compute linear coupling
         call lnrcoup(deth1,qx1(l),qy1(l),dqmin1)
         call lnrcoup(deth2,qx2(l),qy2(l),dqmin2)
c aux, dqminf,betxf,betyf,gammaf
         write(13,*) set_f(l),dqmin1,betxt(1,1),betyt(1,1),gamma,rf
         write(14,*) set_b(l),dqmin2,betxt(m,2),betyt(m,2),gamma,rb
c compute initial first-order dispersion
         do i = 1, 6
           do j = 1, 6
             rp11(j,i) = tr1(j,i,l)
             rp12(j,i) = tr2(j,i,l)
           enddo
         enddo

         call eigen(rp11, 6, q1f, q2f)
         call eigen(rp12, 6, q1b, q2b)

         call disp(tr1(1,1,l), tr1(1,6,l), d1(1,l,0))
         call disp(tr2(1,1,l), tr2(1,6,l), d2(1,l,m))
         d1(5,l,0) = 0.0
         d1(6,l,0) = 1.0
         d2(5,l,m) = 0.0
         d2(6,l,m) = 1.0
c derivative of transfer matrix w.r.t. delta(p)/p
         do 30 i = 1, 6
            aux1(i) = 0.0
            aux2(i) = 0.0
 
            do 20 k = 1, 6
               temp1 = 0.0
               temp2 = 0.0
 
               do 10 j = 1, 6
                  temp1 = temp1 + tt1(i,j,k,l) * d1(j,l,0)
                  temp2 = temp2 + tt2(i,j,k,l) * d2(j,l,m)
 10            continue
 
               aux1(i) = aux1(i) + temp1 * d1(k,l,0)
               rp11(i,k) = 2.0 * temp1
               aux2(i) = aux2(i) + temp2 * d2(k,l,m)
               rp12(i,k) = 2.0 * temp2
 20         continue
 30      continue
 
c compute inital second-order dispersion
         call disp(tr1(1,1,l), aux1, dd1(1,l,0))
         call disp(tr2(1,1,l), aux2, dd2(1,l,m))
         dd1(5,l,0) = 0.0
         dd1(6,l,0) = 0.0
         dd2(5,l,m) = 0.0
         dd2(6,l,m) = 0.0
 
c horizontal motion
         sinmu2 = - tr1(1,2,l) * tr1(2,1,l) -
     +        (tr1(1,1,l) - tr1(2,2,l))**2 / 4.0
         if (sinmu2 .lt. 0.0) sinmu2 = toler
         sinmux = sign(sqrt(sinmu2), tr1(1,2,l))
         qxp1(l) = - (rp11(1,1) + rp11(2,2))/(4.0 * pi * sinmux)
 
         sinmu2 = - tr2(1,2,l) * tr2(2,1,l) -
     +        (tr2(1,1,l) - tr2(2,2,l))**2 / 4.0
         if (sinmu2 .lt. 0.0) sinmu2 = toler
         sinmux = sign(sqrt(sinmu2), tr2(1,2,l))
         qxp2(l) = - (rp12(1,1) + rp12(2,2))/(4.0 * pi * sinmux)
 
c vertical motion
         sinmu2 = - tr1(3,4,l) * tr1(4,3,l) -
     +        (tr1(3,3,l) - tr1(4,4,l))**2 / 4.0
         if (sinmu2 .lt. 0.) sinmu2 = toler
         sinmuy = sign(sqrt(sinmu2), tr1(3,4,l))
         qyp1(l) = - (rp11(3,3) + rp11(4,4))/(4.0 * pi * sinmuy)
 
         sinmu2 = - tr2(3,4,l) * tr2(4,3,l) -
     +        (tr2(3,3,l) - tr2(4,4,l))**2 / 4.0
         if (sinmu2 .lt. 0.) sinmu2 = toler
         sinmuy = sign(sqrt(sinmu2), tr2(3,4,l))
         qyp2(l) = - (rp12(3,3) + rp12(4,4))/(4.0 * pi * sinmuy)

c track dispersion
         if (flag) then
            call trackb(.true., .true., l)
         else
            call track0(.true., .true.)
         endif
 90   continue
      close(13)
      close(14)
c-----------------------------------------------------------------------
      end


      subroutine twisspar(a,mu,beta,alpha) 
      implicit none
c --------------------------------------------------------------------
c Given a 2x2 matrix associates the terms of the transfer map with the twiss
c parameters following the article https://pdfs.semanticscholar.org/4b3f/9b594e44e852e3e2e84ce260b45ccf4c7b3f.pdf
c --------------------------------------------------------------------
      double precision a(2,2)
      double precision cosmu,beta,alpha,sinmu2,sinmu,mu
      double precision toler,pi
      parameter (toler = 1.0d-8)
      parameter (pi = 3.1415926535898d0)

      cosmu = (a(1,1) + a(2,2))/2.0d0
      sinmu2 = -a(1,2)*a(2,1)-0.25d0*(a(1,1)-a(2,2))**2

      if (sinmu2 .lt. 0) sinmu2 = toler

      sinmu = sign(sqrt(sinmu2),a(1,2))

      alpha = (a(1,1) - a(2,2))/(2.0d0*sinmu)
     
      beta = a(1,2)/sinmu
      
      mu = acos(cosmu)/(2*pi) ! Tunes in units of 2*pi

      end


      integer function signdet(det)
      implicit none
c ----------------------------------------------------------------------
c Returns 1 if the input value det is positive and -1 otherwise
c ----------------------------------------------------------------------
      double precision det
      integer signdet     
      signdet = 1
      if (det .le. 0) signdet = -1	
      end function



      subroutine symcheck(a,n,maxtsa)
      implicit none
      double precision a(n,n), toler,det,at(n,n),s(n,n),
     +ats(n,n),atsa(n,n)
      double precision M44DET,M66DET,maxtsa
      integer n,i,j
      parameter (toler = 1.0d-4)

     
      if (mod(n,2) .ne. 0) then
         write(*,*) 'Error- Dimension of matrix',n,', not even'
         stop
      endif
      
      if (n .eq. 2) then
         det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
         s = reshape((/0.,1.,-1.,0./),shape(s))

      elseif (n .eq. 4) then
         det = M44DET(a)
         s = reshape((/0.,1.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,
     +1.,0.,0.,-1.,0./),shape(s))


      else 
         det = M66DET(a)
         s = reshape((/0.,1.,0.,0.,0.,0.,
     +-1.,0.,0.,0.,0.,0.,
     +0.,0.,0.,+1.,0.,0.,
     +0.,0.,-1.,0.,0.,0.,
     +0.,0.,0.,0.,0.,1.,0.,
     +0.,0.,0.,-1.,0./),shape(s)) 
         
      endif
      
      if (abs(1.d0-det) .gt. toler) then
         write(*,*) 'Warning- Determinant of transfer map',
     +det,'different from 1. Symplectic condition not satisfied.'
c         flag = .False.
      endif
      
      at = transpose(a)
      
      call mxmpy(at, s, ats, n, n)
      call mxmpy(ats, a, atsa, n, n)
      maxtsa = 0.0d0 
      do i = 1,n
         do j = 1,n
            if (abs(atsa(i,j) - s(i,j)) .gt. maxtsa) maxtsa =  
     +      abs(atsa(i,j) - s(i,j))
            if (abs(atsa(i,j) - s(i,j)) .gt.toler) then
               
               write(*,*) 'Symplectic condition M^T S M = S',
     +         ', not satisfied'
               if (n .eq. 2) write(*,'(1p,2e14.6)') atsa
               if (n .eq. 4) write(*,'(1p,4e14.6)') atsa 
               if (n .eq. 6) write(*,'(1p,6e14.6)') atsa
               write(*,*) 'Error- Transfer map not symplectic'
c               stop
            endif
         enddo
      enddo
      if (maxtsa.gt. toler)write(*,*)'Maximum difference',maxtsa
      end subroutine


      DOUBLE PRECISION FUNCTION M66DET (A)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6,6), INTENT(IN)  :: A
      DOUBLE PRECISION M66DET
      DOUBLE PRECISION :: A11, A12, A13, A14, A15, A16, A21, 
     +A22, A23, A24, A25, A26, A31, A32, A33, A34, A35, A36, A41, 
     +A42, A43, A44, A45, A46,A51, A52, A53, A54, A55, A56, A61, A62,
     +A63, A64, A65, A66

      A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5);
      A16=A(1,6)
      A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5);
      A26=A(2,6)
      A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5);
      A36=A(3,6)
      A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5);
      A46=A(4,6)
      A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5);
      A56=A(5,6)
      A61=A(6,1); A62=A(6,2); A63=A(6,3); A64=A(6,4); A65=A(6,5);
      A66=A(6,6)

      M66DET = -(A16*A25*A34*A43*A52-A15*A26*A34*A43*A52-
     +A16*A24*A35*A43*A52+A14*A26*A35*A43*A52+A15*A24*A36*A43*A52-
     +A14*A25*A36*A43*A52-A16*A25*A33*A44*A52+A15*A26*A33*A44*A52
     ++A16*A23*A35*A44*A52-A13*A26*A35*A44*A52-A15*A23*A36*A44*A52
     ++A13*A25*A36*A44*A52+A16*A24*A33*A45*A52-A14*A26*A33*A45*A52-
     +A16*A23*A34*A45*A52+A13*A26*A34*A45*A52+A14*A23*A36*A45*A52-
     +A13*A24*A36*A45*A52-A15*A24*A33*A46*A52+A14*A25*A33*A46*A52
     ++A15*A23*A34*A46*A52-A13*A25*A34*A46*A52-A14*A23*A35*A46*A52
     ++A13*A24*A35*A46*A52-A16*A25*A34*A42*A53+A15*A26*A34*A42*A53
     ++A16*A24*A35*A42*A53-A14*A26*A35*A42*A53-A15*A24*A36*A42*A53
     ++A14*A25*A36*A42*A53+A16*A25*A32*A44*A53-A15*A26*A32*A44*A53-
     +A16*A22*A35*A44*A53+A12*A26*A35*A44*A53+A15*A22*A36*A44*A53-
     +A12*A25*A36*A44*A53-A16*A24*A32*A45*A53+A14*A26*A32*A45*A53
     ++A16*A22*A34*A45*A53-A12*A26*A34*A45*A53-A14*A22*A36*A45*A53
     ++A12*A24*A36*A45*A53+A15*A24*A32*A46*A53-A14*A25*A32*A46*A53- 
     +A15*A22*A34*A46*A53+A12*A25*A34*A46*A53+A14*A22*A35*A46*A53- 
     +A12*A24*A35*A46*A53+A16*A25*A33*A42*A54-A15*A26*A33*A42*A54-
     +A16*A23*A35*A42*A54+A13*A26*A35*A42*A54+A15*A23*A36*A42*A54-
     +A13*A25*A36*A42*A54-A16*A25*A32*A43*A54+A15*A26*A32*A43*A54
     ++A16*A22*A35*A43*A54-A12*A26*A35*A43*A54-A15*A22*A36*A43*A54
     ++A12*A25*A36*A43*A54+A16*A23*A32*A45*A54-A13*A26*A32*A45*A54-
     +A16*A22*A33*A45*A54+A12*A26*A33*A45*A54+A13*A22*A36*A45*A54-
     +A12*A23*A36*A45*A54-A15*A23*A32*A46*A54+A13*A25*A32*A46*A54
     ++A15*A22*A33*A46*A54-A12*A25*A33*A46*A54-A13*A22*A35*A46*A54 
     ++A12*A23*A35*A46*A54-A16*A24*A33*A42*A55+A14*A26*A33*A42*A55
     ++A16*A23*A34*A42*A55-A13*A26*A34*A42*A55-A14*A23*A36*A42*A55
     ++A13*A24*A36*A42*A55+A16*A24*A32*A43*A55-A14*A26*A32*A43*A55-
     +A16*A22*A34*A43*A55+A12*A26*A34*A43*A55+A14*A22*A36*A43*A55-
     +A12*A24*A36*A43*A55-A16*A23*A32*A44*A55+A13*A26*A32*A44*A55
     ++A16*A22*A33*A44*A55-A12*A26*A33*A44*A55-A13*A22*A36*A44*A55
     ++A12*A23*A36*A44*A55+A14*A23*A32*A46*A55-A13*A24*A32*A46*A55-
     +A14*A22*A33*A46*A55+A12*A24*A33*A46*A55+A13*A22*A34*A46*A55-
     +A12*A23*A34*A46*A55+A15*A24*A33*A42*A56-A14*A25*A33*A42*A56-
     +A15*A23*A34*A42*A56+A13*A25*A34*A42*A56+A14*A23*A35*A42*A56-
     +A13*A24*A35*A42*A56-A15*A24*A32*A43*A56+A14*A25*A32*A43*A56
     ++A15*A22*A34*A43*A56-A12*A25*A34*A43*A56-A14*A22*A35*A43*A56
     ++A12*A24*A35*A43*A56+A15*A23*A32*A44*A56-A13*A25*A32*A44*A56-
     +A15*A22*A33*A44*A56+A12*A25*A33*A44*A56+A13*A22*A35*A44*A56-
     +A12*A23*A35*A44*A56-A14*A23*A32*A45*A56+A13*A24*A32*A45*A56
     ++A14*A22*A33*A45*A56-A12*A24*A33*A45*A56-A13*A22*A34*A45*A56
     ++A12*A23*A34*A45*A56)*A61+(A16*A25*A34*A43*A51-
     +A15*A26*A34*A43*A51-A16*A24*A35*A43*A51+A14*A26*A35*A43*A51
     ++A15*A24*A36*A43*A51-A14*A25*A36*A43*A51-A16*A25*A33*A44*A51
     ++A15*A26*A33*A44*A51+A16*A23*A35*A44*A51-A13*A26*A35*A44*A51-
     +A15*A23*A36*A44*A51+A13*A25*A36*A44*A51+A16*A24*A33*A45*A51-
     +A14*A26*A33*A45*A51-A16*A23*A34*A45*A51+A13*A26*A34*A45*A51
     ++A14*A23*A36*A45*A51-A13*A24*A36*A45*A51-A15*A24*A33*A46*A51
     ++A14*A25*A33*A46*A51+A15*A23*A34*A46*A51-A13*A25*A34*A46*A51-
     +A14*A23*A35*A46*A51+A13*A24*A35*A46*A51-A16*A25*A34*A41*A53
     ++A15*A26*A34*A41*A53+A16*A24*A35*A41*A53-A14*A26*A35*A41*A53-
     +A15*A24*A36*A41*A53+A14*A25*A36*A41*A53+A16*A25*A31*A44*A53-
     +A15*A26*A31*A44*A53-A16*A21*A35*A44*A53+A11*A26*A35*A44*A53
     ++A15*A21*A36*A44*A53-A11*A25*A36*A44*A53-A16*A24*A31*A45*A53
     ++A14*A26*A31*A45*A53+A16*A21*A34*A45*A53-A11*A26*A34*A45*A53-
     +A14*A21*A36*A45*A53+A11*A24*A36*A45*A53+A15*A24*A31*A46*A53-
     +A14*A25*A31*A46*A53-A15*A21*A34*A46*A53+A11*A25*A34*A46*A53
     ++A14*A21*A35*A46*A53-A11*A24*A35*A46*A53+A16*A25*A33*A41*A54-
     +A15*A26*A33*A41*A54-A16*A23*A35*A41*A54+A13*A26*A35*A41*A54
     ++A15*A23*A36*A41*A54-A13*A25*A36*A41*A54-A16*A25*A31*A43*A54
     ++A15*A26*A31*A43*A54+A16*A21*A35*A43*A54-A11*A26*A35*A43*A54-
     +A15*A21*A36*A43*A54+A11*A25*A36*A43*A54+A16*A23*A31*A45*A54-
     +A13*A26*A31*A45*A54-A16*A21*A33*A45*A54+A11*A26*A33*A45*A54
     ++A13*A21*A36*A45*A54-A11*A23*A36*A45*A54-A15*A23*A31*A46*A54
     ++A13*A25*A31*A46*A54+A15*A21*A33*A46*A54-A11*A25*A33*A46*A54-
     +A13*A21*A35*A46*A54+A11*A23*A35*A46*A54-A16*A24*A33*A41*A55
     ++A14*A26*A33*A41*A55+A16*A23*A34*A41*A55-A13*A26*A34*A41*A55-
     +A14*A23*A36*A41*A55+A13*A24*A36*A41*A55+A16*A24*A31*A43*A55- 
     +A14*A26*A31*A43*A55-A16*A21*A34*A43*A55+A11*A26*A34*A43*A55
     ++A14*A21*A36*A43*A55-A11*A24*A36*A43*A55-A16*A23*A31*A44*A55
     ++A13*A26*A31*A44*A55+A16*A21*A33*A44*A55-A11*A26*A33*A44*A55-
     +A13*A21*A36*A44*A55+A11*A23*A36*A44*A55+A14*A23*A31*A46*A55-
     +A13*A24*A31*A46*A55-A14*A21*A33*A46*A55+A11*A24*A33*A46*A55
     ++A13*A21*A34*A46*A55-A11*A23*A34*A46*A55+A15*A24*A33*A41*A56-
     +A14*A25*A33*A41*A56-A15*A23*A34*A41*A56+A13*A25*A34*A41*A56 
     ++A14*A23*A35*A41*A56-A13*A24*A35*A41*A56-A15*A24*A31*A43*A56
     ++A14*A25*A31*A43*A56+A15*A21*A34*A43*A56-A11*A25*A34*A43*A56-
     +A14*A21*A35*A43*A56+A11*A24*A35*A43*A56+A15*A23*A31*A44*A56-
     +A13*A25*A31*A44*A56-A15*A21*A33*A44*A56+A11*A25*A33*A44*A56
     ++A13*A21*A35*A44*A56-A11*A23*A35*A44*A56-A14*A23*A31*A45*A56
     ++A13*A24*A31*A45*A56+A14*A21*A33*A45*A56-A11*A24*A33*A45*A56-
     +A13*A21*A34*A45*A56+A11*A23*A34*A45*A56)*A62-
     +(A16*A25*A34*A42*A51-A15*A26*A34*A42*A51-A16*A24*A35*A42*A51
     ++A14*A26*A35*A42*A51+A15*A24*A36*A42*A51-A14*A25*A36*A42*A51-
     +A16*A25*A32*A44*A51+A15*A26*A32*A44*A51+A16*A22*A35*A44*A51-
     +A12*A26*A35*A44*A51-A15*A22*A36*A44*A51+A12*A25*A36*A44*A51
     ++A16*A24*A32*A45*A51-A14*A26*A32*A45*A51-A16*A22*A34*A45*A51
     ++A12*A26*A34*A45*A51+A14*A22*A36*A45*A51-A12*A24*A36*A45*A51-
     +A15*A24*A32*A46*A51+A14*A25*A32*A46*A51+A15*A22*A34*A46*A51-
     +A12*A25*A34*A46*A51-A14*A22*A35*A46*A51+A12*A24*A35*A46*A51-
     +A16*A25*A34*A41*A52+A15*A26*A34*A41*A52+A16*A24*A35*A41*A52-
     +A14*A26*A35*A41*A52-A15*A24*A36*A41*A52+A14*A25*A36*A41*A52
     ++A16*A25*A31*A44*A52-A15*A26*A31*A44*A52-A16*A21*A35*A44*A52
     ++A11*A26*A35*A44*A52+A15*A21*A36*A44*A52-A11*A25*A36*A44*A52-
     +A16*A24*A31*A45*A52+A14*A26*A31*A45*A52+A16*A21*A34*A45*A52-
     +A11*A26*A34*A45*A52-A14*A21*A36*A45*A52+A11*A24*A36*A45*A52 
     ++A15*A24*A31*A46*A52-A14*A25*A31*A46*A52-A15*A21*A34*A46*A52
     ++A11*A25*A34*A46*A52+A14*A21*A35*A46*A52-A11*A24*A35*A46*A52
     ++A16*A25*A32*A41*A54-A15*A26*A32*A41*A54-A16*A22*A35*A41*A54
     ++A12*A26*A35*A41*A54+A15*A22*A36*A41*A54-A12*A25*A36*A41*A54-
     +A16*A25*A31*A42*A54+A15*A26*A31*A42*A54+A16*A21*A35*A42*A54-
     +A11*A26*A35*A42*A54-A15*A21*A36*A42*A54+A11*A25*A36*A42*A54
     ++A16*A22*A31*A45*A54-A12*A26*A31*A45*A54-A16*A21*A32*A45*A54
     ++A11*A26*A32*A45*A54+A12*A21*A36*A45*A54-A11*A22*A36*A45*A54-
     +A15*A22*A31*A46*A54+A12*A25*A31*A46*A54+A15*A21*A32*A46*A54-
     +A11*A25*A32*A46*A54-A12*A21*A35*A46*A54+A11*A22*A35*A46*A54-
     +A16*A24*A32*A41*A55+A14*A26*A32*A41*A55+A16*A22*A34*A41*A55-
     +A12*A26*A34*A41*A55-A14*A22*A36*A41*A55+A12*A24*A36*A41*A55
     ++A16*A24*A31*A42*A55-A14*A26*A31*A42*A55-A16*A21*A34*A42*A55
     ++A11*A26*A34*A42*A55+A14*A21*A36*A42*A55-A11*A24*A36*A42*A55-
     +A16*A22*A31*A44*A55+A12*A26*A31*A44*A55+A16*A21*A32*A44*A55-
     +A11*A26*A32*A44*A55-A12*A21*A36*A44*A55+A11*A22*A36*A44*A55
     ++A14*A22*A31*A46*A55-A12*A24*A31*A46*A55-A14*A21*A32*A46*A55
     ++A11*A24*A32*A46*A55+A12*A21*A34*A46*A55-A11*A22*A34*A46*A55
     ++A15*A24*A32*A41*A56-A14*A25*A32*A41*A56-A15*A22*A34*A41*A56
     ++A12*A25*A34*A41*A56+A14*A22*A35*A41*A56-A12*A24*A35*A41*A56-
     +A15*A24*A31*A42*A56+A14*A25*A31*A42*A56+A15*A21*A34*A42*A56-
     +A11*A25*A34*A42*A56-A14*A21*A35*A42*A56+A11*A24*A35*A42*A56
     ++A15*A22*A31*A44*A56-A12*A25*A31*A44*A56-A15*A21*A32*A44*A56
     ++A11*A25*A32*A44*A56+A12*A21*A35*A44*A56-A11*A22*A35*A44*A56-
     +A14*A22*A31*A45*A56+A12*A24*A31*A45*A56+A14*A21*A32*A45*A56-
     +A11*A24*A32*A45*A56-A12*A21*A34*A45*A56+A11*A22*A34*A45*A56)*
     +A63+(A16*A25*A33*A42*A51-A15*A26*A33*A42*A51-A16*A23*A35*A42*
     +A51+A13*A26*A35*A42*A51+A15*A23*A36*A42*A51-A13*A25*A36*A42*
     +A51-A16*A25*A32*A43*A51+A15*A26*A32*A43*A51+A16*A22*A35*A43*
     +A51-A12*A26*A35*A43*A51-A15*A22*A36*A43*A51+A12*A25*A36*A43*
     +A51+A16*A23*A32*A45*A51-A13*A26*A32*A45*A51-A16*A22*A33*A45*
     +A51+A12*A26*A33*A45*A51+A13*A22*A36*A45*A51-A12*A23*A36*A45*
     +A51-A15*A23*A32*A46*A51+A13*A25*A32*A46*A51+A15*A22*A33*A46*
     +A51-A12*A25*A33*A46*A51-A13*A22*A35*A46*A51+A12*A23*A35*A46*
     +A51-A16*A25*A33*A41*A52+A15*A26*A33*A41*A52+A16*A23*A35*A41*
     +A52-A13*A26*A35*A41*A52-A15*A23*A36*A41*A52+A13*A25*A36*A41*
     +A52+A16*A25*A31*A43*A52-A15*A26*A31*A43*A52-A16*A21*A35*A43*
     +A52+A11*A26*A35*A43*A52+A15*A21*A36*A43*A52-A11*A25*A36*A43*
     +A52-A16*A23*A31*A45*A52+A13*A26*A31*A45*A52+A16*A21*A33*A45*
     +A52-A11*A26*A33*A45*A52-A13*A21*A36*A45*A52+A11*A23*A36*A45*
     +A52+A15*A23*A31*A46*A52-A13*A25*A31*A46*A52-A15*A21*A33*A46*
     +A52+A11*A25*A33*A46*A52+A13*A21*A35*A46*A52-A11*A23*A35*A46*
     +A52+A16*A25*A32*A41*A53-A15*A26*A32*A41*A53-A16*A22*A35*A41*
     +A53+A12*A26*A35*A41*A53+A15*A22*A36*A41*A53-A12*A25*A36*A41*
     +A53-A16*A25*A31*A42*A53+A15*A26*A31*A42*A53+A16*A21*A35*A42*
     +A53-A11*A26*A35*A42*A53-A15*A21*A36*A42*A53+A11*A25*A36*A42*
     +A53+A16*A22*A31*A45*A53-A12*A26*A31*A45*A53-
     +A16*A21*A32*A45*A53+A11*A26*A32*A45*A53+A12*A21*A36*A45*A53-
     +A11*A22*A36*A45*A53-A15*A22*A31*A46*A53+A12*A25*A31*A46*A53
     ++A15*A21*A32*A46*A53-A11*A25*A32*A46*A53-A12*A21*A35*A46*A53
     ++A11*A22*A35*A46*A53-A16*A23*A32*A41*A55+A13*A26*A32*A41*A55
     ++A16*A22*A33*A41*A55-A12*A26*A33*A41*A55-A13*A22*A36*A41*A55
     ++A12*A23*A36*A41*A55+A16*A23*A31*A42*A55-A13*A26*A31*A42*A55-
     +A16*A21*A33*A42*A55+A11*A26*A33*A42*A55+A13*A21*A36*A42*A55-
     +A11*A23*A36*A42*A55-A16*A22*A31*A43*A55+A12*A26*A31*A43*A55
     ++A16*A21*A32*A43*A55-A11*A26*A32*A43*A55-A12*A21*A36*A43*A55
     ++A11*A22*A36*A43*A55+A13*A22*A31*A46*A55-A12*A23*A31*A46*A55-
     +A13*A21*A32*A46*A55+A11*A23*A32*A46*A55+A12*A21*A33*A46*A55-
     +A11*A22*A33*A46*A55+A15*A23*A32*A41*A56-A13*A25*A32*A41*A56-
     +A15*A22*A33*A41*A56+A12*A25*A33*A41*A56+A13*A22*A35*A41*A56-
     +A12*A23*A35*A41*A56-A15*A23*A31*A42*A56+A13*A25*A31*A42*A56
     ++A15*A21*A33*A42*A56-A11*A25*A33*A42*A56-A13*A21*A35*A42*A56
     ++A11*A23*A35*A42*A56+A15*A22*A31*A43*A56-A12*A25*A31*A43*A56-
     +A15*A21*A32*A43*A56+A11*A25*A32*A43*A56+A12*A21*A35*A43*A56-
     +A11*A22*A35*A43*A56-A13*A22*A31*A45*A56+A12*A23*A31*A45*A56
     ++A13*A21*A32*A45*A56-A11*A23*A32*A45*A56-A12*A21*A33*A45*A56
     ++A11*A22*A33*A45*A56)*A64-(A16*A24*A33*A42*A51-
     +A14*A26*A33*A42*A51-A16*A23*A34*A42*A51+A13*A26*A34*A42*A51
     ++A14*A23*A36*A42*A51-A13*A24*A36*A42*A51-A16*A24*A32*A43*A51
     ++A14*A26*A32*A43*A51+A16*A22*A34*A43*A51-A12*A26*A34*A43*A51-
     +A14*A22*A36*A43*A51+A12*A24*A36*A43*A51+A16*A23*A32*A44*A51-
     +A13*A26*A32*A44*A51-A16*A22*A33*A44*A51+A12*A26*A33*A44*A51
     ++A13*A22*A36*A44*A51-A12*A23*A36*A44*A51-A14*A23*A32*A46*A51
     ++A13*A24*A32*A46*A51+A14*A22*A33*A46*A51-A12*A24*A33*A46*A51-
     +A13*A22*A34*A46*A51+A12*A23*A34*A46*A51-A16*A24*A33*A41*A52
     ++A14*A26*A33*A41*A52+A16*A23*A34*A41*A52-A13*A26*A34*A41*A52-
     +A14*A23*A36*A41*A52+A13*A24*A36*A41*A52+A16*A24*A31*A43*A52-
     +A14*A26*A31*A43*A52-A16*A21*A34*A43*A52+A11*A26*A34*A43*A52
     ++A14*A21*A36*A43*A52-A11*A24*A36*A43*A52-A16*A23*A31*A44*A52
     ++A13*A26*A31*A44*A52+A16*A21*A33*A44*A52-A11*A26*A33*A44*A52-
     +A13*A21*A36*A44*A52+A11*A23*A36*A44*A52+A14*A23*A31*A46*A52-
     +A13*A24*A31*A46*A52-A14*A21*A33*A46*A52+A11*A24*A33*A46*A52
     ++A13*A21*A34*A46*A52-A11*A23*A34*A46*A52+A16*A24*A32*A41*A53-
     +A14*A26*A32*A41*A53-A16*A22*A34*A41*A53+A12*A26*A34*A41*A53
     ++A14*A22*A36*A41*A53-A12*A24*A36*A41*A53-A16*A24*A31*A42*A53
     ++A14*A26*A31*A42*A53+A16*A21*A34*A42*A53-A11*A26*A34*A42*A53-
     +A14*A21*A36*A42*A53+A11*A24*A36*A42*A53+A16*A22*A31*A44*A53-
     +A12*A26*A31*A44*A53-A16*A21*A32*A44*A53+A11*A26*A32*A44*A53
     ++A12*A21*A36*A44*A53-A11*A22*A36*A44*A53-A14*A22*A31*A46*A53
     ++A12*A24*A31*A46*A53+A14*A21*A32*A46*A53-A11*A24*A32*A46*A53-
     +A12*A21*A34*A46*A53+A11*A22*A34*A46*A53-A16*A23*A32*A41*A54
     ++A13*A26*A32*A41*A54+A16*A22*A33*A41*A54-A12*A26*A33*A41*A54-
     +A13*A22*A36*A41*A54+A12*A23*A36*A41*A54+A16*A23*A31*A42*A54-
     +A13*A26*A31*A42*A54-A16*A21*A33*A42*A54+A11*A26*A33*A42*A54
     ++A13*A21*A36*A42*A54-A11*A23*A36*A42*A54-A16*A22*A31*A43*A54
     ++A12*A26*A31*A43*A54+A16*A21*A32*A43*A54-A11*A26*A32*A43*A54-
     +A12*A21*A36*A43*A54+A11*A22*A36*A43*A54+A13*A22*A31*A46*A54-
     +A12*A23*A31*A46*A54-A13*A21*A32*A46*A54+A11*A23*A32*A46*A54
     ++A12*A21*A33*A46*A54-A11*A22*A33*A46*A54+A14*A23*A32*A41*A56-
     +A13*A24*A32*A41*A56-A14*A22*A33*A41*A56+A12*A24*A33*A41*A56
     ++A13*A22*A34*A41*A56-A12*A23*A34*A41*A56-A14*A23*A31*A42*A56
     ++A13*A24*A31*A42*A56+A14*A21*A33*A42*A56-A11*A24*A33*A42*A56-
     +A13*A21*A34*A42*A56+A11*A23*A34*A42*A56+A14*A22*A31*A43*A56-
     +A12*A24*A31*A43*A56-A14*A21*A32*A43*A56+A11*A24*A32*A43*A56
     ++A12*A21*A34*A43*A56-A11*A22*A34*A43*A56-A13*A22*A31*A44*A56
     ++A12*A23*A31*A44*A56+A13*A21*A32*A44*A56-A11*A23*A32*A44*A56-
     +A12*A21*A33*A44*A56+A11*A22*A33*A44*A56)*A65
     ++(A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51
     ++A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-
     +A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-
     +A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51
     ++A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51
     ++A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-
     +A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-
     +A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-
     +A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-
     +A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52
     ++A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52
     ++A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-
     +A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-
     +A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52
     ++A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52
     ++A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52
     ++A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53
     ++A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-
     +A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-
     +A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53
     ++A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53
     ++A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-
     +A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-
     +A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-
     +A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-
     +A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54
     ++A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54
     ++A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-
     +A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-
     +A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54
     ++A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54
     ++A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54
     ++A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55
     ++A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-
     +A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-
     +A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55
     ++A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55
     ++A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-
     +A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-
     +A11*A23*A32*A44*A55-A12*A21*A33*A44*A55
     ++A11*A22*A33*A44*A55)*A66

      RETURN

      END FUNCTION M66DET


      DOUBLE PRECISION FUNCTION M44DET(A)
      
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN)  :: A
      DOUBLE PRECISION M44DET
      
      M44DET = (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*
     +(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-
     +A(3,3)*A(4,2)))- A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*
     +A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*
     +(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+ A(1,3)*(A(2,1)*(A(3,2)*
     +A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*
     +A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))- A(1,4)*
     +(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)
     +-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
      
      END FUNCTION M44DET




      subroutine eigenmad(inim,n,e,f,deth,r,aux2)
c-----------------------------------------------------------------------
c Given an initial symplectic matrix inim of dimension 4x4 finds out the 
c block diagonal matrix following the procedure of Edwards and Teng and 
c returns the two 2x2 block matrixs a and b. 
c-----------------------------------------------------------------------
      implicit none
      double precision inim(n,n),e(n/2,n/2),f(n/2,n/2)
      integer n,k,l
      double precision A(2,2),B(2,2),C(2,2),D(2,2),
     +BBAR(2,2),h(2,2),amind(2,2),ra(2,2),rm(4,4),
     +rmbar(4,4),rmi(4,4),test(4,4),br(2,2),rb(2,2),
     +r(2,2),rbar(2,2)
      double precision aux,aux2

      double precision toler,deth,phi,dtr,detra,zero
      double precision pi,factor,delta,gfact
      integer signdet,itmax,it,i,j,g,m

      parameter (toler = 1.0d-4)
      parameter (pi = 3.1415926535898d0)
      parameter (zero  = 0.0d0)

      if (n .ne. 4) then
         write(*,*) 'Error- The dimension ',n,
     +   'of the input matrix is not 4'
         stop
      endif
      
      call symcheck(inim,4,aux) 

      call matrzero(A,2,2,1)
      call matrzero(B,2,2,1)
      call matrzero(C,2,2,1)
      call matrzero(D,2,2,1)
      call matrzero(BBAR,2,2,1)
        

      A(1:2,1:2) = inim(1:2,1:2)
      D(1:2,1:2) = inim(3:4,3:4)
      B(1:2,1:2) = inim(1:2,3:4)
      C(1:2,1:2) = inim(3:4,1:2)

c Defining the simplectic conjugate matrix of nod following "https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.2.074001"
      call matbar(2,B,BBAR)

      do k = 1,2
         do l = 1,2
            h(k,l) = C(k,l) + BBAR(k,l) 
            amind(k,l) = A(k,l) - D(k,l)
         enddo
      enddo
         
      deth = h(1,1)*h(2,2) - h(2,1)*h(1,2)
      dtr = 0.5d0*(amind(1,1) + amind(2,2))
      
      Delta = deth + dtr**2

      if (Delta .lt. zero) then
         write(*,*) 'Unstable beam due to coupling' 
         stop
      else
         if (Delta .eq. zero) then
            write(*,*) 'Delta = 0'
            stop
            do k = 1,2
               do l=1,2
                 if (k .eq. l) ra(k,l) = 1.0d0
                 if (k .ne. l) ra(k,l) = 0.0d0
               enddo
            enddo
           
         else
            factor = -(dtr + sign(sqrt(Delta),dtr))
            ra = h/factor
         endif
      endif
 
      detra = ra(1,1)*ra(2,2) - ra(2,1)*ra(1,2)
      gfact = 1.0d0/sqrt(1+detra)
      
c Define the two boxes E and F
      call mxmpy(B, Ra,BR, 2, 2)
      e = A - BR	
      call mxmpy(Ra, B,RB, 2, 2)
      f = D + RB

c Check that the diagonalization in boxes has been correctly done
c diag(M) = Rm^-1*M*Rm = g**2 [I -Rbar] * [A B] * [I Rbar] = [E 0]
c                             [R   I  ]   [C D]   [-R  I ]   [0 F]
      r = ra

      call matbar(2,r,rbar)
      call matrzero(rm,4,4,1)
      do i = 1,4
         rm(i,i) = 1.0d0
      enddo
      rm(3:4,1:2) = r(1:2,1:2)
      rm(1:2,3:4) = -rbar(1:2,1:2)
      call matbar(4,rm,rmbar)
      call mxmpy(rm,inim,rmi,4,4)
      call mxmpy(rmi,rmbar,test,4,4)
      test = gfact**2*test
      


c Check that the resulting matrix is block diagonal  
      aux2 = 0    
      do k = 1,2
         do l = 3,4
            aux2 = max(aux2,test(k,l),test(l,k))
       	    if ((test(k,l) .gt. toler) .or. (test(l,k) .gt. toler))then
               write(*,*) 'Off diagonal block matrix different',
     +         ' from 0'
                     write(*,*) 'test'
                     write(*,'(1p,4e12.4)') test
                     write(*,*) 'e'
                     write(*,'(1p,2e12.4)') e
                     write(*,*) 'f' 
                     write(*,'(1p,2e12.4)') f
c                     stop
            endif
         enddo
      enddo
 
      end



      subroutine matbar(n,a,abar)
      implicit none
      !----------------------------------------------------------------------*
      ! Returns the symplectic conjugate of a symplectic matrix ABAR = -SMtS *
      ! Input:                                                               *
      ! 	a: Squared symplectic matrix                                 *
      !         n: dimension of a, where n is even                           *
      ! output:                                                              *
      ! 	abar: symplectic conjugate of a                              *
      !----------------------------------------------------------------------*
      double precision :: a(n,n), abar(n,n),aux(n,n),aux2(n,n)
      double precision :: s(n,n)
      integer n

      if (n .eq. 2) then
         s = reshape((/0.,1.,-1.,0./),shape(s))
         
      elseif (n .eq. 4) then
         s = reshape((/0.,1.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,
     +1.,0.,0.,-1.,0./),shape(s))
          
      else 
         
         s = reshape((/0.,1.,0.,0.,0.,0.,
     +-1.,0.,0.,0.,0.,0.,
     +0.,0.,0.,+1.,0.,0.,
     +0.,0.,-1.,0.,0.,0.,
     +0.,0.,0.,0.,0.,0.,1.,
     +0.,0.,0.,0.,-1.,0./),shape(s)) 

      endif

      call mxmpy(s, transpose(a),aux, n, n)
      call mxmpy(aux,s,aux2,n,n)
      abar = -aux2
      
      end subroutine matbar



      subroutine lnrcoup(deth,qx,qy,dqmin)
      implicit none
      double precision deth,qx,qy,dqmin
      double precision pi
      parameter (pi = 3.1415926535898d0)

      dqmin = sqrt(abs(deth))/(pi*(sin(2*pi*qx) + sin(2*pi*qy)))
      end subroutine lnrcoup


      subroutine getmaskmap(name,mapname,nname,nmap,mask)
      implicit none
      !----------------------------------------------------------------------*
      ! Creates a mask of length nmap (nmap > nname) signaling the order of  *
      ! appearance of an element of name in mapname. If the element is just  *
      ! present on mapname, mask presents a 0.                               *
      ! Input:                                                               *
      ! 	name: array of dimension nname                               *
      ! 	mapname: array of dimension nmap                             *
      ! Output:                                                              *
      ! 	mask:array of dimension nmap                                 *
      !----------------------------------------------------------------------* 
     
        
      integer nname,nmap,maxlength
      parameter (maxlength = 9)
      character*(maxlength) name(nname),mapname(*),name2(nname)
      integer mask(nmap)
      integer nampos,lastnb
      integer i,j
      character*10 sname
      
      
      call dzero(mask,nmap) ! Initialize mask

      do i = 1, nname
         if (len(name(i)) .le. maxlength) then
            name2(i) = name(i)(1:lastnb(name(i)))
         else
            name2(i) = name(i)(1:maxlength)
         endif
      enddo

      do i = 1, nmap
         sname = mapname(i)(1:lastnb(mapname(i)))
         j = nampos(sname,name2,nname)
         if (j .gt. 0) mask(i) = 1
      enddo
      return

      end subroutine
      
      
      integer function rotate(ipos)
      implicit none
      !----------------------------------------------------------------------*
      ! Chooses as reference the leftmost parasitic interaction of the first *
      ! interaction point appearing in the maps.This way we can match the    *
      ! interaction number assocated to the optical functions and the maps   *
      ! without the need of always introducing IP5.                          * 
      !----------------------------------------------------------------------*
       
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer ipit,ncoll
      integer actlist
      integer ipos(8)
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/

      write(*,*) 'Reference: IP',actlist(1)
      rotate = ipos(actlist(1)) - npar(actlist(1))

      end function

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function eqvoptmap(nopt)
      implicit none
      integer nopt,l,i

      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision clight,ech,emass,erad,pi,pmass,prad
c-----------------------------------------------------------------------
c electron:
c   classical radius [m]:
      parameter         (erad   = 2.817 940 92 d-15)
c   rest mass [GeV]:
      parameter         (emass  = 0.510 999 06 d-03)
 
c proton:
c   classical radius [m]:
      parameter         (prad   = 1.534 698 57 d-18)
c   rest mass [GeV]:
      parameter         (pmass  = 0.938 272 31 d+00)
 
c elementary charge:
      parameter         (ech    = 1.602 189 2  d-19)
 
c velocity of light:
      parameter         (clight = 2.997 924 58 d+08)
c pi:
      parameter         (pi     = 3.1415926535898d0)
      integer ipit,ncoll
      integer actlist
      integer ipos(8)
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/ 
      l = 0
      do i = 1, nmaps
         if (mapmask(i) .ne. 0) l = l + 1
         if (l .eq. nopt) eqvoptmap = i
      enddo
      return
      end function

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prextra(fhead, header,bd,z,index,scale)
      implicit none
      integer i,index,k,l,bd,lastnb
      real z
      double precision scale,zdif,zline,zmax,zmin,zsum,muline
c-----------------------------------------------------------------------
c print value table per bunch and all extra elements for one beam
c the values are scaled by "scale"
c-----------------------------------------------------------------------
 
      integer mdim,mlocal,mpit,mvary,melm,mmaps, msect
      integer orbout, mulist, mucoll, msep, maxsequ, mcnam, ustart,
     +mtrack, max_list, lumilist, iunit
c-----------------------------------------------------------------------
      integer mbuck, mbunch, mdslt, mpar, mcol
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
c maximum number of collisions per pit:
      parameter (mlocal = 2 * mpar + 1)
c maximum number of phase space dimensions
      parameter (mdim = 4)
c maximum number of variables
      parameter (mvary = mdim * mbunch * 2)
c maximum number of observed bunches during tracking
      parameter (max_list = 10)
c input/output units
      parameter (iunit=11, orbout = 22, mulist = 23, mucoll = 24,
     +msep = 25, lumilist=26,mtrack = 30)
      parameter (ustart = 50)
c various array sizes etc.
      parameter (maxsequ = 20000, mcnam = 16, msect = 259)
      double precision zero, one, two, three, ten, ten3m, ten9m, toler
      double precision half
      double precision ten3p,ten6p
      parameter (zero  = 0.0d0)
      parameter (one   = 1.0d0)
      parameter (two   = 2.0d0)
      parameter (three = 3.0d0)
      parameter (ten   = 10.d0)
      parameter (half = 0.5d0)
      parameter (ten3m = 1.0d-3, ten9m = 1.0d-9)
      parameter (toler = 1.0d-8)
      parameter (ten3p = 1.0d3, ten6p = 1.0d6)
      double precision bcurr1,bcurr2,qx1,qx2,qxp1,qxp2,qy1,qy2,
     +qyp1,qyp2, q11,q12,q21,q22
c-----------------------------------------------------------------------
c external code for bunches
      common /buncha/ code(0:mbunch)
      character*4 code
      save /buncha/
c number of particles per bunch
      common /bunchf/ bcurr1(mbunch), bcurr2(mbunch),
     +     qx1(mbunch), qy1(mbunch), qx2(mbunch), qy2(mbunch),
     +     qxp1(mbunch), qyp1(mbunch), qxp2(mbunch), qyp2(mbunch),
     +     q11(mbunch), q21(mbunch), q12(mbunch), q22(mbunch)
      save /bunchf/
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      double precision arad,bcurr,circum,deltap,epsx0,epsy0,frev,gamma,
     +gev,partno, ampx, ampy, sigb, sigem, tmass, tradius, xisign,
     +xifact,hofact,ippos,root2,lumicnt, lumiav,lumifact
c-----------------------------------------------------------------------
c global counters
      common /globa/ title, type, date, hour, timew
      character title*80, type*16, date*10, hour*10, timew*8
      save /globa/
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/
      common /globf/ epsx0, epsy0, deltap, gev, bcurr,ampx(2),ampy(2),
     +sigb, sigem, gamma, arad, partno, frev, circum, tmass, tradius,
     +xisign,xifact,hofact,root2,lumicnt,lumiav,lumifact,ippos(8)
      save /globf/
      common /globc/ seq_name(2),  para_names(mcol)
      save /globc/
      character*(mcnam) seq_name,  para_names
c     number of pits:                     npit
c     total number of bunches:            nbunch
c     collision points per pit:           nlocal
c     number or interaction points:       ninter
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk
      double precision ztr
c   ztr trajectory
      common / mtcommd / ztr(6,mbunch,2)
      real orb_amp
c   initial orbit amplitude ((x=1,y=2),bunch,ring)
      common / mtcommr / orb_amp(2,mbunch,2)
      save /mtcomm/, /mtcommd/, /mtcommr/
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
      integer ipit,ncoll
      double precision si
c-----------------------------------------------------------------------
c pit azimuths and association of collision points with pits
      common /pitc/ pitnam(mpit)
      save /pitc/
      character*4 pitnam
      common /pitf/ si(mpit)
      save /pitf/
      integer actlist
      common /piti/ ncoll(mpit), ipit(mpit),actlist(mpit)
      save /piti/
      common / flagsi/ debug, c_tunes, beamc_f, nturns, outbcnt,
     +                 outpos, outnorm, emitt_f, outblist(max_list)
      save /flagsi/
      integer debug, c_tunes, beamc_f, nturns, outbcnt, outpos, outnorm,
     +emitt_f, outblist
      common /flagsl/ bcfile, w_coll, w_frequ, w_equ, w_set,
     +w_alt, c_orbit, f_coll, f_second, w_detail, all_write, lumi_hist
      logical bcfile, c_orbit, f_second, w_detail, all_write, lumi_hist
      logical f_coll, w_coll, w_frequ, w_equ, w_set, w_alt
      save /flagsl/ 
c------------------------------------------------------------------------
      double precision qx,qy,q
      common /globtun/ qx,qy
      save/ globtun/

      double precision kick
      common/ kck/ kick(2,mbunch,mcol,2)
      save /kck/
      
      double precision orbdev
      double precision betaline, beta
      double precision kicktemp
      double precision templist(nbunch)
      double precision diflist(mcol)
      double precision temp(2,mbunch,2)
      common /tempc/ temp
      save /tempc/
c-----------------------------------------------------------------------
      character*(*) fhead, header
      character * 8  dir
      character * 1000 temp1, temp2, temp3
      dimension z(6,mbunch,0:mmaps+1)
      dimension zline(mbunch), zmax(mbunch), zmin(mbunch), zsum(mbunch)
      dimension zdif(mbunch),muline(mbunch),betaline(mbunch),
     +beta(mbunch)

      character*120 filename, chst,filenm
      common /filesch/ filenm,chst
      save / filesch/ 
      integer nstart,filestart,j
      character*4 kstr
      double precision difmatr(1,nbunch,nbunch),
     +tunematr(2,nbunch,nbunch)
      double precision chromatr(2,nbunch,nbunch)
      double precision enj,enl,betanew,sigmaj,sigmal,sigmapj,sigmapl,xj,
     +xl,xpj,xpl,ql,qj,qpl,qpj,r,phi
      character*3 part
      double precision dmatmax,rms
      double precision difarray(nbunch*(nbunch-1)/2)
      integer diflen,exptemp 
      real(kind = 8) :: cum,cumn,cumtemp,maxtol
      logical cumflag

      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/   
      integer exstart,nextra
      integer iextra(nmaps)
      
c-------------------------------------------------------------------------
       exstart = 250
       if (bd .eq. 1)  then
        dir = 'forward'
        temp1 = fhead(:lastnb(fhead)) // '_f'

      else
        dir = 'backward'
        temp1 = fhead(:lastnb(fhead)) // '_b'
      endif
     
      nextra  =  0
      do i = 0, nmaps-1
         if (mapmask(i+1) .eq.0 ) then
            nextra = nextra + 1
          
            iextra(nextra) = i  
         endif
      enddo
      write (mulist, 910) header, dir
      write (mulist, 920)
      write (mulist, 930) (' ', mname(iextra(i)), i = 1, nextra), '  '
      write (mulist, 920) 
      
      do i = 1, nextra
         zmax(i) =scale*z(index,1,iextra(i))
         zmin(i) = zmax(i)
         zsum(i) = 0.0
      enddo
      
      do i = 1,nextra
         temp2 =  temp1(:lastnb(temp1)) // '.' // mname(iextra(i))
         open(exstart+i, file = temp2)
      enddo

      do l = 1, nbunch
         do i = 1, nextra
            zline(i) = scale*z(index,l,iextra(i))
            zmax(i) = max(zmax(i),zline(i))
            zmin(i) = min(zmin(i),zline(i))
            zsum(i) = zsum(i) + zline(i)
         enddo 
         
         if (bd .eq. 1)  then
           k = set_f(l)
         else
           k = set_b(l)
         endif
         
         if (l .eq. 1 .or. l .eq. nbunch)  then
           write (mulist, '(a5, i5, 16('' '',a,f12.6),''   '')')
     +     code(l), k, (' ', zline(i), i = 1, nextra), '  '
         endif
         do i = 1, nextra
 
           write(exstart+i, *) k, zline(i)
         enddo
      enddo 
      
         
       
      do 40 i = 1, nextra
         close(exstart+i)
         zsum(i) = zsum(i) / nbunch
         zdif(i) = zmax(i) - zmin(i)
 40   continue

      write (mulist, 920)
      write (mulist, 940) 'mx  ', (' ', zmax(i), i = 1, npit), '  '
      write (mulist, 940) 'mn  ', (' ', zmin(i), i = 1, npit), '  '
      write (mulist, 940) 'av  ', (' ', zsum(i), i = 1, npit), '  '
      write (mulist, 940) 'df  ', (' ', zdif(i), i = 1, npit), '  '
      write (mulist, 920)
c-----------------------------------------------------------------------
 910  format(' '/a,' for ',a,' beam:')
 920  format(' ')
 930  format('bunch bucket',16(' ',a,5x,a3,4x))
 940  format(a5,16(' ',a,f12.6),'   ')
c-----------------------------------------------------------------------
      end
      


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pravdisp(fhead,scale)
      implicit none
c-----------------------------------------------------------------------
c Print value table of the average dispersion of all bunches along the 
c ring for beam 1 and 2 once the beam-beam interactions have been 
c introduced.
c The values are scaled by "scale"
c-----------------------------------------------------------------------
      integer mpar,mpit,mcol,melm,mmaps,max_list 
      integer mbuck, mbunch, mdslt
      integer mcnam
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck)
c max. number of parasitic crossings on each side of IP
      parameter (mpar = 50)
c maximum number of pits
      parameter (mpit = 4)
c maximum of h.o. + parasitic
      parameter (mcol = 2 * mpar * mpit + mpit)
c maximum number of extra elements per one side interaction point
      parameter (melm = 10)
c maximum number of maps
      parameter(mmaps = mcol + 2*mpit*melm)
      parameter(mcnam = 16)
      parameter (max_list = 10)
      character*(*) fhead
      double precision scale 	
      double precision betx,bety,delta,dx,dy,s,epsx,epsy,x,xmu,y,ymu,
     +eiv1, eiv2, orb0_1, orb0_2, alfx, alfy
c-----------------------------------------------------------------------
c description of interaction points
      common /optica/ name(mcol,2)
      save /optica/
      character*(mcnam)    name
      common / optici / occur(mcol,2)
      save /optici/
      integer occur
      common /opticf/ delta, s(mcol,2),
     +   x(mcol,2), dx(mcol,2), betx(mcol,2), xmu(mcol,2),
     +   epsx(mbunch,2), y(mcol,2), dy(mcol,2), bety(mcol,2),
     +   ymu(mcol,2), epsy(mbunch,2), alfx(mcol,2), alfy(mcol,2),
     +   eiv1(6,6,max_list), eiv2(6,6,max_list),
     +   orb0_1(6,max_list), orb0_2(6,max_list)
      save /opticf/
            real z1,z1a,z1b,z2,z2a,z2b,d1,d2,dd1,dd2
c-----------------------------------------------------------------------
c phase space coordinates at start and end of system and in collisions
c dispersion per bunch and interaction point
      common /corbit/  z1(6,mbunch,0:mmaps+1), z2(6,mbunch,0:mmaps+1),
     +               z1a(6,mbunch,mmaps), z1b(6,mbunch,mmaps),
     +               z2a(6,mbunch,mmaps), z2b(6,mbunch,mmaps),
     +               dd1(6,mbunch,0:mmaps+1), dd2(6,mbunch,0:mmaps+1),
     +               d1(6,mbunch,0:mmaps+1), d2(6,mbunch,0:mmaps+1) 
       
      integer nbunch,ninter,npar(8),nlocal,npit,iseed,iact,c_turn,
     +amp_bunch, amp_fac, b2_off,n_parasit
      common /globi/ npit, nbunch, nlocal, ninter, npar, c_turn, b2_off,
     +iseed,amp_bunch,amp_fac,n_parasit,iact(8)
      save /globi/  
 
      common /extel/ nmaps, extraelem,mapmask,mname
      integer nmaps,mapmask(mmaps)
      logical extraelem
      character*12 mname(mmaps)
      save /extel/  

      character * 1000 temp
      integer i,j,k,outunit,lastnb
      double precision tempx1,tempx2,tempy1,tempy2

      outunit = 713       
      temp = fhead(:lastnb(fhead))  
      open(outunit, file = temp)

      j = 0
 	
      do i = 1, nmaps
         if (mapmask(i) .gt.0 ) then
            tempx1 = 0
            tempy1 = 0
            tempx2 = 0
            tempy2 = 0
            j = j + 1
            do k = 1, nbunch
               tempx1 = tempx1 + d1(1,k,i)
               tempy1 = tempy1 + d1(3,k,i)
               tempx2 = tempx2 + d2(1,k,i)
               tempy2 = tempy2 + d2(3,k,i)  
            enddo  
            tempx1 = tempx1/nbunch
            tempy1 = tempy1/nbunch
            tempx2 = tempx2/nbunch
            tempy2 = tempy2/nbunch

            write(outunit,*) mname(i), s(j,1), tempx1,tempy1,tempx2,tempy2 
         endif
      enddo
      close(outunit)
      end 

      integer function nstart()
      implicit none
c-----------------------------------------------------------------------
c Returns the number of the first regular bunch, i.e., the first bunch of  
c the filling scheme symmetric for both beams 1 and 2.
c Output:
c	nstart	--number of the first regular bunch
c----------------------------------------------------------------------- 
      integer mbuck, mbunch, mdslt
c number of bunch slots around the machine
      parameter (mbuck = 3564, mbunch = 3000, mdslt = 2 * mbuck) 
      common / mtcomm /
     +colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +        ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +        set_f(mbuck), set_b(mbuck), tcount(mbuck,2),
     +        collsk(2,0:mbuck-1)
      integer colpnt, ntotal_f, ctotal_f, nset_f, equl_f, ordl_f,
     +                ntotal_b, ctotal_b, nset_b, equl_b, ordl_b,
     +                set_f, set_b, tcount, collsk 
      integer nstart,dif,cnteq,i,top,j
      cnteq = 0
      do i = 0, mbuck-1
         dif = abs(collsk(2,i) - collsk(1,i))
         if (dif .eq. 0) cnteq = cnteq + 1
         if (dif .ne. 0)  cnteq = 0 
         if (dif .eq. 0 .and. collsk(1,i) .eq. 1 .and. 
     +    cnteq .gt. 10) then
            top = i
            nstart = 0
            do j = 0,top
              nstart = nstart + collsk(1,j)
            enddo 
            return
         endif   
      enddo
      end function