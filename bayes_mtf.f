! program bayes_mtf
!
! Copyright 2022-2026 Bernard Delley
!  Licensed under the Apache License, Version 2.0
!
! Author: B. Delley 2022-2026
!
!  B. Delley and Y.L. Delley, "Bayesian approach to slanted edge modulation transfer function (MTF) measurements",
!  J. Opt. Soc. Am. A 43, xx  (2026)
!
! bayes_mtf analyzes slanted edge image files of linear P2 type GrayMap files (non de-Bayered) and with black level 0.
! a pgm maximum value up to 16383 is interpreted as a command to extract Bayer green from RGGB pixels.
! a higher maximum value is a command to consider the file as monochrome image.
!
! simple use: (parameters may be default or overridden by supplied by first # comment line in pgm file)
!   bayes_mtf pix.pgm
!  call arguments on first pgm comment line will superseed plain defaults, used for most example on github site

! typical use with two principal parameters specified: lens aperture N_f and and sensor specific cgain
!   bayes_mtf pix.pgm  1.4  3.92

! full use of specific input (overriding presets), as shown in bayes_mtf output, for example
! bayes_mtf  pix_sim_f4_angle17p6_A43.pgm  4.0  3.92 0.0 0.000 0.000  4.345  537.  64 16 0 2 1  64  0.87  4.0


! input argument values with default setting,    see our JOSAA article for discussion of parameters
!     file            ! filename of pgm type P2 ascii image file of slanted edge, see example files
!     anum = 5.6      ! N_f, actual lens stop value, sets frequency limit for MTF
!     cgain = 3.92    ! conversion factor DN->electrons: fundamental Poisson noise, sigma_el=sqrt(N_el)
!     f0 = 0.         ! sigma_0   read noise
!     fp = 0.         ! c2        pixel response non-uniformity
!     fh = 0.         ! sigma_e   (blade) slanted edge roughness
!     pix  = 4.345    ! pixel pitch (in mm internally), input in [mum], needed for correct c/mm frequencies
!     alam = 537.     ! nominal light wavelength in [nm]
!     np   = 64       ! n_p  ROI normal size, base harmonic frequency = 1/n_p
!     naqs = 16       ! m_Y  'resolution' parameter for spline envelope model Y,
!                     !    -1 -> full harmonic model y
!     nes  = 3        ! 1-2 illumination gradient normal, gradient parallel parameters
!                     ! 3-n  +edge_curvature splines 3:order 3, 1 internal free parameter ~ "curvature"
!     ity  = 5        ! 5: is standard: with 3 iterations for Stage_3 yr,yi and Yr,Yi
!                     !    2: do only real part, stop with Stage_2
!     mprbay = 1      ! printlevel
!     nquer = 0       ! -> = np  ROI parallel size, target value
!     apixm = 0       ! pixel aperture for reference OTF. when apixm = 0, then no OTF_deviation_magnified.eps 
                      ! apixm is used in Stage 1, to generate optimal Phi, Dark and Illum_jump to start Stage 2
                      ! pixel aperture is not a parameter in Stage 2,3, 
                      ! an ad-hoc reference OTF gives a hint the sensor pparameter apix
!     anumm = 0       ! -> anum, N_f of reference OTF

! special call for basic cropping and black level removal mode for pgm type P2 files. 
!  creates pix_crop.pgm file and exits:
!   bayes_mtf pix.pgm  crop  iw1 iw2 ih1 ih2 iblack

! special calls for MTF map use and debugging (hindered input params can be passed in via pgm comment-input line)
!  these calls have special defaults anum=8 naqs=12 and ity=10, a convergence threshold limits the number of iterations
!   bayes_mtf pix.pgm map
!   bayes_mtf pix.pgm skipto ktrps kedgs
!   map produces a g_map.eps file showing numbered trapezoids with color coded total system MTFhNyq at each edge.
!   MTFhNyq is MTF at half Nyqvist frequency for (10) Miller direction. 
!   color code is as in: B. Delley and F. van den Bergh,
!   "Fast full-field modulation transfer function analysis for photographic lens quality assessment",
!   Applied Optics. 60 (8). 2197-2206 (2021)

      module image_data
        parameter( mn=65536, nmm=-0)
        parameter( ithsq = 32, itrpe=32 )
        parameter( mprtra=0 ) ! 0,1,2,3,4,5,9

        integer, allocatable :: ipix(:,:), ipxy(:,:), jpix(:,:) !,idelt(:,:)
        integer iw,ih,imdn, lenf
        integer iw1,iw2,ih1,ih2, ktrps, kedgs
        integer ktrpm, monochrom, nnull
        real    gain,rn,darkthres
        character*2 ptype
        character*128 filen
        logical keep_roi, skipto
      end module image_data

      module esf_data
        parameter( ns=32)
        parameter( mprplo=9 ) ! 4,5 mprbay 1,3,5,9 

        parameter( thed = 16. )               ! trapezoid minimum parameter
        parameter( throutly=3.0 )  ! outlyer threshold 3 sigma of pixel noise according hypothesis

! Z8 raw 8280x5520 36x24 4.3478 4.3478   8256x5504  35.9x23.9 4.3484 4.3423 
! Z7     8288x5520       4.3436 4.3478   8256x5504  35.9x23.9  avg= 4.3454
! D850   8288x5520        avg= 4.3457    8256x5504  35.9x23.9                MY_AVG_NOW 4.345  for all

        integer np, na, n1, n2, nm, iroin, iroip, naqs, nes, nquer, nes1
        integer ncd, npoly, npolym, nparm1, ity, mprbay
        real    pix, bpix, alam, cgain, f0, fh, fp
        real    anum, dcut, fcut, anuf, roip
        real    y00 ,esfh, xm50, xhn
        real    one, pi, pi2, pi4, tpi, pii, tpin, pre, pi2deg, epsma
        real, allocatable :: si(:),ci(:)              ! sine cosine tables on oversampled grid
        real, allocatable :: circ(:),papt(:),dcirc(:) ! mtf circular aperture, pixel aperture, derivative circ
        real, allocatable :: yr(:), yi(:)             ! 1 d OTF real and imag part
        real, allocatable :: esf(:),csf(:)            ! ESF, LSF
        real, allocatable :: dsf(:,:,:), fsf(:,:,:)   ! derivatives LSF and ESF
        real, allocatable :: xk2(:,:)                 ! data store for plots later
        real, allocatable :: edge(:,:,:)              ! edge data
        real, allocatable :: vlp(:), pl(:)            ! edge curvature spline parameters now
        real, allocatable :: xes(:), yes(:), y2es(:,:)! edge spline 
        real, allocatable :: cpmm(:),cpix(:)    
        real, allocatable :: afn(:),bfn(:),cfn(:),dfn(:)
        real, allocatable :: yrfm(:)
        real    phi,pos0,fil(4),ve(2),vm(2),vn(2)
        real    phi0,fil0(2),vm0(2),ve0(2),vn0(2)
!       real    ymap(0:63),yrfmap(0:63)           !umapper
        logical debug
      end module esf_data

      module fit_data
        logical, parameter :: umapper=.false.
        logical, parameter :: skip1 = .false.
        integer nwh,npar,nparm,ndat,npar1,npar2,npar3,ihys,nsing
        real, allocatable :: udes(:,:), wdes(:,:), vrk(:,:), cvm(:,:)
        real, allocatable :: rmsi(:), brh(:), awr(:), sngv(:), rv1(:) 
        real, allocatable :: sig(:), bh(:), vdes(:,:), crh(:),drh(:)
        real, allocatable :: rk2(:,:),xk1(:,:) ,erh(:)
        real, allocatable :: xq(:),yqr(:),yqi(:)
        real, allocatable :: bqfn(:,:)
        real     sumry(5)
        real     dyr1,dyrn,dyi1,dyin,scalan,scale1,scalen 
        integer, allocatable :: ke1(:)
        integer  isumry(2)
        integer  nqm,nqm0,nqm1,nqm2
        integer  npx,  npk,   npg,   nph,   npc   ! need to be set correctly for postprob ! 4  5  5->4 6->5 7->6
        logical  uxno, uanum, ugrad, ugnor, ucurv ! need to be set correctly for postprob  
        logical  uimag, udark       ! uimag use with splines, uqspl=.true.
        logical  onlyi,uphi,usigd          ! phi now always used "3",  onlyi   not used now
        logical  uqspl,unats,domap
      end module fit_data

      module plot_data_compare
        real    anumm,phim,apixm
      end module plot_data_compare

      module plot_data
        parameter( xlen= 150., ylen= 100.)
        character*4 xfmt,yfmt
        dimension xtick(21),ytick(21),xval(21),yval(21)
        dimension da(5),db(5)
        integer nx,ny
        logical  dopos
      end module plot_data

      program bayes_mtf
! bayes mtf  main routine

      use esf_data
      use image_data
      use fit_data
      use plot_data_compare
      use plot_data , only : dopos
      character*20 word
      logical edgeok,exst,do_out_sfr

      write(6,'(/2a)')'Bayesian Approach to Slanted Edge',
     I ' OTF and MTF Measurements, version 260701'
      write(6,'(/a)')'Cite work using this program as:'
      write(6,'(a)')
     I '  B. Delley and Y.L. Delley, J. Opt. Soc. Am. A 43, xx  (2026)'
      write(6,'(a)')

      call cpu_time(t0)      !   f95 intrinsic function
      call get_input

      call atrpz

      call preptable

      do_out_sfr=.false.
      if(uqspl) then
         call setqsplin
         if(np.eq.64) do_out_sfr=.true.
         if(do_out_sfr) open(36,file='edge_sfr_values.txt')
      else
!       write(6,'(/a,i2,2i4,a,f6.1)')'#_OTF_model_y [',
        write(6,'(/a,i2,2i4,a,f6.1)')'OTF_model_y [',
     I  nes,np,2*iroip,']'
      endif

!       write(6,'(/2a)')'      trapezoid     edge_pos      len',
        write(6,'(/2a)')'    trapezoid     edge_pos      len',
     I  '     phi     psi   2*roip'

        do ktrp= 1,ktrpm 
         if(ktrps.eq.0 .or. (ktrps.gt.0 .and. ktrp.eq.ktrps)) then
          do kedg= 1,4
           if(kedgs.eq.0 .or. (kedgs.gt.0 .and. kedg.eq.kedgs)) then

            call ck_edge( ktrp, kedg, edgeok )

            anor =  edge(3,kedg,ktrp)  ! length of edge
!           roip = min( anor*0.4, float(iroip))         ! needed back in for chart
            roip = min( anor*0.5-6, float(iroip))       ! now just 5 spare pixels to each side for chart
!           roip = float(iroip)                         ! out to frame OK for blade

!            if(mprtra.gt.1) write(6,'(a)')'one-line-summary:'
             phd = phi*pi2deg
             phis = min( abs(phd), abs(phd-90), abs(phd-180), 
     I              abs(phd-270), abs(phd-360))
             
!      write(6,'(a,L2,i4,i2,6f8.2,2x,5f9.4,2i5)')'#_edge',edgeok
       write(6,'(a,L2,i4,i2,6f8.2,2x,5f9.4,2i5)')'edge',edgeok
     I  ,ktrp,kedg,(edge(i,kedg,ktrp),i=1,5),2*roip
!    I ,vm,edge(3,kedg,ktrp),phd,phis,edge(5,kedg,ktrp) !,sumry,isumry

            if(edgeok) then
              call analy2( ktrp, kedg )
              if(mprplo.gt.4) call mkplot0( ktrp, kedg )
              call get_mtf
              phd=phi*pi2deg
              if(do_out_sfr) call out_sfr( ktrp, kedg, phd )
              edge(8,kedg,ktrp) = xhn
            endif

           endif ! kedgs
          enddo
         endif ! ktrps
        enddo

        if(domap) call mkplotm

        close(36)
      call cpu_time(tt)      
      write(6,'(a,f12.1,a,2f20.9)')'CPU time',tt-t0,' s'

      write(6,'(a)')'done'
      call exit(0)

      end

      subroutine get_input

! Author: Bernard Delley 2026

! input:
!  pix file  ascii pgm graphics file type P2 with slanted edge raw image  
!    im max value < 16384 : extract green from Bayer RGGB
!    im max value > 16383 : interpret as monochrome image (used in example as 8x pixelshift green)

! Nikon camera file uncomressed or loss-less-compressed raw   pre-processing for input here
! dcraw -D -4 -c DSC_????.NEF | pnmtoplainpnm > full_size_file_P2.pgm
!  the full_size_file needs to be cropped and black level subtracted (400 for D500,D850, 1008 for Z7,Z8,Z9 etc)

! args
      use esf_data
      use image_data
      use fit_data
      use plot_data_compare
      use plot_data , only : dopos
      character*10 word
      character*30 fileg
      logical exst

! set default input parameter values (cgain Z8 (Z9?), pix alam D850,Z7,Z8,Z9) 
      anum = 5.6      ! N_f 
      cgain = 3.92    ! conversion factor DN->electrons: fundamental Poisson noise, sigma_el=sqrt(N_el) 
      f0 = 0.         ! sigma_0   read noise
      fp = 0.         ! c2        pixel response non-uniformity
      fh = 0.         ! sigma_e   (blade) slanted edge roughness
      pix  = 4.345    ! pixelpitch (in mm internally), input in [mum]
      alam = 537.     ! nominal light wavelength in [nm]
      np   = 64       ! n_p  ROI normal size, base harmonic frequency = 1/n_p
      naqs = 16       ! m_Y  'resolution' parameter for spline envelope model Y, 
                      !    -1 -> full harmonic model y
      nes  = 3        ! 1-2 illumination gradient normal, gradient parallel parameters 
                      ! 3-n  +edge_curvature splines 3:order 3, 1 internal free parameter ~ "curvature"  
      ity  = 5        ! 5: is standard 3 iterations for Stage_3 yr,yi and Yr,Yi   
                      !    2: do only real part, stop with Stage_2
      mprbay = 1      ! printlevel
      nquer = 0       ! -> = np  ROI parallel target value
      apixm = 0       ! pixel aperture for reference OTF, no OTF_deviation_magnified plot, when apixm = 0
      anumm = 0       ! -> anum, N_f of reference OTF

      uqspl = .true.
      unats = .true.


      call setcns
      dopos = .false.

      call getarg(1,filen)   ! arg filename 
      inquire(file=filen,exist=exst)
      if(exst) then
      open(15,file=filen)
      else
        write(6,'(9a)')'Error: file ',filen(1:lenf),' does not exist'
        stop 'error image file missing'
      endif

! clean up previous graphical output, so one never gets fooled by leftover stuff when a job has failed
      fileg='g_ROI_orientation.eps'
      inquire(file=fileg,exist=exst)
      if(exst) then
      open(95,file=fileg)
      close(95,status='delete')
      endif
      fileg='g_MTF.eps'
      inquire(file=fileg,exist=exst)
      if(exst) then
      open(95,file=fileg)
      close(95,status='delete')
      endif
      fileg='g_MTF_OTF.eps'
      inquire(file=fileg,exist=exst)
      if(exst) then
      open(95,file=fileg)
      close(95,status='delete')
      endif
      fileg='g_OTF_deviation_magnified.eps'
      inquire(file=fileg,exist=exst)
      if(exst) then
      open(95,file=fileg)
      close(95,status='delete')
      endif
      fileg='g_edge_function.eps'
      inquire(file=fileg,exist=exst)
      if(exst) then
      open(95,file=fileg)
      close(95,status='delete')
      endif
! for the other file overwriting is good enough: 
! g_ESF_bright_side.eps  g_ESF_dark_side.eps  g_ESF.eps  g_ESF_magnified.eps  g_LSF.eps g_ROI_map_outliers.eps

      lenf = len_trim(filen)
      call read_image    ! P2 pgm linear ASCII file, including its supplied parameter values, if any

      if(skipto) then ! deal with special branch
      call getarg(3,word)    ! arg    ktrps
      read(word,*,iostat=ierr) inp
      if(inp.gt.0 .and. ierr.eq.0) then
        ktrps = inp
      else
       write(6,'(a,2i5,3a)')
     I 'Error: skipto special branch requires ktrps>0',ktrps,ierr
     I ,'>>',word,'<<'
        k=k+1
      endif

      call getarg(4,word)    ! arg    kedgs
      read(word,*,iostat=ierr) inp
      if(inp.gt.0 .and. ierr.eq.0) then
        kedgs = inp
      else
       write(6,'(a,i3)')
     I 'Error: skipto special branch requires kedgs>0',kedgs
        k=k+1
      endif

      else ! std input, .not.skipto

      k=0
      call getarg(2,word)    ! arg   anum
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0.5) then
       write(6,'(a,f6.1,i3)')'Error anum<0.5, not possible for N_f!!',arg
     I ,ierr
        k=k+1
        anum=arg
      else
        anum=arg
      endif

      call getarg(3,word)    ! arg   cgain
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0.0) then
      write(6,'(a,f8.2)')'Error cgain<0, zero electrons!!',arg
        k=k+1
        cgain=arg
      elseif(arg.ne.0) then
        cgain=arg
      endif

      call getarg(4,word)    ! arg   f0
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0) then
       if(f0.eq.0) then
      write(6,'(a,f8.3)')'Error f0<0, sigma_0 not possible!!',arg
        k=k+1
        f0=arg
       else
        write(6,'(a)')'Warning: f0 reset to f0=0' 
        f0=0 
       endif
      elseif(arg.gt.0) then
        f0=arg
      endif

      endif ! skipto branch specialty
      gain = 1/cgain

      call getarg(5,word)    ! arg   fp
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0) then
       if(fp.eq.0) then
        write(6,'(a,f8.4)')'Error fp<0, c_2 not possible!!',arg
        k=k+1
        fp=arg
       else
        write(6,'(a)')'Warning: fp reset to fp=0' 
        fp=0 
       endif
      elseif(arg.gt.0) then
        fp=arg
      endif

      call getarg(6,word)    ! arg   fh
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0) then
       if(fh.eq.0) then
      write(6,'(a,f8.4)')'Error fh<0, sigma_e not possible!!',arg
        k=k+1
        fh=arg
       else
        write(6,'(a)')'Warning: fh reset to fh=0' 
        fh=0 
       endif
      elseif(arg.ne.0) then
        fh=arg
      endif

      call getarg(7,word)    ! arg   pix
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0.0) then
      write(6,'(a,f8.4)')'Error pix<=0, pixel pitch not possible!',arg
     I ,arg
        k=k+1
        pix=arg
      else
        pix=arg
      endif

      call getarg(8,word)    ! arg   alam
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0.0) then
      write(6,'(a,f8.4)')'Error alam<0, lambda not possible!!',arg
        k=k+1
        alam=arg
      else
        alam=arg
      endif

      call getarg(9,word)    ! inp   np
      read(word,*,iostat=ierr) inp 
      if(inp.eq.0 .or. ierr.ne.0) then ! OK keep as preset
      elseif(inp.lt.16 .or. np.lt.16) then
        write(6,'(a,i4)')'Error np is set <16, not allowed !',inp
        k=k+1
        np=inp
      else
        np=inp
      endif

      call getarg(10,word)   ! inp   naqs
      read(word,*,iostat=ierr) inp 
!     write(6,'(3a,3i5)')'ck_naqs >>',word,'<<',inp,naqs,ierr
      if(inp.eq.0 .or. ierr.ne.0) then ! OK keep as preset
      elseif(inp.lt.0) then
        write(6,'(a)')'Note: asking for full harmonic OTF parameters'
        uqspl = .false.
        naqs=inp       
      elseif(inp.lt.10) then
        write(6,'(a,i4)')'Error naqs<10, not allowed !',inp
        k=k+1
      elseif(inp .gt. min( np/2, 24) )then
       write(6,'(a,i4,a,i4)')'Error naqs=',inp,' is set too large max='
     I ,min(np/2,24)
        k=k+1
      else
        naqs=inp
      endif

      call getarg(11,word)   ! inp   nes
      read(word,*,iostat=ierr) inp 
!     write(6,'(3a,3i5)')'ck_nes  >>',word,'<<',inp,nes,ierr
      if(inp.eq.0 .or. ierr.ne.0) then ! OK keep as preset
      elseif(inp.lt.0) then
       if(nes.eq.0) then
        write(6,'(a,i4)')'Error nes<0, number of correction terms',inp
        k=k+1
       else ! force nes=0
        write(6,'(a)')'Warning: nes reset to nes=0' 
        nes=0
       endif
      else
        nes=inp
      endif

      call getarg(12,word)   ! inp   ity  ! model parameter, ity > 2 -> yi
      read(word,*,iostat=ierr) inp 
!     write(6,'(3a,3i5)')'ck_ity  >>',word,'<<',inp,ity,ierr
      if(inp.eq.0 .or. ierr.ne.0) then ! OK keep as preset
      elseif(inp.lt.-2) then
        write(6,'(a)')'Warning inp<-2, no stage 2,3 ',inp
        ity=inp
      elseif(inp.gt.5) then
        write(6,'(a)')'Warning inp>5, should not be needed,',
     I ' seek cause of convergence problem',inp
        ity=inp
      else
        ity=inp
      endif

      call getarg(13,word)   ! inp  mprbay   ! model parameter
      read(word,*,iostat=ierr) inp 
!     write(6,'(3a,3i5)')'ck_mpr  >>',word,'<<',inp,mpr,ierr
      if(inp.eq.0 .or. ierr.ne.0) then ! OK keep as preset
      else
        mprbay=inp
      endif

      call getarg(14,word)   ! inp   nquer
      read(word,*,iostat=ierr) inp 
!     write(6,'(3a,3i5)')'ck_nquer>>',word,'<<',inp,nquer,ierr
      if(inp.eq.0 .or. ierr.ne.0) then ! OK keep as preset
! usually reset to nquer = np -> iroip = iroin
      elseif(inp.lt.np/4) then
        write(6,'(a)')'Error nquer<np/4, not allowed, ignored', inp
      else
        nquer=inp
      endif

      call getarg(15,word)   ! arg   apixm    ! model parameter
      read(word,*,iostat=ierr) arg
!     write(6,'(3a,2f9.3,i5)')'ck_apixm >>',word,'<<',arg,apixm,ierr
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0.0) then
       if(apixm.eq.0) then
        write(6,'(a,f7.3)')'Error apixm<0, pixel aperture not possible'
     I ,arg
        k=k+1
        apixm=arg
       else
        write(6,'(a)')'Warning: no Reference Model due to apix<=0' 
        apixm=0 
       endif
      else
        apixm=arg
      endif

      call getarg(16,word)   ! arg   anumm   ! model parameter
      read(word,*,iostat=ierr) arg
      if(arg.eq.0.or.ierr.ne.0) then ! 
      elseif(arg.lt.0.5) then
        write(6,'(a,f6.1)')'Error anumm<0.5, set = anum',arg
        k=k+1
        anumm = arg
      else
        anumm=arg
      endif

      na = np*ns
      n1 = 1
      n2 = na
      nm = na/2
      iroin = np/2
      if(nquer.lt.16) nquer=np
      iroip = nquer/2
      nes1 = max(2,nes)
!     npolym = max( 2,nes)  ! internal spline points = nes-2
      npolym = nes1
      nparm1 = 5 + npolym
      if(anumm.lt.0.5) anumm=anum
      phim = phi

      if(skipto) then
      write(6,'(/2a)')'detailed call arguments in force:',
     I'        |> accessible from call with skipto option'
      elseif(domap) then
      write(6,'(/2a)')'detailed call arguments in force:',
     I' (anum input remains accessible through pgm image comment'
      else
      write(6,'(/a)')'detailed call arguments in force:'
      endif
      write(6,'(2a,f7.2,f8.3,f6.2,2f8.4,f7.3,f6.0,6i4,f6.2,f7.2)')
     I       'bayes_mtf ',filen(1:lenf),
     I        anum, cgain, f0, fp, fh, pix, alam,
     I        np, naqs, nes, ity, mprbay, nquer,
     I        apixm, anumm

      if(k.gt.0) then
        write(6,'(a,i3)')'Exit',k
        call exit(k)
      endif
      pix=pix*0.001  ! from now on in [mm]

      write(6,'(/2(a,f8.3),a,f8.5)')'Base Input: Lambda',alam,
     I '   N_f',anum,'   pixel_pitch [mum]',1000*pix

      write(6,'(2(a,f8.3),2(a,f8.5))')'Noise model: cgain',cgain,
     I '   f0',f0,'   fp',fp,'   fh',fh

      write(6,'(a,9i10)')'ROI: 2*roip,2*roin,np,ns',2*iroip,2*iroin
     I ,np,ns

      if(naqs.le.0) then
        write(6,'(a)')'Note: asking for full harmonic OTF parameters'
        uqspl = .false.
      endif

      if(apixm.le.0) then
       write(6,'(5(a,f8.3))')'No Reference model was defined (apix=0)'
      endif

        allocate(si(n1:n2),ci(n1:n2),circ(n1:n2),papt(n1:n2),stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc si,ci,circ,papt',n1,n2
        allocate( dcirc(n1:n2), yr(n1:n2), yi(n1:n2), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc dcirc,yr,yi',n1,n2
        allocate( esf(-nm:nm), csf(-nm:nm), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc esf,csf',n1,n2
        allocate( xk2(-nm:nm,5), vlp(npolym), pl(npolym), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc xk2,vlp,pl',n1,n2
        allocate( xes(nes1), yes(nes1), y2es(nes1,nes1-1), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc xes...',n1,n2
        allocate( cpmm(0:nm),cpix(0:nm), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc cpmm...',n1,n2
        allocate( afn(0:nm), bfn(0:nm), cfn(0:nm), dfn(0:nm),stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc afn...',n1,n2
        allocate( yrfm(0:nm), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc yrfm...',n1,n2

      return
      end

      subroutine read_image

! Author: Bernard Delley 2026
! read in image  *pgm    type P2
!  special branch 'crop iw1 iw2 ih1 ih2 iblack' for cropped image: pix_crop.pgm
!  special branch 'map' -> mprbay=0  anum=8 naqs=12 ity=10 (with conv threshold)
!  special branch 'skipto ktrps kedgs' to pick a specific edge from a chart for detailed analysis

      use image_data
      use esf_data
      use fit_data, only : umapper,domap
      use plot_data_compare

       dimension ihisto(4,nmm:mn),icum(4,nmm:mn)  ! automatic arrays
       real    ainp(7),binp(2)
       integer iinp(6)
       character*128 line,line1
       character*10 word
       logical crop,search,dotr,found

!     if(np.eq.64 .and. umapper) then
!     write(line,'(9a)')'mtf_mapper ',filen(1:lenf),
!    I '  . -q  --bayer green --single-roi --esf-model loess'
!     write(6,'(/9a)')'ck_line ',line
!     call system(line) 
!     write(line,'(9a)')'rm edge_mtf_values.txt edge_line_deviation.txt'
!     call system(line)
!     endif 

           if(mprbay.ge.5) then
      write(6,'(/a,f6.1,f9.3,f6.2,2f8.4,f7.3,f6.0,6i4,f6.2,f6.1)')
     I       'parameters default',
     I        anum, cgain, f0, fp, fh, pix, alam,
     I        np, naqs, nes, ity, mprbay, nquer, 
     I        apixm, anumm
            write(6,'(a)')
            endif

      crop=.false.
      domap=.false.
      skipto=.false.
      call getarg(2,word)   ! arg
      if(word.eq.'crop') crop=.true.
      if(word.eq.'map') domap=.true.
      if(word.eq.'skipto') skipto=.true.
      if(domap) then
        mprbay=0
        naqs = 12
        anum = 8.0
        ity = 10
      endif
      if(skipto) then
        naqs = 12
        anum = 8.0
        ity = 10
      endif

! insert to read input values associated with image file, first comment line

        read(15,*) ptype
        if(mprtra.gt.0) write(6,'(/2a)')   'image_ptype ',ptype
        read(15,'(a)') line
        line1 = '#'

        k=1
        do while(line(1:1).eq.'#')
         if(k.eq.1) write(6,'(2a)')'image file: ',filen(1:lenf)
     I  ,' comments:'
         write(6,'(a)') line
         if(k.eq.1) line1=line
         if(.not.crop) then
          if(k.eq.1) then
            l1=0
            do i=1,7
              read(line(2:),*,iostat=ierr) (ainp(j),j=1,i)
              if(ierr.eq.0) l1=l1+1
            enddo
            if(l1.ge.1 .and. ainp(1).gt.0.5) anum  = ainp(1)
            if(l1.ge.2 .and. ainp(2).gt.0.)  cgain = ainp(2)
            if(l1.ge.3 .and. ainp(3).gt.0.)  f0    = ainp(3)
            if(l1.ge.4 .and. ainp(4).gt.0.)  fp    = ainp(4)
            if(l1.ge.5 .and. ainp(5).gt.0.)  fh    = ainp(5)
            if(l1.ge.6 .and. ainp(6).gt.0.)  pix   = ainp(6)
            if(l1.ge.7 .and. ainp(7).gt.0.)  alam  = ainp(7)
            if(ierr.eq.0) then
             l1=0
             do i=1,6
              read(line(2:),*,iostat=ierr)ainp,(iinp(j),j=1,i)
              if(ierr.eq.0) l1=l1+1
             enddo
             if(l1.ge.1 .and. iinp(1).gt.0) np     = iinp(1)
             if(l1.ge.2 .and. iinp(2).ne.0) naqs   = iinp(2)
             if(l1.ge.3 .and. iinp(3).ge.0) nes    = iinp(3)
             if(l1.ge.4 .and. iinp(4).gt.0) ity    = iinp(4)
             if(l1.ge.5 .and. iinp(5).gt.0) mprbay = iinp(5)
             if(l1.ge.6 .and. iinp(6).gt.0) nquer  = iinp(6)
             if(ierr.eq.0) then
              read(line(2:),*,iostat=ierr)ainp,iinp,binp(1)
              if(ierr.eq.0 .and. binp(1).gt.0) apixm =binp(1)
              read(line(2:),*,iostat=ierr)ainp,iinp,binp(1),binp(2)
              if(ierr.eq.0 .and. binp(2).gt.0) anumm =binp(2)
             endif
            endif
           if(mprbay.ge.5) then
      write(6,'(/a,f6.1,f9.3,f6.2,2f8.4,f7.3,f6.0,6i4,f6.2,f6.1)')
     I       'parameters default and  from image file',
     I        anum, cgain, f0, fp, fh, pix, alam,
     I        np, naqs, nes, ity, mprbay, nquer,
     I        apixm, anumm
            write(6,'(a)')
           endif
          endif
         endif
          read(15,'(a)') line
          k=k+1
        enddo
        read(line,*) iw,ih
      
        read(15,*) imdn      ! maximum value  parameter for digital number
        monochrom=0
        if(imdn.gt.16383) monochrom=1
        iblack=0

        if(mprtra.gt.0) write(6,'(a,2i6,i8)')'ck_iw_ih_im',iw,ih,imdn
        if(mprtra.gt.0) write(6,'(a,2i5,i6,i15)')

! end insert for image file input parameters

!                 patch to crop out a custom pix file
      if(crop) then
        call getarg(3,word)  
        read(word,*,iostat=ierr) iw1
        iw1 = max(0,iw1/2)
        iw1 = 2*iw1 + 1
        if(ierr.ne.0) iw1=-1
        call getarg(4,word)  
        read(word,*,iostat=ierr) iw2
        iw2 = max(1,iw2/2)
        iw2 = min( 2*iw2, iw) 
        if(ierr.ne.0) iw1=-1

        call getarg(5,word)  
        read(word,*,iostat=ierr) ih1
        ih1 = max(0,ih1/2)
        ih1 = 2*ih1 + 1
        if(ierr.ne.0) iw1=-1
        call getarg(6,word)  
        read(word,*,iostat=ierr) ih2
        ih2 = max(1,ih2/2)
        ih2 = min( 2*ih2, ih) 
        if(ierr.ne.0) iw1=-1
        iw3=iw2-iw1+1
        ih3=ih2-ih1+1

        iblack=0
        call getarg(7,word)  
        read(word,*,iostat=ierr) iblack
        if(ierr.ne.0) iblack=0

        write(6,'(/3a,3(i7,i5))')
     I  'bayes_mtf ',filen(1:lenf),'  crop',iw1,iw2,ih1,ih2,iblack

        if( iw1.le.0 .or. iw3.lt.32 .or. ih3.lt.32 ) then
          write(6,'(/a,2(i8,i6))')'Error with crop definition'
     I   ,iw1,iw2,ih1,ih2
          call exit(2)
        endif

!       write(6,'(a,9i6)')'ck_iw,ih',iw1,iw2,iw3,ih1,ih2,ih3
        allocate( ipix(iw,ih) , jpix(iw3,ih3), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc ipix',istat,iw,ih
        read(15,*) ipix
        close(15)

        i1=0
        do i=ih1,ih2
          i1=i1+1
          j1=0
          do j=iw1,iw2
            j1=j1+1
            jpix(j1,i1) = ipix(j,i)-iblack
            jpix(j1,i1) = ipix(j,i)-iblack
          enddo
        enddo

        open(26,file='pix_crop.pgm')

        write(26,'(a)')ptype
        write(26,'(a)')line1
        write(26,'(2i5,i8,3i5)') iw3,ih3  !,iw1,iw2,ih1,ih2
        write(26,'(i5)') imdn
        write(26,'(15i6)') jpix
        if(iw3.gt.32 .or. ih3.gt.32) then
          write(6,'(/a,2i5)')
     I    'Cutout file pixc0.pgm prepared successfully: iw,ih:',iw3,ih3
        close(26)
!       else
!         write(6,'(/a,2i5,a)')'Cutout failed, iw,ih:'
!    I    ,iw3,ih3,' are way too small to use!'
        endif
        close(16)
        write(6,'(a)') 'cropped image pix_crop.pgm done'
        call exit(0)

      else   ! normal full read in
        allocate( ipix(iw,ih) , stat=istat)
        if(istat.ne.0) write(6,*)'Error allocating ipix',istat,iw,ih
        read(15,*) ipix
        close(15)
      endif
      return
      end

      subroutine atrpz

! find trapezoid edge parameters
! Author: Bernard Delley 2022

      use image_data
      use esf_data

       parameter( kptrmx=2500)
       dimension ihisto(4,nmm:mn),icum(4,nmm:mn)  ! automatic arrays
       integer iqub(4,4),iqt(4)
       real rms(2,4),dis(4),rnb(9)
       dimension ncnt(2,4)
       integer, allocatable :: itrpz(:,:,:),jtrpz(:,:,:)
       character*128 line
       logical belo,search,dotr,found

        do j=1,4
          do i=nmm,mn
            ihisto(j,i)=0
          enddo
        enddo

!       if(xfmt.eq.' ') then
         ib1=1
         ib2=2
         ib3=3
         ib4=4
!        xfmt='rggb'
!       endif

       ipixmin=99999
       ipixmax=0
       do i=1,ih
        do j=1,iw
          ipix(j,i) = max(ipix(j,i),0)
          iraw = ipix(j  ,i  )   ! - iblack  
          ipixmin = min( ipixmin, iraw )
          ipixmax = max( ipixmax, iraw )
          if(iraw.lt.nmm) then
            write(6,'(a,i6,a,2i6,a)')'Error: negative pixel value',iraw,
     I '   at',j,i
            stop 'error'
          endif
!      write(6,'(2i5,3i8,i24)')i,j,iraw,ipixmin,ipixmax
          ii=mod(i+1,2)
          jj=mod(j+1,2)
          kk=2*ii+jj+1
          ihisto(kk,iraw) = ihisto(kk,iraw)+1
        enddo
      enddo

      do i=1,4
        icum(i,nmm) = ihisto(i,nmm) 
        do j=nmm+1,mn
          icum(i,j) = icum(i,j-1) + ihisto(i,j)
        enddo
      enddo

!     do j=nmm,mn
!       write(16,'(9i10)')(ihisto(i,j),i=1,4),j,(icum(i,j),i=1,4)
!     enddo

      if(mprtra.gt.0 .or. skipto)
     I write(6,'(/a,2i12)')
     I'lowest highest Octile + Quartile DN values RGGB',ipixmin,ipixmax
      iq0 = ih*iw/4
      iqt(1) = iq0/8
      iqt(2) = iq0/4
      iqt(3) = iq0*3/4
      iqt(4) = iq0*7/8
      do i=1,4
       iq1=iqt(i)
      do j=nmm,mn-1
       do ii=1,4
        belo = icum(ii,j).le.iq1
         if(belo) iqub(ii,i)=j+1  ! mod 240110bd
       enddo
      enddo
      if(mprtra.gt.0 .or. skipto) 
     Iwrite(6,'(a,5i10)')'iqub',(iqub(ii,i),ii=1,4),iqt(i)
      enddo

!     istr =   0.25*(  iqub(2,1)  +iqub(3,1) +iqub(2,3)+iqub(3,3))
!     istr =  0.125*(3*iqub(2,1)+3*iqub(3,1) +iqub(2,3)+iqub(3,3))
!     istr =    0.1*(4*iqub(2,1)+4*iqub(3,1) +iqub(2,3)+iqub(3,3))
      istr = 0.0625*(7*iqub(2,1)+7*iqub(3,1) +iqub(2,3)+iqub(3,3))
!     istr =   0.05*(9*iqub(2,1)+9*iqub(3,1) +iqub(2,3)+iqub(3,3))
      
      if(mprtra.gt.0 .or. skipto) 
     Iwrite(6,'(a,9i10)')'global strike for trapezoids',istr
      ntrp = 0

      allocate( itrpz(2,10,kptrmx), stat=istat)

      do ktrp=1,kptrmx 
       do k1=1,2
        do k2=1,8
         itrpz(k1,k2,ktrp)=0
        enddo
       enddo 
         itrpz(2,1,ktrp) = ih+1
         itrpz(1,2,ktrp) = iw+1
         itrpz(2,5,ktrp) = ih+1
         itrpz(1,6,ktrp) = iw+1
      enddo

      jth=8
      ith=4
      ktrpm=0
      ktrp=1
      do i=1,ih
        jl=0
        jr=0
        do j=1,iw
          iraw = ipix(j  ,i  ) ! - iblack
          ii=mod(i+1,2)
          jj=mod(j+1,2)
          kk=2*ii+jj+1
          dotr=.false.
      if(kk.eq.2 .or. kk.eq.3 .and. iraw.gt.0) then
        if( iraw.lt.istr) then                        ! black spot found
          if(jl.eq.0) jl = j
          jr = j
!         if(j.eq.iw) dotr=.true.
          if(j.ge.iw-1) dotr=.true.
!     if(i.gt.1 .and. i.le.153 .and. mprtra.gt.7) 
!    I  write(6,'(a,5i5,i6,2x,2L2,2(4(i8,i5),2x))')'ck_B',ktrp,jl,jr, ! Black
!    I    j,i,iraw,found,dotr, ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,8)
        else
          if(jr.gt.0) dotr=.true.  ! white after trapez
!     if(i.gt.1 .and. i.le.153 .and. mprtra.gt.7) 
!    I  write(6,'(a,5i5,i6,2x,2L2,2(4(i8,i5),2x))')'ck_W',ktrp,jl,jr, ! White
!    I    j,i,iraw,found,dotr, ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,8)
        endif
        if(dotr) then
          found=.false.
          search = .true.
          ktrp=0
          do while(search)
            ktrp=ktrp+1
            if(ktrp.gt.ktrpm) then
              search=.false.
              ktrpm = ktrp
              if(ktrp.gt.kptrmx) stop 'ktrp exceeded'
!             if(ktrp.gt.73) goto 73
            endif
!           if(jr.ge.itrpz(1,2,ktrp)-1 .and. jl.le.itrpz(1,4,ktrp)) then  
            if(jr.ge.itrpz(1,2,ktrp)-1 .and.jl.le.itrpz(1,4,ktrp)+1)then ! 250728bd 
              if(i.le.itrpz(2,3,ktrp)+2) then
                found=.true.
                search=.false.
              endif
            endif
          enddo
         if(found) then
            if(jl.lt.itrpz(1,2,ktrp) ) then
               itrpz(1,2,ktrp) = jl
               itrpz(2,2,ktrp) = i
            endif 
            if(jl.le.itrpz(1,6,ktrp) ) then
               itrpz(1,6,ktrp) = jl
               itrpz(2,6,ktrp) = i
            endif 
!           if(jr.gt.itrpz(1,4,ktrp) ) then
!              itrpz(1,5,ktrp) = jr         ! backup for when 1 is right, but not as high as 2
!              itrpz(2,5,ktrp) = i
!           endif
            if(jr.ge.itrpz(1,4,ktrp) ) then
               itrpz(1,4,ktrp) = jr         
               itrpz(2,4,ktrp) = i
            endif
            if(jr.gt.itrpz(1,8,ktrp) ) then
               itrpz(1,8,ktrp) = jr         
               itrpz(2,8,ktrp) = i
            endif
!           if(jr.gt.jl ) then
!              itrpz(1,8,ktrp) = jr         ! backup for 4, when not most right
!              itrpz(2,8,ktrp) = i
!           endif
!           if(jl.le.itrpz(1,3,ktrp) ) then
!             itrpz(1,7,ktrp) = jl          ! backup for 3 when left, but not lowest
!             itrpz(2,7,ktrp) = i
!           endif
            itrpz(1,3,ktrp) = jl
            itrpz(2,3,ktrp) = i
            itrpz(1,7,ktrp) = jr
            itrpz(2,7,ktrp) = i
            itrpz(1,9,ktrp) = itrpz(1,9,ktrp) + (jl+jr)*(1+jr-jl) ! cgx
            itrpz(2,9,ktrp) = itrpz(2,9,ktrp) + i*(1+jr-jl)       ! cgy
            itrpz(1,10,ktrp) = itrpz(1,10,ktrp) + (1+jr-jl)       ! w,n
         else  ! .not.found  -> a new trapez
            itrpz(1,1,ktrp) = jr
            itrpz(2,1,ktrp) = i
            itrpz(1,5,ktrp) = jl
            itrpz(2,5,ktrp) = i
               itrpz(1,2,ktrp) = jl
               itrpz(2,2,ktrp) = i
               itrpz(1,6,ktrp) = jl         ! backup for alt 2, when  not most left
               itrpz(2,6,ktrp) = i
            itrpz(2,3,ktrp) = i
            itrpz(1,3,ktrp) = jl
            itrpz(2,7,ktrp) = i
            itrpz(1,7,ktrp) = jr
               itrpz(1,4,ktrp) = jr
               itrpz(2,4,ktrp) = i
               itrpz(1,8,ktrp) = jr
               itrpz(2,8,ktrp) = i
            itrpz(1,9,ktrp) = itrpz(1,9,ktrp) + (jl+jr)*(1+jr-jl) ! cgx
            itrpz(2,9,ktrp) = itrpz(2,9,ktrp) + i*(1+jr-jl)       ! cgy
            itrpz(1,10,ktrp) = itrpz(1,10,ktrp) + (1+jr-jl)       ! w,n
         endif
      if(mprtra.gt.5) 
     I  write(6,'(a,5i5,8x,2L2,2(4(i8,i5),2x))')'ck_F',ktrp,jl,jr,j,i,
     I    found,dotr, ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,8)
            jl = 0
            jr = 0
        endif
      endif
        enddo
      enddo
  73  continue
!     itrpe
      if(mprtra.gt.0) write(6,'(a,i5)')'ck_number_of_trapezoids',ktrpm
      call flush(6)

      do ktrp=1,ktrpm  ! new checkup/fix method 250727bd, check for bad cases with backup option
        found=.false.
      
        if(itrpz(1,1,ktrp) .le. itrpz(1,2,ktrp)+itrpe 
     I .and. itrpz(2,1,ktrp) .eq. itrpz(2,2,ktrp) ) then
          found=.true.
          icase=1
        endif
        if(itrpz(1,2,ktrp) .eq. itrpz(1,3,ktrp) 
     I .and. itrpz(2,2,ktrp) .eq. itrpz(2,3,ktrp) ) then
          found=.true.
          icase=2
        endif
        if(itrpz(1,4,ktrp) .le. itrpz(1,3,ktrp)+itrpe 
     I .and. itrpz(2,4,ktrp) .eq. itrpz(2,3,ktrp) ) then
          found=.true.
          icase=3
        endif
        if(itrpz(1,1,ktrp) .eq. itrpz(1,4,ktrp) 
     I .and. itrpz(2,1,ktrp) .eq. itrpz(2,4,ktrp) ) then
          found=.true.
          icase=4
        endif
        i11 = itrpz(1,2,ktrp)-itrpz(1,1,ktrp)
        i12 = itrpz(1,3,ktrp)-itrpz(1,2,ktrp)
        i13 = itrpz(1,4,ktrp)-itrpz(1,3,ktrp)
        i14 = itrpz(1,1,ktrp)-itrpz(1,4,ktrp)
        i21 = itrpz(2,2,ktrp)-itrpz(2,1,ktrp)
        i22 = itrpz(2,3,ktrp)-itrpz(2,2,ktrp)
        i23 = itrpz(2,4,ktrp)-itrpz(2,3,ktrp)
        i24 = itrpz(2,1,ktrp)-itrpz(2,4,ktrp)
        i123 = i11*i22 -i21*i12
        i341 = i13*i24 -i23*i14
        i15 = itrpz(1,6,ktrp)-itrpz(1,5,ktrp)
        i16 = itrpz(1,7,ktrp)-itrpz(1,6,ktrp)
        i17 = itrpz(1,8,ktrp)-itrpz(1,7,ktrp)
        i18 = itrpz(1,5,ktrp)-itrpz(1,8,ktrp)
        i25 = itrpz(2,6,ktrp)-itrpz(2,5,ktrp)
        i26 = itrpz(2,7,ktrp)-itrpz(2,6,ktrp)
        i27 = itrpz(2,8,ktrp)-itrpz(2,7,ktrp)
        i28 = itrpz(2,5,ktrp)-itrpz(2,8,ktrp)
        i567 = i15*i26 -i25*i16
        i785 = i17*i28 -i27*i18
!     if(skipto .and. ktrp.eq.ktrps) 
!    I  write(6,'(a,i5,4i10,i2,5x,2(4(i8,i5),4x))')'ck_rc',ktrp,
!    I  i123,i567,i341,i785
        icross = i123-i567 +i341-i785
      itrpz(1,9,ktrp) = itrpz(1,9,ktrp)/(itrpz(1,10,ktrp)*2)  ! 2 is a Bayer factor
      itrpz(2,9,ktrp) = itrpz(2,9,ktrp)/itrpz(1,10,ktrp)
      itrpz(2,10,ktrp) = (i123+i341)*1000./itrpz(1,10,ktrp)   ! expected quotient -2*1000
!     if(mprtra.gt.3 .or. (skipto .and. ktrp.eq.ktrps))
      if(mprtra.gt.3)
     I  write(6,'(a,i5,L2,i2,5x,3(4(i8,i5),4x))')'ck_r0',ktrp
     I ,found,icase,((itrpz(k1,k2,ktrp),k1=1,2),k2=1,10)

        if(found) then !  mprtra.gt.4) 
!       if(mprtra.gt.0)
!    I  write(6,'(a,i5,9x,2(4(i8,i5),4x))')'ck_r1',ktrp,
!    I    ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,8)
          itrpz(1,icase,ktrp) = itrpz(1,icase+4,ktrp)
          itrpz(2,icase,ktrp) = itrpz(2,icase+4,ktrp)
        if(mprtra.gt.2 .or. (skipto .and. ktrp.eq.ktrps))
     !  write(6,'(a,2i5,4x,4(i8,i5))')'ck_R1',ktrp,icase,
     I    ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,4)
        endif

        if(icross.gt.0) then
          itrpz(2,10,ktrp) = (i567+i785)*1000./itrpz(1,10,ktrp)    
          if(mprtra.gt.2 .or. (skipto .and. ktrp.eq.ktrps))
     I    write(6,'(a,i5,i9,3(4(i8,i5),4x))')'ck_RC',ktrp,icross,
     I    ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,10)
          do k2=1,4
            do k1=1,2
              itrpz(k1,k2,ktrp) = itrpz(k1,k2+4,ktrp)
            enddo
          enddo
        endif

        j = 3*np/2 
        if(abs(itrpz(2,10,ktrp)) .lt. 1850    ! relative area spanned by corners is too small
     I   .and. itrpz(1,10,ktrp) .gt. 2500     ! and area is greater than min allowed size
     I   .and. itrpz(1,9,ktrp) .gt.j .and. itrpz(1,9,ktrp) .lt. iw-j
     I   .and. itrpz(2,9,ktrp) .gt.j .and. itrpz(2,9,ktrp) .lt. ih-j
     I  ) then
          aside = sqrt( float(itrpz(1,10,ktrp)) )
          jv = itrpz(1,9,ktrp) + 0.5*aside
          iv = itrpz(2,9,ktrp) - 0.5*aside
          if( mod( jv+iv, 2) .eq. 0 ) iv=iv-1
          itrpz(1,1,ktrp) = jv
          itrpz(2,1,ktrp) = iv
          jv = itrpz(1,9,ktrp) - 0.5*aside
          if( mod( jv+iv, 2) .eq. 0 ) iv=iv+1
          itrpz(1,2,ktrp) = jv
          itrpz(2,2,ktrp) = iv
          iv = itrpz(2,9,ktrp) + 0.5*aside
          if( mod( jv+iv, 2) .eq. 0 ) iv=iv+1
          itrpz(1,3,ktrp) = jv
          itrpz(2,3,ktrp) = iv
          jv = itrpz(1,9,ktrp) + 0.5*aside
          if( mod( jv+iv, 2) .eq. 0 ) iv=iv-1
          itrpz(1,4,ktrp) = jv
          itrpz(2,4,ktrp) = iv
        i11 = itrpz(1,2,ktrp)-itrpz(1,1,ktrp)
        i12 = itrpz(1,3,ktrp)-itrpz(1,2,ktrp)
        i13 = itrpz(1,4,ktrp)-itrpz(1,3,ktrp)
        i14 = itrpz(1,1,ktrp)-itrpz(1,4,ktrp)
        i21 = itrpz(2,2,ktrp)-itrpz(2,1,ktrp)
        i22 = itrpz(2,3,ktrp)-itrpz(2,2,ktrp)
        i23 = itrpz(2,4,ktrp)-itrpz(2,3,ktrp)
        i24 = itrpz(2,1,ktrp)-itrpz(2,4,ktrp)
        i123 = i11*i22 -i21*i12
        i341 = i13*i24 -i23*i14
        itrpz(2,10,ktrp) = (i123+i341)*1000./itrpz(1,10,ktrp)   ! expected quotient -2*1000
      if(mprtra.gt.2 .or. (skipto .and. ktrp.eq.ktrps))
     I  write(6,'(a,i5,i4,5x,3(4(i8,i5),4x))')'ck_RR',ktrp
     I ,nint(aside),((itrpz(k1,k2,ktrp),k1=1,2),k2=1,10)
        endif

      enddo

! check now edge length, suppress misformed trapezoids with small edge length
      L0 = +1
      L1=iw-L0
      L2=ih-L0
      L0 = L0+1
      min1=100
      min2=100
      dmin=100
!     if(ktrpm.gt.1) then
!      write(6,'(a,i5)')'ck_k>?',ktrpm
       kk=0
       do ktrp=1,ktrpm
          do i=1,4
            j=mod(i,4)+1
            dis(i) = sqrt( float(
     I        (itrpz(1,j,ktrp) - itrpz(1,i,ktrp))**2 +
     I        (itrpz(2,j,ktrp) - itrpz(2,i,ktrp))**2 ))
          enddo
            dis4 = min(dis(1),dis(2),dis(3),dis(4))
        if(    itrpz(2,1,ktrp).gt.L0 .and. itrpz(1,2,ktrp).gt.L0
     I   .and. itrpz(2,3,ktrp).lt.L2 .and. itrpz(1,4,ktrp).lt.L1
     I   .and.  dis4.gt.ithsq .or. ktrpm.eq.1) then
!    I   .and.  dis4.gt.1 .or. ktrpm.eq.1) then
          kk=kk+1
          do k2=1,10 ! 4   ! 260705 bd all info is moved to kk
           do k1=1,2
            itrpz(k1,k2,kk) = itrpz(k1,k2,ktrp)         ! kept and moved to kk 
           enddo
          enddo
          minv = min( itrpz(2,2,kk)-itrpz(2,1,kk), 
     I                itrpz(2,4,kk)-itrpz(2,1,kk), 
     I                itrpz(2,3,kk)-itrpz(2,2,kk), 
     I                itrpz(2,3,kk)-itrpz(2,4,kk)) 
          minh = min( itrpz(1,1,kk)-itrpz(1,2,kk), 
     I                itrpz(1,3,kk)-itrpz(1,2,kk), 
     I                itrpz(1,4,kk)-itrpz(1,1,kk), 
     I                itrpz(1,4,kk)-itrpz(1,3,kk)) 
      if(mprtra.gt.1) 
     I    write(6,'(a,2i5,6(i8,i5),4x,5f8.1,3i5)')'ck_T',ktrp,kk,
     I    ((itrpz(k1,k2,kk),k1=1,2),k2=1,4),
     I    ((itrpz(k1,k2,kk),k1=1,2),k2=9,10),dis,dis4,minh,minv,ithsq

        else
         if(mprtra.gt.3) 
     I    write(6,'(a,i5,5x,6(i8,i5),4x,5f8.1)')'ck_R',ktrp,
     I    ((itrpz(k1,k2,ktrp),k1=1,2),k2=1,4),
     I    ((itrpz(k1,k2,ktrp),k1=1,2),k2=9,10),dis,dis4
        endif
       enddo
       ktrpm=kk

!     allocate( edge(5,4,ktrpm), stat=istat)
      allocate( edge(8,4,ktrpm), stat=istat)
      if(istat.ne.0) then
         write(6,*)'Error alloc edge',istat
         stop 'alloc'
      endif
      xh = iw/2 +0.5
      yh = ih/2 +0.5
      eox = 1
      eoy = 0
      do kk=1,ktrpm
        cgx=0
        cgy=0 
        do k2=1,4
          k = mod(k2,4) +1
          eoz = eox
          eox = eoy
          eoy = -eoz
          enx = +( itrpz(2,k2,kk) -itrpz(2,k,kk))    ! vn in bright direction
          eny = -( itrpz(1,k2,kk) -itrpz(1,k,kk))
          anor = sqrt( enx**2 + eny**2)
          emx = 0.5*(itrpz(1,k2,kk) +itrpz(1,k,kk))
          emy = 0.5*(itrpz(2,k2,kk) +itrpz(2,k,kk))
          ecx = emx - xh
          ecy = emy - yh
          bnor = sqrt( ecx**2 + ecy**2)
          edge(1,k2,kk) = emx
          edge(2,k2,kk) = emy
          edge(3,k2,kk) = anor
          sphi = max( -1.0, min( 1.0, (enx*eoy - eny*eox)/anor ))
          phi0 = asin( sphi ) + (k2-1)*pi2
          spsi = max( -1.0, min( 1.0, (enx*ecy - eny*ecx)/(anor*bnor)))
          psi0 = asin( abs(spsi) )
          edge(4,k2,kk) = phi0*pi2deg 
          edge(5,k2,kk) = psi0*pi2deg
          cgx=cgx+emx
          cgy=cgy+emy
          edge(6,k2,kk) = itrpz(1,k2,kk)
          edge(7,k2,kk) = itrpz(2,k2,kk)
        enddo
        igx=(cgx+2)/4
        igy=(cgy+2)/4
      if(mprtra.gt.2) 
     I  write(6,'(a,i4,2i5,2x,4(2f8.1,f7.1,2f6.1,2x),4f8.1)')'ck_t',kk,
     I igx,igy,((edge(k,k2,kk),k=1,5),k2=1,4)
      enddo

      if(mprtra.gt.1) 
     Iwrite(6,'(a,3i5,f7.1)')'ck_ktrp,min1,2,d',ktrpm,min1,min2,dmin

      deallocate( itrpz )
!     if(i.gt.0) stop 'stop-here'
      if(ktrpm.le.0) stop 'ktrpm=0'

      return
      end

      subroutine ck_edge( ktrp, kedg, edgeok )

! Author: Bernard Delley 2022

      use image_data
      use esf_data
      use fit_data
      logical edgeok

      edgeok=.false.
      if( edge(1,kedg,ktrpm) .lt. thed .or. 
     I    edge(1,kedg,ktrpm) .gt. iw-thed+1 .or.
     I    edge(2,kedg,ktrpm) .lt. thed .or.
     I    edge(2,kedg,ktrpm) .gt. ih-thed+1 ) then
       if(mprtra.gt.1)
     I write(6,'(a)')'ck_edge_too_close_to frame'
         do i=1,n2
           yr(i) = 0
           yi(i) = 0
         enddo
         do i=1,5
           sumry(i) = 0
         enddo
         do i=1,2
           isumry(i) = 0
         enddo
      else              ! edge is OK
         edgeok=.true.
      endif
      return
      end

      subroutine setcns

! Author: Bernard Delley 1990 +revs

      use esf_data
      one = 1
      pi4 = atan(one)
      pi2 = 2*pi4
      pi  = 2*pi2
      tpi = 2*pi
      pii = one/pi
      pre = 2/pi
      pi2deg = 45/pi4
      return
      end

      subroutine preptable
! Author: Bernard Delley 
      use esf_data
      use fit_data
      use image_data

!     character*20 word
! prepare si, ci, circ and initial papt tables
! initialize OTF   yr = circ*papt

! ncc index of cutoff frequency for circular aperture, highest non zero circ(ncc)

      fac = 1./ns

      bpix = one/(np*pix)
      bp = one/np
      d0 = one/anum
      dcut = anum*alam*1.e-6
      fcut = 1/dcut    ! c/mm
!     fcut = pix*fcut  ! c/p

!     write(6,'(14x,9a12)')'tpi*i','ci','si','fr','ff','circ','fp'
!     do i=0,np
      do i=0,nm
        cpmm(i)=i*bpix
        cpix(i)=i*bp
      enddo
      tpin =   2*pi/na
      do i=n1,n2
        ip=i
        if(i.gt.nm) ip=i-na
        si(i)=sin(tpin*i)
        ci(i)=cos(tpin*i)
        circ(i) = 0
        freq = ip*bpix
        ff = min( one, abs( freq/fcut) )
        sq = sqrt(one -ff*ff)
        circ(i) = pre*( acos(ff) - ff*sq )   ! circular aperture MTF
        dcirc(i) = -2*pre*sq
      enddo

        fmax = nm*bpix
        ncc = fcut/bpix
        ncd = min(ncc,nm)
!       write(6,'(a,i6,2f9.1,a,f10.6,f12.6)')
!    I 'aperture cutoff ncc fcut,fmax[c/mm]' ,ncc, fcut,fmax,
!    I '   fcut_pix [c/p]',pix*fcut,fcut/bpix
        write(6,'(a,i6,f9.1,f11.5,f10.1)')
     I 'aperture cutoff ncc fcut[c/mm] fcut[c/p] fmax[c/mm]',ncc,fcut,
     I  pix*fcut,fmax
        if(fcut.gt.fmax) then
         write(6,'(2a,2i5,3x,a,f10.3)')'Warning, oversampling ',
     I  'insufficient  for this fast aperture',ncc,ncd,
     I  'proposed oversampling factor:',ns*fcut/fmax  !(fcut/fmax + 1.)
          if(ncd.le.1) stop 'ncd error'
        endif

      allocate( dsf(ncd,-nm:nm,2), fsf(ncd,-nm:nm,2), stat=istat)
        if(istat.ne.0) then
          write(6,*)'Error alloc dsf,fsf',istat,ncd,nm
          stop 'alloc'
        endif


        nwh = (1+monochrom)*2.2*iroip*iroin +nm  !iw*ih
        npoly = 0
!       nparm = 8+2*ncd
        nparm = 4 + npolym + 2*ncd

        istata=0
      if(mprbay.gt.4) then
        write(6,'(a,9i9)')'ck_A_nwh nparm npolym ncd',
     I  nwh,nparm,npolym,ncd,nes
      endif
        allocate( udes(nwh,nparm), wdes(nwh,nparm), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc udes wdes  ',istat
     I  ,nwh,nparm
        istata=max(istat,istata)
        allocate( rmsi(nwh), brh(nwh), ipxy(2,nwh), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc rmsi brh   ',istat,nwh
        istata=max(istat,istata)
        allocate( crh(nwh), drh(nwh), erh(nwh), ke1(nwh), stat=istat)
        allocate( rk2(2,nwh), xk1(nwh,5), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc rk2,xk1,ke1',istat,nwh
        istata=max(istat,istata)
        allocate( vrk(nparm,nparm), sngv(nparm), awr(nparm), stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc vrk sngv   ',istat
     I  ,nparm
        istata=max(istat,istata)
        allocate( cvm(nparm,nparm), rv1(nparm),  stat=istat)
        if(istat.ne.0) write(6,*)'Error alloc cvm rv1    ',istat
     I  ,nparm
        istata=max(istat,istata)
        if(istata.ne.0) stop 'Error alloc'

      return
      end

      subroutine init_y00
! Author: Bernard Delley 2026
      use esf_data

!     esfh = 0
      esfh = 4*dcut/(pi*pi*pi*np*pix) ! 0.129006

      y00 = 1 -2*esfh
!     write(6,'(/a,2i6,f8.2,9f10.6)')'ck_y00',-nm-i,i,anum,y00,esfh
!    I ,dcut/pix
      return
      end

      subroutine getesfd
! Author: Bernard Delley 2023
      use esf_data
      use fit_data
!  calculate  line spread fun csf , 
!  and edge spread function esf 
!  and first derivatives to yr yi

      if(mprbay.ge.3) then
      if(ihys.ge.2) then
        write(6,'(a,i4,a,i4)')'recording hypothesis for',npar,' params,
     I   ihys=',ihys
      endif
      endif
      csfmin= 00.
      bp = 1.0/np
      bp2 = bp*2
      bs = 1.0/ns
!     y00 =1
      do i=-nm,nm               ! goes up monotonically
        ip = i
        if(ip.le.0 ) ip=i+na      ! 1 <= ip <= na
        csf(i) = y00*bp
        esf(i) = y00*(i+nm)*bp*bs
        xkyrfm = esf(i)
        do j=1,ncd    ! only nc1->ncd  and nc2-ncd+1->nc2 are needed
          tpij = pii/j
          ii = mod( ip*j-1 , na) +1
          csf(i) = csf(i) + ci(ii)*yr(j)*bp2    ! yr must be even for csf to be real
!         csf(i) = csf(i) + si(ii)*yi(j)*bp2    ! yi must be odd
          csf(i) = csf(i) - si(ii)*yi(j)*bp2    ! yi must be odd  ! 260127bd
          dsf(j,i,1) = ci(ii)*bp2
!         dsf(j,i,2) = si(ii)*bp2
          dsf(j,i,2) =-si(ii)*bp2
!         esf(i) = esf(i) + tpij*(si(ii)*yr(j) -ci(ii)*yi(j))
          esf(i) = esf(i) + tpij*(si(ii)*yr(j) +ci(ii)*yi(j))
          xkyrfm = xkyrfm + tpij* si(ii)*yrfm(j)
          fsf(j,i,1) = tpij*si(ii)
!         fsf(j,i,2) = -tpij*ci(ii)
          fsf(j,i,2) = +tpij*ci(ii)
        enddo
        if(csf(i).lt.csfmin .and. i.gt.0) then
          csfmin = csf(i)  !min( csf(i), csfmin ) 
          isumry(1) = i
        endif
        if(mprplo.gt.1) then
         if(ihys.eq.1) then
          write(6,'(a,i5,9f12.5,9f14.2)')'ck_xk',i,
     I    esf(i),csf(i)
         endif
        if(ihys.eq.2) then     ! perepare plot data
          xk2(i,2) = fil(1)+esf(i)*fil(2)
          xk2(i,1) = i*bs
!         write(6,'(a,i5,9f12.5,9f14.2)')'ck_xk',i,xk2(i,1),xk2(i,2),
!    I    esf(i),csf(i),yr(ip),ci(ip),si(ip) !,fil(1),fil(2)
        endif
        if(ihys.eq.3 .or. ihys.eq.4) then     ! prepare plot data
          xk2(i,ihys) = fil(1)+esf(i)*fil(2)
!         write(6,'(a,i5,9f14.2)')'ck_xk+5',i,xk2(i,1),xk2(i,ihys),
!    I    xk2(i,5)
!    I    esf(i),csf(i),yr(ip),yi(ip)
        endif
        if(ihys.eq.5) then     ! prepare plot data for initial yr=yrfm
          xk2(i,5) = fil(1)+xkyrfm*fil(2)
!         if(i.eq.0) write(6,'(a)')'ck_xkyrfm'
        endif
        endif
      enddo
      sumry(5) = csfmin/csf(0)

      return
      end

      subroutine analy2( ktrp, kedg )
! Author: Bernard Delley 2026

!  analyze image near edge, initialize edge fit data

      use image_data
      use esf_data
      dimension ihisto(4,nmm:mn),icum(4,nmm:mn)  ! automatic arrays
      integer iqub(4,4),iqt(4)
      real avg(2,4),rms(2,4),dis(4)
      integer ncnt(2,4)
      logical belo

      keep_roi=.false.

        do j=1,4
          do i=nmm,mn
            ihisto(j,i)=0
            icum(j,i)=0
          enddo
        enddo

        vm(1) = edge(1,kedg,ktrp)
        vm(2) = edge(2,kedg,ktrp)
!       anor =  edge(3,kedg,ktrp)  ! length of edge
!       roip = min( anor*0.4, float(iroip))
        jroip = roip

        phi   = edge(4,kedg,ktrp)/pi2deg
        ve(1) = cos(phi)
        ve(2) = -sin(phi)
!       ve(3) = 0
        vn(1) = -ve(2)
        vn(2) = ve(1)
        do j=1,npolym
          vlp(j)=0
        enddo

!       iw1 = max( 1, nint(vm(1)) -jroip -iroin )
!       iw2 = min(iw, nint(vm(1)) +jroip +iroin )
!       ih1 = max( 1, nint(vm(2)) -jroip -iroin )
!       ih2 = min(ih, nint(vm(2)) +jroip +iroin )
      jroip = sqrt( roip*roip + iroin*iroin ) +3
       iw1 = max( 1, nint(vm(1)) -jroip )
       iw2 = min(iw, nint(vm(1)) +jroip )
       ih1 = max( 1, nint(vm(2)) -jroip )
       ih2 = min(ih, nint(vm(2)) +jroip )

      if(mprtra.gt.4) then
      write(6,'(a,2i4,4i5,3f9.1,9f9.5)')'ck_a',ktrp,kedg,
     I iw1,iw2,ih1,ih2, vm, edge(4,kedg,ktrp), phi, ve, vn
      endif

        ij = 0
        do i=ih1,ih2
          do j=iw1,iw2
           iraw = ipix(j  ,i  ) ! - iblack
           ii=mod(i+1,2)
           jj=mod(j+1,2)
           kk=2*ii+jj+1
           dxp = j-vm(1)
           dyp = i-vm(2)
!          if(kk.eq.2 .or. kk.eq.3) then !  green
            xpa = dxp*ve(1) + dyp*ve(2)
            xno = dxp*vn(1) + dyp*vn(2)
!           xno = -xno   ! change sign back to original

            if( abs(xpa).lt.roip .and. abs(xno).lt.iroin ) then 
            if(kk.eq.2 .or. kk.eq.3) then !  green
              ij = ij + 1
!     write(6,'(a,2i5,i7,i3,2f7.1)')'ck_e',i,j,ij,kk,xpa,xno
            endif
              ihisto(kk,iraw) = ihisto(kk,iraw)+1
            endif
!          endif
          enddo
        enddo

      do i=1,4
        icum(i,nmm) = ihisto(i,nmm) 
        do j=nmm+1,mn
          icum(i,j) = icum(i,j-1) + ihisto(i,j)
        enddo
      enddo
      if(mprbay.gt.8) then
        write(6,'(/a)')'R-GG-B histograms and cumulative'
        do j=nmm+1,mn
           if( min(ihisto(2,j), ihisto(3,j)) .gt. 0) write(6,
     I '(2i10,i6,2i10)')ihisto(2,j), ihisto(3,j),j,icum(2,j),icum(3,j)
        enddo
      endif
      if(mprbay.gt.4) then
      write(6,'(/a)')'lowest highest quartile octile DN values RGGB'
      write(6,'(a,9i10)')'ncum',(icum(ii,mn),ii=1,4)
      endif

!     iq0 = (ncum(2)+ncum(3))/2
      iq0 = (ij+1)/2
      iqt(1) = iq0/8
      iqt(2) = iq0/4
      iqt(3) = iq0*3/4
      iqt(4) = iq0*7/8
      do i=1,4
       iq1=iqt(i)
      do j=nmm,mn-1
       do ii=1,4
        belo = icum(ii,j).le.iq1
         if(belo) iqub(ii,i)=j+1  ! mod 240109bd
       enddo
      enddo
      if(mprbay.gt.4) then
      write(6,'(a,9i10)')'iqub',(iqub(ii,i),ii=1,4),iqt(i)
      endif
      enddo
  
      do i=1,2
        do ii=1,4
          ncnt(i,ii) = 0
          avg(i,ii) = 0
        enddo
      enddo

!     do j=nmm,mn-1
!       do ii=1,4
!         if(j.le.iqub(ii,1)) then
!          i=1
!          belo=.true.
!         elseif(j.ge.iqub(ii,4)) then
!          i=2
!          belo=.true.
!         else
!          belo=.false.
!         endif
!         if(belo) then
!         avg(i,ii) = avg(i,ii) + ihisto(ii,j)*j
!         ncnt(i,ii) = ncnt(i,ii) + ihisto(ii,j)
!         endif
!       enddo
!     enddo
!     if(mprbay.gt.4) then
!     write(6,'(/a)')'dark and bright area analysis'
!     endif
!     do i=1,2
!       do ii=1,4
!         avg(i,ii) = avg(i,ii)/ncnt(i,ii)
!       enddo
!     if(mprbay.gt.4) then
!     write(6,'(a,4i10)')'nc  ',(ncnt(i,ii),ii=1,4)
!     write(6,'(a,4f10.1)')'avg ',(avg(i,ii),ii=1,4)
!     endif
!     enddo

!     fil(1) = ( avg(1,2) +avg(1,3) )*0.5   ! dark illum level
      fil(1) = ( iqub(2,1) +iqub(3,1) ) *0.5
!     darkthres = avg(1,2) +avg(1,3)        ! 2x starting fil(1) as threshold for fil(1) optimization
      darkthres = 2* fil(1)
      darkthres = 16383
!     fil(2) = ( avg(2,2) +avg(2,3) )*0.5   ! bright illum level
      fil(2) = ( iqub(2,4) +iqub(3,4) ) *0.5
      fil(2) = fil(2)-fil(1)
      fil(3) = 0                            ! x illum gradient
      fil(4) = 0                            ! y illum gradient

      if(mprbay.gt.4) then
      write(6,'(a,9f10.1)')'f illum    ',fil
      write(6,'(a,9f10.3)')'edge pos   ',vm
      write(6,'(a,9f10.5)')'edge vect  ',ve(1),ve(2)
      write(6,'(a,9f10.5)')'edge normal',vn
      endif

      return 
      end

      subroutine update1( damp, damp1, damp2 )
! update 3-8 parameter edge model
! Author: Bernard Delley 2025
      use image_data
      use esf_data
      use fit_data
      use plot_data_compare , only : apixm
      dimension alim(5)  
!     data alim / 0.05, 0.05, 1.0, 1.0, 1.0, 5.0, 0.1  , 5.0 /
!     data alim / 0.5 , 0.05, 0.01, 1.0, 5.0, 1.0, 1.0  , 5.0 / ! rev order 250708bd  250713bd-phi->3
!     data alim / 0.05, 0.01, 0.3 , 1.0, 5.0, 1.0, 1.0  , 5.0 / ! rev balance awr(1)
!     data alim / 0.05, 0.05, 0.001, 1.0, 5.0, 1.0, 1.0  , 5.0 / ! rev balance awr(3) 250804bd
      data alim / 0.40, 0.05, 0.03 , 1.0, 5.0 / ! rev balance awr(1) fil1 awr 1->3 251119bd

      if(mprbay.gt.6) then
      write(6,'(a,i2,2(3x,8f10.5))')'coefN',npar1,(awr(i),i=1,npar1)
      endif
      do i=1,5  ! npar1
        awr0 = awr(i)
        awr(i) = sign( min( alim(i), abs(awr(i)) ), awr(i))    ! limiting corrections in Newton step
       if(mprbay.gt.1) then
        if(awr0.ne.awr(i)) write(6,'(a,i4,2f10.5)')'ck_L',i,awr0,awr(i)
       endif
      enddo

      if(mprbay.gt.6) then
      write(6,'(a,2(5x,8f10.5))')'coefL',(awr(i),i=1,npar1)
      endif

      if(npar1.gt.0) then
       if(udark) then
       fil(1) =  fil(1)*(1 + scale1*damp*awr(1))             ! 1 illum  
!      fil(1) =  fil(1)    + scale1*awr(1)              ! 1 illum  
       endif
       fil(2) = fil(2)*(1 + damp*awr(2))              ! 2 illum  

       if(uphi) then
         phi = phi + awr(3)*damp2*scalan                    ! 7->3 
         ve(1) = cos(phi)
         ve(2) = -sin(phi)
         vn(1) = sin(phi)  !-ve(2)
         vn(2) = cos(phi)
       endif
       if(mprbay.gt.5) then
         write(6,'(a,f9.3,9f10.5)')'slanted edge angle  degrees,+rad'
     I , phi*pi2deg,phi, awr(4) !,awr(8)
         write(6,'(a,9f10.5)')'edge normal',vn
       endif
      endif

!     if(npar1.ge.4) then 
       if( uxno ) then
        pos0 =  pos0 +        awr(npx)*damp1                  ! 4 edge pos
        vm(1) = vm(1) + vn(1)*awr(npx)*damp1
        vm(2) = vm(2) + vn(2)*awr(npx)*damp1
!       pos0 =  pos0 + awr(4)*damp1                    ! 4 edge pos
!       vm(1) = vm(1) + vn(1)*awr(4)*damp1
!       vm(2) = vm(2) + vn(2)*awr(4)*damp1
       endif
     
       if( uanum ) then
         anufn = anuf + damp1*awr(npk)                ! 6->5 effective f-ratio as fit parameter "npk"
         if(anufn.lt.0.9) then
           anuf = 0.5*( anuf + 0.9 )
         else
           anuf = anufn
         endif
!        dcutf = anuf*alam*1.e-6
         if(mprbay.gt.5) then
          write(6,'(a,9f10.5)')'anuf      ',anuf,awr(npk)
         endif
       endif 

       if( ugrad ) then
         fil(3) = fil(3) + damp1*awr(npg)               ! 4->6 illum grad parallel
!      write(6,'(a,i3,9f10.5)')'ck_fil3',npg,fil(3),awr(npg),damp1
       endif
       if( ugnor ) then
         fil(4) = fil(4) + damp1*awr(nph)               ! 5->7 illum grad normal
       endif
!      write(6,'(a,l2,i3,9f10.5)')'ck_fil4',ugnor,nph,fil(4),awr(nph)
!    I ,damp1
       if(mprbay.gt.5) then
       write(6,'(a,9f10.2)')'illum param',fil(3),awr(npg)
       write(6,'(a,9f10.3)')'edge pos   ',vm,pos0,awr(4)
       endif

       if( ucurv ) then
       jpp = npc
       do jp=2,npoly
         vlp(jp) = vlp(jp) + awr(jpp)
         jpp = jpp +1
       enddo
      endif

      if(uanum) then    ! use the effective aperture model in  phase 1
        phis = min( abs(phi), abs(phi-pi2), abs(phi-pi), 
     I   abs(phi-3*pi2), abs(phi-2*pi))
        sth=sin(phis)
        cth=cos(phis)

       if(mprbay.gt.5) then
      write(6,'(a,19f12.5)')'ck_phi',phi*pi2deg,phis*pi2deg,phis,sth,cth 
       endif

         bp = 1.0/np
         dcutf = anuf*alam*1.e-6
!        yrfn(0) = 1
!        do i=1,nm
!          yrfn(i) = 0
!        enddo
         do i=1,ncd
          x = i*bp*apixm       ! phase 1 uses model input apix=apixm     
          papt(i) =     sinc(pi*x*sth)*sinc(pi*x*cth)
          ff = i*bpix*dcutf
          if(ff.ge.one) then
           yr(i) = 0
           yi(i) = 0
           dcirc(i)=0
           circi=0
          else
           sq = sqrt(one -ff*ff)
           circi = pre*( acos(ff) - ff*sqrt(one -ff*ff) )   ! circular aperture MTF
           dcirc(i) = -2*pre*sq

           yr(i) = circi  *papt(i)
!          yrfn(i)=yr(i)
           yi(i) = 0
          endif
!     write(6,'(a,i5,9f12.6)')'ck_d',i,yr(i),circi,papt(i),dcirc(i),ff
         enddo
      endif  ! uanum
!     endif  ! npar1.gt.4

         return
         end

      subroutine update2
! update yr edge model
! Author: Bernard Delley 2025
!     use image_data
      use esf_data
      use fit_data

!     if(onlyi) then
!       do j=1,ncd
!        yi(j) = yi(j) + awr(j+npar1)
!       enddo
!     else
        do j=1,ncd
!        y0 = yr(j)
!        yr(j) = yr(j) + awr(j+npar1)*damp1
         yr(j) = yr(j) + awr(j+npar1)
        enddo
        j1=0
        do j=ncd+1,npar
         j1=j1+1
         yi(j1) = yi(j1) + awr(j+npar1)
        enddo
!     endif

      return
      end

      subroutine init_yr
! init yr edge model, and awr ,ihys
! Author: Bernard Delley 2023
      use image_data
      use esf_data
      use fit_data
      use plot_data_compare , only : apixm

      phd = phi*pi2deg
        phis = min( abs(phd), abs(phd-90), abs(phd-180), abs(phd-270),
     I   abs(phd-360))
        sth=sin(phis/pi2deg)
        cth=cos(phis/pi2deg)
      if(mprbay.gt.5) then
       write(6,'(a,9f12.5)')'ck_phI',phi*pi2deg,phis,phis/pi2deg,sth,cth
      endif

      pos0 = 0
      anuf = anum
      dcut = anum*alam*1.e-6

         bp = 1.0/np
!        yrfn(0)= 1
         do i=1,ncd
          x = i*bp*apixm
          papt(i) =     sinc(pi*x*sth)*sinc(pi*x*cth)
          ff = i*bpix*dcut
          if(ff.ge.one) then
           yr(i) = 0
           yi(i) = 0
           dcirc(i)=0
          else
           sq = sqrt(one -ff*ff)
           circi = pre*( acos(ff) - ff*sqrt(one -ff*ff) )   ! circular aperture MTF
           dcirc(i) = -2*pre*sq

           yr(i) = circi    *papt(i)
!          yrfn(i)=yr(i)
           yi(i) = 0
!         write(6,'(a,i5,9f12.5)')'ck_d',i,yr(i),circ(i),papt(i) !,dcirc(i),dcirc(i)*ff
          endif
         enddo
!     write(6,'(a,f8.3)')'ck_yr Init',apixm

        do i=1,nparm
          awr(i) = 0
        enddo
        ihys = 0
      return
      end

      subroutine init_yrfm
! init yrfm edge model
! Author: Bernard Delley 2024
      use image_data
      use esf_data
      use fit_data
      use plot_data_compare
      character*20 word

        phim= phi*pi2deg
        phd = phi*pi2deg 
      if(apixm.gt.0) then
       write(6,'(/a,5(f8.3,a))')'Reference model: apixm',apixm,
     I '  anumm',anumm,'  phim',phim
      endif
!       phd = phi  ! phim
        phis = min( abs(phd), abs(phd-90), abs(phd-180), abs(phd-270),
     I   abs(phd-360))
        sth=sin(phis/pi2deg)
        cth=cos(phis/pi2deg)
!     write(6,'(5(a,f8.3))')'ck_yrfm anumm',anumm,'  phim',phim
!    I ,'  apixm',apixm
!     write(6,'(a,9f12.5)')'ck_phI',phim,phis,phis/pi2deg,sth,cth

         dct = anumm*alam*1.e-6
         bp = 1.0/np
         yrfm(0)= 1
         do i=1,ncd
          x = i*bp*apixm
          papti =     sinc(pi*x*sth)*sinc(pi*x*cth)
          ff = i*bpix*dct
          if(ff.ge.one) then
           yrfm(i) = 0
          else
           sq = sqrt(one -ff*ff)
           circi = pre*( acos(ff) - ff*sqrt(one -ff*ff) )   ! circular aperture MTF
           yrfm(i) = circi*papti
!         write(6,'(a,i5,9f12.5)')'ck_d',i,yrfm(i),circi,papti
          endif
         enddo

!     if(umapper) then
!        yrfmap(0)= 1
!        bpm = 1.0/64.0
!        bpmx=bpm/pix
!        do i=1,63
!         x = i*bpm*apixm
!         papti =     sinc(pi*x*sth)*sinc(pi*x*cth)
!         ff = i*bpmx*dct
!         if(ff.ge.one) then
!          yrfmap(i) = 0
!         else
!          sq = sqrt(one -ff*ff)
!          circi = pre*( acos(ff) - ff*sqrt(one -ff*ff) )   ! circular aperture MTF
!          yrfmap(i) = circi*papti
c         write(6,'(a,i5,9f12.5)')'ck_d',i,yrfm(i),circi,papti
!         endif
!        enddo
!     endif

      return
      end

      subroutine get_mtf ! (ktrp, kedg )
! perform the refinement iterations
! Author: Bernard Delley 2024
      use image_data
      use esf_data
      use fit_data
      use plot_data_compare , only : apixm
      character*51 fmt
        damp = 1.0  ! mixing parameters, 1: Newton no damping 
        damp1= 1.0  
        damp2= 1.0 
        scalan = 0.1        ! 0.159  rev 250804bd for rev balance phi awr(3)
        scale1 = 1.         ! option rev balance for fil(1)  awr(1)
        scalen = 0.003      ! option rev balance for fil(4) normal gradient
      write(fmt,'(a)') 
     I '(/a,i4,i5,i7,2f10.2, 5f10.5,a,f10.3,a,2f9.3,a,f6.3)'
      onlyi = .false.
      uimag = .false.

      npar =  5                   ! yet intentionally uphi=.false. thus not used for 0'th it
      npar1 = npar 
      npar2 = npar1
      npar3 = npar
      npoly = 0
      uphi=.false.    
      udark=.true.
! fil(1) fil(2) phi-3 xno-npx=4 anum-npk=5 (grad-npg=5 gnor-nph-6 curv-npc-7....)
      uxno=.true.
      npx=4
      uanum=.true.
      npk = 5
      ugnor=.false.
      nph = 5
      ugrad=.false.
      npg = 5       ! to be overwritten while npk=5
      ucurv=.false.
      npc = 6       ! to be overwritten while npk=5
      usigd=.false.
      jtn=0
!     keep start values in case stage 1 fails
      do i=1,2
        vm0(i) = vm(i)
        ve0(i) = ve(i)
        vn0(i) = vn(i)
        fil0(i) = fil(i)
      enddo

      call init_yr
      call init_y00
      call init_yrfm
     
      if(mprbay.gt.0) then
      write(6,'(/3a)') 'Stage 1    M     N   singular      z',
     I  '           y        conv       mx_coef_corr  phi_deg'
     I ,'    dark  j-bright'
      endif


      ihys=5
      call getesfd
      ihys=0
        keep_roi=.false.
        nnull=0
!     call walltime(wtime0,'other',2)

       if(mprbay.gt.2) then
        write(6,fmt)'Coef',jtn,npar
     I ,ndat,fil(1),fil(2),phi,pos0,anuf,fil(4),fil(3)    
     I ,'   phi_deg',phi*pi2deg
        write(6,'(20f10.5)') (yr(i),i=1,ncd)
        write(6,'(a,20x,2a)')'chisq','chisq    (chisq-n)/2 ',
     I '    zfun        yfun       yfun/n     Qgam        Conv'
       endif

      call postprob( yfun )


        call svdfitb(udes,brh,ndat,npar,nwh,nparm,sngv,awr,vrk,cvm,
     I  wdes,rv1,zfun,yfun,qgam,conv,jtn)

      if(skip1) then
        write(6,'(a)')'Skipping Stage1'
      else

      keep_roi=.true.
        ihys=0
      npar = 5
      npar1 = npar
      npar2 = npar1
      npar3 = npar
      uphi=.true.

      do itr=1,2
      conv = 1
      do while( conv.gt.0.001 )
        jtn=jtn+1
        call update1( damp, damp1, damp2 )
!     phi = phim/pi2deg
        call getesfd

        if(mprbay.gt.2 .and. .not.keep_roi) then
          write(6,'(a)')'Resetting ROI'
        endif
        if(.not.keep_roi) nnull=0

       if(mprbay.gt.2) then
        write(6,fmt)'Coef',jtn,npar
     I ,ndat,fil(1),fil(2),phi,pos0,anuf,fil(4),fil(3)    
     I ,'   phi_deg',phi*pi2deg
        write(6,'(20f10.5)') (yr(i),i=1,ncd)
       endif

        call postprob( yfun )

        call svdfitb(udes,brh,ndat,npar,nwh,nparm,sngv,awr,vrk,cvm,
     I  wdes,rv1,zfun,yfun,qgam,conv,jtn)

        convr=conv
        if(jtn.eq.20) conv=0
!       if(jtn.eq.200) conv=0
        keep_roi=.true.
      enddo  ! conv
        keep_roi=.false.
        jtn=0
        pos0=0
        ihys=2
      enddo
      endif ! skip1
        keep_roi=.true.
        sumry(4)=jtn+conv
        ihys=0

        call mkwipe
      
      if(mprbay.gt.0) then
      write(6,'(a)')'Stage 2'
      call flush(6)
      endif

! fil(1) fil(2) phi-3 xno-npx=4 gno-nph=5 grad-npg=6 curv-npc=7....
!     uxno=.true.    ! remains
!     npx=4
         uanum=.false.
!        npoly = nes-1 !  1 is minimum here for including npg: nes=2, for curv use npoly=2 and higher
         npoly = 0
         if(nes.gt.0) call mkesfn
!        if(npoly.ge.1) ugrad=.true.
         if(nes.ge.1) ugnor=.true.
         nph = 5
         if(nes.ge.2) ugrad=.true.
         npg = 6       
         if(nes.ge.3) ucurv=.true.
         npoly = nes-1 ! for curv use nes=3 and higher
         npc = 7
         npar1 = 4 + nes
         npar2 = npar1
       if(uqspl) then
         nqm0 = nqm 
         if(unats) nqm0 = nqm-2
         npar = npar1 + nqm0
         call mkqfn
       else
         npar = npar1 + ncd
       endif
         npar3 = npar
         write(fmt(21:22),'(i2)') 3 + max(2,nes)
         npk = npar1+1       ! value to be overwritten

      do itn=-2,ity ! 2 5 
        if(itn.gt.3 .and. conv.lt.0.001) goto 222
        call update1( damp, damp1, damp2 ) ! used for npar1  1...4
        if(uqspl) then
          call update3
        else
          call update2
        endif

!       ihys=0
!       if(itn.eq. 1 ) ihys=2
        if(itn.eq. 2 ) ihys=3
        if(itn.ge. 5 ) ihys=4

        call getesfd

        if(itn .eq. 3 ) then
      if(mprbay.gt.0) then
      write(6,'(a)')'Stage 3'
      endif
! fil(1) fil(2) phi-3           gno-nph=5 grad-npg=6 curv-npc=7....
!         udark=.false.
          uxno=.false.             
!     write(6,'(a)')'ck_Warning uxno kept true'
          ugnor=.false. ! should it remain, no
          nph = 4     
!         ugrad=      ! remains 
          npg = 4       
!         ucurv=      ! remains
          npc = 5
          npar1 = 2 + max(1,nes)
          npar2 = npar1
          uimag = .true.
         if(uqspl) then
          npar = npar1 + nqm0*2
         else
          npar = npar1 + ncd*2
         endif
        endif

       if(mprbay.gt.2) then
        write(6,fmt)'Coef',itn,npar
     I ,ndat,fil(1),fil(2),phi,pos0,anuf,fil(4),fil(3)
     I ,(vlp(jp),jp=2,npoly)
     I ,'   phi_deg',phi*pi2deg
       if(uqspl) then
        if(unats) then
         write(6,'(20f10.5)') (yqr(i),i=1,nqm)
         if(uimag)
     I   write(6,'(20f10.5)') (yqi(i),i=1,nqm)
        else 
         write(6,'(20f10.5)') (yqr(i),i=1,nqm0),dyr1,dyrn
        endif
         write(6,'(a)')'-------------------------------'
       endif
         write(6,'(20f10.5)') (yr(i),i=1,ncd)
         if(uimag)
     I   write(6,'(20f10.5)') (yi(i),i=1,ncd)
       endif

!       if(mprplo.gt.4. and. itn.eq.ity) call mkploto
        call postprob( yfun )
!       if(mprplo.gt.4. and. itn.eq.ity) call mkplote

        call svdfitb(udes,brh,ndat,npar,nwh,nparm,sngv,awr,vrk,cvm,
     I  wdes,rv1,zfun,yfun,qgam,conv,itn)
        if(itn.eq.1) sumry(1)=zfun
        if(itn.eq.2) sumry(2)=zfun
        if(itn.eq.4) sumry(3)=zfun
        if(itn.eq.4) sumry(4)=qgam

        ihys=0
       enddo ! itn

 222  continue 
        if(mprplo.gt.4) call mkplotou

        if(uqspl) then
          call update3
        else
          call update2
        endif
        if(itn.eq. 2 ) then
         ihys=3
         call getesfd
!       elseif(itn.eq. 3 ) then
        elseif(itn.gt. 3 ) then
         ihys=4
         call getesfd
        endif

      if(mprbay.gt.0) then
      if(uqspl) then
        write(6,'(/a,2i4)')'inferred envelope model Y parameters:'
!    I ,nqm0,nqm !<<>>
      else
        write(6,'(/a)')'inferred full harmonic model y parameters:'
      endif
        rmsa = fil(1)*sqrt(cvm(1,1))
        rmsb = fil(2)*sqrt(cvm(2,2))
        if(uxno.and.uanum) then
          write(6,'(a,2L3,9x,(20f10.5))')'RMS a',uanum,uxno,
     I  (sqrt(cvm(i,i)),i=1,npar1)
!    I  (rmsa,rmsb,sqrt(cvm(i,i)),i=3,npar1)
        elseif(uxno) then
          write(6,'(a,2L3,9x,(4f10.5,    20f10.5))')'RMS b',uanum,uxno,
     I  (sqrt(cvm(i,i)),i=1,npar1)
!    I  (rmsa,rmsb,sqrt(cvm(i,i)),i=3,npar1)
        else
          write(6,'(a,2L3,9x,(3f10.5,20x,20f10.5))')'RMS c',uanum,uxno,
     I  (sqrt(cvm(i,i)),i=1,npar1)
!    I  (rmsa,rmsb,sqrt(cvm(i,i)),i=3,npar1)
        endif
        if(npar.gt.npar1) then
         if(nqm.gt.1) then
           write(6,'(20f10.5)') (rv1(i),i=1,nqm0)
           if(uimag) write(6,'(20f10.5)') (rv1(i),i=nqm0+1,2*nqm0)
         else
           write(6,'(20f10.5)') (rv1(i),i=1,ncd)
          if(npar.gt.npar1+ncd)
     I     write(6,'(20f10.5)') (rv1(i),i=ncd+1,2*ncd)
         endif
        endif
      endif
        i=1
        do while( yr(i)**2+yi(i)**2 .gt.0.25 )
          i=i+1
        enddo
        y2 = sqrt(yr(i)**2+yi(i)**2)
        i=i-1
        y1 = sqrt(yr(i)**2+yi(i)**2)
        p = (0.5-y1)/(y2-y1)
        xm50 = cpix(i) +p/np
        i=np/4
        xhn = sqrt( yr(i)**2 + yi(i)**2 )

        write(fmt(21:22),'(i2)') 2 + max(2,nes)
        write(6,fmt)'COEF',itn-1,npar
!    I ,ndat,fil(1),fil(2),phi,pos0,anuf,fil(4),fil(3)
     I ,ndat,fil(1),fil(2),phi,pos0, fil(4),fil(3)
     I ,(vlp(jp),jp=2,npoly) 
     I ,'   phi_deg',phi*pi2deg,'   edge_pos',vm,'   MTFhNyq',xhn
       
      if(mprbay.gt.0) then
       if(uqspl) then
        if(unats) then
         write(6,'(20f10.5)') (yqr(i),i=2,nqm-1)
         if(uimag)
     I   write(6,'(20f10.5)') (yqi(i),i=2,nqm-1)
        else
         write(6,'(20f10.5)') (yqr(i),i=1,nqm0),dyr1,dyrn
        endif
         write(6,'(a)')'------------------------------OTF y:'
       endif
         write(6,'(20f10.5)') (yr(i),i=1,ncd)
         if(uimag)
     I   write(6,'(20f10.5)') (yi(i),i=1,ncd)

       if(mprplo.gt.0) call mkplotf
       if(mprplo.gt.1) call mkploth
       call mkplotc
       if(ucurv) call mkplotz
      endif
       if(nnull.gt.0) then
         write(6,'(/a,i5,a)')'Warning: ROI contains',nnull,
     I  ' pixel values <=0 , results are not trust worthy!'
         write(6,'(2a)')
     I  ' a smaller ROI may avoid the problematic sub-region ',
     I  '(try a reduced np value)'
       endif
       if(conv.gt.0.03) then
        write(6,'(/a,f10.6,a)')'Warning: OTF not converged',conv,
     I  ', results are wrong!'
        xm50=0
       elseif(conv.gt.0.001 ) then
        write(6,'(/a,f10.6,a)')'Warning: OTF not converged',conv,
     I  ' results are not trust worthy or wrong!'
        xm50=0
       elseif(mprbay.gt.0) then
        if(apixm.gt.0) then
         i=1
         do while( yrfm(i) .gt.0.5 )
          i=i+1
         enddo
         y2 = yrfm(i)
         i=i-1
         y1 = yrfm(i)
         p = (0.5-y1)/(y2-y1)
         xmm = cpix(i) +p/np
         write(6,'(/a,9(f6.3,a))')'MTFhNyq=',xhn,'   MTF50=',xm50,
     I ' c/p  MTF50ref=',xmm,' c/p'
        else
         write(6,'(/a,9(f6.3,a))')'MTFhNyq=',xhn,'   MTF50=',xm50,
     I ' c/p'
        endif
       endif
       call checkgap

       return
       end

      subroutine postprob( yfun )
!  posterior probability given the observation, div by marginal likelihood
! Author: Bernard Delley 2024
      use image_data
      use esf_data
      use fit_data
      use plot_data
      logical activ,record,chk1
      double precision yfuns

      yfuns = 0
      i=roip
      if(uqspl)then
        do i=1,npar
          do j=1,ij
            udes(j,i) = 0
          enddo
        enddo
      endif
!     nnull=0

      bs = 1.0/ns
      ij = 0   !  data point
      if(keep_roi) ij=1
      k = 0
      jroip = sqrt( roip*roip + iroin*iroin ) +3
       iw1 = max( 1, nint(vm(1)) -jroip )
       iw2 = min(iw, nint(vm(1)) +jroip )
       ih1 = max( 1, nint(vm(2)) -jroip )
       ih2 = min(ih, nint(vm(2)) +jroip )
      do i=ih1,ih2
        do j=iw1,iw2
          raw = ipix(j  ,i  ) ! - iblack
          ii=mod(i+1,2)
          jj=mod(j+1,2)
          kk=2*ii+jj+1
          dxp = j-vm(1)
          dyp = i-vm(2)
          if(kk.eq.2 .or. kk.eq.3 .or. monochrom.eq.1) then !  green
            xpa = dxp*ve(1) + dyp*ve(2)
            xno = dxp*vn(1) + dyp*vn(2)
            xno = -xno   ! change sign back to original
            activ=.false.
            if(keep_roi) then
              if(i.eq.ipxy(1,ij) .and. j.eq.ipxy(2,ij)) activ=.true.
!      if(ij.eq.1) write(6,'(a,l2,9i9)')'ck_BaY',keep_roi,ij,
!    I j,ipxy(2,ij),i,ipxy(1,ij),iw1,iw2,ih1,ih2
            else
              if( abs(xpa).lt. roip .and. abs(xno).lt.iroin ) then  !  i,j  is in ROI
                activ=.true.
                ij = ij + 1
                ipxy(1,ij) = i
                ipxy(2,ij) = j
!      if(ij.eq.1) write(6,'(a,l2,9i9)')'ck_BAY',keep_roi,ij,
!    I j,ipxy(2,ij),i,ipxy(1,ij),iw1,iw2,ih1,ih2
                if(raw.le.0.) then
                  nnull=nnull+1
!                 write(6,'(a,3i5)')'ck_W',i,j,nnull
                endif
              endif
            endif
            if(activ) then
              edi = xno
              if(npoly.ge.2) then   ! polynom correction along edge position 251204bd
                xpar = xpa/roip
!               xpar = xpa/iroip
!               call Legen1(xpar,pl,npoly)
!               call Fourier(xpar,pl,npoly)
                do jp=2,npoly
                  yes(jp) = 1
                  call splint(xes,yes,y2es(1,jp),nes,xpa,pl(jp)) ! do spline basisfns at xpar, for derivatve later
                  yes(jp) = 0
                  edi = edi + vlp(jp)*pl(jp)
                enddo
!               write(6,'(a,i5,99f8.3)')'ck_pl',ij,xpa,(pl(jp),jp=2,npoly)
              endif
              edins = max(-float(nm), min( float(nm), edi*ns) ) +nm       ! shifted
              iedi = edins                               ! rounded down
              edip = edins - iedi                        ! partial
              iedi = max( -nm, min( nm-1, iedi-nm))      ! shifted back to centered

              edin  = edi*ns
              iedj = iedi + 1

              csfi0=0
             edip0 = edip
             edip = max( 0., edip)
                pq = (1-edip)*edip*bs
                c0 = csf(iedi)*bs +esf(iedi) -esf(iedj)
                c1 = -csf(iedj)*bs -esf(iedi) +esf(iedj)
                esfi0 = (1-edip)*esf(iedi) +edip*esf(iedj)
                esfi  = (1-edip)*(esf(iedi) +pq*c0)
     I               +  edip*(esf(iedj) +pq*c1)
                csfi0 = (1-edip)*csf(iedi) +edip*csf(iedj)

      hypo =fil(1) +esfi*fil(2)*(1 +0.001*xpa*fil(3) +scalen*xno*fil(4)) 

              hypo0 = fil(1) +esfi*fil(2)
              hypod = hypo - hypo0       ! hypod contains illumination gradient correction for raw in g4-g7
!             hypod = esfi*fil(2)*( +0.001*( xpa*fil(3) ))

              dev = hypo-raw  ! x,yt
!             rms2 = raw*gain !   approach a) my not so good standard before 250803bd
!             rms2 = hypo*gain !  approach b,c)   starting with l8331l
      rms2 = hypo*gain + (hypo*fp)**2 + f0**2 + (csfi0*fil(2)*fh)**2  ! approach b with fp f0 fh 251021bd, 260102bd
      rms2 = max( 1.0, rms2 )
         
      postlog  = 0.5*( dev*dev/rms2 -1)
              yfuns = yfuns + postlog

              rmsd = sqrt(rms2)
            if(ihys.gt.1) then
              rk2(1,ij) = raw-rmsd  -hypod
              rk2(2,ij) = raw+rmsd  -hypod
              xk1(ij,1)=edi
              xk1(ij,2)=raw -hypod 
              xk1(ij,3)=hypo
              xk1(ij,4)=xpa
              xk1(ij,5)=rmsd !csfi0*fil(2)*fh
            endif
              rmsi(ij) = 1/rmsd   ! sqrt(rms2)
              brh(ij) = (raw-hypo)*rmsi(ij)
!             if(usigd) rmsi(ij) = 1/rmsd + brh(ij)*rmsi(ij)*gain   ! with sigma derivative, without csfi0 term

            if(npar1.gt.0) then
!             udes(ij,1) =        rmsi(ij) *scale1       ! d/dfil1      illum
              udes(ij,1) = fil(1)*rmsi(ij) *scale1       ! d/dfil1      illum
              udes(ij,2) = fil(2)*esfi*rmsi(ij)          ! d/dfil2
              udes(ij,3) = -csfi0*fil(2)*rmsi(ij)*xpa*scalan ! d/dphi 
!             udes(ij,3) = -csfi0*fil(2)*rmsi(ij)*xpa        ! d/dphi 0.1 for spectrum balance with damp2=0.1 or 0.159
              udes(ij,npx) = csfi0*fil(2)*rmsi(ij)  ! d/xno   -> vm1,2  edge pos
              udes(ij,nph) = fil(2)*esfi*rmsi(ij)*xno*scalen ! d/dfil4   illum/brightness grad normal to edge "nph" 
              udes(ij,npg) = fil(2)*esfi*rmsi(ij)*xpa*0.001 ! d/dfil3   illum/brightness grad parallel to edge "npg"
              jpp = npc      ! npc corresponds to first edge curvature parameter (jp=2) 251204bd
              do jp = 2,npoly
                udes(ij,jpp) = csfi0*fil(2)*rmsi(ij)*pl(jp)   ! d/dvlp(jp)  edge curve pl(jp)
                jpp = jpp + 1
              enddo
             if(uanum) then
              udes(ij,npk) = 0                              ! aperture model parameter      "npk" >=4
              do jp=1,ncd
                u1 = fil(2)*((1-edip)*fsf(jp,iedi,1)              ! desf/dyr(jp)
     I                         + edip*fsf(jp,iedj,1))*rmsi(ij)
                u2 = papt(jp) *dcirc(jp) *jp*bpix
                udes(ij,npk) = udes(ij,npk) + u1*u2
              enddo
                udes(ij,npk) = udes(ij,npk) *alam*1.e-6
             endif
            endif  ! npar1 > 0

            if(uqspl) then  ! spline yqr yqi models
             if(npar.gt.npar1) then
               u0 = fil(2)*rmsi(ij)
               call udesqsplin( iedi, iedj, ij, u0, edip)
             endif
!           elseif(npar.ge.npar1+ncd) then   !   straight yr,yi models
!              u0 = fil(2)*rmsi(ij)
!              call udesyryi( iedi,iedj, ij, u0, edip)
            else
              jpp=npar1 !0
             if(npar.ge.npar1+ncd) then
              if(.not.onlyi) then
              do jp=1,ncd
                jpp=jpp+1
                udes(ij,jpp) = fil(2)*(
     I          (1-edip)*   fsf(jp,iedi,1)
!    I             +pq*(dsf(jp,iedi,1) +fsf(jp,iedi,1) -fsf(jp,iedj,1)))    ! when ! remove ( before fsf
     I           + edip *   fsf(jp,iedj,1)
!    I             +pq*(-dsf(jp,iedj,1)-fsf(jp,iedi,1) +fsf(jp,iedj,1)))
     I          ) *rmsi(ij)
!      write(6,'(a,4i5,f12.2,9f15.9)')'ck_u',ij,jp,iedi,iedj,
!    I fil(2)*rmsi(ij), ((1-edip)*fsf(jp,iedi,1) + edip*fsf(jp,iedj,1))
!    I ,1-edip,fsf(jp,iedi,1),edip,fsf(jp,iedj,1)
              enddo
              endif
             if(npar.gt.npar1+ncd .or. onlyi) then
              do jp=1,ncd
                jpp=jpp+1
                udes(ij,jpp) = fil(2)*(
     I          (1-edip)*   fsf(jp,iedi,2)
!    I             +pq*(dsf(jp,iedi,2) +fsf(jp,iedi,2) -fsf(jp,iedj,2)))    ! when ! remove ( before fsf
     I           + edip *   fsf(jp,iedj,2)
!    I             +pq*(-dsf(jp,iedj,2)-fsf(jp,iedi,2) +fsf(jp,iedj,2)))
     I          ) *rmsi(ij)
              enddo
             endif ! npar.ge.npar1+2*ncd
             endif ! npar.ge.npar1+ncd
            endif ! uqspl
      if(mprbay.gt.5) then
!     if(postlog.gt.36) then
      write(6,'(a,2i5,i7,i3,i6,2f7.1,f9.4,5f12.1,9f12.3)')'ck_E',
     I  i,j,ij
     I ,kk,iedi ,xpa,xno,esfi,hypo,raw,dev,postlog,yfun
!    I ,(udes(ij,jp),jp=ncd+1,ncd+5)
      endif
!     if(mprplo.gt.4) then   ! check raw intensity in relative coord of edge: reasonable edge guess?
!     if(dopos) then
!     dd = diam*sqrt(raw/(fil(1)+fil(2)))
!     if(brh(ij).gt.throutly) then
!       dd = max( dd, diamin*2)
!       call symbrgb(1, xno, -xpa, dd, 1.0, 0.2, 0.2)
!     elseif(brh(ij).lt.-throutly) then
!       dd = max( dd, diamin*2.5)
!       call symbrgb(1, xno, -xpa, dd, 0.0, 0.7, 1.0)
!     else
!       dd = max( dd, diamin)
!       call symbrgb(1, xno, -xpa, dd, 1.0, 1.0, 1.0)
!     endif
!     endif  ! dopos
!     endif

            if(keep_roi) ij=ij+1    ! prepare for next in this mode
            endif  ! ROI
          endif ! green
        enddo
        if(.not.keep_roi) ndat = ij
      enddo

!     write(6,'(a,L2,9i9)')'ck_ndat',keep_roi,ndat,ij-1
      if(mprbay.gt.5) then
       write(6,'(29x,2a7,a9,9a12)')'xpa','xno','esf','hypo',
     I 'raw','dev','post','yfun'
      endif

      if(ij.gt.nwh) then
        write(6,'(a,i10)')'Error, insufficient allocation nwh',nwh
        stop 'error_nwh'
      endif

      yfun = yfuns

      return
      end

      function sinc(x)

      if(x.eq.0) then
        sinc = 1.
      else
        sinc = sin(x)/x
      endif
      return
      end

      subroutine out_sfr( ktrp, kedg, phis )
! Author: Bernard Delley 2024
      use esf_data
      use image_data
      real sfr(0:63)
      ycurr=0
      sfr(0) = 1
      if(nm.lt.63) then
        write(6,'(a,3i5)')'Warning, nm too small for out_sfr np,ns,nm'
     I ,np,ns,nm
        return
      endif
!     nmax=min( 63, nm)
!     if(np.eq.64) then
       sfmx=0
       do i=1,63
        sfr(i) = min( 1., sqrt( yr(i)**2 + yi(i)**2 ) )
!       sfmx = max( sfmx, sfr(i))
       enddo
!     endif
!     if(sfmx.gt.2.01) then
!       do i=1,nmax
!         sfr(i) = 0
!       enddo
!     endif
      write(36,'(i3,3f10.3,f7.2,64f7.4)')
     I  ktrp,vm,phis,edge(5,kedg,ktrp),sfr
      call flush(36)
      return
      end

      subroutine mkplotm
! Author: Bernard Delley 2026
!  plot map of tetrahedrons with detection number and MTF @ Nyq/2 color codes
      use image_data
      use esf_data
      use plot_data
      real  rgb(3)
      character*5 word
      character*9 wcoord(4)
      xmin = 1
      xmax = iw
      ymin = -ih
      ymax = -1

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)
      do i=1,ny
        yval(i) = -ytick(i)
      enddo
      xfmt='i4'
      yfmt='i4'
      xle = (xmax-xmin)  
      yle = (ymax-ymin)  
      scy = 160./yle
      xle = xle*scy
      yle = yle*scy
!     write(6,'(a,9f9.3)')'ck_scy',xle,yle,scy
      if(xle.gt.240) then
        scx = 240./xle
        xle = xle*scx
        yle = yle*scx
!     write(6,'(a,9f9.3)')'ck_scx',xle,yle,scx
      endif
      dd = 0.35*yle/24.

       call frameb(xmin,xmax, xle, nx, xtick, xtick, xfmt,
     I             ymin,ymax, yle, ny, ytick, yval, yfmt,
     I 1.5, 0.5, 'X','Y', 'g_map.eps')

        da(1) = xmin
        db(1) = ymin
        da(2) = xmax
        db(2) = db(1)
        da(3) = da(2)
        db(3) = ymax
        da(4) = xmin
        db(4) = ymax
!     call areargb(4,da,db, 0.8,0.8 ,0.8)        ! background

       write(98,'(a)')'TFnt'
      sat=1
      val=0.75
      do ktrp=1,ktrpm
        do k=1,4
          da(k) = edge(6,k,ktrp)
          db(k) = -edge(7,k,ktrp)
          if(ih.lt.400) then
            write(wcoord(k),'(ai3,a,i3,a)')
     I     '(',nint(edge(7,k,ktrp)),'-',nint(edge(6,k,ktrp)),')'
          endif 
        enddo
        call areargb(4,da,db, 0.5,0.5 ,0.5)      ! trapezoid
        cgx=0
        cgy=0
        do k=1,4
          if(ih.lt.400) then
            write(98,'(a)')'1 0 0 sc'
            dx=0
            dy=0
            call PSshow( da(k),db(k),dx,dy,wcoord(k))
          endif
         cgx = cgx +edge(1,k,ktrp)
         cgy = cgy -edge(2,k,ktrp)
         hue=280*edge(8,k,ktrp)
         call hsl2rgb( rgb(1), rgb(2), rgb(3), hue, sat, val)        
         call symbrgb(1, edge(1,k,ktrp), -edge(2,k,ktrp), dd,
     I   rgb(1), rgb(2), rgb(3))
!        write(6,'(a,i4,i2,7f8.1,2f8.3)')'ck_e',ktrp,k,da(k),db(k),
!    I   edge(1,k,ktrp),edge(2,k,ktrp),cgx/k,cgy/k,hue,edge(8,k,ktrp)
!    I   ,xmax
        enddo
         cgx=cgx*0.25
         cgy=cgy*0.25
         write(word,'(a,i3,a)')'(',ktrp,')'
       write(98,'(a)')'0 0 0 sc'
         dx=-6
         dy=-1
         call PSshow( cgx,cgy,dx,dy,word)
      enddo
       
      call endfrm( nx, xtick, ny, ytick)

      return
      end

      subroutine mkplot0( ktrp, kedg )
! Author: Bernard Delley 2024
! plot ROI in image
      use image_data
      use esf_data
      use plot_data
      use fit_data, only : domap

      if(domap) return
      j1 = max(  1, nint(vm(1))-255 )
      j2 = min( iw, nint(vm(1))+255 )
      i1 = max(  1, nint(vm(2))-255 )
      i2 = min( ih, nint(vm(2))+255 )
      xmin = j1
      xmax = j2
      ymin = -i2
      ymax = -i1
       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)
      do i=1,ny
        yval(i) = -ytick(i)
      enddo
      xfmt='i4'
      yfmt='i4'
!     xle = (xmax-xmin)  
!     yle = (ymax-ymin)  
!     gle = max(xle,yle)
!     xle = xle*150./gle
!     yle = yle*150./gle

      xle = (xmax-xmin)
      yle = (ymax-ymin)
      scy = 160./yle
      xle = xle*scy
      yle = yle*scy
!     write(6,'(a,9f9.3)')'ck_scy',xle,yle,scy
      if(xle.gt.240) then
        scx = 240./xle
        xle = xle*scx
        yle = yle*scx
!     write(6,'(a,9f9.3)')'ck_scx',xle,yle,scx
      endif


       call frameb(xmin,xmax, xle, nx, xtick, xtick, xfmt,
     I             ymin,ymax, yle, ny, ytick, yval, yfmt,
     I 1.5, 0.5, 'X','Y', 'g_ROI_orientation.eps')

        da(1) = xmin
        db(1) = ymin
        da(2) = xmax
        db(2) = db(1)
        da(3) = da(2)
        db(3) = ymax
        da(4) = xmin
        db(4) = ymax
      call areargb(4,da,db, 0.3,0.3 ,0.3)          ! background

      if( max(i2,j2) .gt. 250) then
        ifac=7
      elseif( max(i2,j2) .gt. 100) then
        ifac=3
        ifac=1
      else
        ifac=1
      endif
      do i=i1,i2,ifac
        do j=j1,j2,ifac
          ii=mod(i+1,2)
          jj=mod(j+1,2)
          kk=2*ii+jj+1
          if(monochrom.eq.1) then ! full lattice
            raw = ipix(j  ,i  ) ! - iblack
            dd = 0.9*ifac*sqrt(raw/(fil(1)+fil(2)))*xle/(j2-j1)
            da1 = +j
            db1 = -i
            call symbrgb(1, da1, db1, dd, 0.9, 0.9, 0.9) ! white dots
          elseif(kk.eq.2 .or. kk.eq.3) then 
            raw = ipix(j  ,i  ) ! - iblack
            dd = 1.2*ifac*sqrt(raw/(fil(1)+fil(2)))*xle/(j2-j1)
            da1 = +j
            db1 = -i
            call symbrgb(1, da1, db1, dd, 0.9, 0.9, 0.9) ! white dots
          endif
        enddo
      enddo
       anor =  edge(3,kedg,ktrp)  ! length of edge
        da(1) = +vm(1) - 0.5*anor*ve(1)
        db(1) = -vm(2) + 0.5*anor*ve(2)
        da(2) = +vm(1) + 0.5*anor*ve(1)
        db(2) = -vm(2) - 0.5*anor*ve(2)
      call linrgb(2,da,db, 0.0, 0.7, 0.0)
        da(1) = +vm(1) 
        db(1) = -vm(2) 
        da(2) = +vm(1) + 50*vn(1)
        db(2) = -vm(2) - 50*vn(2)
!       write(6,'(a,2f8.1,2f9.4)')'ck_v',vm,vn,ve(1),ve(2)
!     call linrgb(2,da,db, 0.0, 0.7, 0.0)
      call symbrgb(1, da, db, 3., 0.8, 0.0, 0.8)  ! pink dot
        da(1) = +vm(1) -ve(1) *roip -vn(1)*iroin
        db(1) = -vm(2) +ve(2) *roip +vn(2)*iroin
        da(2) = +vm(1) -ve(1) *roip +vn(1)*iroin
        db(2) = -vm(2) +ve(2) *roip -vn(2)*iroin
        da(3) = +vm(1) +ve(1) *roip +vn(1)*iroin
        db(3) = -vm(2) -ve(2) *roip -vn(2)*iroin
        da(4) = +vm(1) +ve(1) *roip -vn(1)*iroin
        db(4) = -vm(2) -ve(2) *roip +vn(2)*iroin
        da(5) = da(1)
        db(5) = db(1)
      call linrgb(5,da,db, 0.9, 0.0, 0.0)    ! ROI

      call endfrm( nx, xtick, ny, ytick)

      return
      end

      subroutine mkplotf
! Author: Bernard Delley 2024
! plot OTF real and imaginary part, +MTF
      use esf_data
      use plot_data
      use fit_data
      use plot_data_compare

      xmin = 0 ! -cpmm(1) ! 0
        nmax = ncd
      xmax = fcut*1.03
      ymin = -0.15 ! -0.1
      ymax =  1.00

      if(mprbay.gt.5) then
       write(6,'(a,i6,2f10.3)')'ck_max_plot_freq',np,cpmm(np)
      endif

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)
       yfmt='f3.1'

       call frameb(xmin,xmax, xlen, nx, xtick, xtick, xfmt,
     I             ymin,ymax, ylen, ny, ytick, ytick, yfmt,
!    I 1.5, 0.3, 'c/mm','MTF, OTF', 'g_MTF_OTF.eps')
     I 1.5, 0.3, 'c/mm','     MTF,)sr( OTF', 'g_MTF_OTF.eps')

      da(1) = cpmm(np/2)  ! Nyqvist
      da(2) = da(1)
      db(1) = ymin
      db(2) = ymin*0.95 + ymax*0.05
      write(98,'(a)')'flw'
      call linrgb( 2, da, db, 0.5, 0.5, 0.5)     ! Nyqvist mark cpmm(np/2)
      da(1) = cpmm(np)  ! 2*Nyqvist = 1c/pix
      da(2) = da(1)
      db(2) = ymin*0.96 + ymax*0.04
      call linrgb( 2, da, db, 0.5, 0.5, 0.5)     ! 1 c/pix
      da(1) = cpmm(2*np)  ! 2c/pix
      da(2) = da(1)
      db(2) = ymin*0.96 + ymax*0.04
      call linrgb( 2, da, db, 0.5, 0.5, 0.5)     ! 2 c/pix
      da(1) = fcut
      da(2) = da(1)
      db(1) = -0.05
      db(2) = 0.05
      call linrgb( 2, da, db, 0.0, 0.0, 1.0)     ! fcut mark
      da(1)=xmin
      da(2)=xmax
      db(1)=0
      db(2)=0
      write(98,'(a)')'slw'
      call linrgb( 2, da, db, 0.0, 0.0, 0.0)     ! zero line
      write(98,'(a)')'flw'

      nplo=nmax+1

      if(apixm.gt.0) then
      afn(0)=1
      do i=1,nmax
        afn(i)=yrfm(i)
!       bfn(i)=0
      enddo
      call linrgb(nplo, cpmm, afn, 0.7, 0.7, 0.7)  !   yrfm : given model
      endif

      fcorr=1
      cfn(0)=1
      dfn(0)=1
      do i=1,nmax
        bfn(i)=yr(i)*fcorr
        afn(i)=yr(i)*fcorr
        dfn(i)=sqrt( yr(i)**2 +yi(i)**2 ) !*sign(1.,yr(i)) *fcorr
        cfn(i)=yi(i)*fcorr
      enddo
!     call symbrgb(nplo, cpmm, cfn, 1.0, 1.0, 0.0, 0.0)    ! MTF            Bayesian inferred   
      rap = min(1.0,46./np)
        ram = rap*0.7
      if(uimag)
     Icall symbrgb(nplo, cpmm, cfn, rap, 0.0, 0.7, 0.0)    ! imag part OTF  Bayesian inferred
      call symbrgb(nplo, cpmm, afn, rap, 0.9, 0.0, 0.0)    ! real part OTF  Bayesian inferred
      call symbrgb(nplo, cpmm, dfn, rap, 0.0, 0.0, 0.0)    ! abs OTF  Bayesian inferred

!     if(umapper) then
C generate like: mtf_mapper  pix_.ppm . -q  --bayer green  --esf-model loess --single-roi
!     open(35,file='edge_sfr_values.txt',iostat=ierr)
!     if(ierr.eq.0) then
!       read(35,*)(ymap(i),i=1,5),(ymap(i),i=0,63)
c       call symbrgb(64,cpmm,ymap,ram, 0.5, 0.7, 1.0)
!     endif
!     endif

      call endfrm( nx, xtick, ny, ytick)
      if(mprbay.gt.5) then
      write(6,*)'end making plotf-1',ncd,np
      endif

!     if(nqm.gt.0 .and. unats) then ! unats MTF plot with error bar
      xmin= 0 ! -0.05
      xmax = xmax*pix ! cpix(np)
      if(ncd.gt.np)then
       if(ncd*1.03.le.nm) then
        xmax = cpix(ncd)*1.03
       else
        xmax = cpix(nm)
       endif
      endif
      da(1) = cpix(np/2)
      da(2) = da(1)

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call frameb(xmin,xmax, ylen, nx, xtick, xtick, xfmt,
     I             ymin,ymax, ylen, ny, ytick, ytick, yfmt,
     I 1.5, 0.3, 'c/p','MTF','g_MTF.eps')

      da(1) = xmin
      da(2) = xmax
      db(1) = 0
      db(2) = 0
      call linrgb( 2, da, db, 0.5, 0.5, 0.5)     ! zero line

         call symbrgb(nplo, cpix, dfn, ram, 0.0, 0.0, 0.0)
      write(98,'(a)')'flw'
      da(1) = fcut*pix
      da(2) = da(1)
      db(1) = -0.1
      db(2) =  0.1
      call linrgb( 2, da, db, 0.0, 0.0, 1.0)     ! fcut mark
!     da(1) = xmin
!     da(2) = xmax
!     db(1) = 0
!     db(2) = 0
!     call linrgb( 2, da, db, 1.0, 0.0, 0.0)     ! zero line
      write(98,'(a)')'slw'
      call endfrm( nx, xtick, ny, ytick)
      if(mprbay.gt.5) then
      write(6,*)'end making plotf-mtf'
      endif
      
      if(apixm.gt.0) then
      ymin = -0.051 !-0.1 ! -0.05
      ymax =  0.051 ! 0.1 ! 0.05
       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)

       call frameb(xmin,xmax, xlen, nx, xtick, xtick, xfmt,
     I             ymin,ymax, ylen, ny, ytick, ytick, yfmt,
     I 1.5, 0.3, 'c/pix','magnified difference SFR',
     I 'g_OTF_deviation_magnified.eps')

      do i=1,nqm
      da(1) = xq(i)/np
      da(2) = da(1)
      db(1) = ymin
      db(2) = ymax
      call linrgb( 2, da, db, 0.5, 0.5, 0.5)     ! xq
      enddo
      da(1) = fcut*pix
      da(2) = da(1)
      db(1) = -0.01
      db(2) = 0.01
      write(98,'(a)')'flw'
      call linrgb( 2, da, db, 0.0, 0.0, 1.0)     ! fcut mark
      write(98,'(a)')'slw'

      afn(0)=0
      bfn(0)=0
      cfn(0)=0
      if(nqm.gt.0) then  ! show error for spline envelope
       if(unats) then
        if(uimag) then
         do i=1,nqm0
          afn(i-1) = xq(i+1)/np
          cfn(i-1) = -rv1(i+nqm0)
          bfn(i-1) = +rv1(i+nqm0)
         enddo   
!        call area_bc(nqm0,afn,bfn,cfn, 0.5, 0.5, 0.6)
        endif
         do i=1,nqm0
          afn(i-1) = xq(i+1)/np
          cfn(i-1) = -rv1(i)
          bfn(i-1) = +rv1(i)
         enddo   
         call area_bc(nqm0,afn,bfn,cfn, 0.8, 0.8, 0.9)
       else
        do i=1,nqm
          afn(i-1) = xq(i)/np
          cfn(i-1) = -rv1(i)
          bfn(i-1) = +rv1(i)
        enddo   
        call area_bc(nqm,afn,bfn,cfn, 0.8, 0.75, 0.8)
       endif
      else
       do i=1,nmax ! min( ncd, np)
        cfn(i) = -rv1(i)  
        bfn(i) = +rv1(i)  
!       write(6,'(a,i5,f8.2,9f15.5)')'ck_Y',i,cpmm(i),udes(i,1)
       enddo
       call area_bc(nplo,cpix,bfn,cfn, 0.8, 0.8, 0.8)
      endif ! nqm splines

      da(1) = fcut*pix
      da(2) = da(1)
      db(1) = -0.01
      db(2) = 0.01
      write(98,'(a)')'flw'
      call linrgb( 2, da, db, 0.0, 0.0, 1.0)     ! fcut mark
      write(98,'(a)')'slw'

      if(mprbay.gt.2) then
      write(6,'(a,5i5,2f9.2)')
     I'mkplotf MTF MTF-MTF_ref',np,nm,ncd,nmax,nplo,ncd*bpix,xmax
      endif

      afn(0)=0
!     if(ierr.eq.0 .and. umapper) then  ! mtf_mapper for comparison
!         ymap(0) = 0
!       do i=1,63
!         afn(i)=i/64.
!         ymap(i) = ymap(i) -yrfmap(i)
!       enddo
!       call linrgb(64, afn,ymap, 0.3, 0.3, 1.0)
!     endif
      cfn(0)=0
!     ymap(0)=0
      do i=1,nmax
!       afn(i) = dfn(i) - yrfm(i)
        afn(i) = dfn(i) - abs(yrfm(i))
!       ymap(i) = ymap(i) - yrfm(i)
        cfn(i) = yr(i) - yrfm(i)
        bfn(i) = yi(i)
      enddo
c     call symbrgb(nplo, cpix, afn, rap, 0.0, 0.0, 0.0)
      call symbrgb(nplo, cpix, cfn, ram, 0.9, 0.0, 0.0)
      if(uimag)
     Icall symbrgb(nplo, cpix, bfn, ram, 0.0, 0.7, 0.0)

!     call linrgb( 2, da, db, 0.0, 0.0, 1.0)     ! 

      call endfrm( nx, xtick, ny, ytick)
      endif  ! apixm>0

      if(mprbay.gt.5) then
      write(6,*)'end making plotf-2',ncd,np
      endif

      return
      end

      subroutine checkgap
! Author: Bernard Delley 2025
      use esf_data
      use plot_data
      use fit_data

      dmax=0
      x0=drh(1)
      do i=2,ndat
        x1 = drh(i)    ! already ordered, positions
!       if(x1.gt.-5 .and. x1.lt.+5) then
          if( x1-x0 .gt. dmax) then
            dmax = x1-x0
            xg=x1
          endif
!         dmax = max(dmax, x1-x0)
!       endif
        x0=x1
      enddo

      if(mprbay.gt.0 .or. 0.5/dmax .lt. fcut*pix) then
       write(6,'(/a,f8.5,a,f6.1,5(a,f8.3))')
     I 'Sampling gap',dmax,'p   at',xg,'   Safe_lim_freq',0.5/dmax,
     I 'c/p   N_f lim_freq',fcut*pix,'c/p' !   angle',phi*pi2deg,' deg'
      endif

      if(0.5/dmax .lt. fcut*pix) then
       write(6,'(/2a)')'WARNING: Safe_limit_frequency too low!'
     I ,'   This can be due to different causes:'
       write(6,'(2a)')'- the slant angle is problematic, see ',
     I'supplementary Table 1 and inspect g_ESF_magnified.eps graphic'
       write(6,'(2a)')'- ROI extends beyond image limits, ',
     I' inspect g_ROI_orientation.eps graphic'
!     if(abs(xg).gt.8.) write(6,'(2a)')'  As the sampling gap is far',
!    I' out, the Result may be OK if g_OTF_magnified.eps is OK'
       write(6,'(a,f6.1)')'- ROI nquer too small ',2*roip
!    I ' for the N_f value input, set from prior knowledge'
!      write(6,'(2a)')'inspect g_ESF_magnified.eps graphic'
!      if( anum*fcut*pix*dmax*2 .lt. 7 ) then
!       write(6,'(a,f7.2,a)')'by setting N_f to:',
!    I  anum*fcut*pix*dmax*2,
!    I' will truncate high frequencies and limit N_f lim_freq by Prior'
!      else
!       write(6,'(2a)')'or measure with a slightly different angle',
!    I' further away from problematic angle'
!      endif      

      endif
        write(6,'(a)')
       
      return
      end

      subroutine mkploth
! Author: Bernard Delley 2024, revs 2026 BD
!  plot ESF (H* MAP)
      use esf_data
      use plot_data
      use fit_data
      character*21 gfile

      xmax =  1 !  0.2  ! 0 !np+1
      xmin = -iroin-1 ! -np/2-2 
!     rms = sqrt(fil(2)/cgain)
!     ymin = max( -0.,fil(1)-2*rms)
      hypo = fil(1)+esf(-nm)*fil(2)
      rms = sqrt(max( (hypo/cgain + (hypo*fp)**2 +f0**2), hypo*0.1, 1.))
      db(1)= hypo
      ymin = max( 0., db(1) -6*rms)
!     ymax = fil(1)+3*rms
      ymax = ymin + 44*rms
      if(mprbay.gt.5) then
      write(6,'(a,2i8,9f9.2)')'ploth Minmax',ndat,nm,xmin,xmax,ymin,ymax
     I ,rms,hypo
      endif
      hypo =fil(1)+esf(nm)*fil(2)
      rms = sqrt(max( (hypo/cgain + (hypo*fp)**2 +f0**2), hypo*0.1 ))
!     write(6,'(a,9e12.3)')'ck_,min,max',rms,ymin,ymax

        db(3)=db(1)

      do iplo=1,4

      if(iplo.eq.1) then
        gfile='g_ESF_dark_side.eps'
        len=19
      elseif(iplo.eq.2) then
        gfile='g_ESF.eps'
        len=9
      elseif(iplo.eq.3) then
        gfile='g_ESF_magnified.eps'
        len=19
      else
        gfile='g_ESF_bright_side.eps'
        len=21
      endif

      if(mprbay.gt.5) then
      write(6,'(a,2i8,5f9.1)')'ploth minmax',ndat,nm,xmin,xmax,ymin,ymax
     I ,rms
      endif

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)

       call frameb(xmin,xmax, xlen, nx, xtick, xtick, xfmt,
     I             ymin,ymax, ylen, ny, ytick, ytick, yfmt,
     I 1.5, 0.3,'distance from edge [pixels]','data & hypothesis',
     I  gfile(1:len) )

      jj=nint(0.0027*(ndat-npar+nsing))
      if(iplo.eq.2) write(6,'(/a,i6,a,i4,a,i3,a,9x,a,i4,a)')
     I 'N',ndat,'  M*',npar-nsing,' <',jj,
     I ' >  position      sig   x_nor       raw*       hypo        xpa' 
     I,'expected number of 3 sigma outliers <',jj,' >'
      jj=0
      rsy = 0.4
      if(iplo.eq.3) rsy = 1.2
      rsyo = 1.5*rsy
      jj=0
      do ii=1,ndat
        i=ke1(ii)       ! make it ordered in position from edge
        da(1) = xk1(i,1)
!       da(2) = xk1(i,1)
        if( abs(brh(i)) .le. throutly) then  
           da(2)=xk1(i,2)
        call symbrgb(1,da(1), da(2), rsy, 0.6, 0.6, 0.6)    ! values
        else  ! mark outliers
           da(2)=xk1(i,2)
        call symbrgb(1,da(1), da(2), rsyo, 0.9, 0., 0.)
          if(iplo.eq.2) then
           jj=jj+1
          if(jj.le.20 .or. jj.eq.jjj) then
           write(6,'(a,i7,i7,i6,f8.2,f8.2,9f11.2)')'outlier_pixel',
     I     jj,i,ii,brh(i),(xk1(i,j),j=1,4)
          endif
          ijj=i
          iijj=ii
          endif
        endif
      enddo
      if(jj.gt.20) then
           i=ijj
           ii=iijj
           write(6,'(a,i7,i7,i6,f8.2,f8.2,9f11.2)')'outlier_pixel',
     I     jj,i,iijj,brh(i),(xk1(ijj,j),j=1,4)
      endif

      da(1)=xmin+0.1
      da(2)=xmax-0.1
      da(3)=da(1)

      if(iplo.eq.1) then ! add rms bar to the left
        da(3) = -iroin-0.2 
        da(4) = da(3)
!       rms1 = sqrt(fil(1)/cgain)
        rms1 = sqrt(db(1)/cgain)
        db(4) = db(1)-rms1
        db(5) = db(1)+rms1
        call linrgb(2, da(3), db(4), 0.7, 0.0, 0.7)   
        rms1 = sqrt(rms1*rms1 + f0**2)
        da(3) = -iroin-0.4 
        da(4) = da(3)
        db(4) = db(1)-rms1
        db(5) = db(1)+rms1
        call linrgb(2, da(3), db(4), 0.7, 0.0, 0.7)   
        da(3) = 0
        da(4) = da(3)
        db(4) = db(1)
        db(5) = db(1) -esfh*fil(2)
        call linrgb(2, da(3), db(4), 0.0, 0.5, 1.0)   
!       write(6,'(a,f10.5)')'ck_mkploth_esfh',esfh*fil(2)
      endif

!     if(iplo.eq.2 .or. iplo.eq.3) then
      if(iplo.eq.2 ) then
        db(2) = db(1)
        call linrgb(2, da, db, 0.0, 0.0, 0.0)      ! fil levels 
        db(4) = crh(1)
        db(5) = crh(1)
        call linrgb(2, da, db(4), 0.0, 0.0, 0.0)      ! crh "0" reference 
        call linrgb(ndat,drh,crh, 0.0, 0.6, 0.0)
        call linrgb(ndat,drh,erh, 0.9, 0.0, 1.0)
        db(3)=fil(1)+esf(nm)*fil(2)
      endif
        db(4)=db(3)
      if(iplo.ne.3)
     Icall linrgb(2, da, db(3), 0.0, 0.0, 0.0)      ! fil levels 

!     if(ihys.ge.5)
!     if( iplo.ne.3 .and. iplo.ne.2 ) then
!     call linrgb(na, xk2(-nm,1), xk2(-nm,2), 0.0, 0.0, 1.0)  ! ihys=2     stage1 model
!     call linrgb(na, xk2(-nm,1), xk2(-nm,5), 0.8, 0.0, 0.8)  ! plot yrfm, stage1 starting model
!     endif
      call linrgb(na, xk2(-nm,1), xk2(-nm,3), 0.0, 0.5, 0.0)  ! ihys=3     stage2 final model
      if(ihys.ge.4)
     Icall linrgb(na, xk2(-nm,1), xk2(-nm,4), 0.9, 0.0, 0.0)  !            stage3 final model

      if(iplo.eq.4) then ! add rms bar to the right
        da(1) = iroin+0.2 
        da(2) = da(1)
        rms1  = sqrt(db(3)/cgain)
        rms = sqrt(rms1*rms1 + (db(3)*fp)**2 )
!       write(6,'(a,9f15.3)')'ck_h',db(3),rms1,rms
        db(4) = db(3)-rms
        db(5) = db(3)+rms
        call linrgb(2, da, db(4), 0.7, 0.0, 0.7)   
        da(1) = iroin+0.4 
        da(2) = da(1)
        db(4) = db(3)-rms1
        db(5) = db(3)+rms1
        call linrgb(2, da, db(4), 0.7, 0.0, 0.7)   
      endif

      call endfrm( nx, xtick, ny, ytick)

      if(iplo.eq.1) then ! prepare for g5.ps
!      xmin = -np/2-2 !-iroin-2  !-1.01
       xmax = iroin+1 ! np/2+2 ! iroin+2
!      ymin = 0.5*fil(1)  ! 0.50*fil(2)
!      ymax = 1.10*fil(2)
       ymin = fil(1) -3*rms
       ymax = fil(1)+fil(2) +5*rms
      elseif(iplo.eq.2) then ! prepare for new  g6.ps
       xmin = -1.4 !3
       xmax = +1.4 !3
      elseif(iplo.eq.3) then ! prepare for new  g7.ps ~ former g6.ps
       xmin = -1
       xmax = iroin+1
       ymin = (fil1)+fil(2) -4*rms
       ymax = fil(1)+fil(2) +4*rms
!     else
!      xmin = -1 ! -0.2  !-1.01
      endif

      enddo ! iplo

      if(mprbay.gt.5) then
      write(6,'(a,4f8.2)')'end mkploth'
      endif

      return
      end

      subroutine mkplotc
! Author: Bernard Delley 2025
!  plot csf     LSF
      use esf_data
      use plot_data
      use fit_data

      xmin = -iroin -1
      xmax =  iroin +1
      ymin=-0.0011
      do i=-nm,nm
        ymin=min(ymin,csf(i))
      enddo
      ymax =  csf(0)*1.05
      ymin = -0.25  !2*ymin
      ymax = 0.5    !-ymin
      

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)

       call frameb(xmin,xmax, xlen, nx, xtick, xtick, xfmt,
     I             ymin,ymax, ylen, ny, ytick, ytick, yfmt,
     I 1.5,0.3,'distance from edge [pixels]','line spread function csf'
     I ,'g_LSF.eps')

      call linrgb(na, xk2(-nm,1), csf(-nm), 0.0, 0.2, 1.0)

      call endfrm( nx, xtick, ny, ytick)

      return
      end

      subroutine mkplotz
! show fitted blade edge
! Author: Bernard Delley 2025
      use esf_data
      use plot_data

!       write(6,'(a,25x,99f10.5)')'ck_z',(vlp(j),j=2,npoly)
      xmin = -0.03
      xmax =  0.03
      ymin = -iroip-1
      ymax =  iroip+1
      do i=2,nes-1
        yes(i) = vlp(i)
        xmin = min(xmin,vlp(i))
        xmax = max(xmax,vlp(i))
      enddo
      xmin = xmin*1.05
      xmax = xmax*1.05
      dye = 1.e30
      call spline(xes,yes,nes,dye,dye,y2es,dfn)     !  dfn <- u

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)

       call frameb(xmin,xmax, 100., nx, xtick, xtick, xfmt,
     I             ymin,ymax, 130., ny, ytick, ytick, yfmt,
     I 1.5, 0.3, 'distance pixels','distance pixels'
     I , 'g_edge_function.eps')

      do i=2,nes-1
        db(1) = xes(i)
        db(2) = xes(i)
        da(1) = xmin
        da(2) = xmax
        call linrgb(2,da,db,0.5,0.5,0.5)
      enddo
      dx=2*roip/128.
      xpa = -roip
      do i=0,128
        dfn(i)=xpa
        call splint(xes,yes,y2es,nes,xpa,afn(i))
        xpa = xpa + dx
      enddo
!     do i=0,128
!       dfn(i)=xpa
!       xpar = xpa/roip
!       afn(i) = 0
!       call Legen1(xpar,pl,npoly)
!       call Fourier(xpar,pl,npoly)
!       do jp=2,npoly
!          afn(i) = afn(i) + vlp(jp)*pl(jp)
!       enddo
!       write(6,'(a,i5,99f10.5)')'ck_z',i, dfn(i),afn(i)
!    I ,(pl(j),j=2,npoly)
!       xpa = xpa + dx
!     enddo
      call linrgb( 129, afn, dfn, 0.9, 0.0, 0.0)

      call endfrm( nx, xtick, ny, ytick)

      return
      end

      subroutine mkplotou
! image document outliers
! Author: Bernard Delley 2024
      use image_data , only : monochrom 
      use esf_data
      use plot_data
      use fit_data

      if(domap) return
      dopos = .true.

!     diam = 160.*sqrt(2.0-monochrom)/(2*max(iroip,iroin))
!     write(6,'(a,2i5,f10.5)')'ck_diam',iroip,iroin,diam
!     diamin = 0.25*diam*roip/32.
!     diamin = min( 0.85*diam, diamin)
!dd = 0.9*ifac*sqrt(raw/(fil(1)+fil(2)))*xle/(j2-j1)

      if(mprbay.gt.5) then
      write(6,'(a,4f8.2)')'begin mkplotou'
      endif
      ymin = -iroip-1
      ymax =  iroip+1
      xmin = -iroin-2
      xmax =  iroin+2
      xl1 = max( xmax, ymax)

      xle = (xmax-xmin)
      yle = (ymax-ymin)
      scy = 160./yle
      xle = xle*scy
      yle = yle*scy
!     write(6,'(a,9f9.3)')'ck_scy',xle,yle,scy
      if(xle.gt.240) then
        scx = 240./xle
        xle = xle*scx
        yle = yle*scx
!     write(6,'(a,9f9.3)')'ck_scx',xle,yle,scx
      endif
      if(monochrom.eq.1) then
        diam = 0.9*xle/(xmax-xmin)
!       write(6,'(a,f9.3)')'ck_A',diam
      else
        diam = 1.2*xle/(xmax-xmin)
!       write(6,'(a,f9.3)')'ck_a',diam
      endif
      diamin = 0.1 *diam*roip/32.
      diamin = min( 0.80*diam, diamin)

!     write(6,'(a,i5,4f9.2)')'ploto min,max',ndat,xmin,xmax,ymin,ymax

       call autotic( xtick, nx, xfmt, xmax, xmin)
       call autotic( ytick, ny, yfmt, ymax, ymin)
       do i=1,ny
         yval(i) = -ytick(i)
       enddo


       call frameb(xmin,xmax, xle, nx, xtick, xtick, xfmt,
     I             ymin,ymax, yle, ny, ytick, yval, yfmt,
     I 1.5, 0.3, 'distance pixels' ,'distance pixels',
     I 'g_ROI_map_outliers.eps')

        da(1) = xmin
        db(1) = ymin
        da(2) = xmax
        db(2) = db(1)
        da(3) = da(2)
        db(3) = ymax
        da(4) = xmin
        db(4) = ymax
      call areargb(4,da,db, 0.40, 0.40, 0.40)  ! black background

      do ij=1,ndat
        da(1) = xk1(ij,1)
       dd = diam*sqrt(xk1(ij,3)/(fil(1)+fil(2)))
       if(brh(ij).gt.throutly) then
        dd = 0.75*diam !max( dd, diamin*2)
        call symbrgb(1, xk1(ij,1), -xk1(ij,4) , dd, 1.0, 0.2, 0.2)
       elseif(brh(ij).lt.-throutly) then
        dd = 0.8*diam !max( dd, diamin*2.5)
        call symbrgb(1, xk1(ij,1), -xk1(ij,4) , dd, 0.0, 0.7, 1.0)
       else
        dd = max( dd, diamin)
        call symbrgb(1, xk1(ij,1), -xk1(ij,4) , dd, 1.0, 1.0, 1.0)
       endif
      enddo

!     return
!     end

!     subroutine mkplote
! Author: Bernard Delley 2024
!     use esf_data
!     use plot_data
!     use fit_data, only : domap

!     if(domap) return
      dopos = .false.
      call endfrm( nx, xtick, ny, ytick)
      if(mprbay.gt.5) then
      write(6,*)'end making ploto'
      endif

      return
      end

      subroutine indexx(n, arr, indx)

c     from numerical recipes Press, Flannery, Teukolsky and Vetterling
c     Cambridge Univ Press 1986, pg 233
C* Revision history: bugfix 952710 BD

!     implicit double precision (a-h,o-z)
      integer  indx(n)
      dimension   arr(n)

      do 11 j=1,n
       indx(j) = j
  11  continue

      L = n/2+1
      ir = n

  10  continue
      if(L.gt.1) then
        L = L-1
        indxt = indx(L)
        q = arr(indxt)
      else
        indxt = indx(ir)
        q = arr(indxt)
        indx(ir) = indx(1)
        ir = ir-1
        if(ir.le.1) then
          indx(1) = indxt
          return
        end if
      end if
      i = L
      j = L+L

  20  continue
      if(j.le.ir) then
        if(j.lt.ir) then
          if(arr(indx(j)).le.arr(indx(j+1))) j = j+1
        end if
        if(q.le.arr(indx(j))) then     !
          indx(i) = indx(j)
          i = j
          j = j+j
        else
          j = ir+1
        end if
        go to 20
      end if
      indx(i) = indxt
      go to 10

      end

      subroutine svdcmp(a,m,n,mp,np,w,v,rv1)

c     from numerical recipes Press, Flannery, Teukolsky and Vetterling
c     Cambridge Univ Press 1986, pg 233
c     mod b.d.

! rev 2023bd
! a   m,n matrix  m>n
! m
! n
! mp  leading dimension of a matrix  mp >= m
! np  leading dimension of v matrix  np >= n
! w   singular values vector
! v   n,n orthogonal matrix
! rv1 here internal work array

      dimension a(mp,np),w(np),v(np,np)
      dimension rv1(n)  ! work array
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if (i.le.m) then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if (scale.ne.0.0) then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            if (i.ne.n) then
              do 15 j=l,n
                s=0.0
                do 13 k=i,m
                  s=s+a(k,i)*a(k,j)
13              continue
                f=s/h
                do 14 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
14              continue
15            continue
            endif
            do 16 k= i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if ((i.le.m).and.(i.ne.n)) then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if (scale.ne.0.0) then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            if (i.ne.m) then
              do 23 j=l,m
                s=0.0
                do 21 k=l,n
                  s=s+a(j,k)*a(i,k)
21              continue
                do 22 k=l,n
                  a(j,k)=a(j,k)+s*rv1(k)
22              continue
23            continue
            endif
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if (i.lt.n) then
          if (g.ne.0.0) then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=n,1,-1
        l=i+1
        g=w(i)
        if (i.lt.n) then
          do 33 j=l,n
            a(i,j)=0.0
33        continue
        endif
        if (g.ne.0.0) then
          g=1.0/g
          if (i.ne.n) then
            do 36 j=l,n
              s=0.0
              do 34 k=l,m
                s=s+a(k,i)*a(k,j)
34            continue
              f=(s/a(i,i))*g
              do 35 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
35            continue
36          continue
          endif
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
!       do 48 its=1,30
        do 48 its=1,40
          do 41 l=k,1,-1
            nm=l-1
            if ((abs(rv1(l))+anorm).eq.anorm)  go to 2
!           if(nm.gt.0) then  ! w(nm=0) is not reached when l=1   231025bd except when uij=infty
            if ((abs(w(nm))+anorm).eq.anorm)  go to 1
!           else
!        write(6,'(a,3i5,2f15.6)')'Warning w(nm)',nm,l,k,its,w(nm),anorm
!           endif
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            if ((abs(f)+anorm).ne.anorm) then
              g=w(i)
              h=sqrt(f*f+g*g)
              w(i)=h
              h=1.0/h
              c= (g*h)
              s=-(f*h)
              do 42 j=1,m
                y=a(j,nm)
                z=a(j,i)
                a(j,nm)=(y*c)+(z*s)
                a(j,i)=-(y*s)+(z*c)
42            continue
            endif
43        continue
2         z=w(k)
          if (l.eq.k) then
            if (z.lt.0.0) then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            go to 3
          endif
!         if (its.eq.30) then
          if (its.eq.50) then
!          write(6,'(a)') 'no convergence in 30 iterations'
           write(6,'(a)') 'no convergence in 50 iterations'
           stop 'error no convergence'
          endif
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=sqrt(f*f+1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=sqrt(f*f+h*h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 nm=1,n
              x=v(nm,j)
              z=v(nm,i)
              v(nm,j)= (x*c)+(z*s)
              v(nm,i)=-(x*s)+(z*c)
45          continue
            z=sqrt(f*f+h*h)
            w(j)=z
            if (z.ne.0.0) then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 nm=1,m
              y=a(nm,j)
              z=a(nm,i)
              a(nm,j)= (y*c)+(z*s)
              a(nm,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      end

      subroutine svbksb(u,w,v,m,n,mp,np,b,tmp,x)

c     from numerical recipes Press, Flannery, Teukolsky and Vetterling
c     Cambridge Univ Press 1986, pg 233
c     mod B.D. 

      dimension u(mp,np),w(np),v(np,np),b(mp),x(np)
      dimension tmp(np)  ! work array
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      end

      subroutine svdvar(v,ma,np,w,cvm,ncvm,wti)

c     from numerical recipes Press, Flannery, Teukolsky and Vetterling
c     Cambridge Univ Press 1986, pg 233
c     mod B.D.

      dimension v(np,np),w(np),cvm(ncvm,ncvm)
      dimension wti(ma) ! work array
      do 11 i=1,ma
        wti(i)=0.
        if(w(i).ne.0.) wti(i)=1./(w(i)*w(i))
11    continue
      do 14 i=1,ma
        do 13 j=1,i
          sum=0.
          do 12 k=1,ma
            sum=sum+v(i,k)*v(j,k)*wti(k)
12        continue
          cvm(i,j)=sum
          cvm(j,i)=sum
13      continue
14    continue
      return
      end

      subroutine svdfitb( udes,brh,ndat,npar0,nwh,np,sngv,awr,vrk,cvm,
     I  wdes,rv1,zfun,yfun,qgam,conv,iter)

! do SVD fit with checking output
! Author: Bernard Delley 2024

      use fit_data , only : npar1,npar2,npar3,crh,drh,xk1, erh,
     I         ke1,nqm,nqm1,nqm0,uqspl,uimag,uxno,uanum
      use esf_data , only : mprbay,ncd,fil,yr, phi, pi2deg,fil
!     parameter( tol=1.e-5 )
!     parameter( tol=6.e-5 )
!     parameter( tol=1.e-4 )
!     parameter( tol=8.e-4 )
      parameter( tol=3.e-4 )
!     parameter( tol=1.e-3 )
!     parameter( tol=1.2e-3 ) !too high: 2.e-3:3 ! better 1.e-3:0 no jmp
      dimension udes(nwh,np),wdes(nwh,np),vrk(np,np),cvm(np,np)
      dimension brh(nwh),sngv(np),awr(np),rv1(np),av(5),jv(5),sv(5)
      double precision chisq,csum
        npar=npar0
        do i=1,npar
          do j=1,ndat
            wdes(j,i) = udes(j,i)
          enddo
        enddo
        call svdcmp( udes, ndat, npar, nwh, np, sngv, vrk, rv1)
        do j=1,npar
            amx=0
            do i=1,npar
              if(abs(vrk(i,j)).gt.amx) then
                 imx=i
                 amx=abs( vrk(i,j) )
              endif
            enddo
!         write(6,'(a,2i4,0p,2f10.5,1p,e10.2)')'ck_V',imx,j,amx,
!    I  vrk(3,j),sngv(j) 
        enddo
!      endif
        wmax=1.e-30
!       do i=1,5
!         av(i) = 0
!       enddo
        ami=1.e6
        do j=1,npar
          if(sngv(j).gt.wmax) wmax=sngv(j)
          if(sngv(j).lt.ami) then
            ami = sngv(j)
            jmi= j
          endif
        enddo
        thresh=tol*wmax
        ami=1.e6
        do j=1,npar
           if(sngv(j).lt.ami .and. sngv(j).ge.thresh) then
             ami=sngv(j)
             jmi=j
           endif
        enddo
        k=0
        do j=1,npar
            amx=0
            do i=1,npar
              if(abs(vrk(i,j)).gt.amx) then
                 imx=i
                 amx=abs( vrk(i,j) )
              endif
            enddo
         if(sngv(j).lt.thresh .or. j.eq.jmi) then
           if(mprbay.gt.5) then
            write(6,'(a,2i4,0p,f10.5,1p,2e10.2)')'ck_SNG'
     I     ,j,imx,amx,sngv(j),sngv(j)/wmax
            write(6,'(20f10.3)')(vrk(i,j),i=1,npar1)
            if(nqm.gt.1) then
              write(6,'(20f10.3)') (vrk(i,j),i=npar1+1,npar1+nqm0)
              if(uimag)
     I        write(6,'(20f10.3)') (vrk(i,j),i=npar1+nqm0+1,npar)
            else
              write(6,'(20f10.3)') (vrk(i,j),i=npar1+1,npar1+ncd)
              if(npar.gt.npar1+ncd)
     I        write(6,'(20f10.3)') (vrk(i,j),i=npar1+ncd+1,npar)
            endif
           endif
         else
           if(mprbay.gt.5) then
            write(6,'(a,2i4,0p,f10.5,1p,2e10.2)')'ck_SNG'
     I     ,j,imx,amx,sngv(j),sngv(j)/wmax  
           endif
         endif
         if(sngv(j).lt.thresh) then
           sngv(j)=0
           k=k+1
         endif
        enddo
       if(mprbay.gt.4) then
       write(6,'(a,3i4,1p,4e10.2,0p,i4,1p,e10.2/(20e10.2))')'sing',
     I  npar,k,jmi,wmax,ami,ami/wmax,tol
       endif

        call svbksb( udes, sngv, vrk, ndat, npar, nwh, np, brh, rv1,awr)
        chisq=0
        call indexx( ndat, xk1(1,1), ke1)
        do i=1,ndat
          io=ke1(i)
          sum = 0
          do j=1,npar
            sum = sum + awr(j)*wdes(io,j)     ! contains rmsi^1 factor
          enddo
          udes(io,1) = sum                  ! overwritten for fit function output
          cterm = (brh(io)-sum)**2 
          chisq = chisq + cterm               ! brh also contains rmsi^1 factor
          crh(i) = chisq
        enddo
!       zfun = (chisq-ndat)*0.5/ndat
        zfun = (chisq-ndat+npar)*0.5/(ndat-npar)

        afac = 10000. !  ! normalize zfun=1 -> crh(io) to saturate at 10000 for plot g5.ps
        if(fil(2).gt.20000) afac=20000.
        if(fil(2).gt.40000) afac=40000.
        dfac = 0.5/(ndat-npar)

!       call indexx( ndat, xk1(1,1), ke1)
        csum = +npar      
        crh00 = csum*dfac*afac
        crhn = crh(ndat)
        crhf = float(ndat-npar)/ndat
        do i=1,ndat
          io = ke1(i)
          erh(i) = xk1(io,5)*100  +0.5*afac
          drh(i) = xk1(io ,1)
          crh(i) = ((crh(i)-i*crhf)*dfac   +0.5)*afac
        enddo 
        an2 = (ndat-npar)*0.5
        ch2 = chisq*0.5
        conv=0
        conw=0
        do i=1,npar
!        if(i.eq.3 .and. .not.uphi) then
!        else
          conv = conv + awr(i)**2
          if(abs(awr(i)).gt.abs(conw)) then
            iconv=i
            conw = awr(i)
          endif
!        endif
        enddo
        conv=sqrt(conv)
       yfun = yfun + npar*0.5    ! 251030bd    now yfun, yfun/(ndat-npar) should only differ by rounding from zfun on conv

      if(mprbay.gt.0) then
        write(6,'(a,i3,i5,i7,i7,2f12.3,2f12.6,i5,f9.3,f9.1,f9.0)')'conv'
     I ,iter,npar,ndat,k,zfun,yfun/(ndat-npar),conv,conw,iconv 
     I ,phi*pi2deg,fil(1),fil(2)
      endif
        nsing = k
       if(mprbay.gt.2) then
        write(6,'(a,2L2,9i5)')'ck_uxno',uxno,uanum,npar1,nes
        if(uxno.and.uanum) then
          write(6,'(a,i4,11x,(20f10.5))')'coefn',iter,(awr(i),i=1,npar1)
        elseif(uxno) then
          write(6,'(a,i4,11x,(4f10.5,10x,20f10.5))')'coefn',iter,
     I    (awr(i),i=1,npar1)
        else
          write(6,'(a,i4,11x,(3f10.5,30x,20f10.5))')'coefn',iter,
     I    (awr(i),i=1,npar1)
        endif

         if(nqm.gt.1) then
           write(6,'(20f10.5)') (awr(i),i=npar1+1,npar1+nqm0)
         if(uimag)
     I     write(6,'(20f10.5)') (awr(i),i=npar1+nqm0+1,npar)
         else
           write(6,'(20f10.5)') (awr(i),i=npar1+1,npar1+ncd)
          if(npar.gt.npar1+ncd)
     I     write(6,'(20f10.5)') (awr(i),i=npar1+ncd+1,npar)
         endif
       endif

       call svdvar(vrk,npar,np,sngv,cvm,np,rv1)
!      write(6,'(a,9f10.5)')'ck_RMS',(sqrt(cvm(i,i)),i=1,npar1)
        if(npar.gt.npar1) then
         do i=npar1+1,npar ! ncd
           rv1(i-npar1)= sqrt(cvm(i,i))
         enddo
        endif
       if(mprbay.gt.5) then
        write(6,'(a,1p,20e10.2)')'variances ',(cvm(i,i),i=1,npar2)
        write(6,'(1p,20e10.2)')(cvm(i,i),i=npar2+1,npar)
       endif
!     write(6,'(a,9i5)')'ck_nqm0',npar1,npar2,npar,nqm0
       if(mprbay.gt.4) then
!       if(nqm.gt.1 .and. uimag) then  ! check crosscorrelations
!         write(6,'(a)')'Variance_r Covariance_ri Variance_i, *10^10'
!         fx=1.e10
!         write(6,'(20f10.0)')(fx*cvm(i,i),i=npar1+1,npar1+nqm0)
!         write(6,'(20f10.0)')(fx*cvm(i+nqm0,i),i=npar1+1,npar1+nqm0)
!         write(6,'(20f10.0)')(fx*cvm(i,i),i=npar1+1+nqm0,npar1+2*nqm0)
!       endif


       endif

       return
       end

      subroutine spline(x,y,n,yp1,ypn,y2,u)

c     from numerical recipes Press, Flannery, Teukolsky and Vetterling
c     Cambridge Univ Press 1986, pg 233
c     mod B.D.
 !    workspace u

      dimension x(n),y(n),y2(n),u(n)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end

      subroutine splint(xa,ya,y2a,n,x,y)

c     from numerical recipes Press, Flannery, Teukolsky and Vetterling
c     Cambridge Univ Press 1986, pg 233
c     mod B.D.

      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
!     if (h.eq.0.) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

      subroutine setqsplin
!     use mod_qsplines
! Author: Bernard Delley 2025
      use fit_data
      use esf_data , only : fcut,pix,bpix,np,ncd,naqs,nes,iroip,roip
      logical debug

!     intention have xq(1)=2 and xq(16) ~ 64 : find shifted exponential grid
! use program trymap to find starting aqs for a changed intention
!     xq1 =    1.     ! 1   ! 2  ! 4
      xq1 =    1.*64./np
!     n =      24 ! 0.25  [256] 
!     n =      20 ! 0.5   [128,256]
!     n =      16 ! 1      [64,128]
!     n =      12  ! 1     [32,64]
      n = naqs
      debug=.false.
!     debug=.true.
   1  continue
      if(debug) write(6,'(5x,9a)')
     I 'i    daqs         aqs            bqs            fra',
     I '           dfrada         exp(aqs)'
      xqn = 64 ! at 1 c/p ~ fcut(f/8)
      rxq = xqn/xq1
      aqs = 42. /(n*n*xq1) 
      do i=1,9
        fra = exp(aqs*n) -rxq*exp(aqs) +rxq -1
        dfrada = n*exp(aqs*n) -rxq*exp(aqs)
        daqs = -fra/dfrada
        bqs = xq1/(exp(aqs)-1)
      if(debug) then
      write(6,'(i6,1p,e10.2,0p,19f15.7)')i,daqs,aqs,bqs,
     I fra,dfrada,exp(aqs*n),exp(aqs) !,rxq
      endif
        aqs = max( aqs +daqs , aqs*0.8)
      enddo
      if(debug) stop 'debug aqs'
      if(abs(daqs).gt.1.e-6) then
        debug=.true.
        goto 1
      endif

!  find dimension nqm of spline parametrization for current fcut
      i=1
      xqi=0
       xqn = 64*fcut*pix
!      xqn = 2*np*fcut*pix
      do while( xqi.lt. xqn )
        i=i+1
        xqi= bqs*(exp(aqs*i)-1)
      enddo
      nqm = i+1

!  adjust grid parameters to hit fcut precisely with end point
      n = i
      rxq = xqn/xq1
      do i=1,6
        fra = exp(aqs*n) -rxq*exp(aqs) +rxq -1
        dfrada = n*exp(aqs*n) -rxq*exp(aqs)
        daqs = -fra/dfrada
        bqs = xq1/(exp(aqs)-1)
!     write(6,'(2i3,1p,e10.2,0p19f15.5)')n,i,daqs,aqs,bqs,
!    I fra,dfrada,exp(aqs*n),exp(aqs),xqn
      aqs = aqs +daqs
      enddo

! adapt bqs if np different from 64
      if(np.ne.64) then
      bqs = np*bqs/64.
      endif     
    
      nqm1 = nqm+1
      nqm2 = nqm1 + 1
      allocate( xq(nqm), yqr( nqm2), yqi( nqm),
     I  bqfn(ncd,nqm2), stat=istat)
       if(istat.ne.0) write(6,*)'Error alloc qspline ',istat,nqm

      do i=1,nqm
        xq(i)= bqs*(exp(aqs*(i-1))-1)
        yqr(i)=0
        if(mprbay.gt.4) 
     I  write(6,'(a,i5,9f15.5)')'ck_xq',i,xq(i),xq(i)/np
      enddo
!     write(6,'(/a,i2,2i3,2i4,a,f6.1)')'#_OTF_model_Y [1-',
      write(6,'(/a,i2,2i3,2i4,a,f6.1)')'OTF_model_Y [1-',
     I naqs,nqm-2,nes,np,2*iroip,']'

      return
      end

      subroutine udesyryi( iedi,iedj, ij, u0, edip)
! Author: Bernard Delley
      use fit_data
      use esf_data

!     u1 = fil(2)*((1-edip)*fsf(jp,iedi,1)                         ! desf/dyr(jp)
!    I              + edip *fsf(jp,iedj,1))*rmsi(ij)
!     u0 = fil(2)*rmsi(ij)
!     dhypo(ij)/dyr(j) = u1

      imdo = 1
      if(uimag) imdo = 2


       jp = npar1
       do imd=1,imdo 
      do j=1,ncd
        u1 = u0*((1-edip)*fsf(j,iedi,imd) + edip*fsf(j,iedj,imd))   ! dhypo/dyr /dyi
        udes(ij,jp+j) = u1
       write(6,'(a,5i4,f12.2,9f15.9)')'ck_u',imd,ij,j,iedi,iedj,
     I u0, ((1-edip)*fsf(j,iedi,1) + edip*fsf(j,iedj,1))
     I ,1-edip,fsf(j,iedi,1),edip,fsf(j,iedj,1)
      enddo
        jp = npar1 +ncd
       enddo

      return
      end

      subroutine udesqsplin( iedi,iedj, ij, u0, edip)
! Author: Bernard Delley
      use fit_data
      use esf_data

!     u1 = fil(2)*((1-edip)*fsf(jp,iedi,1)                         ! desf/dyr(jp)
!    I              + edip *fsf(jp,iedj,1))*rmsi(ij)
!     u0 = fil(2)*rmsi(ij)
!     dhypo(iedi,iedj,xno)/dyr(j) = u1
!     dyr(j)/dqr(im...iq) = bqfn(j,im...iq)  else 0

      imdo = 1
      if(uimag) imdo = 2

       npa11 = npar1   
       do imd=1,imdo 
        do j=1,ncd

         u1 = u0*((1-edip)*fsf(j,iedi,imd) + edip*fsf(j,iedj,imd))   ! dhypo/dyr(j) /dyi
         do i=1,nqm0
          udes(ij,npa11+i) = udes(ij,npa11+i) + u1*bqfn(j,i+1)       ! * dyr(j)/dyqr(i)
         enddo

        enddo
        npa11 = npa11 +nqm0
       enddo
!         write(6,'(a,i5,99f9.4)')'ck_u',ij,(udes(ij,i),i=1,npar)
      return
      end

      subroutine update3
! update yqr and yr edge model, and imag parts
! onlyi not supported, obsolete now
! Author: Bernard Delley
      use esf_data
      use fit_data

        dye=1.e30
      if(unats) then
        do j=1,nqm0 
          yqr(j+1) = yqr(j+1) + awr(j+npar1)
        enddo
        call spline(xq,yqr,nqm,dye,dye,vrk,cvm)
      else   ! unats -------------------------------------------------
        dyr1 = dyr1 + awr(npar1+1)
        dyrn = dyrn + awr(npar1+nqm)
        do j=2,nqm-1 
          yqr(j) = yqr(j) + awr(j+npar1)
        enddo
!     write(6,'(a,(t5,20f10.5))')'ck_C',(yqr(i),i=1,nqm)
        call spline(xq,yqr,nqm,dyr1,dyrn,vrk,cvm)
      endif  ! unats -------------------------------------------------
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yr(j))
        enddo

      if(uimag) then
       if(unats) then
        do j=1,nqm0 
          yqi(j+1) = yqi(j+1) + awr(j+npar1+nqm0)
        enddo
        yqi(1) = 0
        yqi(nqm) = 0
        call spline(xq,yqi,nqm,dye,dye,vrk,cvm)
       else   ! unats ------------------------------------------------
        dyi1 = dyi1 + awr(npar1+nqm+1)
        dyin = dyin + awr(npar)
        yqi(1) = 0
        yqi(nqm) = 0
        j1=1
        do j=npar1+nqm+2,npar-1
         yqi(j1) = yqi(j1) + awr(j+npar1+nqm)
         j1=j1+1
        enddo
        call spline(xq,yqi,nqm,dyi1,dyin,vrk,cvm)
       endif  ! unats ------------------------------------------------
        do j=1,ncd
         xj=j
         call splint(xq,yqi,vrk,nqm,xj,yi(j))
        enddo
      endif  ! uimag -------------------------------------------------

      return
      end

      subroutine mkesfn
! spline function for edge curvature, endpoints fixed at 0, 
! fit parameters: nes-2  inner points 
! Author: Bernard Delley 2025
      use esf_data
      use fit_data

      dye = 1.e30
      xpa = -roip
      dxpa = 2*roip/(nes-1)
      do i=1,nes
        xes(i) = xpa
        yes(i) = 0
        xpa = xpa + dxpa
      enddo
      do i=2,nes-1        ! prepare y2es for each edge spline basis fn
         yes(i) = 1
         call spline(xes,yes,nes,dye,dye,y2es(1,i),cvm)     !  cvm <- u
         yes(i) = 0
!        write(6,'(a,i3,99f10.5)')'ckes',i,(y2es(j,i),j=1,nes)
      enddo
      return
      end

      subroutine mkqfn
! prepare spline envelope functions and derivative bqfn
! Author: Bernard Delley 2025
      use esf_data
      use fit_data

      dye = 1.e30

!     write(6,'(a,9i5)')'ck_mkqfn',nqm,nqm0,nqm1,nqm2
      call flush(6)
      do i=1,nqm2
        do j=1,ncd
          udes(j,i) = 0
          bqfn(j,i) = 0
        enddo
        yqr(i) = 0
      enddo

      nqm2 = nqm+2
         call spline(xq,yqr,nqm,1.,dye,vrk,cvm)     ! vrk <- y2     cvm <- u
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yj)
         bqfn(j,1) = yj    ! basisfunctions for parameters 1...nqm at j=1...ncd
        enddo
      do i=2,nqm-1
         yqr(i) = 1
         call spline(xq,yqr,nqm,dye,dye,vrk,cvm)
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yj)
         bqfn(j,i) = yj
        enddo
         yqr(i) = 0
      enddo
         call spline(xq,yqr,nqm,dye,1.,vrk,cvm)
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yj)
         bqfn(j,nqm) = yj
         brh(j) = yr(j)
        enddo
        i=1
         yqr(i) = 1
         call spline(xq,yqr,nqm,dye,dye,vrk,cvm)
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yj)
         bqfn(j,nqm1) = yj
        enddo
         yqr(i) = 0
        i=nqm
         yqr(i) = 1
         call spline(xq,yqr,nqm,dye,dye,vrk,cvm)
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yj)
         bqfn(j,nqm2) = yj
        enddo
         yqr(i) = 0

      

      conv = 1
      iter = -9
      yio=1.
      io=0
      do i=1,nqm
!       yqr(i) = 1 - xq(i)/xq(nqm)
        do while(io+0.9.lt.xq(i))
          io=io+1
          yio = yr(io)
        enddo
        yqr(i) = yio
        yqi(i) = 0 
!       write(6,'(a,2i5,19f10.5)')'ck_mkqfn',i,io,xq(i),yqr(i)
      enddo
      dyr1 = yqr(1)
      dyrn = yqr(nqm)
      yqr(1) = 1      !  will remain fixed
      yqr(nqm) = 0    !  will remain fixed
!     dyi1 = 0
!     dyin = 0
!       call spline(xq,yqr,nqm,dyr1,dyrn,vrk,cvm)
        call spline(xq,yqr,nqm,dye,dye,vrk,cvm)
        do j=1,ncd
         xj=j
         call splint(xq,yqr,vrk,nqm,xj,yr(j))
        enddo
      do i=1,nparm ! nqm2
        awr(i) = 0
      enddo

      return
      end

      subroutine mkwipe
! Author: Bernard Delley 2026
      use esf_data
      use fit_data

      do i=1,ncd
        yr(i) = 1 - i/float(ncd)
      enddo
      return
      end
