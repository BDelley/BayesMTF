! Postscript graphics by Fortran calls
!
! Copyright 1987-2026 Bernard Delley
!  Licensed under the Apache License, Version 2.0
!
      module bgraph
        integer iout, ixlow, iylow, ixhigh, iyhigh
        real scalex,scaley,xmin,xmax,ymin,ymax,deltax,deltay
        real xmaxs,xmins,ymaxs,ymins
      end module bgraph

      subroutine frameb(xmn,xmx,width, nxtick,xtick,xtval,xform,
     I                  ymn,ymx,height,nytick,ytick,ytval,yform,
     |                  framewidth,linewidth,tith,titv,gfile)

c      Tektronix version: B. Delley jan 87
c      Postscript version: A. Kessi apr 94   bd nov/dec 94, apr98_P
c
c    revisions:
c  001101 bd
c  260617 bd, .eps , gfile
c...
      use bgraph
       character*(*) xform,yform,tith,titv,gfile
      character*11 string  ! max format i9 or f9.x + ()
        character*13 fmt
      integer nxtick,nytick,i,inttemp, itlen
      real cmx,cmy      ! pt/mm
!     real deltax,deltay      ! lower left edge of frame 
c                               in page coordinates
!     real scaley,scalex      ! data -> page coordinates
!     real xmin,xmax,ymin,ymax      ! in page coordinates
      real xmn,xmx,ymn,ymx      ! in data coordinates
      real framewidth,linewidth      ! in points
      real width,height      ! in mm
      real xtick(nxtick),ytick(nytick) ! list of ticks (position)
      real xtval(nxtick),ytval(nytick) ! list of tick label values
      real curx,cury
      real dx,dy,dxy,dyy,wchar      ! adjust labels
      integer lxint,lyint      ! number of letters in labels
      real ticklenx,tickleny,ticklx,tickly
!     common/bbox/ ixlow,iylow
      data cmx,cmy / 2.8346457, 2.8346457 /
      data dxy,dyy,wchar / -24., -0., 10.0 /
      data ticklx,tickly / 5., 5. /

      if(height.lt.1. .or. height.gt.330) goto 90
      if(width.lt.1. .or. width.gt.330) goto 90
      deltay = 750. - cmy*height
      deltay = 60.
  5      format(1x,i1)
      read(xform,5,err=91,end=91) lxint
      dx = - wchar*lxint*0.5
      read(yform,5,err=92,end=92) lyint
      dy = - wchar*lyint - wchar*0.5   !2

!     deltax = 70. + wchar*lyint + wchar
      deltax = 30. + wchar*lyint + wchar
      xmin=xmn      ! store parameters in common block
      ymin=ymn
      xmax=xmx
      ymax=ymx
      scalex=cmx*width/(xmax-xmin)
      scaley=cmy*height/(ymax-ymin)
      ticklenx=ticklx/scalex
      tickleny=tickly/scaley

        ixlow = 11. ! 44. ! deltax + dy-20.
        iylow = 11. ! deltay + dxy-30.
        call PSprolog(18,framewidth,linewidth,gfile)  ! font scale 18

        itlen = len(titv)
        k=0
        do j=1,itlen-4
          if(titv(j:j+1).eq.")s") k=k+20
        enddo
        curx = 25.     ! 260617bd
        cury = deltay + (ymax-ymin)*0.5*scaley - wchar*itlen*0.5 +k
        write(iout,'(2f8.2,4a)')curx,cury,' m  gs 90 rotate (',titv,
     I  ') show gr'

        itlen = len(tith)
        curx = deltax + (xmax-xmin)*0.5*scalex - wchar*itlen*0.5
        if(xform.eq.' ') then
          cury = deltay + dxy
          write(iout,'(2f8.2,4a)')curx,cury,' m (',tith,') show'
        else
          cury = deltay - 44. 
          write(iout,'(2f8.2,4a)')curx,cury,' m (',tith,') show'
        endif

      do i=1,nytick      ! y ticks
        cury=ytick(i)
        if((cury .ge. ymin) .and. (cury .le. ymax)) then
          curx=xmax-ticklenx
          call PSmoveto(xmax,cury)
          call PSlineto(curx,cury)
          curx=xmin+ticklenx
          call PSmoveto(xmin,cury)
          call PSlineto(curx,cury)

          if( yform .eq. ' ') then

          elseif (yform(1:1).eq.'f'.or.yform(1:1).eq.'F' ) then
              write(fmt,'(a,a,a)')'(a,',yform,',a)'
              write(string,fmt,err=12) '(',ytval(i),')'
              call PSshow(xmin,cury,dy,dyy,string)   

           elseif( yform(1:1).eq.'e'.or.yform(1:1).eq.'E' ) then
              write(fmt,'(a,a,a)')'(a,1p,',yform,',a)'
              write(string,fmt,err=12) '(',ytval(i),')'
              call PSshow(xmin,cury,dy,dyy,string)  

           elseif (yform(1:1).eq.'i'.or.yform(1:1).eq.'I') then
              inttemp=nint(ytval(i))
              write(fmt,'(a,a,a)')'(a,',yform,',a)'
              write(string,fmt,err=12) '(',inttemp,')'
              call PSshow(xmin,cury,dy,dyy,string)     
          end if
12          continue
        end if
      enddo


      do i=1,nxtick            ! x ticks
        curx=xtick(i)
        if ((curx .ge. xmin) .and. (curx .le. xmax)) then
          cury=ymax-tickleny
          call PSmoveto(curx,ymax)
          call PSlineto(curx,cury)
          cury=ymin+tickleny
          call PSmoveto(curx,ymin)
          call PSlineto(curx,cury)

          if( xform .eq. ' ') then

          elseif ((xform(1:1) .eq. 'f') 
     I       .or. (xform(1:1) .eq. 'F')) then
             write(fmt,'(a,a,a)')'(a,',xform,',a)'
             write(string,fmt,err=22) '(',xtval(i),')'
             call PSshow(curx,ymin,dx,dxy,string)

          elseif ((xform(1:1) .eq. 'i') 
     I       .or. (xform(1:1) .eq. 'I')) then
             inttemp=nint(xtval(i))
             write(fmt,'(a,a,a)')'(a,',xform,',a)'
             write(string,fmt,err=22) '(',inttemp,')'
             call PSshow(curx,ymin,dx,dxy,string)

           elseif( xform(1:1).eq.'e'.or.xform(1:1).eq.'E' ) then
             write(fmt,'(a,a,a)')'(a,1p,',xform,',a)'
             write(string,fmt,err=22) '(',xtval(i),')'
             call PSshow(curx,ymin,dx,dxy,string)

          end if
22        continue
        end if
      enddo

      write(iout,'(a)')'st'

c      draw frame - has to be after ticks and labels (clipping)
c      no PSgsave/PSgrestore, since clipping has to be enabled 
c       after this
      call PSmoveto(xmin,ymin)
      call PSlineto(xmax,ymin)
      call PSlineto(xmax,ymax)
      call PSlineto(xmin,ymax)
      write(iout,'(a)')'closepath gs st gr clip slw'
      write(iout,'(a)')'%closepath st slw'

      return

  90  continue
      write(6,*)'Error with frameb width height',width,height
      write(6,*)'argument values',
     I                  xmn,xmx,width, nxtick,xtick,xtval,xform,
     I                  ymn,ymx,height,nytick,ytick,ytval,yform,
     I                  framewidth,linewidth,tith,titv
      stop 'error'
  91  continue
      write(6,*)'Error with frameb-xform:',xform
      stop 'error'
  92  continue
      write(6,*)'Error with frameb-yform:',yform
      stop 'error'
      end

	subroutine PSprolog(isiz,framewidth,linewidth,gfile)
c	ak apr 94/bd nov 94, 2000bd 20171207bd
c       eps gfile 260617 bd
c...
	use bgraph
        character*(*) gfile
        real framewidth,linewidth
	integer isiz               ! unit # , script siz

c...
	iout=98
        if(gfile.ne.'') then
	  open( iout, file=gfile, access='SEQUENTIAL',
     I	     status='UNKNOWN' )
        else
	  open(iout,file='graph.ps',access='SEQUENTIAL',
     I	     status='UNKNOWN')
        endif
c	write Postscript prolog to file
          write(iout,'(a)') '%!PS-Adobe-2.0 EPSF-2.0'
          ixhigh = deltax + scalex*(xmax-xmin) + 2
          iyhigh = deltay + scaley*(ymax-ymin) + 2
          write(iout,'(a,4i5,a)') '%%BoundingBox:',
     I       ixlow, iylow, ixhigh, iyhigh
          write(iout,'(a,4f7.1)') '%%HiResBoundingBox:',
     I       float(ixlow), float(iylow), float(ixhigh), float(iyhigh)
	write(iout,'(a)')'%%EndComments','%%BeginProlog','save'
     I  ,'countdictstack','mark','newpath','/showpage {} def'
     I  ,'/setpagedevice {pop} def','%%EndProlog','%%Page 1 1'
cwrite(iout,'(a)') '/mm'
cwrite(iout,'(a)') '   { 2.8346457 mul }'
cwrite(iout,'(a)') '   def'
cwrite(iout,'(a)') '/inch'
cwrite(iout,'(a)') '   { 72 mul }'
cwrite(iout,'(a)') '   def'
cwrite(iout,'(a)') '   def'
	write(iout,'(a)') '/c { 0 360 arc } def % x y r'
        write(iout,'(a)') '/f {fill} def'
	write(iout,'(a)') '/l {lineto} def'
	write(iout,'(a)') '/m {moveto} def'
        write(iout,'(a)') '/n {newpath} def'
        write(iout,'(a)') '/r {rmoveto} def'
        write(iout,'(a)') '/rl {rlineto} def'
!       write(iout,'(a)') '/s {show} def'
        write(iout,'(a)') '/t {translate} def'
        write(iout,'(a)') '/sc {setrgbcolor} def'
        write(iout,'(a)') '/st {stroke} def'
        write(iout,'(a)') '/sf {sc f} def'
        write(iout,'(a)') '/ls {l stroke} def'
        write(iout,'(a)') '/lf {l f} def'
        write(iout,'(a)') '/gs {gsave} def'
        write(iout,'(a)') '/gr {grestore} def'
        write(iout,'(a,a)') '/c0  { 1.20 c .0 .0 .0 sf} def % x y '
        write(iout,'(a,a)') '/cr  { c 0.8 0.0 0.0 sf} def % x y r'
        write(iout,'(a,a)') '/cb  { c 0.0 0.2 1.0 sf} def % x y r'
        write(iout,'(a)') '/lw {setlinewidth} def' 
        write(iout,'(a,f5.2,a)') '/flw {',framewidth,' lw} def' 
        write(iout,'(a,f5.2,a)') '/slw {',linewidth,' lw} def' 
c	write(iout,'(a,i5,a)') '%%/LabelFont /Times-Roman findfont'
        write(iout,'(a,i5,a)') '/LbFnt {/Helvetica findfont'
     I  ,isiz,' scalefont setfont} def'
        write(iout,'(a,i5,a)') '/SFnt {/Helvetica findfont'
     I  ,12,' scalefont setfont} def'
        write(iout,'(a,i5,a)') '/TFnt {/Helvetica findfont'
     I  ,6,' scalefont setfont} def'
        write(iout,'(a,i5,a)') '/BFnt {/Helvetica findfont'
     I  ,24,' scalefont setfont} def'
        write(iout,'(a,i5,a)') '/GFnt {/Symbol findfont'
     I  ,isiz,' scalefont setfont} def'
        write(iout,'(2a)') '/sr { show 1 0 0 sc } def'  
     I  ,'% modify next part of title like )sr('
        write(iout,'(a)') '/sg { show 0 0.8 0 sc } def'
        write(iout,'(a)') '/sb { show 0 0 1 sc } def'
        write(iout,'(a)') '/sS { show /GFnt } def'
        write(iout,'(a)') '/sL { show /LbFnt } def'
	write(iout,'(/a)') 'LbFnt flw'
        write(iout,'(a)')' gs'

      xmaxs = deltax+scalex*(xmax-xmin)
      ymaxs = deltay+scaley*(ymax-ymin)
      xmins = deltax
      ymins = deltay

	return
	end

        subroutine endfrm(nxtick,xtick,nytick,ytick)

! rev 260617 bd eps
        use bgraph
        implicit real(a-h,o-z)
        real xtick(*),ytick(*)
      data ticklx,tickly / 5., 5. /

      ticklenx=ticklx/scalex
      tickleny=tickly/scaley

      write(iout,'(a)')' gr flw 0 0 0 sc n'
      do i=1,nytick      ! y ticks
        cury=ytick(i)
        if((cury .ge. ymin) .and. (cury .le. ymax)) then
          curx=xmax-ticklenx
          call PSmoveto(xmax,cury)
          call PSlineto(curx,cury)
          curx=xmin+ticklenx
          call PSmoveto(xmin,cury)
          call PSlineto(curx,cury)
        end if
      enddo
      do i=1,nxtick            ! x ticks
        curx=xtick(i)
        if ((curx .ge. xmin) .and. (curx .le. xmax)) then
          cury=ymax-tickleny
          call PSmoveto(curx,ymax)
          call PSlineto(curx,cury)
          cury=ymin+tickleny
          call PSmoveto(curx,ymin)
          call PSlineto(curx,cury)
        end if
      enddo
      write(iout,'(a)')'st'

c redraw frame
      call PSmoveto(xmin,ymin)
      call PSlineto(xmax,ymin)
      call PSlineto(xmax,ymax)
      call PSlineto(xmin,ymax)
      write(iout,'(a)')'closepath st','showpage','%%Trailer'
     I ,'cleartomark','countdictstack','exch sub { end } repeat'
     I ,'restore','%%EOF'

      close(iout)

      return
      end

	subroutine PSmoveto(x,y)
c...
	use bgraph
        integer mx,my
        real x,y,xp,yp, cut
        data mx,my, cut /0, 0, 999.99/
c...
        xp = deltax+scalex*(x-xmin)
        yp = deltay+scaley*(y-ymin)

        if(xp.gt.-200 .and. xp.lt.cut .and. 
     I     yp.gt.-200 .and. yp.lt.cut) then
c its ok
            write(iout,'(2f8.2,a,3f8.2)') xp,yp,' m' !,xo,yo,thr
        else
            write(iout,'(a,2f9.2,a)')'%',xp,yp,' m % out of bounds' !,xo,yo,thr
          mx = mx+1
!         if(mx.le.5) write(16,'(a,2f15.2,1p,2e11.2)')
!    I    ' move outside frame',xp,yp,x,y
          if(x.lt.xmin) xp = deltax
          if(y.lt.ymin) yp = deltay
          if(x.gt.xmax) xp = deltax+scalex*(xmax-xmin)
          if(y.gt.ymax) yp = deltay+scaley*(ymax-ymin)
            write(iout,'(2f8.2,a,3f8.2)') xp,yp,' m % out of bounds'
        endif

	return
	end

	subroutine PSlineto(x,y)
c...
	use bgraph
        integer mx,my
        real x,y,xp,yp,xo,yo,thr, cut
        save xo,yo
        data mx,my, cut, thr /0, 0, 999.99, .5/
        data xo,yo / -3., -3. /
c...
        xp = deltax+scalex*(x-xmin)
        yp = deltay+scaley*(y-ymin)

        if(xp.gt.-200 .and. xp.lt.cut .and. 
     I     yp.gt.-200 .and. yp.lt.cut) then
c its ok
          if(abs(xp-xo).gt.thr .or. abs(yp-yo).gt.thr) then
            write(iout,'(2f8.2,a,3f8.2)') xp,yp,' l' !,xo,yo,thr
            xo = xp
            yo = yp
          endif
        else
          mx = mx+1
!         if(mx.le.5) write(16,'(a,2f15.2,1p,2e11.2)')
!    I    ' line outside frame',xp,yp,x,y
        endif

	return
	end

        subroutine PSshow(x,y,dx,dy,string)

! rev 260617 bd
c
c       Open a new graphics environment and print string at (x,y)
c       to file iout.
c       This allows to print text even while a path is not completed.
c...
        use bgraph
        character*(*) string
        real x,y,dx,dy,xx,yy

        xx = deltax+dx+scalex*(x-xmin)
        yy = deltay+dy+scaley*(y-ymin)
        if( yy.lt.iyhigh .and. yy.ge.iyhigh-13) yy = iyhigh-13     ! shift it in
        if( xx.lt.ixhigh .and. xx-2*dx.gt.ixhigh) then
           xx = ixhigh+2*dx ! shift it in
        endif

        write(iout,'(2f8.2,4a)') xx,yy,' m ',string,' show'
        return
        end

c	collection of graphics routines for Postscript
c---
	subroutine linrgb(n,x,y,cr,cg,cb)

c	B. Delley dec 2000
c       rev dec 2005 bd, 20260617 bd
c...
        use bgraph
	integer n,i,mx
        character*1 ml
	dimension x(n),y(n)
        logical start

        start = .true.

        call psline( n, x, y, start)
        write(iout,'(a,3f5.2,a)')' ',cr,cg,cb,' sc st'
        
	return
	end

        subroutine psline(n,x,y,start)

c 051230 bd revision: elementary line in ps with checks
c 251126 bd revision,simplification with new fnd routine after bay_mtf_v6 -> g6.ps failure
c rev 260617 bd  bgraph module

c input:
c  n       number of points on line
c  x,y     arrays
c  start   logical to start line rather than continue
c output:
c  start   false= drawing actually has started 

      use bgraph
      integer i,i1,i2,i3,n,j
      real x(*), y(*), x1,x2,y1,y2,xp,yp, xc,yc, xi,yi,thrd
     I  , cut,x3,y3,xb,yb,val(4)

      logical start, draw, drawnc, go_out, come_in, drawn,drawp
      character*2 sml
      data thrd / 0.3 /
      save xp,yp, xi,yi
      save x1,y1,x3,y3
      save j
      save go_out

      if(n.ge.1) then
        i1=1
        i2=n
        i3=1
      elseif(-n.gt.1) then
        i1=-n
        i2=1
        i3=-1
      else
        write(6,'(a,i9)')'Error: psline requires at least n=2, n=',n
        stop 'psline'
      endif
      j = -1

      if(start) then
        sml = ' m'
        write(iout,'(a,3i5)')'n %',i1,i2,i3
        x1 = x(i1)
        y1 = y(i1)
        j = 0
        go_out =.false.
      else
        sml = ' l'
        write(iout,'(a,3i5)')'% reverse',i1,i2,i3
        start = .true.
      endif
      draw = .false.
      drawp = .false.

      do i=i1,i2,i3

         drawn = .false.
         drawnc = .false.
         x2 = x(i)
         y2 = y(i)
!        write(iout,'(a,2L2,9(f10.2,f8.2))')'%ck_0',draw,drawp,
!    I   x2,y2,x1,y1,xp,yp,xi,yi
           call fnd_border( xc,yc, x2,y2, x1,y1, xp,yp,
     I     draw,drawp,xi,yi,sml,j)
!        write(iout,'(a,2L2,9(f10.2,f8.2))')'%ck_f',draw,drawp,
!    I   xc,yc,xi,yi,xp,yp
         if(start) then
           xc = deltax+scalex*(x2-xmin)
           yc = deltay+scaley*(y2-ymin)
           if(x2.ge.xmin .and. x2.le.xmax .and.
     I        y2.ge.ymin .and. y2.le.ymax) then
             write(iout,'(2f8.2,2a,3f8.2)') xc,yc,sml
             sml = ' l'
             j = 1
             start=.false.
             drawn=.true.
             drawp=.false.
           endif
             xi = xc
             yi = yc
         elseif(draw) then
           j = j + 1
           if(j.gt.2) then
               call collinear( xc,yc, xi,yi, xp,yp, drawnc, thrd,val)
               if(drawnc) then
                 write(iout,'(2f8.2,2a,i6,9(f10.2,f8.2))') xi,yi,sml
!    I ,' %T',j,xc,yc,xp,yp,val
                 xp = xi
                 yp = yi
                 xi = xc    
                 yi = yc
                 drawp=.true.  ! may need to be drawn
               else
!                write(iout,'(2f8.2,2a,i6,9(f10.2,f8.2))') xi,yi,sml
!    I ,' %F',j,xc,yc,xp,yp,val
               endif
           elseif(j.eq.1) then
           write(iout,'(2f8.2,2a,i3)') xc,yc,sml   ,' %',j
             xp = xi
             yp = yi
             xi = xc
             yi = yc
             drawn=.true.
             drawp=.true.
           elseif(j.eq.2) then
             xp = xi
             yp = yi
             xi = xc
             yi = yc
!     write(iout,'(a,2L2,9(f10.2,f8.2))')'%ck_2',draw,drawn,xi,yi,xp,yp
             drawp=.true.
           endif
         endif
        x1 = x2
        y1 = y2
      enddo
!     write(iout,'(a,3L2,2f8.2,2a)')'%ck_?',draw,drawnc,drawn,xi,yi,sml 
      if(.not.drawn) then
        if(draw) then
          write(iout,'(2f8.2,2a)') xc,yc,sml ,' %c'
        elseif(j.gt.1) then
          write(iout,'(2f8.2,2a)') xi,yi,sml ,' %i'
        endif
      endif
      return
      end

      subroutine collinear( xc,yc, xi,yi, xp,yp, drawnc, thrd,val)

c output
c  drawnc  true when non-collinear

c  input
c xp,yp, xi,yi

      use bgraph
      real xa,ya, xb,yb, xc,yc, xi,yi, xp,yp, zz, z0, thrd, val(4)
     I ,xj,yj
      logical drawnc

      xj=xi
      yj=yi
      xa = xc - xi
      ya = yc - yi
      xb = xp - xi
      yb = yp - yi
      xx = xa*xa + ya*ya
      yy = xb*xb + yb*yb
      xx = sqrt(xx)
      yy = sqrt(yy)
      xy = min( xx, yy)
      z0 = xa*yb - ya*xb
      zz = (xa*yb - ya*xb)/max( xy, 1.e-6)
        drawnc = .true.
      if(abs(zz).lt.thrd) then
        drawnc = .false.
        xi = xc
        yi = yc
      else
        drawnc = .true.
      endif
      val(1)=z0
      val(2)=zz
      val(3)=xx
      val(4)=yy
!     write(iout,'(a,L6,9(f10.2,f8.2))')'%ck_col',drawnc,
!    I xc,yc,xi,yi,xp,yp,xa,ya,xb,yb,z0,zz,xx,yy
!    I xi,yi,xc,yc,xj,yj,xp,yp,val
      return
      end

      subroutine fnd_border( xc,yc, x2,y2, x1,y1, xp,yp,
     I  draw,drawp,xi,yi,sml,j)

! input
!  x2,y2 current orig coord

      use bgraph
      real x2,y2, xb, yb, xc, yc, xp,yp, x1,y1
      real xi,yi,xo,yo,xj,yj
      character*2 sml
      logical draw, drawp, drx, dry, go_out, come_in

      xb = deltax+scalex*(x2-xmin)   ! current  point
      yb = deltay+scaley*(y2-ymin)
      xc = deltax+scalex*(x1-xmin)   ! previous point
      yc = deltay+scaley*(y1-ymin)
      xo = xb
      yo = yb
      xj = xc
      yj = yc
!     write(iout,'(a,2L2,9(f9.2,f8.2))') 
!    I '%xo,yo,xj,yj',draw,drawp,xo,yo,xj,yj

      draw = .false.
      drx = .false.
      dry = .false.
      go_out = .false.
      come_in = .false.
      inxc=0
      inyc=0
      inxb=0
      inxb=0

      if( xc.gt.xmins .and. xc.lt.xmaxs) then    ! previous inside domain
       if( xb.lt.xmins) then ! go out
        yb = yc + (yb-yc)*(xmins-xc)/(xb-xc)
        xb = xmins
        go_out = .true.
       elseif( xb.gt.xmaxs) then  !  then go out
        yb = yc + (yb-yc)*(xmaxs-xc)/(xb-xc)
        xb = xmaxs
        go_out = .true.
       endif
      else
        drawp=.false.  ! if out, do not draw it in any case
      endif

      if( yc.gt.ymins .and. yc.lt.ymaxs) then
       if( yb.lt.ymins) then  ! go out
        xb = xc + (xb-xc)*(ymins-yc)/(yb-yc)
        yb = ymins
        go_out = .true.
       elseif( yb.gt.ymaxs) then
        xb = xc + (xb-xc)*(ymaxs-yc)/(yb-yc)
        yb = ymaxs
        go_out = .true.
       endif
      else
        drawp=.false.  ! if out, do not draw it in any case
      endif

      if( xb.ge.xmins .and. xb.le.xmaxs) then   ! current in domain
       if( xc.lt.xmins) then !               go in
        yc = yc + (yb-yc)*(xmins-xc)/(xb-xc)
        xc = xmins
        come_in = .true.
       elseif( xc.gt.xmaxs) then  !     then go in
        yc = yc + (yb-yc)*(xmaxs-xc)/(xb-xc)
        xc = xmaxs
        come_in = .true.
       endif
       drx=.true.
      endif

      if( yb.ge.ymins .and. yb.le.ymaxs) then
       if( yc.lt.ymins) then  ! go in
        xc = xc + (xb-xc)*(ymins-yc)/(yb-yc)
        yc = ymins
        come_in = .true.
       elseif( yc.gt.ymaxs) then
        xc = xc + (xb-xc)*(ymaxs-yc)/(yb-yc)
        yc = ymaxs
        come_in = .true.
       endif
       dry=.true.
      endif

      if(come_in) then
      if(drx.and.dry) then
!       sml = ' m'
        write(iout,'(2f8.2,2a,9(f9.2,f8.2))') xc,yc,sml,' %fndi'
     I   ,xj,yj,xo,yo
        sml = ' l'
          xi=xc
          yi=xc
       if(.not.go_out) then
!!      write(iout,'(2f8.2,2a,9(f9.2,f8.2))') xb,yb,sml,' %fnde' ! needed for g0.ps log9945_g ...A52
!    I   ,xj,yj,xo,yo
       endif
        j=1                 ! does not work for case log9945_g g0.ps
      endif
!       come_in = .false.   ! 260107bd commented out
        draw=.true.         ! inserted 260107bd ??
      endif
      if(drx.and.dry) then
      if(go_out) then
        if(drawp) then
          write(iout,'(2f8.2,a,9(f9.2,f8.2))') xj,yj,' l %Fndo-prev'
          drawp=.false.
        endif
        write(iout,'(2f8.2,a,9(f9.2,f8.2))') xb,yb,' l %fndo'
     I   ,xo,yo  ,xj,yj
        j=0
      else
        draw=.true.
      endif
      endif

      if((come_in.and.drx.and.dry).or.go_out) then  
!       already done case
!       write(iout,'(a,4(f9.2,f8.2))') ' %<done>'
!    I   ,xj,yj,xo,yo,xc,yc,xb,yb
      else                !  any line crossing domain ? 260107bd
       drx = .false.
       dry = .false.
       come_in = .false.
       go_out = .false.
       inxc=0
       inyc=0
       inxb=0
       inyb=0
       if(xb.lt.xmins) then
        inxb=-1
       elseif(xb.gt.xmaxs) then
        inxb=+1
       endif
       if(yb.lt.ymins) then
        inyb=-1
       elseif(yb.gt.ymaxs) then
        inyb=+1
       endif
       if(xc.lt.xmins) then
        inxc=-1
       elseif(xc.gt.xmaxs) then
        inxc=+1
       endif
       if(yc.lt.ymins) then
        inyc=-1
       elseif(yc.gt.ymaxs) then
        inyc=+1
       endif
       if(abs(inxb+inxc).gt.1 .or. abs(inyb+inyc).gt.1) then
!        line sure outside
!       write(iout,'(a,2(f9.2,f8.2),9i3)') ' %<=>'
!    I   ,xj,yj,xo,yo,inxc,inyc,inxb,inyb
       elseif(inxb.eq.0.and.inyb.eq.0.and.inxc.eq.0.and.inyc.eq.0) then
!        line sure inside
       else
!       write(iout,'(a,2(f9.2,f8.2),9i3)') ' %<->'
!    I   ,xj,yj,xo,yo,inxc,inyc,inxb,inyb
       if( xj.lt.xmins) then !               go in
        yc = yj + (yo-yj)*(xmins-xj)/(xo-xj)
        xc = xmins
        if(yc.gt.ymins.and.yc.lt.ymaxs) come_in = .true.
       elseif( xj.gt.xmaxs) then  !     then go in
        yc = yj + (yo-yj)*(xmaxs-xj)/(xo-xj)
        xc = xmaxs
        if(yc.gt.ymins.and.yc.lt.ymaxs) come_in = .true.
       endif
!      drx=.true.
       if( yj.lt.ymins .and. .not.come_in) then  ! go in
        xc = xj + (xb-xj)*(ymins-yj)/(yb-yj)
        yc = ymins
        if(xc.gt.xmins.and.xc.lt.xmaxs) come_in = .true.
       elseif( yj.gt.ymaxs .and. .not.come_in) then
        xc = xj + (xb-xj)*(ymaxs-yj)/(yb-yj)
        yc = ymaxs
        if(xc.gt.xmins.and.xc.lt.xmaxs) come_in = .true.
       endif
!      dry=.true.
       if(come_in) then
        write(iout,'(2f8.2,2a,9(f9.2,f8.2))') xc,yc,sml,' %->in'
     I   ,xj,yj,xo,yo
         sml=' l'
       endif
       if( xo.lt.xmins) then ! go out
        yb = yj + (yo-yj)*(xmins-xj)/(xo-xj)
        xb = xmins
        if(yb.gt.ymins.and.yb.lt.ymaxs) go_out = .true.
       elseif( xo.gt.xmaxs) then  !  then go out
        yb = yj + (yo-yj)*(xmaxs-xj)/(xo-xj)
        xb = xmaxs
        if(yb.gt.ymins.and.yb.lt.ymaxs) go_out = .true.
       endif
       if( yo.lt.ymins) then  ! go out
        xb = xj + (xo-xj)*(ymins-yj)/(yo-yj)
        yb = ymins
        if(xb.gt.xmins.and.xb.lt.xmaxs) go_out = .true.
       elseif( yo.gt.ymaxs) then
        xb = xj + (xo-xj)*(ymaxs-yj)/(yo-yj)
        yb = ymaxs
        if(xb.gt.xmins.and.xb.lt.xmaxs) go_out = .true.
       endif
       if(go_out) then
        write(iout,'(2f8.2,2a,9(f9.2,f8.2))') xb,yb,sml,' %->out'
     I   ,xj,yj,xo,yo
       endif
      endif
      endif
        xc=xb      ! assign current point perhaps put on border to xc for output
        yc=yb

!     if(xb.le.xmaxs .and. yb.le.ymaxs .and.
!    I   xb.ge.xmins .and. yb.ge.ymins ) then
!      write(iout,'(a,5L2,9(f10.2,f8.2))')'%ck_fnd',
!    I draw,drx,dry,come_in,go_out,
!    I xc,yc, xo,yo, xj,yj
!     endif

      return
      end

	subroutine symbrgb(n,x,y,diam,cr,cg,cb)

c       colored symbols  20001220 bd
! rev 2606176 bd

        use bgraph
	integer n,i
	real cmx,cmy
	real x,y,diam,r,  cr,cg,cb
	data cmx,cmy/ 2.8346457, 2.8346457 /
	dimension x(*),y(*)
        data cut / 751. /
c
	if (diam .gt. 0.) then		! in mm
	  r=cmy*diam*0.5
	else				! in pts
	  r=-diam*0.5
	end if

        do i=1,n
          xp = deltax+scalex*(x(i)-xmin)
          yp = deltay+scaley*(y(i)-ymin)

          if(xp .gt. xmins .and. xp .lt. xmaxs .and.
     I       yp .gt. ymins .and. yp .lt. ymaxs)      then

              write(iout,'(a,2f8.2,f7.2,a,3f5.2,a)')
     I        'n',xp,yp,r,' c',cr,cg,cb,' sf'
          endif
        enddo
	return
	end

c	collection of graphics routines for Postscript
c---
	subroutine areargb(n,x,y,cr,cg,cb)

c	B. Delley Jan 2001

        use bgraph
	integer isymb,n,i
	dimension x(n),y(n),x2(2),y2(2)
        logical start

        if(n.le.1) then
          write(6,*)'warning Tabulation length !',n
          if(n.le.0)return
        endif
        if(n.gt.19000) then
          write(6,*)'warning Tabulation Length !',n
        endif

        start = .true.

        write(iout,'(a)')'% areargb new'
        call psline( n, x, y, start)
        write(iout,'(3f6.3,a)')cr,cg,cb,' sf'

	return
	end

        subroutine area_bc(n,x,y,z,cr,cg,cb)

c       B.Delley Dec 2005

        use bgraph
        integer isymb,n,i
        common/wsave/xsave1,ysave1,xsave2,ysave2
        dimension x(n),y(n),z(n)
        logical start

        if(n.le.1) then
          write(6,*)'warning Tabulation length !',n
          if(n.le.0)return
        endif
        if(n.gt.19000) then
          write(6,*)'warning tabulation Length !',n
        endif

        start = .true.

        write(iout,'(a)')'% area_bc'
        call psline( n, x, y, start)
        call psline(-n, x, z, start)
        write(iout,'(3f5.2,a)')cr,cg,cb,' sf'

        return
        end

        subroutine autotic(tick,ntick,fmt,amax,amin)

        dimension tick(10)
        character*4 fmt

        fs = 1.d7
c       if(amax.gt.999) then
        do i=-20,20
          ts = 10.0**i
          if(amax-amin.gt.3*ts) then  !  3 =< nt < 6
            fs=ts
          endif
          if(amax-amin.gt.6*ts) then  !  3 =< nt < 7
            fs=2*ts
          endif
          if(amax-amin.gt.15*ts) then !  3 =< nt < 6
            fs=5*ts
          endif
        enddo

        i1 = nint(amin/fs + 0.0)
        i2 = nint(amax/fs - 0.51)
        k = 0
        do i=i1,i2
          k = k + 1
          tick(k) = i * fs
c         write(6,*)i,k,tick(k)
        enddo
        ntick = k

        gs = 9.5
        if(amin.lt.0) then
          ld = 2
          do i=2,8
            if( amax .gt. gs) ld = max( ld, i)
            if(-amin .gt. gs) ld = max( ld, i+1)
c       write(6,'(a,i5,3f15.6)')'ld-',ld,amin,amax,gs
            gs = gs * 10
          enddo
        else
          ld = 1
          do i=2,8
            if( amax .gt. gs) ld = max( ld, i)
            if( amin .gt. gs) ld = max( ld, i)
            gs = gs * 10
          enddo
        endif
          gs = 0.95
          md = 0
          do i=1,6
            if( fs .lt. gs)  md = i
c       write(6,'(a,i5,2f15.6)')'md ',md,fs,gs
            gs= 0.1 * gs
          enddo

        if(md.eq.0 .and. ld.lt.8) then
          write(fmt,'(a,i1)')    'i',ld
        elseif(ld+md+1.le.8) then
          write(fmt,'(a,i1,a,i1)')'f',ld+md+1,'.',md
        else
          write(fmt,'(a)')'e8.1'
        endif
!       write(6,'(9a)')'ck_fmt >>',fmt,'<<'

c       write(6,*)'ck_tick ',fmt,'<<',amax,amin,ntick,fs,i1,i2,ls,lr
        ts = amax+amin
        if(ts.gt.-1.e30 .and. ts.lt.1.e30) then
        else
        write(6,*)'CK_tick ',fmt,'<<',amax,amin,ntick,fs
        endif
        if(amax.le.amin) then
          write(6,'(a,1p,2e12.4)')'Error-autotic amin, amax',amin,amax
          stop
        endif

        return
        end

      subroutine hsl2rgb(r,g,b, h0,s0, gl0)

c Author: B.Delley 2008

c 2008 bd after http://en.wikipedia.org/wiki/HSV_color_space
c 2010 bd modulo+limiter

c input
c  h    hue 0-360. color circle, 0-270 rainbow red-violet, 300=magenta
c  s    saturation 0-1
c  gl   lightness 0-1

c output
c   rgb values

      h = mod( h0, 360.)
      s = max(0.,min(s0,1.))
      gl= max(0.,min(gl0,1.))

      if(gl.lt.0.5) then
        q = gl*(1+s)
      else
        q = gl+s-gl*s
      endif
      p = 2*gl-q
      hn = h/360
      ihn = hn
      hn = hn - ihn
      thrd = 1./3.

      tr = hn + thrd
      if(tr.gt.1) tr = tr-1
      tg = hn
      tb = hn - thrd
      if(tb.lt.0) tb = tb+1

      b1 = 1./6.
      b2 = 0.5
      b3 = 2./3.
      if(tr.lt.b1) then
        r = p+((q-p)*6*tr)
      elseif(tr.lt.b2) then
        r = q
      elseif(tr.lt.b3) then
        r = p+((q-p)*6*(b3-tr))
      else
        r = p
      endif

      if(tg.lt.b1) then
        g = p+((q-p)*6*tg)
      elseif(tg.lt.b2) then
        g = q
      elseif(tg.lt.b3) then
        g = p+((q-p)*6*(b3-tg))
      else
        g = p
      endif

      if(tb.lt.b1) then
        b = p+((q-p)*6*tb)
      elseif(tb.lt.b2) then
        b = q
      elseif(tb.lt.b3) then
        b = p+((q-p)*6*(b3-tb))
      else
        b = p
      endif

c     write(6,'(20f8.2)')h,s,gl,p,q,tr,tg,tb,r,g,b

      return
      end
