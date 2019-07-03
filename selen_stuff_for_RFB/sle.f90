!
! This is program "SLE.F90" 
!
!       GS 04-11-2008 "Intel port..."
! 		Modified by GS september 2008 for horizontal movements 
! 		Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! 		*** Reviewed GS & FC November 2009 - Porting under gfortran 
! 		Also reviewed in April 2010 for the implementation of degree 1 
! 		(remember the mode coupling issue...) 
! 		Revised May 26 GS for new ice SH  
!       Revised Jun 16, 2011 DM for dynamic memory allocation
!       Revised Jun 2011 DM for parallel execution
!
!       Modified RH March 2018 - Adapted for rotational feedback
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! ============================================
! Solves the SLE by the pseudo-spectral method 
! ============================================
!
! Input files:
!       - px-table.dat 
!	- shof.dat 
!	- shice.dat
!	- eb.dat
!	- ebu.dat
!       - ebv.dat  <<< new... 
! 	- sh.bin
!
! Output files:
!	- shu.bin 
!	- shn.bin 
!	- shs.bin 
!       - shz.bin 
!       - shv.bin <<< new... 
!
!
! INCLUDE "harmonics.f90"
 PROGRAM SLE 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER ijunk
 INTEGER :: npix, na
 CHARACTER*22, PARAMETER :: FILES='shs.bin', FILEU='shu.bin', & 
 			    FILEN='shn.bin', FILEZ='shz.bin', & 
			    FILEV='shv.bin'
 CHARACTER*12 HEADER 
 INTEGER I, J, K, L, P, IS, LJ, MJ, LI, DOM, IND
! INTEGER LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP) 
! REAL*4 ALF(JMAX,NANCH), LONP(NP), LATP(NP), X(NP,0:NN)  
 INTEGER, ALLOCATABLE :: LL(:), MM(:), DM(:), ANC(:), WET(:)
 REAL*8, ALLOCATABLE :: ALF(:,:), LONP(:), LATP(:), X(:,:)  
! REAL*4 BETAS(0:LMAX,0:NN), ES(0:LMAX)
! REAL*4 BETAU(0:LMAX,0:NN), EU(0:LMAX)
! REAL*4 BETAN(0:LMAX,0:NN), EN(0:LMAX)
! REAL*4 BETAV(0:LMAX,0:NN), EV(0:LMAX)      
! REAL*4 RESH, IMSH, AAVV(0:NN), BBVV(0:NN)
 REAL*8, ALLOCATABLE :: BETAS(:,:), ES(:)
 REAL*8, ALLOCATABLE :: BETAU(:,:), EU(:)
 REAL*8, ALLOCATABLE :: BETAN(:,:), EN(:)
 REAL*8, ALLOCATABLE ::  BETAV(:,:), EV(:)      
 REAL*8 RESH, IMSH
 REAL*8, ALLOCATABLE :: AAVV(:), BBVV(:) 
! COMPLEX*16 LONG_TABLE(0:LMAX,NP), OC(JMAX) 
! COMPLEX*16 IIII(JMAX,0:NN), ZE(JMAX,0:NN), SE(JMAX,0:NN)
! COMPLEX*16 AAAA(JMAX,0:NN), AAAA_MOD(JMAX,0:NN)
! COMPLEX*16 BBBB(JMAX,0:NN), BBBB_MOD(JMAX,0:NN) 
! COMPLEX*16 HHHH(JMAX,0:NN), KKKK(JMAX,0:NN) 
! COMPLEX*16 S(JMAX,0:NN), N(JMAX,0:NN), U(JMAX,0:NN), V(JMAX,0:NN) 
! COMPLEX*16 Z(JMAX,0:NN,0:SMAX)
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:)
 COMPLEX*16, ALLOCATABLE :: IIII(:,:), ZE(:,:), SE(:,:)
 COMPLEX*16, ALLOCATABLE :: AAAA(:,:), AAAA_MOD(:,:)
 COMPLEX*16, ALLOCATABLE :: BBBB(:,:), BBBB_MOD(:,:) 
 COMPLEX*16, ALLOCATABLE :: HHHH(:,:), KKKK(:,:) 
 COMPLEX*16, ALLOCATABLE :: S(:,:), N(:,:), U(:,:), V(:,:)
 COMPLEX*16, ALLOCATABLE :: Z(:,:,:)
 REAL*8 RHOE, RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW 
 INTEGER TSTART, T0, T1, TMAX
!
!
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
! CHARACTER*3, PARAMETER :: RFB_OPT='y'
 INTEGER, ALLOCATABLE :: JJ(:)
 REAL*8, ALLOCATABLE :: BETASr(:,:), ESr(:)
 REAL*8, ALLOCATABLE :: BETAUr(:,:), EUr(:)
 REAL*8, ALLOCATABLE :: BETANr(:,:), ENr(:)
 REAL*8, ALLOCATABLE :: BETAVr(:,:), EVr(:)
 REAL*8, ALLOCATABLE :: ROT(:,:,:)
 COMPLEX*16, ALLOCATABLE :: IIIIr(:,:), Zr(:,:,:)
! ========================================================
! ========================================================
!
!
! --- Determining the number of pixels
!
  call system_clock(COUNT_MAX=tmax)
  call system_clock(tstart)
  t0 = tstart
  open(1,file='anchor.tmp',status='old')
  read(1,*) na
  close(1)
  Write(*,*) "    - Found ", na, " anchor pixels in file px-lat.dat"
  npix=np
!
!
! ========================================================
! ========================================================
! ========================================================
!
! --- Allocate memory space
!
 ALLOCATE( LL(JMAX), MM(JMAX), DM(JMAX), ANC(NPIX), WET(NPIX) )
 ALLOCATE( ALF(JMAX,NA), LONP(NPIX), LATP(NPIX), X(NPIX,0:NN) )
 ALLOCATE( LONG_TABLE(0:LMAX,NPIX), OC(JMAX) )
 ALLOCATE( IIII(JMAX,0:NN), ZE(JMAX,0:NN), SE(JMAX,0:NN) )
 ALLOCATE( AAAA(JMAX,0:NN), AAAA_MOD(JMAX,0:NN) )
 ALLOCATE( BBBB(JMAX,0:NN), BBBB_MOD(JMAX,0:NN) )
 ALLOCATE( HHHH(JMAX,0:NN), KKKK(JMAX,0:NN) )
 ALLOCATE( Z(JMAX,0:NN,0:SMAX) )
 ALLOCATE( S(JMAX,0:NN), N(JMAX,0:NN), U(JMAX,0:NN), V(JMAX,0:NN) )
 ALLOCATE( BETAS(0:LMAX,0:NN), ES(0:LMAX) )
 ALLOCATE( BETAU(0:LMAX,0:NN), EU(0:LMAX) )
 ALLOCATE( BETAN(0:LMAX,0:NN), EN(0:LMAX) )
 ALLOCATE( BETAV(0:LMAX,0:NN), EV(0:LMAX) )     
 ALLOCATE( AAVV(0:NN), BBVV(0:NN) ) 
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
 ALLOCATE( JJ(0:3) )
 AllOCATE( BETASr(0:2,0:NN), ESr(0:2) )
 AllOCATE( BETAUr(0:2,0:NN), EUr(0:2) )
 AllOCATE( BETANr(0:2,0:NN), ENr(0:2) )
 AllOCATE( BETAVr(0:2,0:NN), EVr(0:2) )
 ALLOCATE( ROT(0:3,0:NN,0:NN) )
 ALLOCATE( IIIIr(0:3,0:NN), Zr(0:3,0:NN,0:SMAX) )
! ========================================================
! ========================================================

!
!
! ********************************************************
! ********************************************************
! ********************************************************
!
! --- Extracting the average Earth density from 
!     the TABOO or ALMA log files... GS July 09
!
 IF(CDE.GE. 0) IND=1 
 IF(CDE.EQ.-1) IND=2 
 CALL AVERAGE_EARTH_DENSITY(IND, RHOE)
 RHOI_O_RHOE_X3 = 3.*RHOI/RHOE 
 RHOW_O_RHOE_X3 = 3.*RHOW/RHOE 
 RHOI_O_RHOW    =    RHOI/RHOW   
!
!
! --- Pre-computing 'l' and 'm' corresponding to degree 'J'
!
	do j=1, jmax 
		mm(j)=mj(j) 
		ll(j)=lj(j)
		dm(j)=2-dom(j)
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
! --- Pre-computing degree 'J' for RFB components:
!     (l,m) = (0,0), (2,-1), (2,0), (2,1)  
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
                jj(ll(j)+mm(j)) = j
        endif
! ========================================================
! ========================================================
	enddo	
!    if(rfb_opt=='y') Write(*,*) '    - JJ = ',JJ ! test test test!
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a22,f8.3,a3)') '    - Pre-Computing (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a22,f8.3,a3)') '    - Pre-Computing (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
! --- Reading the ALFs table from <<sh.bin>>
! 
 	Write(*,*) '    - Reading the ALFs from file sh.bin'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	Close(3) 
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
!
! --- Examining the pixels table & extracting information
!
	Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a12)')header
	Enddo
	Do i=1, np 
		Read (1,*) lonp(i), latp(i), anc(i), k, wet(i) 	
	Enddo
	close(1)
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a29,f8.3,a3)') '    - Reading px-table.dat (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a29,f8.3,a3)') '    - Reading px-table.dat (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
!
! --- Reading the SH OF coefficients from shof.dat 
!
 	Write(*,*) '    - Reading the SH OF coeff. from shof.dat'
 	open(3,file='shof.dat',status='unknown')
 	do j=1, jmax   
		read(3,*) k, resh, imsh 
                oc(j)=cmplx(resh, imsh)	
 	enddo
 	close(3)
!
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
!
! --- Reading the <<I>> array 
	IIII(:,:)=0. 
	If(imode/=5) then 
 		Write(*,*) "    - Reading array 'I'"
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)
!
        call system_clock(t1)
        if(t1<t0) then
            Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
        else
            Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
        endif
        t0=t1
    Endif
!
!
!
! --- Reading the E and BETA arrays 
 	Write(*,*) "    - Reading arrays 'E' and 'BETA'"
	open(1,file='ebs.dat', status='unknown')    
	open(2,file='ebu.dat',status='unknown')    
	open(3,file='ebn.dat',status='unknown')    
	open(4,file='ebv.dat',status='unknown')    

	do li=0, lmax
		read(1,*) l, Es(l), (betas(l,k), k=0,nn) 
		read(2,*) l, Eu(l), (betau(l,k), k=0,nn) 
		read(3,*) l, En(l), (betan(l,k), k=0,nn) 
		read(4,*) l, Ev(l), (betav(l,k), k=0,nn) 
	enddo
	close(4); close(3) ; close(2) ; close(1) 
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
!
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
 if(rfb_opt=='y') then
!
! --- Reading the Er and BETAr arrays 
 	Write(*,*) "    - Reading arrays 'Er' and 'BETAr'"
	open(1,file='ebsr.dat', status='unknown')    
	open(2,file='ebur.dat',status='unknown')    
	open(3,file='ebnr.dat',status='unknown')    
	open(4,file='ebvr.dat',status='unknown')    

	do li=0, 2
		read(1,*) l, Esr(l), (betasr(l,k), k=0,nn) 
		read(2,*) l, Eur(l), (betaur(l,k), k=0,nn) 
		read(3,*) l, Enr(l), (betanr(l,k), k=0,nn) 
		read(4,*) l, Evr(l), (betavr(l,k), k=0,nn) 
	enddo
	close(4); close(3) ; close(2) ; close(1) 
!     write(*,*)"ESr=",ESr 
!     write(*,*)"EUr=",EUr 
!     write(*,*)"ENr=",ENr 
!     write(*,*)"EVr=",EVr 
!     write(*,*)"BETASr=",BETASr 
!     write(*,*)"BETAUr=",BETAUr 
!     write(*,*)"BETANr=",BETANr 
!     write(*,*)"BETAVr=",BETAVr 
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
!
! --- Reading the ROTATIONAL FEEDBACK matrix 
 	Write(*,*) "    - Reading matrix 'ROT'"
    open(1,file='rfb.dat',status='unknown')
!  
    do li=1, 4*(NN+1)
      read(1,*) j, k, (ROT(j,k,l), l=0, nn)
    enddo
!
    close(1)
!
!      write(*,*)"rot0=",rot(0,:,:) 
!      write(*,*)"rot1=",rot(1,:,:) 
!      write(*,*)"rot2=",rot(2,:,:) 
!      write(*,*)"rot3=",rot(3,:,:) 
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
    endif
    t0=t1
!
!
! --- Computing rotational Load contributors for degrees 'l = 0, 2'
    IIIIr(0,:) = matmul(ROT(0,:,:), (IIII(jj(0),:) - IIII(jj(2),:)/sqrt(5.))) ! l=0, m=0
    IIIIr(1,:) = 0. ! IIII(jj(1),:)                                           ! l=2, m=-1
    IIIIr(2,:) = matmul(ROT(2,:,:), (IIII(jj(0),:) - IIII(jj(2),:)/sqrt(5.))) ! l=2, m=0
    IIIIr(3,:) = matmul(ROT(3,:,:), IIII(jj(3),:))                            ! l=2, m=1
    IIIIr(:,:) = RHOI * IIIIr(:,:)
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a44,f8.3,a3)') "    - Computing rotational potential 'Ir' (",(t1-t0+tmax)/1000.," s)"
    else
        Write(*,'(a44,f8.3,a3)') "    - Computing rotational potential 'Ir' (",(t1-t0)/1000.," s)"
    endif
    t0=t1
!
 endif
! ========================================================
! ========================================================
!
!
! --- Computing the eustatic Z array...
	ze(:,:) = 0.  		
	do k=0,nn	
		ze(:,k) = - rhoi_o_rhow*(iiii(1,k)/oc(1))*oc(:)
	enddo
!
!
! --- Computing the eustatic S array...
	se(:,:) = 0.
	se(1,:) = - rhoi_o_rhow*(iiii(1,:)/oc(1)) 
!
!
! --- Computing the A array...
	aaaa(:,:)=0.
	do j=1, jmax 
 	      do k=0,NN
              aaaa(j,k) = ES(ll(j))*IIII(j,k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIII(j,p)            )*BETAS(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIII(j,p)-IIII(j,p-1))*BETAS(ll(j),k-p)  
 	      enddo 
	      aaaa(j,k)= RHOI_o_RHOE_x3*aaaa(j,k)
	      enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Sa: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              aaaa(j,k) = aaaa(j,k) + ESr(ll(j))*IIIIr(ll(j)+mm(j),k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)                       )* &
                                   BETASr(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-IIIIr(ll(j)+mm(j),p-1))* &
                                   BETASr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
!
! --- Computing the ocean average of A 
	do k=0, NN 
		aavv(k)=0.
		do j=1, jmax 
		aavv(k) = aavv(k) + & 
		          dm(j)*real(oc(j)*conjg(aaaa(j,k)))/oc(1)  
		enddo
	enddo 
!
!
! --- Computing the modified ocean average of A 
	aaaa_mod(:,:) = aaaa(:,:)
	aaaa_mod(1,:) = aaaa(1,:)-aavv(:) 
!
! print *,'Task ',my_rank,' starting loop 1'
!
! --- Computing the R-array...
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP     SHARED(X,ALF,ANC,DM,AAAA_MOD,LONG_TABLE,MM,NPIX) &
!$OMP     SCHEDULE(GUIDED)
	do i=1, npix
	  do k=0, NN
	    x(i,k)=0.			
	    do j=1, jmax  
           x(i,k) = x(i,k) + & 
			 ALF(j,anc(i))*dm(j)*real(aaaa_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	  enddo
	enddo 
!$OMP END PARALLEL DO
!
! print *,'Task ',my_rank,' starting loop 2'

        hhhh = 0.
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,HHHH,X,ALF,ANC,LONG_TABLE,MM,NPIX) &
!$OMP       SCHEDULE(GUIDED)
      do j=1,jmax
 	  do i=1,npix
	    if(wet(i)==1) then
	      do k=0, NN
                 hhhh(j,k) = hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i)) 	          
	      end do
	    end if
	  end do
	end do   
!$OMP END PARALLEL DO
!
! It was:
!	do k=0, NN
!	do j=1, jmax 
!	    hhhh(j,k)=0.  
!	    do i=1, np 
!	    if(wet(i)==1) hhhh(j,k) = & 
!			  hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
!	    enddo
!	enddo
!	enddo 
!
!
!
!
	hhhh(:,:)=hhhh(:,:)/float(np) + ze(:,:)
!
!
! --- Initializing the Z and S arrays 
	Z(:,:,0) = ZE(:,:)
	S(:,:)   = SE(:,:) 
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y') then
	    Zr(0,:,0) = matmul(ROT(0,:,:), (Z(jj(0),:,0) - Z(jj(2),:,0)/sqrt(5.))) ! l=0, m=0
	    Zr(1,:,0) = 0. ! Z(jj(1),:,0)                    ! l=2, m=-1
	    Zr(2,:,0) = matmul(ROT(2,:,:), (Z(jj(0),:,0) - Z(jj(2),:,0)/sqrt(5.))) ! l=2, m=0
	    Zr(3,:,0) = matmul(ROT(3,:,:),  Z(jj(3),:,0))                         ! l=2, m=1
	    Zr(:,:,0) = RHOW * Zr(:,:,0)
        endif
! ========================================================
! ========================================================
!
    call system_clock(t1)
    if(t1<t0) then
    Write(*,'(a34,f8.3,a3)') "    - Computing 'A', 'Z0', 'S0' (",(t1-t0+tmax)/1000.," s)"
    else
    Write(*,'(a34,f8.3,a3)') "    - Computing 'A', 'Z0', 'S0' (",(t1-t0)/1000.," s)"
    endif
    t0=t1
!
!
! --- No recursion for the "Explicit approach"
	if(imode==6.or.imode==7.or.imode==3.or.SMAX==0) goto 2000     
!
!
!
! -----------------------
! ---    Recursion    ---
! -----------------------
!
	Write(*,*) "    - Starting the recursion"
!
 	do 1000 is = 1, SMAX     
!
!
        write(*,'(a12,i2,a3,i2)') '     - step ', is, ' of', SMAX
!
!
! --- Computing the 'B' array...
        bbbb(:,:)=0.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,K,P) &
!$OMP     SHARED(BBBB,ES,LL,MM,Z,IS,BETAS,rhow_o_rhoe_x3,ESr,BETASr,ROT,Zr) &
!$OMP     SCHEDULE(GUIDED)
        do j=1, jmax 
	      do k=0,NN
	      bbbb(j,k) = ES(ll(j))*Z(j,k,is-1)        
	      do p=0, k
	      if(p==0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-0.           )*BETAS(ll(j),k-p)
	      if(p/=0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-Z(j,p-1,is-1))*BETAS(ll(j),k-p)		 
	      enddo 
              bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
              enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Sb: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              bbbb(j,k) = bbbb(j,k) + ESr(ll(j))*Zr(ll(j)+mm(j),k,is-1)        
 	      do p=0, k
              if(p==0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,is-1)-0.            )* &
                                   BETASr(ll(j),k-p)
    	      if(p/=0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,is-1)-Zr(ll(j)+mm(j),p-1,is-1))* &
                                   BETASr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!$OMP END PARALLEL DO
!	
!
! --- Computing the ocean-average of array B array
	bbvv(:)=0.
	do k=0, NN 
		bbvv(k)=0.
		do j=1, jmax 
		bbvv(k) = bbvv(k) + & 
		          dm(j)*real(oc(j)*conjg(bbbb(j,k)))/oc(1)
		enddo
	enddo 
!
!
! --- Computing modified 'B' array
	bbbb_mod(:,:)=bbbb(:,:)
	bbbb_mod(1,:)=bbbb(1,:)-bbvv(:) 
!
!
! --- Computing array K...
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,K,J) &
!$OMP    SHARED(X,ALF,ANC,DM,BBBB_MOD,LONG_TABLE,MM,NPIX) SCHEDULE(GUIDED)
	do i=1, npix 
!	
	do k=0, NN
	    x(i,k)=0.			
	    do j=1, jmax  
            x(i,k) = x(i,k) + & 
			    ALF(j,anc(i))*dm(j)*real(bbbb_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	enddo

	enddo 
!$OMP END PARALLEL DO
!
!
!
      kkkk = 0.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,KKKK,X,ALF,ANC,LONG_TABLE,MM,NPIX) &
!$OMP       SCHEDULE(GUIDED)
       do j=1,jmax
 	  do i=1,npix
	    if(wet(i)==1) then
	      do k=0, NN
                 kkkk(j,k) = kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i)) 	          
	      end do
	    end if
	  end do
	end do   
!$OMP END PARALLEL DO
!
!--IT WAS:---------------------------------------------------------------------------
!	do k=0, NN
!	do j=1, jmax 
!	    kkkk(j,k)=0.  
!	    do i=1, np 
!	    if(wet(i)==1) kkkk(j,k) = & 
!			  kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
!	    enddo
!	enddo
!	enddo 
!------------------------------------------------------------------------------------
!
!
	kkkk(:,:)=kkkk(:,:)/float(np) 
!
!
! --- Solving for arrays 'Z' and 'S' 
	Z(:,:,is) = HHHH(:,:) + KKKK(:,:) 
	S(:,:)    = AAAA_MOD(:,:) + SE(:,:) + BBBB_MOD(:,:)
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y') then
	    Zr(0,:,is) = matmul(ROT(0,:,:), (Z(jj(0),:,is) - Z(jj(2),:,is)/sqrt(5.))) ! l=0, m=0
	    Zr(1,:,is) = 0. ! Z(jj(1),:,0)                    ! l=2, m=-1
	    Zr(2,:,is) = matmul(ROT(2,:,:), (Z(jj(0),:,is) - Z(jj(2),:,is)/sqrt(5.))) ! l=2, m=0
	    Zr(3,:,is) = matmul(ROT(3,:,:),  Z(jj(3),:,is))                         ! l=2, m=1
	    Zr(:,:,is) = RHOW * Zr(:,:,is)
        endif
! ========================================================
! ========================================================
!
!
    call system_clock(t1)
    if(t1<t0) then
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0+tmax)/1000.,' s)'
    else
        Write(*,'(a8,f8.3,a3)') '      (',(t1-t0)/1000.,' s)'
    endif
    t0=t1

! ------------------------------
! ---    End of recursion    ---
! ------------------------------
!
1000 CONTINUE      
!
!
!
2000 CONTINUE       
!
!
! --- Eustatic solution: U=0, S=N, V=0  
        if(imode==3.or.smax==0) then 
				U(:,:) = 0.
				S(:,:) = SE(:,:)
				N(:,:) = S(:,:)
				V(:,:) = 0. 
				goto 3000
				endif
!
! --- Array "B" for vertical displacement 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EU(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAU(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAU(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Ub: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              bbbb(j,k) = bbbb(j,k) + EUr(ll(j))*Zr(ll(j)+mm(j),k,SMAX)        
 	      do p=0, k
              if(p==0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,SMAX)-0.            )* &
                                   BETAUr(ll(j),k-p)
    	      if(p/=0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,SMAX)-Zr(ll(j)+mm(j),p-1,SMAX))* &
                                   BETAUr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
! --- Array "A" for vertical displacement 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = EU(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAU(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAU(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Ua: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              aaaa(j,k) = aaaa(j,k) + EUr(ll(j))*IIIIr(ll(j)+mm(j),k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-0.                    )* &
                                   BETAUr(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-IIIIr(ll(j)+mm(j),p-1))* &
                                   BETAUr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
!
! --- Vertical displacement
	U(:,:) = aaaa(:,:) + bbbb(:,:)

!
! --- Array "B" for Geoid heigth   
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EN(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAN(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAN(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Nb: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              bbbb(j,k) = bbbb(j,k) + ENr(ll(j))*Zr(ll(j)+mm(j),k,SMAX)        
 	      do p=0, k
              if(p==0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,SMAX)-0.            )* &
                                   BETANr(ll(j),k-p)
    	      if(p/=0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,SMAX)-Zr(ll(j)+mm(j),p-1,SMAX))* &
                                   BETANr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
! --- Array "A" for Geoid heigth  
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = EN(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAN(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAN(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Na: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              aaaa(j,k) = aaaa(j,k) + ENr(ll(j))*IIIIr(ll(j)+mm(j),k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-0.                    )* &
                                   BETANr(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-IIIIr(ll(j)+mm(j),p-1))* &
                                   BETANr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
! --- Geoid undulations 
	N(:,:) = aaaa(:,:) + bbbb(:,:)
!
! --- Adding a constant to geoid undulations 
	N(1,:) = N(1,:) +  SE(1,:) - AAVV(:) - BBVV(:)  
!
! --- Geoid undulations (previous formulation) 
!       N(:,:) = S(:,:) + U(:,:)
!
!
! --- Array "B" for horizontal displacement 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EV(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAV(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAV(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Vb: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              bbbb(j,k) = bbbb(j,k) + EVr(ll(j))*Zr(ll(j)+mm(j),k,SMAX)        
 	      do p=0, k
              if(p==0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,SMAX)-0.            )* &
                                   BETAVr(ll(j),k-p)
    	      if(p/=0) bbbb(j,k) = bbbb(j,k)-(Zr(ll(j)+mm(j),p,SMAX)-Zr(ll(j)+mm(j),p-1,SMAX))* &
                                   BETAVr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
! --- Array "A" for horizontal displacement 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = EV(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAV(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAV(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
        if(rfb_opt=='y' .and. (ll(j)==0 .or. ll(j)==2) .and. mm(j)/=2) then
!	      write(*,*) "Rotational feedback added on Va: l=",ll(j)," and m=",mm(j),", j=",j
!
 	      do k=0,NN
              aaaa(j,k) = aaaa(j,k) + EVr(ll(j))*IIIIr(ll(j)+mm(j),k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-0.                    )* &
                                   BETAVr(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIIIr(ll(j)+mm(j),p)-IIIIr(ll(j)+mm(j),p-1))* &
                                   BETAVr(ll(j),k-p) 
 	      enddo 
	      enddo
	endif
! ========================================================
! ========================================================
	enddo
!
!
! --- Horizontal displacement
	V(:,:) = aaaa(:,:) + bbbb(:,:)
!
!
3000 CONTINUE 
!
!
!
! --- Writing SH coefficients of sea level variations, S 
 	open(3,file=files,status='unknown',form='unformatted') 
        write(3) S ; close(3)
!
! --- Writing SH coefficients of vertical displacement, U  
 	open(3,file=fileu,status='unknown',form='unformatted') 
        write(3) U ; close(3) 
!
! --- Writing SH coefficients of geoid undulations, N    
 	open(3,file=filen,status='unknown',form='unformatted') 
        write(3) N ; close(3) 
!
! --- Writing SH coefficients of 'reduced' sea level change, Z=OS 
 	open(3,file=filez,status='unknown',form='unformatted') 
        write(3) Z ; close(3)
!
! --- Writing SH coefficients of horizontal displacement, V 
 	open(3,file=filev,status='unknown',form='unformatted') 
        write(3) V ; close(3)
!
	Write(*,*) "    +++ SLE solved <3 +++" 
!
!
!
! ========================================================
! ========================================================
! ========================================================
!
! --- Free up memory space
!
 DEALLOCATE( LL, MM, DM, ANC, WET )
 DEALLOCATE( ALF, LONP, LATP, X )
 DEALLOCATE( LONG_TABLE, OC )
 DEALLOCATE( IIII, ZE, SE )
 DEALLOCATE( AAAA, AAAA_MOD )
 DEALLOCATE( BBBB, BBBB_MOD )
 DEALLOCATE( HHHH, KKKK )
 DEALLOCATE( Z )
 DEALLOCATE( S, N, U, V )
 DEALLOCATE( BETAS, ES )
 DEALLOCATE( BETAU, EU )
 DEALLOCATE( BETAN, EN )
 DEALLOCATE( BETAV, EV )     
 DEALLOCATE( AAVV, BBVV ) 
! ========================================================
!             !NEW! ROTATIONAL  FEEDBACK !NEW!
! ========================================================
 DEALLOCATE( BETASr, ESr )
 DEALLOCATE( BETAUr, EUr )
 DEALLOCATE( BETANr, ENr )
 DEALLOCATE( BETAVr, EVr )     
 DEALLOCATE( IIIIr, Zr )     
 DEALLOCATE( ROT, JJ )     
! ========================================================
! ========================================================
!
    call system_clock(t1)
    if(t1<t0) then
	    Write(*,'(a23,f8.3,a3)') '     - post-recursion (',(t1-t0+tmax)/1000.,' s)' 
!	    Write(*,'(a19,f8.3,a2)') '       total time: ',(t1-tstart+tmax)/1000.,' s' 
    else
	    Write(*,'(a23,f8.3,a3)') '     - post-recursion (',(t1-t0)/1000.,' s)' 
!	    Write(*,'(a19,f8.3,a2)') '       total time: ',(t1-tstart)/1000.,' s' 
    endif
!
!
	End program SLE
!
!
