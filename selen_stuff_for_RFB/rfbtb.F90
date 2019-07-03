!
!
! This is program "RFBTB.F90" 
!
! *** Modified by RH March 2018, for rotational feedback in version 2.9.12 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! =============================================================================
!
! 1.: Computes the arrays BETAr(l,k) and Er(l), the TIDAL rheology--dependent 
! quantities involved in the solution of the SLE by the PS method. Here the 
! program TABOO is used as a subroutine (in tb-utilities.F90) to compute the 
! TIDAL LDCs. The structure is based on the original tb.F90 for the LOADING
! analysis.
! 2.: Computes the matrix RFB(l+m,NN,NN), including all rheology-time-dependent
! quantities for the rotational feedback involved in the solution of the SLE. 
! The computation follows Milne & Mitrovica, 1998 
!	
!
! Input file:
!	- task_1.dat 
!         [built here according to the input parameters of "data.inc"] 
!
! Output files:  
! 	- spectrum.dat
!	- ih.dat, il.dat, & ik.dat 
!	- h.dat,   l.dat, &  k.dat 
!	- ebsr.dat    ! Sea level green function 
!	- ebur.dat    ! Vertical displacement Green function
!       - ebnr.dat    ! Geoid elevation Green function
!       - ebvr.dat    ! Hhorizontal displacement Green function
!       - rfb.dat     ! Rotational potential matrix for the RFB
!
! =============================================================================
!
!

! INCLUDE "harmonics.f90"
 PROGRAM RFBBB
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
! --- Utility
 CHARACTER*1000 :: line
 INTEGER I, J, L, K, M, P, KK, NOUT  
 INTEGER, PARAMETER :: NV_MAX=9             ! max. number of viscoelastic layers 
 INTEGER, PARAMETER :: NROOTS_MAX=4*NV_MAX  ! max. number of roots    
 INTEGER, PARAMETER :: NROOTS = 4*NV
 INTEGER, PARAMETER :: NIMP=100
! REAL*4  DEN
 REAL*8 LTHIC, JUNK 
 CHARACTER*1 CHAR       
 CHARACTER*100 SS(NIMP)
 CHARACTER*200 ROW
!
! --- LOAD LOVE NUMBERS
 REAL*8  KLE2, KLF2                      ! Elastic, Fluidal k (deg=2, load) 
 REAL*8, ALLOCATABLE :: KLV2(:)          ! Residues of k (deg=2, load)
 REAL*8, ALLOCATABLE :: S2(:)            ! Relaxation times (deg=2)
!
! --- TIDE LOVE NUMBERS
 INTEGER, PARAMETER :: L1=1, LMP1=L1+1, L2=2
 REAL*4  VIS_LM, VIS_UM, VIS_TZ, VIS       
 REAL*4, ALLOCATABLE ::  R_K (:,:) 
 REAL*4, ALLOCATABLE ::  R_H (:,:)
 REAL*4, ALLOCATABLE ::  R_L (:,:) 
 REAL*4, ALLOCATABLE ::  TEMP(:,:) 
 REAL*4, ALLOCATABLE ::  H_E(:), K_E(:), L_E(:), K_F(:)
 REAL*8, ALLOCATABLE ::  VSC(:) 
 REAL*8, ALLOCATABLE ::  BETASr(:,:), ESr(:)
 REAL*8, ALLOCATABLE ::  BETAUr(:,:), EUr(:)
 REAL*8, ALLOCATABLE ::  BETANr(:,:), ENr(:)
 REAL*8, ALLOCATABLE ::  BETAVr(:,:), EVr(:)
!
! --- ROTATIONAL FEEDBACK
 REAL*8, ALLOCATABLE :: LAMBDA(:), EE(:)
 REAL*8  EEpart, D2            
 REAL*8, ALLOCATABLE :: D3(:), D4(:)
 REAL*8, ALLOCATABLE :: ROT(:,:,:)
!
 real*8, parameter :: Omega = 2.*PI/(86164.1)       ! (rad/s) mean Earth rotation rate (23h56m4.1s)
 real*8, parameter :: CF = 2.*PI/(433.1 * 86400.)   ! (rad/sec) Chandler Frequency
 real*8, parameter :: Grav = 9.80665                ! (m/s^2) standard gravity (WGS-84)
 real*8, parameter :: Re = 6371009.                 ! (m) mean radius (GRS 1980)
 real*8, parameter :: Re_eq = 6378136.6             ! (m) Earth equatorial radius
 real*8, parameter :: Me = 5.973698968d24           ! (kg) Earth mass
 real*8, parameter :: Ci = 0.3307007*Me*Re_eq*Re_eq ! (kg.m^2) Earth axial moment of inertia
 real*8, parameter :: Ai = 0.3296108*Me*Re_eq*Re_eq ! (kg.m^2) Earth equatorial moment of inertia
!
! ========================================================
! ========================================================
 IF (rfb_opt == 'y') THEN !RFB-mode on?
   write(*,*)"    - Rotational Feedback is switched ON"
! ========================================================
! ========================================================
!
! --- Allocate memory space
!
 ALLOCATE( R_K (0:L2,NROOTS_MAX) ) 
 ALLOCATE( R_H (0:L2,NROOTS_MAX) )
 ALLOCATE( R_L (0:L2,NROOTS_MAX) )
 ALLOCATE( TEMP(0:L2,NROOTS_MAX) )
 ALLOCATE( H_E(0:L2), K_E(0:L2), L_E(0:L2), K_F(0:L2) ) 
 ALLOCATE( VSC(NV) )
 ALLOCATE( BETASr(0:L2,0:NN), ESr(0:L2) )
 ALLOCATE( BETAUr(0:L2,0:NN), EUr(0:L2) )
 ALLOCATE( BETANr(0:L2,0:NN), ENr(0:L2) )
 ALLOCATE( BETAVr(0:L2,0:NN), EVr(0:L2) )
 ALLOCATE( ROT(0:3,0:NN,0:NN) )
 ALLOCATE( D3(0:NN), D4(0:NN) )
 ALLOCATE( KLV2(0:(NROOTS-1)), S2(0:(NROOTS-1)) )
 ALLOCATE( LAMBDA(0:(NROOTS-2)), EE(0:(NROOTS-2)) )
!
! ********************************************************
! ********************************************************
! ********************************************************
!
! --- VISCOUS RELAXATION TIMES (ka^-1): equal for LOAD and TIDES
!
 j=0
 open(5,file='ss.dat',status='unknown')
 do !read correct line 
   read(5,'(a)') line
   if(line(:4)=='   2') then
     read(line,*) l, junk, s2(j)
     j=j+1
   endif
   if(j>(NROOTS-1)) exit
 enddo
 close(5)
 s2(:) = -s2(:)
!
!
! --- LOADING k LOVE NUMBER (degree l=2): from previous TABOO call
!
 open(5,file='kk.dat',status='unknown')
 do
   read(5,'(a)') line
   if(line(:4)=='   2') exit
 enddo
 close(5)
!
 read(line,*) l, kle2, klf2, klv2  
!
!
! --- TIDAL k LOVE NUMBER (degree l=2): calling TABOO again
!
! # Reading the lithospheric thickness & viscosity profile from "visco.dat"
!
 open(12,file=visco_model,status='unknown') 
 read(12,'(a200)') row 
      call scan_string (row, 1, ss, nout)
      call CHAR100_2_REAL(ss(1), LTHIC)	
 do j=1, nv 
     read(12,'(a200)') row
     call scan_string (row, 1, ss, nout)
     call CHAR100_2_REAL(ss(1),vsc(j))
 enddo 
 close(12)
!
!
! # Building the input file <<task_1.dat>>... 
!
 open(1,file='task_1.dat',status='unknown') 
!  
! General settings
 write(1,'(a6)' )'Active' 
 write(1,'(a16)')'Harmonic_Degrees'
 write(1,*) l1, l2   ! min and max degrees  
 write(1,*) '0'      ! verbose mode 
 write(1,*) '0'      ! i_tides 
!
! Earth model 
 write(1,'(a10)')'Make_Model' 
!
 write(1,*)  NV      ! 
 write(1,*)  CDE     ! 
 write(1,*)  LTHIC   ! thickness of the lithosphere (km) 
 write(1,*) '0' 
! 
 do j=1, nv
     write(1,*) vsc(j)
 enddo
!
! Requesting normalized residues 
 write(1,'(a19)')'Normalized_Residues'  
 do k=1,3
     write(1,*) '1' 
 enddo
!
! Requesting LDCs 
 write(1,'(a15)')'El_Fluid_Viscel'  
 do k=1,3
     write(1,*) '1' 
 enddo
 close(1) 
!
! ******************************* 
!
  write(*,*) '    - Calling TABOO for TIDAL degree l=2' 
  call TABOO
!
! *******************************
!
!
!###################################################
! Building the "E" arrays (Elastic Green functions)
!###################################################
!
! --- Reading the elastic TDCs... 
!
   If(imode==1.or.imode==2.or.imode==5)  then 
     OPEN(1,FILE='h.dat',STATUS='old')
     OPEN(2,FILE='k.dat',STATUS='old')
     OPEN(3,FILE='l.dat',STATUS='old')
!
     do J=1, 2
       READ(1,'(a1)') char
       READ(2,'(a1)') char
       READ(3,'(a1)') char
     enddo
!
     do l=L1, L2
      	READ (1,'((i4,1x,24(1x,e20.8)))') p, h_e(l)
      	READ (2,'((i4,1x,24(1x,e20.8)))') p, k_e(l), k_f(l)
      	READ (3,'((i4,1x,24(1x,e20.8)))') p, l_e(l)
     enddo
!
     CLOSE(1)
     CLOSE(2)
     CLOSE(3)
!
     write(*,*) '    - Approximated fluid k Love number (tidal, l=2) = ', k_f(2)
!
   endif !imode == 1, 2 or 5
!
!
! --- Initialization of the "E" arrays: Eustatic solution 
!
   ESr(:)=0.   
   EUr(:)=0.
   ENr(:)=0.
   EVr(:)=0.
!
!
! --- Setting the ELASTIC Green Functions for Degree 1
!
 IF(DEG1==1) THEN 				
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree l=1 are set by the corresponding Love numbers 
!
   l=1
!
!   GSC or Elastic solution:
	If(imode==1.or.imode==2.or.imode==5)  then 
	    EUr(l) = h_e(l) 
        EVr(l) = l_e(l)
	    ESr(l) = ENr(l) - EUr(l)	
	Endif
!
!   Woodward solution: 
	If(imode==4)  then 
		ENr(l) = 1. 
	Endif 	
!
!   den=2.*float(nd)+1. ! NOT FOR THE TIDALS
!
!   EU(nd) = EU(nd)/den ! NOT FOR THE TIDALS
!   EN(nd) = EN(nd)/den ! NOT FOR THE TIDALS
!   EV(nd) = EV(nd)/den ! NOT FOR THE TIDALS
!   ES(nd) = ES(nd)/den ! NOT FOR THE TIDALS
!
 ENDIF 
!
!
! --- Setting the ELASTIC Green Functions for Degree > 1 
! 
 Do l=LMP1, L2
!
!   GSC or Elastic solution:
	If(imode==1.or.imode==2.or.imode==5)  then 
	    ENr(l) =   1.+k_e(l)
	    EUr(l) =      h_e(l)    
	    EVr(l) =      l_e(l) 
    Endif
!       
!   Woodward solution:
	if(imode==4) then 
		ENr(l) = 1.
	endif
 enddo
!
    ESr(LMP1:) = ENr(LMP1:) - EUr(LMP1:) 
!
!   den=	(2.*float(l)+1.) ! NOT FOR THE TIDALS 
!
!   EU(l) = EU(l)/den ! NOT FOR THE TIDALS
!	EN(l) = EN(l)/den ! NOT FOR THE TIDALS
!	EV(l) = EV(l)/den ! NOT FOR THE TIDALS
!	ES(l) = ES(l)/den ! NOT FOR THE TIDALS
!
!
!############################################################
! Building the "BETA" arrays (Visco-Elastic Green functions)
!############################################################
!
! --- Reading the normalized residues and relaxation times
   if(imode==1.or.imode==5) then
     open(1,file='ih.dat',status='old')
     open(2,file='ik.dat',status='old')
     open(3,file='il.dat',status='old')
     open(4,file='spectrum.dat',status='old')
!
     do k=1,2
           read(1,'(a1)') char
           read(2,'(a1)') char
           read(3,'(a1)') char	 
     enddo
     do k=1,7
	   read(4,'(a1)') char
     enddo
!
     r_h(:,:) =0.
     r_k(:,:) =0.
     r_l(:,:) =0. 
     temp(:,:)=0. 
!
     do l = L1, L2
           do j=1, 4  
           	   read(j, '(a1)') char
	       enddo
!     
	       do m=1, nroots
	           read  (1, '(i4,1x,e15.7)') p,  r_h (l,m)
	           read  (2, '(i4,1x,e15.7)') p,  r_k (l,m)
	           read  (3, '(i4,1x,e15.7)') p,  r_l (l,m)
	           read  (4, '(i4,1x,5(e15.7,1x))') kk, junk, junk, junk, temp(l,m)
!
               r_h(l,m)=-r_h(l,m)      		! **  I have (h/s) in memory  ** !
	           r_k(l,m)=-r_k(l,m)      		! **  I have (k/s) in memory  ** !
	           r_l(l,m)=-r_l(l,m)      		! **  I have (k/s) in memory  ** !
!
               temp(l,m) = temp(l,m)/1000.    	! ** relaxation times in kyrs ** !
           enddo
     enddo
     close(1); close(2); close(3); close(4)  
   endif !imode==1 or 5
!
!
!
! Initialization of the "BETA" arrays
   BETASr(:,:)=0. 
   BETAUr(:,:)=0. 
   BETANr(:,:)=0.
   BETAVr(:,:)=0.
!
!
! --------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree 1 ======
! --------------------------------------------------------------------
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree l=1 are set by the corresponding Love numbers 
!
   if(imode==1.or.imode==5) then
     IF(DEG1==1) THEN 
!
       do k=0, NN 
         l=1 
!         den = 2.*float(l)+1. ! NOT FOR TIDALS 
!
         do m=1, nroots
! 
           BETAUr(l,k) = BETAUr(l,k) + & 
                        (r_h(l,m))* & 
                        (1. - exp(-k*delta/temp(l,m))) !/den NOT FOR TIDALS
!
           BETAVr(l,k) = BETAVr(l,k) + & 
                        (r_l(l,m))* & 
                        (1. - exp(-k*delta/temp(l,m))) !/den NOT FOR TIDALS
!
         enddo
       enddo
       BETASr(l,:) = BETANr(l,:) - BETAUr(l,:)		
!
     ENDIF
!
!
! ----------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree > 1 ======
! ----------------------------------------------------------------------
!
     do k=0, NN 
       do l=lmp1, L2 
!         den = 2.*float(l)+1. ! NOT FOR TIDALS  
         do m=1, nroots
!
			BETANr(l,k) = BETANr(l,k) + & 
			             (r_k(l,m))* & 
			             (1.0 - exp(-k*delta/temp(l,m))) !/den NOT FOR TIDALS
! 
			BETAUr(l,k) = BETAUr(l,k) + & 
			             (r_h(l,m))* & 
			             (1.0 - exp(-k*delta/temp(l,m))) !/den NOT FOR TIDALS
!
			BETAVr(l,k) = BETAVr(l,k) + & 
			             (r_l(l,m))* & 
			             (1.0 - exp(-k*delta/temp(l,m))) !/den NOT FOR TIDALS
!				    
         enddo
       enddo
     enddo
     BETASr(LMP1:,:) = BETANr(LMP1:,:) - BETAUr(LMP1:,:)
!
   endif !imode==1 or 5
!
!
! ----------------------------------- 
! ====== Reporting the results ======
! ----------------------------------- 
!
!
 open(1,file='ebsr.dat',status='unknown')
 open(2,file='ebur.dat',status='unknown')
 open(3,file='ebnr.dat',status='unknown')
 open(4,file='ebvr.dat',status='unknown')
!  
 do l=0, L2
 	 write(1,*) l, ESr(l), (betasr(l,k), k=0, nn)
 	 write(2,*) l, EUr(l), (betaur(l,k), k=0, nn)
 	 write(3,*) l, ENr(l), (betanr(l,k), k=0, nn)		 
 	 write(4,*) l, EVr(l), (betavr(l,k), k=0, nn)
 enddo
!
 close(1)
 close(2)
 close(3)
 close(4) 
!
!
! ##############################################################################
! ---- Compute the True Polar Wanderer and rotational potential for the RFB ----
! ##############################################################################
!
! --- Initialize D2, D3(:), D4(:)
     D2 = 0.
     D3(:) = 0.
     D4(:) = 0.
!
   if(imode==1.or.imode==5) then
!
! --- Determine the roots lambda()
       call find_lambda(nroots-1,s2,r_k(2,1:nroots),lambda)
!
       do k=0,(NROOTS-2)
         EE(k) = (1.+klf2) * product(s2(:)-lambda(k)) / lambda(k)
!
         do l=0,(NROOTS-1)
           EEpart = klv2(l)/s2(l)
           if(l>0)        EEpart = EEpart * product(s2(:l-1)-lambda(k))
           if(l<(NROOTS-1)) EEpart = EEpart * product(s2(l+1:)-lambda(k))
           EE(k) = EE(k) + EEpart
         enddo
!
         if(k>0)        EE(k) = EE(k) / product(lambda(:k-1)-lambda(k))
         if(k<(NROOTS-2)) EE(k) = EE(k) / product(lambda(k+1:)-lambda(k))
       enddo
       EE(:) = -EE(:)
!
       D2 = (1.+klf2)*product(s2)/product(lambda) * delta
!
       do k=0,NN
         D3(k) = sum(  EE(:)/lambda(:) * (exp(-lambda(:)*(k-1)*delta) - exp(-lambda(:)*k*delta)))
         D4(k) = sum(klv2(:)/s2(:)     * (exp(-s2(:)    *(k-1)*delta) - exp(-s2(:)    *k*delta)))
       enddo
     endif
!
     ROT(:,:,:)=0
!
     do k=0,NN
       ROT(:,k,k) = 1.+kle2
     do l=0,k-1
       ROT(0,k,l)= D4(k-l)      ! l=0, m=0
!       ROT(1,k,l)= D2 + D3(k-l) ! l=2, m=-1
       ROT(2,k,l)= D4(k-l)      ! l=2, m=0
       ROT(3,k,l)= D2 + D3(k-l) ! l=2, m=1
     enddo
     enddo
!
! --- Computing rotational factors for degrees 'l = 0, 2'
     ROT(0,:,:) = ROT(0,:,:) * (-16./9./Ci)               ! l=0, m=0
!     ROT(1,:,:) = ROT(1,:,:) * sqrt(16./225.)*omega/Ai/cf ! l=2, m=-1
     ROT(2,:,:) = ROT(2,:,:) * 16./9./Ci/sqrt(5.)         ! l=2, m=0
     ROT(3,:,:) = ROT(3,:,:) * sqrt(16./225.)*omega/Ai/cf ! l=2, m=1
     ROT(:,:,:) = ROT(:,:,:) * PI * Re**6 * omega**2
     ROT(:,:,:) = ROT(:,:,:) / Grav
!
! ----------------------------------- 
! ====== Reporting the results ======
! ----------------------------------- 
!
     open(1,file='rfb.dat',status='unknown')
!  
     do j=0, 3
       do k=0, NN
 	     write(1,*) j, k, (ROT(j,k,l), l=0, nn)
       enddo
     enddo
!
     close(1)
!
! ========================================================
! ========================================================
! ========================================================
!
! --- Free up memory space
!
 DEALLOCATE( R_K, R_H, R_L )
 DEALLOCATE( TEMP )
 DEALLOCATE( H_E, K_E, L_E, K_F ) 
 DEALLOCATE( VSC )
 DEALLOCATE( BETASr, ESr )
 DEALLOCATE( BETAUr, EUr )
 DEALLOCATE( BETANr, ENr )
 DEALLOCATE( BETAVr, EVr )
 DEALLOCATE( ROT )
 DEALLOCATE( D3, D4 )
 DEALLOCATE( KLV2, S2 )
 DEALLOCATE( LAMBDA, EE )
!
! ========================================================
! ========================================================
 ELSE ! RFB-mode off?
   write(*,*)"    - Rotational Feedback is switched OFF -> nothing to do here!"
 ENDIF
! ========================================================
! ========================================================
!
!================
  end program RFBBB 
!================
!
!
!##############################################################################
!Subroutine which calculates the polynom Y(x)=c1*(x-x1)+c2*(x-x2)+...+cN*(x-xN)
!for given x-values. The input gives the parameters:
! - DEG: The degree N of the polynom
! - NPTS: Number of Y-values to be calculated/Number of given x-values 
! - X_ARRAY: Array of x-values on which Y(x) is calculated
! - X_OFF: The values x1,...,xN
! - COEF: The coefficients c1,...,cN
!##############################################################################
 SUBROUTINE make_poly(LSG, X_ARRAY, X_OFF, COEF, NPTS, DEG)
     IMPLICIT NONE
!
     INTEGER DEG, NPTS
     REAL*4  COEF(0:DEG)
     REAL*8  X_OFF(0:DEG)
     REAL*8  X_ARRAY(0:NPTS-1)
     REAL*8  LSG(0:NPTS-1)
!
     INTEGER K,L
     REAL*8  PART(0:NPTS-1)
!
     lsg(:) = 0.
     do k=0,(deg)
         part(:) = coef(k)
         if(k>0)          then
             do l=0, k-1
                 part(:) = part(:) * (x_array(:)-x_off(l))
             enddo
         endif
         if(k<(deg)) then
             do l=k+1, deg
                 part(:) = part(:) * (x_array(:)-x_off(l))
             enddo
         endif
         lsg(:) = lsg(:) + part(:)
     enddo
!
 END SUBROUTINE MAKE_POLY
!
!
!##############################################################################
!Recursive subroutine to find the positions where Y(s)=!=0 in a given accuracy
!The search is reported in a logfile, the found roots are recorded and counted
!##############################################################################
 RECURSIVE SUBROUTINE FIND_ROOTS(X_ARRAY, X_OFF, COEF, NPTS, DEG, EPS, ROOTS, ROOT_COUNT, LOGFILE_UNIT)
     IMPLICIT NONE
!
     INTEGER DEG, NPTS, ROOT_COUNT, LOGFILE_UNIT
     REAL*8 EPS
     REAL*8 X_ARRAY(0:NPTS-1)
     REAL*8 X_OFF(0:DEG)
     REAL*4 COEF(0:DEG)
     REAL*8 ROOTS(0:DEG-1)
!
     INTEGER I, J
     REAL*8 F(0:NPTS-1)
     REAL*8 X_PART(0:10)
     REAL*8 INC
!
!      write(*,*) "size x_array",size(x_array)
     call make_poly(F, X_ARRAY, X_OFF, COEF, NPTS, DEG)
!
     do i=1, npts-1
         if(abs(f(i))<eps) then
             write(logfile_unit,*) "potential root! f(i) = ",f(i)," < eps = ",eps
             if(abs(f(i))<=abs(f(i-1)) .and. abs(f(i))<abs(f(i+1))) then
                 if(root_count>deg-1) exit
                 roots(root_count) = x_array(i)
!                 write(*,*) "counter", root_count
                 root_count = root_count + 1
                 write(logfile_unit,*) "Root no. ",root_count," found at x =", x_array(i), "f = ", f(i)
             endif
         elseif( (f(i)/abs(f(i))) /= (f(i-1)/abs(f(i-1))) .and. abs(f(i-1))>eps ) then
             if(root_count>deg-1) exit
             write(logfile_unit,*) "Root found between x = ", x_array(i-1), " and x = ", x_array(i)
             write(logfile_unit,*) "                   f = ", f(i-1), " and f = ", f(i)
             inc = (x_array(i)-x_array(i-1))*1d-1
             do j=0,10
                 x_part(j) = j*inc
             enddo
             x_part(:) = x_part(:) + x_array(i-1)
             if(size(x_array)==size(x_part) .and. sum(x_array(:)-x_part(:))==0) then
                 roots(root_count) = x_array(i)
!                 write(*,*) "counter", root_count
                 root_count = root_count + 1
                 write(logfile_unit,*) "Root no. ",root_count," determined by mx precision at x =", x_array(i), "f = ", f(i)             
             else
                 write(logfile_unit,*) "---> looking deeper..."             
                 call find_roots(x_part, x_off, coef, size(x_part), deg, eps, roots, root_count, logfile_unit)
             endif
         endif
     enddo    
 END SUBROUTINE FIND_ROOTS
!
!
!##############################################################################
!-----   Subroutine which determines the roots LAMBDA of a polynom Y(s)   -----
!The polynom is defined as: Y(s)= c1*(s-s1)+c2*(s-s2)+...+cN*(s-sN)
!the input parametrs are:
! - DEG: The degree N of the polynom
! - S_OFF: The values s1,...,sN
! - COEF: The coefficients c1,...,cN
!##############################################################################
 SUBROUTINE FIND_LAMBDA(DEG, S_OFF, COEF, LAMBDA)
     IMPLICIT NONE
!
     INTEGER DEG
     REAL*4 COEF(0:DEG)
     REAL*8 S_OFF(0:DEG)
     REAL*8  LAMBDA(0:DEG-1)
!     
     INTEGER  I, E, N
     INTEGER  EXP_MIN, EXP_MAX
     INTEGER  ROOT_COUNT, LOGFILE_UNIT
     INTEGER, PARAMETER :: ITV_PTS = 89999
     REAL*8, PARAMETER :: EPS = 1.d-12
     REAL*8, ALLOCATABLE :: S(:)
!
     exp_min = int(log10(minval(s_off)))
     exp_max = int(log10(maxval(s_off))) + 1
     ALLOCATE( S(0:(EXP_MAX-EXP_MIN+1)*(ITV_PTS+1)-1) )
!
     n=0
     do e = exp_min, exp_max
         do i=0, itv_pts
             s(i+n)=(i+10000)*1d-5 * 10.**e
         enddo
         n=n+itv_pts+1
     enddo
!
     root_count = 0
!     
     logfile_unit=96
     open(logfile_unit,file='lambda_roots.log',status='unknown')
     write(logfile_unit,*) "################################################################"
     write(logfile_unit,*) "This is the log file for the root-finding algorithm of RFBTB.F90"
     write(logfile_unit,*) "################################################################"
     write(logfile_unit,*) "----------------------------------------------------------------"
     write(logfile_unit,*) "coef  = ",coef
     write(logfile_unit,*) "s_off = ",s_off
     write(logfile_unit,*) "s_off in [",minval(s_off)," , ",maxval(s_off),"]"
     write(logfile_unit,*) "x     in [",minval(s)," , ",maxval(s),"]"
!
     write(logfile_unit,*) "----------------------------------------------------------------"
     call find_roots(s, s_off, coef, size(s), deg, eps, lambda, root_count, logfile_unit)
     write(logfile_unit,*) "----------------------------------------------------------------"
     write(logfile_unit,*) "ROOTs found:"
     write(logfile_unit,*) "lambda = ", lambda
     write(logfile_unit,*) "----------------------------------------------------------------"
!
     if(root_count>size(lambda)) then
         write(logfile_unit,*) "WARNING: More roots (",root_count,") found then expected (",size(lambda),")!"
         write(*,*)            "WARNING: More roots (",root_count,") found then expected (",size(lambda),")!"
     elseif(root_count<size(lambda)) then
         write(logfile_unit,*) "WARNING: Less roots (",root_count,") found then expected (",size(lambda),")!"
         write(*,*)            "WARNING: Less roots (",root_count,") found then expected (",size(lambda),")!"
     endif
!
     close(logfile_unit)
!
     DEALLOCATE( S )
!
 END SUBROUTINE FIND_LAMBDA
!
!
