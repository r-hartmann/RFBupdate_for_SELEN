!
!
! This is program "TB.F90" 
!
! Last modified GS 04-11-2008 [INTEL port]
! Modified on August 2008, for version 2.7 (horizontal movements) 
! Also touched October 2008 for the PAULSON et al. profile
! ----->>> New as from August 2008: This program also computes arrays 
! Beta^v(l,k) and E^v(l), pertaining to the horizontal movements.
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed GS May 2010 - Porting under g95 
! *** Reviewed GS May 30 2010 - Degree 1 interface 
! *** Reviewed DM Dec 22 2010 - Adjusted REAL sizes for gfortran compatibility
!
! *** Adapted by RH March 2018 - Excluded TABOO code to tb-utilities.F90 
!                                to allow multiple usage of TABOO for the 
!                                rotational feedback analysis. 
!                                This tb.F90 is stil doing the LOADING analysis
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
! Computes the arrays Beta(l,k) and E(l), and Beta^u(l,k) and E^u(l), the only 
! rheology--dependent quantities involved in the solution of the SLE by the PS 
! method.  Here the program TABOO is used as a subroutine to compute the LDCs. 
!	
!
! Input file:
!	- task_1.dat 
!         [built here according to the input parameters of "data.inc"] 
!
! Output files:  
! 	- spectrum.dat
!	- ih.dat, il.dat, & ik.dat (ihh.dat, ill.dat, & ikk.dat
!	- h.dat,   l.dat, &  k.dat ( hh.dat,  ll.dat, &  kk.dat)
!	- ebs.dat    ! Sea level green function 
!	- ebu.dat    ! Vertical displacement Green function
!       - ebn.dat    ! Geoid elevation Green function
!       - ebv.dat    ! Hhorizontal displacement Green function
!
! =============================================================================
!
!

! INCLUDE "harmonics.f90"
 PROGRAM BB
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
 INTEGER, PARAMETER :: L1=1, LMP1=L1+1, L2=LMAX
!
 INTEGER, PARAMETER :: NV_MAX=9             ! max. number of viscoelastic layers 
 INTEGER, PARAMETER :: NROOTS_MAX=4*NV_MAX  ! max. number of roots    
 INTEGER, PARAMETER :: NIMP=100
 INTEGER J, L, K, M, P, KK, NOUT, NROOTS  
 REAL*4 R_K (0:LMAX,NROOTS_MAX), & 
        R_H (0:LMAX,NROOTS_MAX), &
	R_L (0:LMAX,NROOTS_MAX), & 
        TEMP(0:LMAX,NROOTS_MAX), & 
        H_E (0:LMAX), &
	K_E (0:LMAX), &
	L_E (0:LMAX), & 
	VIS_LM, VIS_UM, VIS_TZ, VIS, JUNK       
 REAL*8 VSC(NV)
! 
! Imported from AM.F90:---------
 INTEGER ND
 REAL*4 ES(0:LMAX), & 
        EU(0:LMAX), &
        EN(0:LMAX), &  
        EV(0:LMAX) 
 REAL*4 BETAS(0:LMAX,0:NN), & 
        BETAU(0:LMAX,0:NN), & 
        BETAN(0:LMAX,0:NN), & 
        BETAV(0:LMAX,0:NN)  
 REAL*4  DEN
! ------------------------------
 REAL*8 LTHIC 
 CHARACTER*1 CHAR       
 CHARACTER*100 SS(NIMP)
 CHARACTER*200 ROW
!
! =============================================================================
!
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
 write(1,*) '1'      ! i_loading 
!
! Earth model 
 write(1,'(a10)')'Make_Model' 
!
 NROOTS=4*NV   	     ! Number of Roots = 4*Number of v.e. layers  
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
!
!  
  write(*,*) '    - Calling TABOO' 
! ******************************* 
!
  call TABOO
!
! *******************************

!
!

!
!
! Initialization of the "BETA" arrays
  BETAS(:,:)=0. 
  BETAU(:,:)=0. 
  BETAN(:,:)=0.
  BETAV(:,:)=0.
!
!
! Reading the normalized residues and relaxation times
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
!
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
	      if(imode==2.or.imode==3.or.imode==4) then 
	      					r_h(l,m)=0.
	      					r_k(l,m)=0.
						r_l(l,m)=0.
			   			   endif 
	      if(imode==7)			   then
	      				        r_k(l,m)=0.
						   endif
!	      	      
	      temp(l,m) = temp(l,m)/1000.    	! ** relaxation times in kyrs ** !
!	      
	   enddo
        enddo
	close(1); close(2); close(3); close(4)  
!
!
!
!
!--- beta.f: reading the elastic LDCs... 
!
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
      	READ (2,'((i4,1x,24(1x,e20.8)))') p, k_e(l)
      	READ (3,'((i4,1x,24(1x,e20.8)))') p, l_e(l)
   enddo
!
   CLOSE(1)
   CLOSE(2)
   CLOSE(3) 
!
!
!###################################################
! Building the "E" arrays (Elastic Green functions)
!###################################################
!
! Initialization of the "E" arrays 
   ES(:)=0.   
   EU(:)=0.
   EN(:)=0.
   EV(:)=0.
!
! --------------------------------------------------------------
! ====== Setting the ELASTIC Green Functions for Degree 1 ======
! --------------------------------------------------------------
!
 IF(DEG1==1) THEN 				
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree nd=1 are set by the corresponding Love numbers 
!
   nd=1
!
! ---- GSC or Elastic solution
	If(imode==1.or.imode==2.or.imode==5)  then 
	EU(nd) = h_e(nd) 
        EN(nd) = 0. 
        EV(nd) = l_e(nd)
	ES(nd) = EN(nd) - EU(nd)	
	Endif
!
! ---- Eustatic or solution 
	If(imode==3)  then 
		EU(nd) = 0.
		EN(nd) = 0. 
		EV(nd) = 0. 	
	        ES(nd) = 0.
	Endif 
!
! ---- Woodward solution 
	If(imode==4)  then 
		EU(nd) = 0.
		EN(nd) = 1. 
		EV(nd) = 0. 	
	        ES(nd) = 0.
	Endif 	
!
   den=2.*float(nd)+1. 	
!
   EU(nd) = EU(nd)/den 
   EN(nd) = EN(nd)/den 
   EV(nd) = EV(nd)/den 
   ES(nd) = ES(nd)/den 	
!
 ENDIF 
!
!
! ----------------------------------------------------------------
! ====== Setting the ELASTIC Green Functions for Degree > 1 ======
! ----------------------------------------------------------------
! 
 Do 10 l=LMP1, LMAX
!
	EN(l) =   1.+k_e(l)
	EU(l) =      h_e(l)    
	EV(l) =      l_e(l) 
      	ES(l) = EN(l)-EU(l)    
!       
! ---- Eustatic solution 
	if(imode==3) then 
			EN(l) = 0.
			EU(l) = 0.
			EV(l) = 0.
	                ES(l) = EN(l)-EU(l) 
		     endif 
!
! ---- Woodward solution			
	if(imode==4) then 
			EN(l) = 1.
			EU(l) = 0.
			EV(l) = 0.
	                ES(l) = EN(l)-EU(l) 
		     endif
!
        den=	(2.*float(l)+1.) 
!
        EU(l) = EU(l)/den 
	EN(l) = EN(l)/den 
	EV(l) = EV(l)/den 
	ES(l) = ES(l)/den 
!
10 Continue
!
!
!
!############################################################
! Building the "BETA" arrays (Visco-Elastic Green functions)
!############################################################
!
!
 BETAS(:,:)=0. 
 BETAU(:,:)=0. 
 BETAN(:,:)=0.
 BETAV(:,:)=0.
!
!
! --------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree 1 ======
! --------------------------------------------------------------------
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree nd=1 are set by the corresponding Love numbers 
!
 IF(DEG1==1) THEN 
!
	do k=0, NN 
!
	l=1 
!
                BETAN(l,k)=0. 
                BETAU(l,k)=0. 
                BETAV(l,k)=0. 
!
   		den = 2.*float(l)+1. 
!
        		do m=1, nroots				    
! 
			BETAU(l,k) = BETAU(l,k) + & 
			            (r_h(l,m))* & 
			            (1. - exp(-k*delta/temp(l,m)))/den
!
			BETAV(l,k) = BETAV(l,k) + & 
			            (r_l(l,m))* & 
			            (1. - exp(-k*delta/temp(l,m)))/den
!				    
			enddo
!			
	BETAS(l,k) = BETAN(l,k) - BETAU(l,k)		
!			
	enddo
!	
 ENDIF
!
!
! ----------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree > 1 ======
! ----------------------------------------------------------------------
!
	do k=0, NN 
!
        do l=lmp1, lmax 
!
                BETAU(l,k)=0.0 
                BETAN(l,k)=0.0 
                BETAV(l,k)=0.0 
!
   	        den = 2.*float(l)+1.  
!
        		do m=1, nroots
!
			BETAN(l,k) = BETAN(l,k) + & 
			            (r_k(l,m))* & 
			            (1.0 - exp(-k*delta/temp(l,m)))/ den 
! 
			BETAU(l,k) = BETAU(l,k) + & 
			            (r_h(l,m))* & 
			            (1.0 - exp(-k*delta/temp(l,m)))/ den 
!
			BETAV(l,k) = BETAV(l,k) + & 
			            (r_l(l,m))* & 
			            (1.0 - exp(-k*delta/temp(l,m)))/ den 
!				    
			enddo
!			
	      BETAS(l,k)=BETAN(l,k) - BETAU(l,k)   			
!			
	      enddo
	enddo


!
! ----------------------------------- 
! ====== Reporting the results ======
! ----------------------------------- 
!
!
 open(1,file='ebs.dat',status='unknown')
 open(2,file='ebu.dat',status='unknown')
 open(3,file='ebn.dat',status='unknown')
 open(4,file='ebv.dat',status='unknown')
!  
 do l=0, lmax
 	 write(1,*) l, ES(l), (betas(l,k), k=0, nn)
 	 write(2,*) l, EU(l), (betau(l,k), k=0, nn)
 	 write(3,*) l, EN(l), (betan(l,k), k=0, nn)		 
 	 write(4,*) l, EV(l), (betav(l,k), k=0, nn)
 enddo
!
 close(1)
 close(2)
 close(3)
 close(4) 
!
!
!
!================
  end program BB 
!================
!
!
