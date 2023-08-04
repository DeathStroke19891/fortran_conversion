
c This code calculates the dc conductivity,  density of states, 
c and the site occupancy on the lattice(whether a site is singly 
c or doubly occupied)
c It takes the eigen vector and eigen value as the input along with 
c the band hopping parameters, the chemical potential, the temperature
c and the energy window (tolerance) within which the delta function
c is evaluated in the conductivity term

        PROGRAM Main
	INTEGER, PARAMETER :: L=30
	INTEGER, PARAMETER :: N_z=1
	INTEGER, PARAMETER :: N=(L**2)*N_z
	INTEGER, PARAMETER :: Ntemp=15
	INTEGER, PARAMETER :: totdospts=400
	REAL, PARAMETER    :: U= 0.0
	REAL, PARAMETER    :: tx=1.0
	REAL, PARAMETER    :: ty=1.0
	REAL, PARAMETER    :: tz=1.0

	integer :: LDA, i, j, k, m, itemp, Ndos(totdospts) 
	integer :: idos, icounter 
        complex :: Z( N, N )
	real    :: W( N )
        real    :: temperature(Ntemp), Temp
	real    :: Ef, dcsigma1,dcsigma2,dcsigma3,dcsigma4
	real    ::  dcsigma5,dcsigma6,dcsigma7,dcsigma8,dcsigma9
	real    :: el,eu,em1,em2,bandwidth
	real    :: sigx1, sigx2, sigx3, sigx4, sigx5, 
     1              sigx6, sigx7,sigx8,sigx9,sigx10,om_ref, u1    
	real    :: ebracketdown, ebracketup, estep, gstateenergy  
	real    :: dos(totdospts),ener(totdospts)  
        EXTERNAL sigdc_T

        open(1,file="eigenvec",status="old") 
        open(2,file="energy",status="old") 
        open(3,file="dcsigma",status="unknown") 
        open(4,file="dos-ener",status="unknown") 

        gstateenergy = 0.0

	do i = 1,N 	 
          do j = 1,N 	 
!         i = site index , j = state index
            read(1,*) k, m, Z(i,j) 
          enddo
!            read(2,*)k, W(i) 
            read(2,*)u1,k, W(i) 
            if (i .LE. (N/2)) then
            gstateenergy = gstateenergy + W(i)
            endif
	enddo

        If (u1.NE.U) THEN
         WRITE(*,*) "Error in U, Program Terminated!!"
         STOP
        Endif 

            el = W(1)
            eu = W(N)
            bandwidth = eu - el
            estep = bandwidth/totdospts
            em1 = W(N/2)
            em2 = W((N/2) + 1)
            Ef = (em1 + em2)/2.0

            om_ref = 0.56   !  
!            om_ref = 0.009   ! N = 18*18*9 = 2916, bandwidth ~ 12t, BW/N = 0.004
!************************************************************************************
!            om_ref = 0.05  ! Taken from Pinaki's paper 
                           ! changed to 0.05 on 27th June 2007 
!************************************************************************************

        do icounter = 1,totdospts
           Ndos(icounter) = 0
           ener(icounter) = 0.0
        enddo    

        do idos = 1, totdospts 
         ebracketdown = el + (idos - 1)*estep
         ebracketup = el + idos*estep
         ener(idos) = el + (idos - 0.5)*estep 
          do i = 1,N
          if((ebracketdown.LE.W(i)).and.(ebracketup.GT.W(i))) then
          Ndos(idos) = Ndos(idos) + 1
          endif
          enddo
        enddo
         write(4,*)  gstateenergy/N
        do icounter = 1,totdospts
        dos(icounter) = REAL(Ndos(icounter))/N
        write(4,*) icounter,bandwidth,ener(icounter),Ndos(icounter),
     $   dos(icounter)
        enddo

!    Define the temperature array here
        temperature(1) = .0005 
        temperature(2) = .001 
        temperature(3) = .005 
        temperature(4) = .007 
        temperature(5) = .01 
        temperature(6) = .03 
        temperature(7) = .05 
        temperature(8) = .07 
        temperature(9) = .09 
        temperature(10) = .1 
        temperature(11) = .2 
        temperature(12) = .4 
        temperature(13) = .6 
        temperature(14) = .8 
        temperature(15) = 1.0 

!  Run the loop over temperature 

        write(3,*) N, U, Ef, bandwidth, Temp, om_ref 

!        do itemp = 1,Ntemp
!        do itemp = 1,12
        do itemp = 1,1
        Temp = temperature(itemp)

        write(*,*) "Calling dcsigma subroutine"
	call sigdc_T(Z,W,tx,ty,tz,sigx1,sigx2,sigx3,sigx4,sigx5,
     &     sigx6,sigx7,sigx8,sigx9,sigx10,
     &     dcsigma1,dcsigma2,dcsigma3,dcsigma4,dcsigma5,dcsigma6,
     &     dcsigma7, dcsigma8, dcsigma9,om_ref,Ef,Temp)

!        write(3,*) Temp, sigx2,sigx3,dcsigma1, dcsigma2
        write(3,*) om_ref,dcsigma1
        write(3,*) 2*om_ref,dcsigma2
        write(3,*) 3*om_ref,dcsigma3
        write(3,*) 4*om_ref,dcsigma4
        write(3,*) 5*om_ref,dcsigma5
        write(3,*) 6*om_ref,dcsigma6
        write(3,*) 7*om_ref,dcsigma7
        write(3,*) 8*om_ref,dcsigma8
        write(3,*) 9*om_ref,dcsigma9

        enddo

        END PROGRAM Main 

!**************************************************************

        
! Code for low freq. optical conductivity. 
! For L_system > L_(mean free path) it gives
! a good estimate for dc-conductivity, ideally one needs a '1/N -> 0'
! scaling to get the dc-conductivity.


! INPUTS:  wavefunctions: (Z) , eigenvalues: (W), hopping parameters: 
! (tx, ty, tz) ( although only tx is used here )
! chemical pot. : (Ef), Temperature: (Temp)

! OUTPUTS: \sigma(\omega) at two low \omega values: (dcsigma1,dcsigma2)
!
! This subroutine can be generalised to compute full \sigma(\omega) also.


      subroutine sigdc_T(Z,W,tx,ty,tz,sigmax1,sigmax2,sigmax3,sigmax4,   $
     $   sigmax5,sigmax6,sigmax7,sigmax8,sigmax9,sigmax10,dcsigma1,      $
     $   dcsigma2,dcsigma3,dcsigma4, dcsigma5,dcsigma6,dcsigma7,         $
     $   dcsigma8,dcsigma9,om_ref,Ef,Temp)

        implicit none

        integer nx,ny,nz

        parameter (nx = 22)
        parameter (ny = 22)
        parameter (nz = 5)

        integer N,N_plane,N_line,num(100),mm,nn,ncount,kkk,i,j,k
        integer icounter
        real W(nx*ny*nz)
        real Temp,avsig(1000),Ef
        real dcsigma1,dcsigma2,dcsigma3,dcsigma4
        real dcsigma5,dcsigma6,dcsigma7,dcsigma8,dcsigma9
        real sigmax1,sigmax2,sigmax3,sigmax4,sigmax5
        real sigmax6,sigmax7,sigmax8,sigmax9,sigmax10
        real tx,ty,tz,omega,f1,f2,fermi, x, y
        real Pi,dsigz(10),om(10),om_ref
        complex Z(nx*ny*nz,nx*ny*nz)
        complex sum1,sum2,sum3,sum4,sum5,sum6
*        EXTERNAL FUNCTION fermi
*        EXTERNAL fermi

       open(8,file="countercheck",status="unknown") 

       N = nx*ny*nz
       N_plane = nx*ny
       N_line = nx
       Pi = acos(-1.)

!  ************************************************************************
!  COMPUTATION OF MARTIX ELEMENTS OF CURRENT OPERATOR    !
!  BETWEEN DIFFERENT PAIRS OF STATES OF THE HAMILTONIAN.      !
!*************************************************************************

!       om_ref = 0.01*(1000./real(N))           !  om_ref ~ bandwidth/N

       om(1) = om_ref*1.
       om(2) = om_ref*2.
       om(3) = om_ref*3.
       om(4) = om_ref*4.
       om(5) = om_ref*5.
       om(6) = om_ref*6.
       om(7) = om_ref*7.
       om(8) = om_ref*8.
       om(9) = om_ref*9.
       om(10) = om_ref*10.


       do kkk = 1,10
          omega = om(kkk)
       dsigz(kkk) = 0.0

       icounter = 0

       do mm = 1,N-1
       do nn = mm+1,N      ! look for all possible states with "dE = w"

            sum1 = (0.,0.)
            sum2 = (0.,0.)
            sum3 = (0.,0.)
            sum4 = (0.,0.)
            sum5 = (0.,0.)
            sum6 = (0.,0.)

        if(abs(W(mm)-W(nn)).ge.10.**(-9.)) then
        if(abs(W(mm)-W(nn)).le.omega) then
       
        x = (W(mm)-Ef)/Temp
        y = (W(nn)-Ef)/Temp

        icounter = icounter + 1

!        f1 = fermi(x)
!        f2 = fermi(y)
        f1 = 1.0/(exp(x) + 1.0) 
        f2 = 1.0/(exp(y) + 1.0) 

! j_x
!          do i = 1,nx-1
!          do j = 1,ny
!          do k = 1,nz
! j_z
          do k = 1,nz-1
          do j = 1,ny
          do i = 1,nx

! expression for j_x

!       sum1 = sum1 + (conjg(Z((k-1)*N_plane+(j-1)*N_line+i,mm))*
!     &  Z((k-1)*N_plane+(j-1)*N_line+i+1,nn)) - 
!     &    (conjg(Z((k-1)*N_plane+(j-1)*N_line+i+1,mm))*
!     &  Z((k-1)*N_plane+(j-1)*N_line+i,nn))

! expression for j_z

         sum2 = conjg(Z((k-1)*N_plane+(j-1)*N_line+i,mm))
         sum3 = Z(k*N_plane+(j-1)*N_line+i,nn)
         sum4 = conjg(Z(k*N_plane+(j-1)*N_line+i,mm))
         sum5 = Z((k-1)*N_plane+(j-1)*N_line+i,nn)
         sum6 = (sum2*sum3) - (sum4*sum5) 
         sum1 = sum1 + sum6

!       sum1 = sum1 + (conjg(Z((k-1)*N_plane+(j-1)*N_line+i,mm))*
!     &  Z(k*N_plane+(j-1)*N_line+i,nn)) - 
!     &    (conjg(Z(k*N_plane+(j-1)*N_line+i,mm))*
!     &  Z((k-1)*N_plane+(j-1)*N_line+i,nn))

!       if(((cabs(sum1)**2).GT.0.0001).AND.             
!     &  (((f1-f2)/(W(nn)-W(mm))).GT.0))THEN 
!       write(*,*) mm,nn,i,j,k,cabs(sum1)**2,
!     &  (cabs(sum1)**2.)*(f1-f2)/(W(nn)-W(mm))
!       write(*,*) (f1-f2)/(W(nn)-W(mm))
!        endif
         
!        IF(((i.EQ.nx).AND.(j.EQ.ny).AND.(k.EQ.nz-1))
!     &   .AND.(cabs(sum1)**2.).GT.0.0001)THEN
!        write(*,*) (cabs(sum1)**2.) 
!        ENDIF
!        IF(REAL(sum6).GT.0.001) THEN
!         write(*,*) REAL(sum1), REAL(sum6)
!        ENDIF

          enddo
          enddo
          enddo               !  <J_mn>
!       write(*,*) mm,nn,i,j,k,(cabs(sum1)**2.)


        dsigz(kkk) = dsigz(kkk) + (cabs(sum1)**2.)*(f1-f2)/(W(nn)-W(mm))

!        if(((f1-f2)/(W(nn)-W(mm))).LT. 0.0) then
!        write(*,*) "warning: coeff is -ve", dsigz(kkk)
!        STOP
!        endif 


!        if((cabs(sum1)**2.).GT.0.0000001) then
!        if(((f1-f2)/(W(nn)-W(mm))).GT.0.01) then
!        write(*,*) (f1-f2)/(W(nn)-W(mm)),(cabs(sum1)**2.)
!        write(*,*) (f1-f2)/(W(nn)-W(mm))*(cabs(sum1)**2.)
!        endif 

        endif
        endif


        enddo
        enddo
        write(8,*) Temp, omega, icounter 
        dsigz(kkk) = tz**2.*dsigz(kkk)/real(N)
        write(*,*) dsigz(kkk)
        enddo

        sigmax1 = dsigz(1)
        sigmax2 = dsigz(2)
        sigmax3 = dsigz(3)
        sigmax4 = dsigz(4)
        sigmax5 = dsigz(5)
        sigmax6 = dsigz(6)
        sigmax7 = dsigz(7)
        sigmax8 = dsigz(8)
        sigmax9 = dsigz(9)
        sigmax10 = dsigz(10)

!************************************************************
c    computing sigma(om) using forward differences.
c------------------------------------------------------------
       dcsigma1 = (dsigz(2)-dsigz(1))/(om(2)-om(1))

       dcsigma2 = (dsigz(3)-dsigz(2))/(om(3)-om(2))

       dcsigma3 = (dsigz(4)-dsigz(3))/(om(4)-om(3))

       dcsigma4 = (dsigz(5)-dsigz(4))/(om(5)-om(4))

       dcsigma5 = (dsigz(6)-dsigz(5))/(om(6)-om(5))

       dcsigma6 = (dsigz(7)-dsigz(6))/(om(7)-om(6))

       dcsigma7 = (dsigz(8)-dsigz(7))/(om(8)-om(7))

       dcsigma8 = (dsigz(9)-dsigz(8))/(om(9)-om(8))

       dcsigma9 = (dsigz(10)-dsigz(9))/(om(10)-om(9))

        write(*,*) dcsigma1, dcsigma2, dcsigma3, dcsigma4,   
     $  dcsigma5, dcsigma6, dcsigma7,dcsigma8,dcsigma9

          return
          end

!***************************************************************
c FUNCTION FERMI DEFINED HERE
!****************************************************************
!      REAL FUNCTION fermi( x )
!      IMPLICIT NONE
!      REAL   x
!       fermi(x) = 1.0/(exp(x) + 1.0)
!      return
!      END 
!       END FUNCTION fermi

c******************************************************************
