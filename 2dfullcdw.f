c      Three dimensional Hubbard model:
c        metal-barrier-metal configuration        
        integer LDA,N,LRWORK,LIWORK,L,m
        integer LWORK,INFO
        parameter (L=30,LDA=(L**2)+1, N=(L**2),Mfib=2584)
        parameter (LRWORK=1 + 5*N + 2*N**2)
        parameter (LIWORK=5*N+3+2,LWORK=2*N+N**2+1)
        parameter (ust=1,uf=5.1) 
        parameter (est=0.0,ef=0.0) 
        complex A( LDA, N )
        complex WORK( LWORK )
        integer IWORK( LIWORK )
        real W( N ),Wup(N),Wdn(N)
        real RWORK( LRWORK )
        character*1  JOBZ,UPLO 
	dimension anup(6000),andn(6000)
        dimension ary1(2584),ary2(2584),ary3(2584) 
	dimension site(2584),anup1(6000)
	dimension andn1(6000),cnumber(6000)
	dimension anti(6000)
	open(4,file='profiler21',status='unknown')
	open(1,file='energy21',status='unknown')
        open(2,file='chsp21.dat',status='unknown')
        open(50,file='eigenvec',status='unknown')
        open(8,file='gapcdw',status='unknown')
        JOBZ='V'
        UPLO='U'
	u1=0.0     !initally it was 2.39
c	write(*,*)"u1",u1
        anfil=nfil
        checkup=0.
        checkdn=0.
	prec=0.000001

c	INITIALIZATION OF NUP & NDOWN
c----------------------------------------------------------------------
38	do 10 i=1,L-1,2
	anup(i)=0.+0.3
 	anup(i+1)=0.+0.7  
 	andn(i)=0.+0.7
 	andn(i+1)=0.+0.3
	anup(i+l*l)=(1-anup(i))
 	anup(i+1+l*l)=(1-anup(i+1))
 	andn(i+l*l)=(1-andn(i))
 	andn(i+1+l*l)=(1-andn(i+1))

        do 44 j=1,(L-2)/2
        k=i+2*j*L
        write(*,*) k
        anup(k)=anup(i)
        anup(k+1)=anup(i+1)
        andn(k)=andn(i)
        andn(k+1)=andn(i+1)
        anup(k+l*l)=(1-anup(i))
        anup(k+1+l*l)=(1-anup(i+1))
        andn(k+l*l)=(1-andn(i))
        andn(k+1+l*l)=(1-andn(i+1))
44      continue
10      continue
c--------------------------------------------------------------------        
	do 17 i=L+1,2*L-1,2
 	anup(i)=0.+0.7
 	anup(i+1)=0.+0.3
 	andn(i)=0.+0.3
 	andn(i+1)=0.+0.7
 	anup(i+l*l)=(1-anup(i))
 	anup(i+1+l*l)=(1-anup(i+1))
 	andn(i+l*l)=(1-andn(i))
 	andn(i+1+l*l)=(1-andn(i+1))
        do 19 j=1,(L-2)/2
        k=i+2*j*L
        anup(k)=anup(i)
        anup(k+1)=anup(i+1)
        andn(k)=andn(i)
        andn(k+1)=andn(i+1)
        anup(k+l*l)=(1-anup(i))
        anup(k+1+l*l)=(1-anup(i+1))
        andn(k+l*l)=(1-andn(i))
        andn(k+1+l*l)=(1-andn(i+1))
19      continue
17      continue
c-------------------------------------------------------------------
138	do 510 i=1,L-1,2
	site(i)=0.0   ! initially it was 1.0
 	site(i+1)=0.0  

        do 404 j=1,(L-2)/2
        k=i+2*j*L
        write(*,*) k
        site(k)=site(i)
        site(k+1)=site(i+1)
        site(k+l*l)=site(i+1)
        site(k+1+l*l)=site(i)
404      continue
510      continue
c--------------------------------------------------------------------        
	do 117 i=L+1,2*L-1,2
 	site(i)=-0.0
 	site(i+1)=0.0   ! initiall it was 1.0
        do 109 j=1,(L-2)/2
        k=i+2*j*L
        site(k)=site(i)
        site(k+1)=site(i+1)
        site(k+l*l)=site(i+1)
        site(k+1+l*l)=site(i)
109      continue
117      continue
c-----------------------------------------------------------------------
c----------------------------------------------------------------------
        do 43 i=1,L*L
43      write(*,*)"anup", i,site(i) 
c	MATRIX PREPARATION & DIAGONALISATION FOR UP
25        call amatp(LDA,N,L,A,u1,site,andn)
       call CHEEVD( JOBZ, UPLO, N, A, LDA, Wup, WORK, LWORK,
     $                  RWORK, LRWORK, IWORK, LIWORK, INFO )
 	call order(LDA,N,L,A,anup1,annup)
         write(*,*) 'INFO= ',info
         do 31 k=1,N
31       continue
         sum1=0
         do 30 k=1,N/2
         sum1=sum1+Wup(k)
30       continue
         sum1=sum1/N
         write(*,*) sum1,INFO  
c	MATRIX PREPARATION & DIAGONALISATION FOR DOWN
        call amatp(LDA,N,L,A,u1,site,anup)
       call CHEEVD( JOBZ, UPLO, N, A, LDA, Wdn, WORK, LWORK,
     $                  RWORK, LRWORK, IWORK, LIWORK, INFO )
	call order(LDA,N,L,A,andn1,anndn)
         write(*,*) 'INFO= ',info
         sum2=0
         do 64 k=1,N/2
         sum2=sum2+Wdn(k)
64       continue
         sum2=sum2/N
	nhope=nhope+1
c	CHECK FOR CONVERGENCE	
        checkup1=checkup
        checkdn1=checkdn
	checkup=0.
	checkdn=0.
	do 32 i=1,L*L
 	checkup=checkup+(anup(i)-anup1(i))*(anup(i)-anup1(i))
 	checkdn=checkdn+(andn(i)-andn1(i))*(andn(i)-andn1(i))
c	checkup=checkup+(anup(i)-anup1(i))
c	checkdn=checkdn+(andn(i)-andn1(i))
32	continue
       write(*,*)'ok1',checkup,checkdn,nhope,u1
         write(*,*) "check", abs(checkup-checkup1),
     1   abs(checkdn-checkdn1)
c       if(abs(checkup).lt.prec.and.abs(checkdn).lt.prec) then
c	go to 100
c	else
         if(abs(checkup-checkup1).lt.prec.and.abs(checkdn-checkdn1)
     1  .lt.prec) go to 100
c	if (nhope.gt.nhopeless) go to 110
	do 40 i=1,L*L
	anup(i)=anup1(i)
	andn(i)=andn1(i)
40	continue
        go to 25
c	endif
100     continue
	do 137 i=1, l*l
        anti(i)=abs(anup(i)-andn(i))
        cnumber(i)=anup(i)+andn(i)
        write(4,*) i,anup(i),andn(i),anti(i),cnumber(i)
137     continue
c       write(*,*) fill(1),fill(2),f,amu
	spin=0.
	charge=0.
        bnup=0
        bndn=0
c	write(*,*)"ITER ",nhope,u1
c-----------------------------------------------------------------------
	do 50 i=L*L+1,L*L+L
	spin=spin+((-1)**i)*(anup(i)-andn(i))
	charge=charge+((-1)**i)*(anup(i)+andn(i)-1)
c	write(*,*)i ,anup(i)," ",andn(i)
        bnup=bnup+anup(i)
        bndn=bndn+andn(i)
        do 42 j=1,(L-2)/2
        k=i+2*j*L
	spin=spin+((-1)**k)*(anup(k)-andn(k))
	charge=charge+((-1)**k)*(anup(k)+andn(k)-1)
42      continue
        
50	continue
c--------------------------------------------------------------------------
	do 47 i=L*L+L+1,L*L+2*L
	spin=spin+((-1)**(i+1))*(anup(i)-andn(i))
	charge=charge+((-1)**(i+1))*(anup(i)+andn(i)-1)
c	write(*,*)i ,anup(i)," ",andn(i)
        bnup=bnup+anup(i)
        bndn=bndn+andn(i)
        do 41 j=1,(L-2)/2
        k=i+2*j*L
	spin=spin+((-1)**(k+1))*(anup(k)-andn(k))
	charge=charge+((-1)**(k+1))*(anup(k)+andn(k)-1)
41      continue
47	continue
	spin=abs(spin/L*L)
	charge=abs(charge/L*L)
c--------------------------------------------------------------------------
	shift=0.
c110	write(*,*) "ITERATION FAILS"
        write(2,*) eb,u1,charge,spin,
     2   W(((L*L)/2)+1)-W((L*L)/2)
c110	write(*,*) "ITERATION FAILS"
c120	write(*,*) "ok"
         do 134 k=1,N
         write(1,*) u1,k, Wup(k),Wdn(K)
134       continue
         sum=sum1+sum2
         write(1,*)sum,sum1,sum2 ,Wup((N/2)+1)-Wup(N/2),  
     1    Wdn((N/2)+1)-Wdn(N/2)  
         write(8,*) u1,eb,l,W((N/2)+1)-W(N/2)  
c        eb=eb+est
c  	u1=u1+ust
c        if(u1.gt.uf) then
c        go to 77
c        else
c        go to 38
c        endif
        nhope=0   ! might need to be commented out later
        
      do i=1,N
        do j=1,N
         write(50,*) i,j,A(i,j)
        enddo
      enddo


77	stop
	end
c **********************************************************
c		SUBROUTINE ORDER
c******************************************************	
c*****  Calculates Nup & Ndown
	subroutine order(LDA,N,L,A,an,ann)
	dimension an(6000)
        integer LDA,N,L 
        complex A(LDA,N)
c	common /tag1/ ea,eb,u1,nfil
c	DETERMINATION OF MU
	do 7 i=1,L*L
	an(i)=0.
        ank1=0
	do 9 j=1,(L*L)/2
9	an(i)=an(i)+A(i,j)*conjg(A(i,j))
c	write(*,*)"No.s ",an(i)
7	continue
        ann=0
        do 34 i=1,L*L
 34     ann=ann+an(i)
c       write(*,*) '***',ann,amu
	return
	end
c*******************************************************
c******************************************************
c		SUBROUTINE AMATP
c **********************************************************
c*****  MATRIX prepare kora hoy
        subroutine amatp(LDA,N,L,A,u1,site,an)
        integer LDA,N,L
	complex A(LDA,N)
        real  u1,an(6000),site(6000)
c	common /tag1/ ea,eb,u1,n

	do 5 i=1,L*L
	do 5 j=1,L*L
	A(i,j)=(0.,0)
5	continue

	 do 10 j=1,L
	 do 10 i=1,L-1
         k=i+(j-1)*L
         A(k+1,k)=(-1,0)
	 A(k,k+1)=(-1,0)
10	continue

	do 48 j=1,L-1
	do 48 i=1,L
         k=i+(j-1)*L
	A(k,k+L)=(-1.0,0) 
	A(k+L,k)=(-1.0,0)
48      continue
        
        do 68 k=1,L*L
 	A(k,k)=site(k)+u1*an(k) 
68      continue

        do 37 i=1,L*L
        do 37 j=1,L*L
c	write(2,*) i,j,A(i,j)
37      continue
	return
        end
c---------------------------------------------------------------
c ***********************************************************
c		SUBROUTINE FIBOGEN
c ************************************************************
c****** Fibonacci Generator
	subroutine fibogen(Mfib,ary3,eb)
c	implicit real*8(a-h,o-z)
        integer M
	dimension ary1(2584),ary2(2584),ary3(2584)
c	common /tag1/ ea,eb,u1,nfil
	nfibo1=1
	nfibo2=2
	ary1(1)=ea
	ary2(1)=ea
	ary2(2)=eb
3	nfibo3=nfibo1+nfibo2
	do 5 i=1,nfibo2
5	ary3(i)=ary2(i)	
	do 7 i=1,nfibo1
7	ary3(i+nfibo2)=ary1(i)
	do 9 i=1,nfibo2
9	ary1(i)=ary2(i)	
	do 11 i=1,nfibo3
11	ary2(i)=ary3(i)	
	nfibo1=nfibo2
	nfibo2=nfibo3
	if (nfibo3.lt.Mfib) then
	go to 3
	else
	endif
c	write(*,*) (ary3(i),i=1,Mfib)
	return
	end
c******************************************************
      SUBROUTINE CHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                   LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               RWORK( * ), W( * )
      COMPLEX            A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CHEEVD computes all eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix A.  If eigenvectors are desired, it uses a
*  divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) REAL array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.
*          If N <= 1,                LWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
*          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) REAL array,
*                                         dimension (LRWORK)
*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*
*  LRWORK  (input) INTEGER
*          The dimension of the array RWORK.
*          If N <= 1,                LRWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
*          If JOBZ  = 'V' and N > 1, LRWORK must be at least
*                         1 + 5*N + 2*N**2.
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the RWORK array,
*          returns this value as the first entry of the RWORK array, and
*          no error message related to LRWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2,
     $                   INDWRK, ISCALE, LIOPT, LIWMIN, LLRWK, LLWORK,
     $                   LLWRK2, LOPT, LROPT, LRWMIN, LWMIN
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANHE, SLAMCH
      EXTERNAL           LSAME, CLANHE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHETRD, CLACPY, CLASCL, CSTEDC, CUNMTR, SSCAL,
     $                   SSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( N.LE.1 ) THEN
         LWMIN = 1
         LRWMIN = 1
         LIWMIN = 1
         LOPT = LWMIN
         LROPT = LRWMIN
         LIOPT = LIWMIN
      ELSE
         IF( WANTZ ) THEN
            LWMIN = 2*N + N*N
            LRWMIN = 1 + 5*N + 2*N**2
            LIWMIN = 3 + 5*N
         ELSE
            LWMIN = N + 1
            LRWMIN = N
            LIWMIN = 1
         END IF
         LOPT = LWMIN
         LROPT = LRWMIN
         LIOPT = LIWMIN
      END IF
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -8
      ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LOPT
         RWORK( 1 ) = LROPT
         IWORK( 1 ) = LIOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHEEVD', -INFO )
         RETURN 
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         IF( WANTZ )
     $      A( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = CLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL CLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call CHETRD to reduce Hermitian matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      INDRWK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWORK = LWORK - INDWRK + 1
      LLWRK2 = LWORK - INDWK2 + 1
      LLRWK = LRWORK - INDRWK + 1
      CALL CHETRD( UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
      LOPT = MAX( REAL( LOPT ), REAL( N )+REAL( WORK( INDWRK ) ) )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, first call
*     CSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
*     tridiagonal matrix, then call CUNMTR to multiply it to the
*     Householder transformations represented as Householder vectors in
*     A.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL CSTEDC( 'I', N, W, RWORK( INDE ), WORK( INDWRK ), N,
     $                WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK,
     $                IWORK, LIWORK, INFO )
         CALL CUNMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ),
     $                WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
         CALL CLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
         LOPT = MAX( LOPT, N+N**2+INT( WORK( INDWK2 ) ) )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 ) = LOPT
      RWORK( 1 ) = LROPT
      IWORK( 1 ) = LIOPT
*
      RETURN
*
*     End of CHEEVD
*
      END
      SUBROUTINE CHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  CHETD2 reduces a complex Hermitian matrix A to real symmetric
*  tridiagonal form T by a unitary similarity transformation:
*  Q' * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the unitary
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the unitary matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) REAL array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) REAL array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) COMPLEX array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   HALF = ( 0.5E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      COMPLEX            ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CHEMV, CHER2, CLARFG, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHETD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         A( N, N ) = REAL( A( N, N ) )
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            ALPHA = A( I, I+1 )
            CALL CLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
               CALL CHEMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*CDOTC( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL CAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL CHER2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
            ELSE
               A( I, I ) = REAL( A( I, I ) )
            END IF
            A( I, I+1 ) = E( I )
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        Reduce the lower triangle of A
*
         A( 1, 1 ) = REAL( A( 1, 1 ) )
         DO 20 I = 1, N - 1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            ALPHA = A( I+1, I )
            CALL CLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
               CALL CHEMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*CDOTC( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL CAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL CHER2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
            ELSE
               A( I+1, I+1 ) = REAL( A( I+1, I+1 ) )
            END IF
            A( I+1, I ) = E( I )
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     End of CHETD2
*
      END
      SUBROUTINE CHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CHETRD reduces a complex Hermitian matrix A to real symmetric
*  tridiagonal form T by a unitary similarity transformation:
*  Q**H * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the unitary
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the unitary matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) REAL array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) REAL array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) COMPLEX array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHER2K, CHETD2, CLATRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.
*
         NB = ILAENV( 1, 'CHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHETRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         NX = MAX( NB, ILAENV( 3, 'CHETRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'CHETRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL CLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
*
*           Update the unreduced submatrix A(1:i-1,1:i-1), using an
*           update of the form:  A := A - V*W' - W*V'
*
            CALL CHER2K( UPLO, 'No transpose', I-1, NB, -CONE,
     $                   A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA )
*
*           Copy superdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL CHETD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NX, NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL CLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
*           an update of the form:  A := A - V*W' - W*V'
*
            CALL CHER2K( UPLO, 'No transpose', N-I-NB+1, NB, -CONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
*
*           Copy subdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL CHETD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of CHETRD
*
      END
      SUBROUTINE CLACGV( N, X, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  CLACGV conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) COMPLEX array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )
     $      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
*
*     End of CLACGV
*
      END
      SUBROUTINE CLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
*          is accessed; if UPLO = 'L', only the lower trapezium is
*          accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) COMPLEX array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
*
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of CLACPY
*
      END
      SUBROUTINE CLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, M, N
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), RWORK( * )
      COMPLEX            A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  CLACRM performs a very simple matrix-matrix multiplication:
*           C := A * B,
*  where A is M by N and complex; B is N by N and real;
*  C is M by N and complex.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A and of the matrix C.
*          M >= 0.
*
*  N       (input) INTEGER
*          The number of columns and rows of the matrix B and
*          the number of columns of the matrix C.
*          N >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA, N)
*          A contains the M by N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >=max(1,M).
*
*  B       (input) REAL array, dimension (LDB, N)
*          B contains the N by N matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >=max(1,N).
*
*  C       (input) COMPLEX array, dimension (LDC, N)
*          C contains the M by N matrix C.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >=max(1,N).
*
*  RWORK   (workspace) REAL array, dimension (2*M*N)
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, REAL
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = REAL( A( I, J ) )
   10    CONTINUE
   20 CONTINUE
*
      L = M*N + 1
      CALL SGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO,
     $            RWORK( L ), M )
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE
*
      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = AIMAG( A( I, J ) )
   50    CONTINUE
   60 CONTINUE
      CALL SGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO,
     $            RWORK( L ), M )
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = CMPLX( REAL( C( I, J ) ),
     $                  RWORK( L+( J-1 )*M+I-1 ) )
   70    CONTINUE
   80 CONTINUE
*
      RETURN
*
*     End of CLACRM
*
      END
      COMPLEX FUNCTION CLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      COMPLEX            X, Y
*     ..
*
*  Purpose
*  =======
*
*  CLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX
*  Y       (input) COMPLEX
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
      CALL SLADIV( REAL( X ), AIMAG( X ), REAL( Y ), AIMAG( Y ), ZR,
     $             ZI )
      CLADIV = CMPLX( ZR, ZI )
*
      RETURN
*
*     End of CLADIV
*
      END
      SUBROUTINE CLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK,
     $                   IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDQ, LDQS, N, QSIZ
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            Q( LDQ, * ), QSTORE( LDQS, * )
*     ..
*
*  Purpose
*  =======
*
*  Using the divide and conquer method, CLAED0 computes all eigenvalues
*  of a symmetric tridiagonal matrix which is one diagonal block of
*  those from reducing a dense or band Hermitian matrix and
*  corresponding eigenvectors of the dense or band matrix.
*
*  Arguments
*  =========
*
*  QSIZ   (input) INTEGER
*         The dimension of the unitary matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, the diagonal elements of the tridiagonal matrix.
*         On exit, the eigenvalues in ascending order.
*
*  E      (input/output) REAL array, dimension (N-1)
*         On entry, the off-diagonal elements of the tridiagonal matrix.
*         On exit, E has been destroyed.
*
*  Q      (input/output) COMPLEX array, dimension (LDQ,N)
*         On entry, Q must contain an QSIZ x N matrix whose columns
*         unitarily orthonormal. It is a part of the unitary matrix
*         that reduces the full dense Hermitian matrix to a
*         (reducible) symmetric tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  IWORK  (workspace) INTEGER array,
*         the dimension of IWORK must be at least
*                      6 + 6*N + 5*N*lg N
*                      ( lg( N ) = smallest integer k
*                                  such that 2^k >= N )
*
*  RWORK  (workspace) REAL array,
*                               dimension (1 + 3*N + 2*N*lg N + 3*N**2)
*                        ( lg( N ) = smallest integer k
*                                    such that 2^k >= N )
*
*  QSTORE (workspace) COMPLEX array, dimension (LDQS, N)
*         Used to store parts of
*         the eigenvector matrix when the updating matrix multiplies
*         take place.
*
*  LDQS   (input) INTEGER
*         The leading dimension of the array QSTORE.
*         LDQS >= max(1,N).
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  =====================================================================
*
*  Warning:      N could be as big as QSIZ!
*
*     .. Parameters ..
      REAL               TWO
      PARAMETER          ( TWO = 2.E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM,
     $                   IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM,
     $                   J, K, LGN, LL, MATSIZ, MSD2, SMLSIZ, SMM1,
     $                   SPM1, SPM2, SUBMAT, SUBPBS, TLVLS
      REAL               TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CLACRM, CLAED7, SCOPY, SSTEQR, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     IF( ICOMPQ .LT. 0 .OR. ICOMPQ .GT. 2 ) THEN
*        INFO = -1
*     ELSE IF( ( ICOMPQ .EQ. 1 ) .AND. ( QSIZ .LT. MAX( 0, N ) ) )
*    $        THEN
      IF( QSIZ.LT.MAX( 0, N ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLAED0', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      SMLSIZ = ILAENV( 9, 'CLAED0', ' ', 0, 0, 0, 0 )
*
*     Determine the size and placement of the submatrices, and save in
*     the leading elements of IWORK.
*
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
*
*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
*     using rank-1 modifications (cuts).
*
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
*
      INDXQ = 4*N + 3
*
*     Set up workspaces for eigenvalues only/accumulate new vectors
*     routine
*
      TEMP = LOG( REAL( N ) ) / LOG( TWO )
      LGN = INT( TEMP )
      IF( 2**LGN.LT.N )
     $   LGN = LGN + 1
      IF( 2**LGN.LT.N )
     $   LGN = LGN + 1
      IPRMPT = INDXQ + N + 1
      IPERM = IPRMPT + N*LGN
      IQPTR = IPERM + N*LGN
      IGIVPT = IQPTR + N + 2
      IGIVCL = IGIVPT + N*LGN
*
      IGIVNM = 1
      IQ = IGIVNM + 2*N*LGN
      IWREM = IQ + N**2 + 1
*     Initialize pointers
      DO 50 I = 0, SUBPBS
         IWORK( IPRMPT+I ) = 1
         IWORK( IGIVPT+I ) = 1
   50 CONTINUE
      IWORK( IQPTR ) = 1
*
*     Solve each submatrix eigenproblem at the bottom of the divide and
*     conquer tree.
*
      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         LL = IQ - 1 + IWORK( IQPTR+CURR )
         CALL SSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                RWORK( LL ), MATSIZ, RWORK, INFO )
         CALL CLACRM( QSIZ, MATSIZ, Q( 1, SUBMAT ), LDQ, RWORK( LL ),
     $                MATSIZ, QSTORE( 1, SUBMAT ), LDQS,
     $                RWORK( IWREM ) )
         IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
         CURR = CURR + 1
         IF( INFO.GT.0 ) THEN
            INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
            RETURN
         END IF
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
*
*     Successively merge eigensystems of adjacent submatrices
*     into eigensystem for the corresponding larger matrix.
*
*     while ( SUBPBS > 1 )
*
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
*
*     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
*     into an eigensystem of size MATSIZ.  CLAED7 handles the case
*     when the eigenvectors of a full or band Hermitian matrix (which
*     was reduced to tridiagonal form) are desired.
*
*     I am free to use Q as a valuable working space until Loop 150.
*
            CALL CLAED7( MATSIZ, MSD2, QSIZ, TLVLS, CURLVL, CURPRB,
     $                   D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS,
     $                   E( SUBMAT+MSD2-1 ), IWORK( INDXQ+SUBMAT ),
     $                   RWORK( IQ ), IWORK( IQPTR ), IWORK( IPRMPT ),
     $                   IWORK( IPERM ), IWORK( IGIVPT ),
     $                   IWORK( IGIVCL ), RWORK( IGIVNM ),
     $                   Q( 1, SUBMAT ), RWORK( IWREM ),
     $                   IWORK( SUBPBS+1 ), INFO )
            IF( INFO.GT.0 ) THEN
               INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
               RETURN
            END IF
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
*
*     end while
*
*     Re-merge the eigenvalues/vectors which were deflated at the final
*     merge step.
*
      DO 100 I = 1, N
         J = IWORK( INDXQ+I )
         RWORK( I ) = D( J )
         CALL CCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100 CONTINUE
      CALL SCOPY( N, RWORK, 1, D, 1 )
*
      RETURN
*
*     End of CLAED0
*
      END
      SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
     $                   LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM,
     $                   GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ,
     $                   TLVLS
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
     $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      REAL               D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * )
      COMPLEX            Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CLAED7 computes the updated eigensystem of a diagonal
*  matrix after modification by a rank-one symmetric matrix. This
*  routine is used only for the eigenproblem which requires all
*  eigenvalues and optionally eigenvectors of a dense or banded
*  Hermitian matrix that has been reduced to tridiagonal form.
*
*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)
*
*    where Z = Q'u, u is a vector of length N with ones in the
*    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
*
*     The eigenvectors of the original matrix are stored in Q, and the
*     eigenvalues are in D.  The algorithm consists of three stages:
*
*        The first stage consists of deflating the size of the problem
*        when there are multiple eigenvalues or if there is a zero in
*        the Z vector.  For each such occurence the dimension of the
*        secular equation problem is reduced by one.  This stage is
*        performed by the routine SLAED2.
*
*        The second stage consists of calculating the updated
*        eigenvalues. This is done by finding the roots of the secular
*        equation via the routine SLAED4 (as called by SLAED3).
*        This routine also calculates the eigenvectors of the current
*        problem.
*
*        The final stage consists of computing the updated eigenvectors
*        directly using the updated eigenvalues.  The eigenvectors for
*        the current problem are multiplied with the eigenvectors from
*        the overall problem.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  CUTPNT (input) INTEGER
*         Contains the location of the last eigenvalue in the leading
*         sub-matrix.  min(1,N) <= CUTPNT <= N.
*
*  QSIZ   (input) INTEGER
*         The dimension of the unitary matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N.
*
*  TLVLS  (input) INTEGER
*         The total number of merging levels in the overall divide and
*         conquer tree.
*
*  CURLVL (input) INTEGER
*         The current level in the overall merge routine,
*         0 <= curlvl <= tlvls.
*
*  CURPBM (input) INTEGER
*         The current problem in the current level in the overall
*         merge routine (counting from upper left to lower right).
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, the eigenvalues of the rank-1-perturbed matrix.
*         On exit, the eigenvalues of the repaired matrix.
*
*  Q      (input/output) COMPLEX array, dimension (LDQ,N)
*         On entry, the eigenvectors of the rank-1-perturbed matrix.
*         On exit, the eigenvectors of the repaired tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  RHO    (input) REAL
*         Contains the subdiagonal element used to create the rank-1
*         modification.
*
*  INDXQ  (output) INTEGER array, dimension (N)
*         This contains the permutation which will reintegrate the
*         subproblem just solved back into sorted order,
*         ie. D( INDXQ( I = 1, N ) ) will be in ascending order.
*
*  IWORK  (workspace) INTEGER array, dimension (4*N)
*
*  RWORK  (workspace) REAL array,
*                                 dimension (3*N+2*QSIZ*N)
*
*  WORK   (workspace) COMPLEX array, dimension (QSIZ*N)
*
*  QSTORE (input/output) REAL array, dimension (N**2+1)
*         Stores eigenvectors of submatrices encountered during
*         divide and conquer, packed together. QPTR points to
*         beginning of the submatrices.
*
*  QPTR   (input/output) INTEGER array, dimension (N+2)
*         List of indices pointing to beginning of submatrices stored
*         in QSTORE. The submatrices are numbered starting at the
*         bottom left of the divide and conquer tree, from left to
*         right and bottom to top.
*
*  PRMPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in PERM a
*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
*         indicates the size of the permutation and also the size of
*         the full, non-deflated problem.
*
*  PERM   (input) INTEGER array, dimension (N lg N)
*         Contains the permutations (from deflation and sorting) to be
*         applied to each eigenblock.
*
*  GIVPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in GIVCOL a
*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
*         indicates the number of Givens rotations.
*
*  GIVCOL (input) INTEGER array, dimension (2, N lg N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (input) REAL array, dimension (2, N lg N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            COLTYP, CURR, I, IDLMDA, IND1, IND2, INDX,
     $                   INDXC, INDXP, IQ, IW, IZ, K, N1, N2, PTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACRM, CLAED8, SLAED9, SLAEDA, SLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
*        INFO = -1
*     ELSE IF( N.LT.0 ) THEN
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -2
      ELSE IF( QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLAED7', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are for bookkeeping purposes only.  They are
*     integer pointers which indicate the portion of the workspace
*     used by a particular array in SLAED2 and SLAED3.
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ = IW + N
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
     $             GIVCOL, GIVNUM, QSTORE, QPTR, RWORK( IZ ),
     $             RWORK( IZ+N ), INFO )
*
*     When solving the final problem, we no longer need the stored data,
*     so we will overwrite the data from this level onto the previously
*     used storage space.
*
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
*
*     Sort and Deflate eigenvalues.
*
      CALL CLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, RWORK( IZ ),
     $             RWORK( IDLMDA ), WORK, QSIZ, RWORK( IW ),
     $             IWORK( INDXP ), IWORK( INDX ), INDXQ,
     $             PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ),
     $             GIVCOL( 1, GIVPTR( CURR ) ),
     $             GIVNUM( 1, GIVPTR( CURR ) ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         CALL SLAED9( K, 1, K, N, D, RWORK( IQ ), K, RHO,
     $                RWORK( IDLMDA ), RWORK( IW ),
     $                QSTORE( QPTR( CURR ) ), K, INFO )
         CALL CLACRM( QSIZ, K, WORK, QSIZ, QSTORE( QPTR( CURR ) ), K, Q,
     $                LDQ, RWORK( IQ ) )
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*
*     Prepare the INDXQ sorting premutation.
*
         N1 = K
         N2 = N - K
         IND1 = 1
         IND2 = N
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of CLAED7
*
      END
      SUBROUTINE CLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA,
     $                   Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      REAL               D( * ), DLAMDA( * ), GIVNUM( 2, * ), W( * ),
     $                   Z( * )
      COMPLEX            Q( LDQ, * ), Q2( LDQ2, * )
*     ..
*
*  Purpose
*  =======
*
*  CLAED8 merges the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny element in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  K      (output) INTEGER
*         Contains the number of non-deflated eigenvalues.
*         This is the order of the related secular equation.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  QSIZ   (input) INTEGER
*         The dimension of the unitary matrix used to reduce
*         the dense or band matrix to tridiagonal form.
*         QSIZ >= N if ICOMPQ = 1.
*
*  Q      (input/output) COMPLEX array, dimension (LDQ,N)
*         On entry, Q contains the eigenvectors of the partially solved
*         system which has been previously updated in matrix
*         multiplies with other partially solved eigensystems.
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max( 1, N ).
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.  On exit, D contains the trailing (N-K) updated
*         eigenvalues (those which were deflated) sorted into increasing
*         order.
*
*  RHO    (input/output) REAL
*         Contains the off diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined. RHO is modified during the computation to
*         the value required by SLAED3.
*
*  CUTPNT (input) INTEGER
*         Contains the location of the last eigenvalue in the leading
*         sub-matrix.  MIN(1,N) <= CUTPNT <= N.
*
*  Z      (input) REAL array, dimension (N)
*         On input this vector contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).  The contents of Z are
*         destroyed during the updating process.
*
*  DLAMDA (output) REAL array, dimension (N)
*         Contains a copy of the first K eigenvalues which will be used
*         by SLAED3 to form the secular equation.
*
*  Q2     (output) COMPLEX array, dimension (LDQ2,N)
*         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
*         Contains a copy of the first K eigenvectors which will be used
*         by SLAED7 in a matrix multiply (SGEMM) to update the new
*         eigenvectors.
*
*  LDQ2   (input) INTEGER
*         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ).
*
*  W      (output) REAL array, dimension (N)
*         This will hold the first k values of the final
*         deflation-altered z-vector and will be passed to SLAED3.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         This will contain the permutation used to place deflated
*         values of D at the end of the array. On output INDXP(1:K)
*         points to the nondeflated D-values and INDXP(K+1:N)
*         points to the deflated eigenvalues.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         This will contain the permutation used to sort the contents of
*         D into ascending order.
*
*  INDXQ  (input) INTEGER array, dimension (N)
*         This contains the permutation which separately sorts the two
*         sub-problems in D into ascending order.  Note that elements in
*         the second half of this permutation must first have CUTPNT
*         added to their values in order to be accurate.
*
*  PERM   (output) INTEGER array, dimension (N)
*         Contains the permutations (from deflation and sorting) to be
*         applied to each eigenblock.
*
*  GIVPTR (output) INTEGER
*         Contains the number of Givens rotations which took place in
*         this subproblem.
*
*  GIVCOL (output) INTEGER array, dimension (2, N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (output) REAL array, dimension (2, N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0E0, ZERO = 0.0E0, ONE = 1.0E0,
     $                   TWO = 2.0E0, EIGHT = 8.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      REAL               C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SLAPY2
      EXTERNAL           ISAMAX, SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CLACPY, CSROT, SCOPY, SLAMRG, SSCAL,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -8
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLAED8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL SSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1
*
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL SSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL SLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
*
*     Calculate the allowable deflation tolerance
*
      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = SLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     -- except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         DO 50 J = 1, N
            PERM( J ) = INDXQ( INDX( J ) )
            CALL CCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   50    CONTINUE
         CALL CLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ), LDQ )
         RETURN
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      K = 0
      GIVPTR = 0
      K2 = N + 1
      DO 60 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N )
     $         GO TO 100
         ELSE
            JLAM = J
            GO TO 70
         END IF
   60 CONTINUE
   70 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 90
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( JLAM )
         C = Z( J )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = SLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( J ) = TAU
            Z( JLAM ) = ZERO
*
*           Record the appropriate Givens rotation
*
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            CALL CSROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                  Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   80       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 80
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 70
   90 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  100 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      DO 110 J = 1, N
         JP = INDXP( J )
         DLAMDA( J ) = D( JP )
         PERM( J ) = INDXQ( INDX( JP ) )
         CALL CCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  110 CONTINUE
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      IF( K.LT.N ) THEN
         CALL SCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
         CALL CLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2, Q( 1, K+1 ),
     $                LDQ )
      END IF
*
      RETURN
*
*     End of CLAED8
*
      END
      REAL             FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CLANHE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  complex hermitian matrix A.
*
*  Description
*  ===========
*
*  CLANHE returns the value
*
*     CLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in CLANHE as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          hermitian matrix A is to be referenced.
*          = 'U':  Upper triangular part of A is referenced
*          = 'L':  Lower triangular part of A is referenced
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, CLANHE is
*          set to zero.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The hermitian matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of A contains the upper triangular part
*          of the matrix A, and the strictly lower triangular part of A
*          is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of A contains the lower triangular part of
*          the matrix A, and the strictly upper triangular part of A is
*          not referenced. Note that the imaginary parts of the diagonal
*          elements need not be set and are assumed to be zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) REAL array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
               VALUE = MAX( VALUE, ABS( REAL( A( J, J ) ) ) )
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( REAL( A( J, J ) ) ) )
               DO 30 I = J + 1, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is hermitian).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( REAL( A( J, J ) ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( REAL( A( J, J ) ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL CLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL CLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         DO 130 I = 1, N
            IF( REAL( A( I, I ) ).NE.ZERO ) THEN
               ABSA = ABS( REAL( A( I, I ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      CLANHE = VALUE
      RETURN
*
*     End of CLANHE
*
      END
      SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX            TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CLARF applies a complex elementary reflector H to a complex M-by-N
*  matrix C, from either the left or the right. H is represented in the
*  form
*
*        H = I - tau * v * v'
*
*  where tau is a complex scalar and v is a complex vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
*  tau.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) COMPLEX array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) COMPLEX
*          The value tau in the representation of H.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CGERC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL CGEMV( 'Conjugate transpose', M, N, ONE, C, LDC, V,
     $                  INCV, ZERO, WORK, 1 )
*
*           C := C - v * w'
*
            CALL CGERC( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL CGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL CGERC( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of CLARF
*
      END
      SUBROUTINE CLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  CLARFB applies a complex block reflector H or its transpose H' to a
*  complex M-by-N matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'C': apply H' (Conjugate transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) COMPLEX array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) COMPLEX array, dimension (LDT,K)
*          The triangular K-by-K matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEMM, CLACGV, CTRMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL CCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL CLACGV( N, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL CGEMM( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C( K+1, 1 ), LDC,
     $                        V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL CTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK,
     $                        LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - CONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL CCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL CGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL CTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ),
     $                        LDV, ONE, C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL CCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL CLACGV( N, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL CGEMM( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V, LDV, WORK, LDWORK,
     $                        ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               CONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL CCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL CGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE,
     $                        C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL CCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL CLACGV( N, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL CGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL CTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL CGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - CONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL CCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C( 1, K+1 ), LDC,
     $                        V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL CTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL CGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL CCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL CLACGV( N, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL CGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE, C,
     $                        LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL CGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE, V,
     $                        LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               CONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL CCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL CGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of CLARFB
*
      END
      SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX            ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  CLARFG generates a complex elementary reflector H of order n, such
*  that
*
*        H' * ( alpha ) = ( beta ),   H' * H = I.
*             (   x   )   (   0  )
*
*  where alpha and beta are scalars, with beta real, and x is an
*  (n-1)-element complex vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a complex scalar and v is a complex (n-1)-element
*  vector. Note that H is not hermitian.
*
*  If the elements of x are all zero and alpha is real, then tau = 0
*  and H is taken to be the unit matrix.
*
*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) COMPLEX
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) COMPLEX array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) COMPLEX
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      REAL               SCNRM2, SLAMCH, SLAPY3
      COMPLEX            CLADIV
      EXTERNAL           SCNRM2, SLAMCH, SLAPY3, CLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSCAL, CSSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = SCNRM2( N-1, X, INCX )
      ALPHR = REAL( ALPHA )
      ALPHI = AIMAG( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
*
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL CSSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = SCNRM2( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI )
            BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = CLADIV( CMPLX( ONE ), ALPHA-BETA )
            CALL CSCAL( N-1, ALPHA, X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = CLADIV( CMPLX( ONE ), ALPHA-BETA )
            CALL CSCAL( N-1, ALPHA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of CLARFG
*
      END
      SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  CLARFT forms the triangular factor T of a complex block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) COMPLEX array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) COMPLEX array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) COMPLEX array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX            VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CLACGV, CTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
*
                  CALL CGEMV( 'Conjugate transpose', N-I+1, I-1,
     $                        -TAU( I ), V( I, 1 ), LDV, V( I, I ), 1,
     $                        ZERO, T( 1, I ), 1 )
               ELSE
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
*
                  IF( I.LT.N )
     $               CALL CLACGV( N-I, V( I, I+1 ), LDV )
                  CALL CGEMV( 'No transpose', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
                  IF( I.LT.N )
     $               CALL CLACGV( N-I, V( I, I+1 ), LDV )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL CTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
*
                     CALL CGEMV( 'Conjugate transpose', N-K+I, K-I,
     $                           -TAU( I ), V( 1, I+1 ), LDV, V( 1, I ),
     $                           1, ZERO, T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
*
                     CALL CLACGV( N-K+I-1, V( I, 1 ), LDV )
                     CALL CGEMV( 'No transpose', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     CALL CLACGV( N-K+I-1, V( I, 1 ), LDV )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL CTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of CLARFT
*
      END
      SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL               CFROM, CTO
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CLASCL multiplies the M by N complex matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) REAL
*  CTO     (input) REAL
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,M)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of CLASCL
*
      END
      SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX            ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CLASET initializes a 2-D array A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set. The lower triangle
*                      is unchanged.
*          = 'L':      Lower triangular part is set. The upper triangle
*                      is unchanged.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          On entry, M specifies the number of rows of A.
*
*  N       (input) INTEGER
*          On entry, N specifies the number of columns of A.
*
*  ALPHA   (input) COMPLEX
*          All the offdiagonal array elements are set to ALPHA.
*
*  BETA    (input) COMPLEX
*          All the diagonal array elements are set to BETA.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
*                   A(i,i) = BETA , 1 <= i <= min(m,n)
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of CLASET
*
      END
      SUBROUTINE CLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL               C( * ), S( * )
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n complex matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) REAL arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               CTEMP, STEMP
      COMPLEX            TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of CLASR
*
      END
      SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  CLASSQ returns the values scl and ssq such that
*
*     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
*  assumed to be at least unity and the value of ssq will then satisfy
*
*     1.0 .le. ssq .le. ( sumsq + 2*n ).
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*            i
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  SCALE and SUMSQ are overwritten by scl and ssq respectively.
*
*  The routine makes only one pass through the vector X.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) COMPLEX array, dimension (N)
*          The vector x as described above.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) REAL
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with the value  scl .
*
*  SUMSQ   (input/output) REAL
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with the value  ssq .
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      REAL               TEMP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, REAL
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( REAL( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( REAL( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( AIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( AIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of CLASSQ
*
      END
      SUBROUTINE CLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      REAL               E( * )
      COMPLEX            A( LDA, * ), TAU( * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  CLATRD reduces NB rows and columns of a complex Hermitian matrix A to
*  Hermitian tridiagonal form by a unitary similarity
*  transformation Q' * A * Q, and returns the matrices V and W which are
*  needed to apply the transformation to the unreduced part of A.
*
*  If UPLO = 'U', CLATRD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', CLATRD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  This is an auxiliary routine called by CHETRD.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NB      (input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit:
*          if UPLO = 'U', the last NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements above the diagonal
*            with the array TAU, represent the unitary matrix Q as a
*            product of elementary reflectors;
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements below the diagonal
*            with the array TAU, represent the  unitary matrix Q as a
*            product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  E       (output) REAL array, dimension (N-1)
*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*          elements of the last NB columns of the reduced matrix;
*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
*          the first NB columns of the reduced matrix.
*
*  TAU     (output) COMPLEX array, dimension (N-1)
*          The scalar factors of the elementary reflectors, stored in
*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*          See Further Details.
*
*  W       (output) COMPLEX array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n) H(n-1) . . . H(n-nb+1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
*  and tau in TAU(i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
*  and tau in TAU(i).
*
*  The elements of the vectors v together form the n-by-nb matrix V
*  which is needed, with W, to apply the transformation to the unreduced
*  part of the matrix, using a Hermitian rank-2k update of the form:
*  A := A - V*W' - W*V'.
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5 and nb = 2:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  a   a   a   v4  v5 )              (  d                  )
*    (      a   a   v4  v5 )              (  1   d              )
*    (          a   1   v5 )              (  v1  1   a          )
*    (              d   1  )              (  v1  v2  a   a      )
*    (                  d  )              (  v1  v2  a   a   a  )
*
*  where d denotes a diagonal element of the reduced matrix, a denotes
*  an element of the original matrix that is unchanged, and vi denotes
*  an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE, HALF
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   ONE = ( 1.0E+0, 0.0E+0 ),
     $                   HALF = ( 0.5E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IW
      COMPLEX            ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CGEMV, CHEMV, CLACGV, CLARFG, CSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              Update A(1:i,i)
*
               A( I, I ) = REAL( A( I, I ) )
               CALL CLACGV( N-I, W( I, IW+1 ), LDW )
               CALL CGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL CLACGV( N-I, W( I, IW+1 ), LDW )
               CALL CLACGV( N-I, A( I, I+1 ), LDA )
               CALL CGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL CLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = REAL( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               ALPHA = A( I-1, I )
               CALL CLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = ALPHA
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL CHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL CGEMV( 'Conjugate transpose', I-1, N-I, ONE,
     $                        W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO,
     $                        W( I+1, IW ), 1 )
                  CALL CGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL CGEMV( 'Conjugate transpose', I-1, N-I, ONE,
     $                        A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO,
     $                        W( I+1, IW ), 1 )
                  CALL CGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL CSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*CDOTC( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL CAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        Reduce first NB columns of lower triangle
*
         DO 20 I = 1, NB
*
*           Update A(i:n,i)
*
            A( I, I ) = REAL( A( I, I ) )
            CALL CLACGV( I-1, W( I, 1 ), LDW )
            CALL CGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL CLACGV( I-1, W( I, 1 ), LDW )
            CALL CLACGV( I-1, A( I, 1 ), LDA )
            CALL CGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL CLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = REAL( A( I, I ) )
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               ALPHA = A( I+1, I )
               CALL CLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*
               CALL CHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL CGEMV( 'Conjugate transpose', N-I, I-1, ONE,
     $                     W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO,
     $                     W( 1, I ), 1 )
               CALL CGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL CGEMV( 'Conjugate transpose', N-I, I-1, ONE,
     $                     A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO,
     $                     W( 1, I ), 1 )
               CALL CGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL CSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*CDOTC( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL CAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of CLATRD
*
      END
      SUBROUTINE CSROT( N, CX, INCX, CY, INCY, C, S )
*
*     applies a plane rotation, where the cos and sin (c and s) are real
*     and the vectors cx and cy are complex.
*     jack dongarra, linpack, 3/11/78.
*
*  =====================================================================
*
*     .. Scalar Arguments ..
      INTEGER           INCX, INCY, N
      REAL              C, S
*     ..
*     .. Array Arguments ..
      COMPLEX           CX( * ), CY( * )
*     ..
*     .. Local Scalars ..
      INTEGER           I, IX, IY
      COMPLEX           CTEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*        code for unequal increments or equal increments not equal
*          to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         CTEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - S*CX( IX )
         CX( IX ) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
   20 DO 30 I = 1, N
         CTEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - S*CX( I )
         CX( I ) = CTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE CSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK,
     $                   LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CSTEDC computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*  The eigenvectors of a full or band complex Hermitian matrix can also
*  be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See SLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original Hermitian matrix
*                  also.  On entry, Z contains the unitary matrix used
*                  to reduce the original matrix to tridiagonal form.
*
*  N       (input) INTEGER
*          The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ,N)
*          On entry, if COMPZ = 'V', then Z contains the unitary
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original Hermitian matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If  COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If eigenvectors are desired, then LDZ >= max(1,N).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.
*          If COMPZ = 'V' and N > 1, LWORK must be at least N*N.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) REAL array,
*                                         dimension (LRWORK)
*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*
*  LRWORK  (input) INTEGER
*          The dimension of the array RWORK.
*          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.
*          If COMPZ = 'V' and N > 1, LRWORK must be at least
*                         1 + 3*N + 2*N*lg N + 3*N**2 ,
*                         where lg( N ) = smallest integer k such
*                         that 2**k >= N.
*          If COMPZ = 'I' and N > 1, LRWORK must be at least
*                         1 + 4*N + 2*N**2 .
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the RWORK array,
*          returns this value as the first entry of the RWORK array, and
*          no error message related to LRWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.
*          If COMPZ = 'V' or N > 1,  LIWORK must be at least
*                                    6 + 6*N + 5*N*lg N.
*          If COMPZ = 'I' or N > 1,  LIWORK must be at least
*                                    3 + 5*N .
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            END, I, ICOMPZ, II, J, K, LGN, LIWMIN, LL,
     $                   LRWMIN, LWMIN, M, SMLSIZ, START
      REAL               EPS, ORGNRM, P, TINY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANST
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY, CLACRM, CLAED0, CSTEQR, CSWAP, SLASCL,
     $                   SLASET, SSTEDC, SSTEQR, SSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( N.LE.1 .OR. ICOMPZ.LE.0 ) THEN
         LWMIN = 1
         LIWMIN = 1
         LRWMIN = 1
      ELSE
         LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( ICOMPZ.EQ.1 ) THEN
            LWMIN = N*N
            LRWMIN = 1 + 3*N + 2*N*LGN + 3*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1
            LRWMIN = 1 + 4*N + 2*N**2
            LIWMIN = 3 + 5*N
         END IF
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -8
      ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN 
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN 
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 )
     $      Z( 1, 1 ) = ONE
         RETURN 
      END IF
*
      SMLSIZ = ILAENV( 9, 'CSTEDC', ' ', 0, 0, 0, 0 )
*
*     If the following conditional clause is removed, then the routine
*     will use the Divide and Conquer routine to compute only the
*     eigenvalues, which requires (3N + 3N**2) real workspace and
*     (2 + 5N + 2N lg(N)) integer workspace.
*     Since on many architectures SSTERF is much faster than any other
*     algorithm for finding eigenvalues only, it is used here
*     as the default.
*
*     If COMPZ = 'N', use SSTERF to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL SSTERF( N, D, E, INFO )
         RETURN 
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPZ.EQ.0 ) THEN
            CALL SSTERF( N, D, E, INFO )
            RETURN 
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            CALL CSTEQR( 'I', N, D, E, Z, LDZ, RWORK, INFO )
            RETURN 
         ELSE
            CALL CSTEQR( 'V', N, D, E, Z, LDZ, RWORK, INFO )
            RETURN 
         END IF
      END IF
*
*     If COMPZ = 'I', we simply call SSTEDC instead.
*
      IF( ICOMPZ.EQ.2 ) THEN
         CALL SLASET( 'Full', N, N, ZERO, ONE, RWORK, N )
         LL = N*N + 1
         CALL SSTEDC( 'I', N, D, E, RWORK, N, RWORK( LL ), LRWORK-LL+1,
     $                IWORK, LIWORK, INFO )
         DO 20 J = 1, N
            DO 10 I = 1, N
               Z( I, J ) = RWORK( ( J-1 )*N+I )
   10       CONTINUE
   20    CONTINUE
         RETURN 
      END IF
*
*     From now on, only option left to be handled is COMPZ = 'V',
*     i.e. ICOMPZ = 1.
*
*     Scale.
*
      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO )
     $   RETURN
*
      EPS = SLAMCH( 'Epsilon' )
*
      START = 1
*
*     while ( START <= N )
*
   30 CONTINUE
      IF( START.LE.N ) THEN
*
*     Let END be the position of the next subdiagonal entry such that
*     E( END ) <= TINY or END = N if no such subdiagonal exists.  The
*     matrix identified by the elements between START and END
*     constitutes an independent sub-problem.
*
         END = START
   40    CONTINUE
         IF( END.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( END ) ) )*SQRT( ABS( D( END+1 ) ) )
            IF( ABS( E( END ) ).GT.TINY ) THEN
               END = END + 1
               GO TO 40
            END IF
         END IF
*
*        (Sub) Problem determined.  Compute its size and solve it.
*
         M = END - START + 1
         IF( M.GT.SMLSIZ ) THEN
            INFO = SMLSIZ
*
*           Scale.
*
            ORGNRM = SLANST( 'M', M, D( START ), E( START ) )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M,
     $                   INFO )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ),
     $                   M-1, INFO )
*
            CALL CLAED0( N, M, D( START ), E( START ), Z( 1, START ),
     $                   LDZ, WORK, N, RWORK, IWORK, INFO )
            IF( INFO.GT.0 ) THEN
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) +
     $                MOD( INFO, ( M+1 ) ) + START - 1
               RETURN 
            END IF
*
*           Scale back.
*
            CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M,
     $                   INFO )
*
         ELSE
            CALL SSTEQR( 'I', M, D( START ), E( START ), RWORK, M,
     $                   RWORK( M*M+1 ), INFO )
            CALL CLACRM( N, M, Z( 1, START ), LDZ, RWORK, M, WORK, N,
     $                   RWORK( M*M+1 ) )
            CALL CLACPY( 'A', N, M, WORK, N, Z( 1, START ), LDZ )
            IF( INFO.GT.0 ) THEN
               INFO = START*( N+1 ) + END
               RETURN 
            END IF
         END IF
*
         START = END + 1
         GO TO 30
      END IF
*
*     endwhile
*
*     If the problem split any number of times, then the eigenvalues
*     will not be properly ordered.  Here we permute the eigenvalues
*     (and the associated eigenvectors) into ascending order.
*
      IF( M.NE.N ) THEN
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 60 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 50 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
   50       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL CSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
   60    CONTINUE
      END IF
*
      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of CSTEDC
*
      END
      SUBROUTINE CSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band complex Hermitian matrix can also
*  be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  Hermitian matrix.  On entry, Z must contain the
*                  unitary matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the unitary
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original Hermitian matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) REAL array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is unitarily similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   THREE = 3.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),
     $                   CONE = ( 1.0E0, 0.0E0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      REAL               ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANST, SLAPY2
      EXTERNAL           LSAME, SLAMCH, SLANST, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASET, CLASR, CSWAP, SLAE2, SLAEV2, SLARTG,
     $                   SLASCL, SLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = SLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = SLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL SLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL CLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL SLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL CLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL SLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL CLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL SLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL CLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL SLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL CSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
*
*     End of CSTEQR
*
      END
      SUBROUTINE CUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CUNM2L overwrites the general complex m-by-n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'C', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'C',
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k) . . . H(2) H(1)
*
*  as returned by CGEQLF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'C': apply Q' (Conjugate transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          CGEQLF in the last k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) COMPLEX array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by CGEQLF.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the m-by-n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      COMPLEX            AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNM2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. NOTRAN .OR. .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) or H(i)' is applied to C(1:m-k+i,1:n)
*
            MI = M - K + I
         ELSE
*
*           H(i) or H(i)' is applied to C(1:m,1:n-k+i)
*
            NI = N - K + I
         END IF
*
*        Apply H(i) or H(i)'
*
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = CONJG( TAU( I ) )
         END IF
         AII = A( NQ-K+I, I )
         A( NQ-K+I, I ) = ONE
         CALL CLARF( SIDE, MI, NI, A( 1, I ), 1, TAUI, C, LDC, WORK )
         A( NQ-K+I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of CUNM2L
*
      END
      SUBROUTINE CUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CUNM2R overwrites the general complex m-by-n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'C', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'C',
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by CGEQRF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'C': apply Q' (Conjugate transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          CGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) COMPLEX array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by CGEQRF.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the m-by-n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX            AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) or H(i)' is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) or H(i)' is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i) or H(i)'
*
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = CONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL CLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC,
     $               WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of CUNM2R
*
      END
      SUBROUTINE CUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CUNMQL overwrites the general complex M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'C':      Q**H * C       C * Q**H
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k) . . . H(2) H(1)
*
*  as returned by CGEQLF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**H from the Left;
*          = 'R': apply Q or Q**H from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'C':  Transpose, apply Q**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          CGEQLF in the last k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) COMPLEX array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by CGEQLF.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IINFO, IWS, LDWORK, LWKOPT,
     $                   MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
      COMPLEX            T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARFB, CLARFT, CUNM2L, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'CUNMQL', SIDE // TRANS, M, N, K,
     $             -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNMQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'CUNMQL', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL CUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i+ib-1) . . . H(i+1) H(i)
*
            CALL CLARFT( 'Backward', 'Columnwise', NQ-K+I+IB-1, IB,
     $                   A( 1, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(1:m-k+i+ib-1,1:n)
*
               MI = M - K + I + IB - 1
            ELSE
*
*              H or H' is applied to C(1:m,1:n-k+i+ib-1)
*
               NI = N - K + I + IB - 1
            END IF
*
*           Apply H or H'
*
            CALL CLARFB( SIDE, TRANS, 'Backward', 'Columnwise', MI, NI,
     $                   IB, A( 1, I ), LDA, T, LDT, C, LDC, WORK,
     $                   LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of CUNMQL
*
      END
      SUBROUTINE CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CUNMQR overwrites the general complex M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'C':      Q**H * C       C * Q**H
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by CGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**H from the Left;
*          = 'R': apply Q or Q**H from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'C':  Conjugate transpose, apply Q**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          CGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) COMPLEX array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by CGEQRF.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
      COMPLEX            T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARFB, CLARFT, CUNM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'CUNMQR', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'CUNMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL CUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL CLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL CLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,
     $                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of CUNMQR
*
      END
      SUBROUTINE CUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CUNMTR overwrites the general complex M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'C':      Q**H * C       C * Q**H
*
*  where Q is a complex unitary matrix of order nq, with nq = m if
*  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
*  nq-1 elementary reflectors, as returned by CHETRD:
*
*  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**H from the Left;
*          = 'R': apply Q or Q**H from the Right.
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangle of A contains elementary reflectors
*                 from CHETRD;
*          = 'L': Lower triangle of A contains elementary reflectors
*                 from CHETRD.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'C':  Conjugate transpose, apply Q**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  A       (input) COMPLEX array, dimension
*                               (LDA,M) if SIDE = 'L'
*                               (LDA,N) if SIDE = 'R'
*          The vectors which define the elementary reflectors, as
*          returned by CHETRD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
*
*  TAU     (input) COMPLEX array, dimension
*                               (M-1) if SIDE = 'L'
*                               (N-1) if SIDE = 'R'
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by CHETRD.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >=M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, UPPER
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           ILAENV, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CUNMQL, CUNMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) )
     $          THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'CUNMQL', SIDE // TRANS, M-1, N, M-1,
     $                      -1 )
            ELSE
               NB = ILAENV( 1, 'CUNMQL', SIDE // TRANS, M, N-1, N-1,
     $                      -1 )
            END IF
         ELSE
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, M-1, N, M-1,
     $                      -1 )
            ELSE
               NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, M, N-1, N-1,
     $                      -1 )
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNMTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( LEFT ) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to CHETRD with UPLO = 'U'
*
         CALL CUNMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C,
     $                LDC, WORK, LWORK, IINFO )
      ELSE
*
*        Q was determined by a call to CHETRD with UPLO = 'L'
*
         IF( LEFT ) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL CUNMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU,
     $                C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of CUNMTR
*
      END
      SUBROUTINE SLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  SLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) REAL array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of SLACPY
*
      END
      SUBROUTINE SLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
*     ..
*
*  Purpose
*  =======
*
*  SLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) REAL
*  B       (input) REAL
*  C       (input) REAL
*  D       (input) REAL
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) REAL
*  Q       (output) REAL
*          The scalars p and q in the above expression.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               E, F
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
*
      RETURN
*
*     End of SLADIV
*
      END
      SUBROUTINE SLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B, C, RT1, RT2
*     ..
*
*  Purpose
*  =======
*
*  SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, and RT2
*  is the eigenvalue of smaller absolute value.
*
*  Arguments
*  =========
*
*  A       (input) REAL
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) REAL
*          The (1,2) and (2,1) elements of the 2-by-2 matrix.
*
*  C       (input) REAL
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) REAL
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) REAL
*          The eigenvalue of smaller absolute value.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
*     ..
*     .. Local Scalars ..
      REAL               AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of SLAE2
*
      END
      SUBROUTINE SLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED0 computes all eigenvalues and corresponding eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*
*  Arguments
*  =========
*
*  ICOMPQ  (input) INTEGER
*          = 0:  Compute eigenvalues only.
*          = 1:  Compute eigenvectors of original dense symmetric matrix
*                also.  On entry, Q contains the orthogonal matrix used
*                to reduce the original matrix to tridiagonal form.
*          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
*                matrix.
*
*  QSIZ   (input) INTEGER
*         The dimension of the orthogonal matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, the main diagonal of the tridiagonal matrix.
*         On exit, its eigenvalues.
*
*  E      (input) REAL array, dimension (N-1)
*         The off-diagonal elements of the tridiagonal matrix.
*         On exit, E has been destroyed.
*
*  Q      (input/output) REAL array, dimension (LDQ, N)
*         On entry, Q must contain an N-by-N orthogonal matrix.
*         If ICOMPQ = 0    Q is not referenced.
*         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
*                          orthogonal matrix used to reduce the full
*                          matrix to tridiagonal form corresponding to
*                          the subset of the full matrix which is being
*                          decomposed at this time.
*         If ICOMPQ = 2    On entry, Q will be the identity matrix.
*                          On exit, Q contains the eigenvectors of the
*                          tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  If eigenvectors are
*         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
*
*  QSTORE (workspace) REAL array, dimension (LDQS, N)
*         Referenced only when ICOMPQ = 1.  Used to store parts of
*         the eigenvector matrix when the updating matrix multiplies
*         take place.
*
*  LDQS   (input) INTEGER
*         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
*         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
*
*  WORK   (workspace) REAL array,
*         If ICOMPQ = 0 or 1, the dimension of WORK must be at least
*                     1 + 3*N + 2*N*lg N + 2*N**2
*                     ( lg( N ) = smallest integer k
*                                 such that 2^k >= N )
*         If ICOMPQ = 2, the dimension of WORK must be at least
*                     4*N + N**2.
*
*  IWORK  (workspace) INTEGER array,
*         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
*                        6 + 6*N + 5*N*lg N.
*                        ( lg( N ) = smallest integer k
*                                    such that 2^k >= N )
*         If ICOMPQ = 2, the dimension of IWORK must be at least
*                        3 + 5*N.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.E0, ONE = 1.E0, TWO = 2.E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM,
     $                   IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM,
     $                   J, K, LGN, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1,
     $                   SPM2, SUBMAT, SUBPBS, TLVLS
      REAL               TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMM, SLACPY, SLAED1, SLAED7, SSTEQR,
     $                   XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.2 ) THEN
         INFO = -1
      ELSE IF( ( ICOMPQ.EQ.1 ) .AND. ( QSIZ.LT.MAX( 0, N ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED0', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      SMLSIZ = ILAENV( 9, 'SLAED0', ' ', 0, 0, 0, 0 )
*
*     Determine the size and placement of the submatrices, and save in
*     the leading elements of IWORK.
*
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
*
*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
*     using rank-1 modifications (cuts).
*
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
*
      INDXQ = 4*N + 3
      IF( ICOMPQ.NE.2 ) THEN
*
*        Set up workspaces for eigenvalues only/accumulate new vectors
*        routine
*
         TEMP = LOG( REAL( N ) ) / LOG( TWO )
         LGN = INT( TEMP )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IPRMPT = INDXQ + N + 1
         IPERM = IPRMPT + N*LGN
         IQPTR = IPERM + N*LGN
         IGIVPT = IQPTR + N + 2
         IGIVCL = IGIVPT + N*LGN
*
         IGIVNM = 1
         IQ = IGIVNM + 2*N*LGN
         IWREM = IQ + N**2 + 1
*
*        Initialize pointers
*
         DO 50 I = 0, SUBPBS
            IWORK( IPRMPT+I ) = 1
            IWORK( IGIVPT+I ) = 1
   50    CONTINUE
         IWORK( IQPTR ) = 1
      END IF
*
*     Solve each submatrix eigenproblem at the bottom of the divide and
*     conquer tree.
*
      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         IF( ICOMPQ.EQ.2 ) THEN
            CALL SSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                   Q( SUBMAT, SUBMAT ), LDQ, WORK, INFO )
            IF( INFO.NE.0 )
     $         GO TO 130
         ELSE
            CALL SSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                   WORK( IQ-1+IWORK( IQPTR+CURR ) ), MATSIZ, WORK,
     $                   INFO )
            IF( INFO.NE.0 )
     $         GO TO 130
            IF( ICOMPQ.EQ.1 ) THEN
               CALL SGEMM( 'N', 'N', QSIZ, MATSIZ, MATSIZ, ONE,
     $                     Q( 1, SUBMAT ), LDQ, WORK( IQ-1+IWORK( IQPTR+
     $                     CURR ) ), MATSIZ, ZERO, QSTORE( 1, SUBMAT ),
     $                     LDQS )
            END IF
            IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
            CURR = CURR + 1
         END IF
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
*
*     Successively merge eigensystems of adjacent submatrices
*     into eigensystem for the corresponding larger matrix.
*
*     while ( SUBPBS > 1 )
*
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
*
*     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
*     into an eigensystem of size MATSIZ.
*     SLAED1 is used only for the full eigensystem of a tridiagonal
*     matrix.
*     SLAED7 handles the cases in which eigenvalues only or eigenvalues
*     and eigenvectors of a full symmetric matrix (which was reduced to
*     tridiagonal form) are desired.
*
            IF( ICOMPQ.EQ.2 ) THEN
               CALL SLAED1( MATSIZ, D( SUBMAT ), Q( SUBMAT, SUBMAT ),
     $                      LDQ, IWORK( INDXQ+SUBMAT ),
     $                      E( SUBMAT+MSD2-1 ), MSD2, WORK,
     $                      IWORK( SUBPBS+1 ), INFO )
            ELSE
               CALL SLAED7( ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB,
     $                      D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS,
     $                      IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ),
     $                      MSD2, WORK( IQ ), IWORK( IQPTR ),
     $                      IWORK( IPRMPT ), IWORK( IPERM ),
     $                      IWORK( IGIVPT ), IWORK( IGIVCL ),
     $                      WORK( IGIVNM ), WORK( IWREM ),
     $                      IWORK( SUBPBS+1 ), INFO )
            END IF
            IF( INFO.NE.0 )
     $         GO TO 130
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
*
*     end while
*
*     Re-merge the eigenvalues/vectors which were deflated at the final
*     merge step.
*
      IF( ICOMPQ.EQ.1 ) THEN
         DO 100 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL SCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100    CONTINUE
         CALL SCOPY( N, WORK, 1, D, 1 )
      ELSE IF( ICOMPQ.EQ.2 ) THEN
         DO 110 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL SCOPY( N, Q( 1, J ), 1, WORK( N*I+1 ), 1 )
  110    CONTINUE
         CALL SCOPY( N, WORK, 1, D, 1 )
         CALL SLACPY( 'A', N, N, WORK( N+1 ), N, Q, LDQ )
      ELSE
         DO 120 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
  120    CONTINUE
         CALL SCOPY( N, WORK, 1, D, 1 )
      END IF
      GO TO 140
*
  130 CONTINUE
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
*
  140 CONTINUE
      RETURN
*
*     End of SLAED0
*
      END
      SUBROUTINE SLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, INFO, LDQ, N
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            INDXQ( * ), IWORK( * )
      REAL               D( * ), Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED1 computes the updated eigensystem of a diagonal
*  matrix after modification by a rank-one symmetric matrix.  This
*  routine is used only for the eigenproblem which requires all
*  eigenvalues and eigenvectors of a tridiagonal matrix.  SLAED7 handles
*  the case in which eigenvalues only or eigenvalues and eigenvectors
*  of a full symmetric matrix (which was reduced to tridiagonal form)
*  are desired.
*
*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)
*
*     where Z = Q'u, u is a vector of length N with ones in the
*     CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
*
*     The eigenvectors of the original matrix are stored in Q, and the
*     eigenvalues are in D.  The algorithm consists of three stages:
*
*        The first stage consists of deflating the size of the problem
*        when there are multiple eigenvalues or if there is a zero in
*        the Z vector.  For each such occurence the dimension of the
*        secular equation problem is reduced by one.  This stage is
*        performed by the routine SLAED2.
*
*        The second stage consists of calculating the updated
*        eigenvalues. This is done by finding the roots of the secular
*        equation via the routine SLAED4 (as called by SLAED3).
*        This routine also calculates the eigenvectors of the current
*        problem.
*
*        The final stage consists of computing the updated eigenvectors
*        directly using the updated eigenvalues.  The eigenvectors for
*        the current problem are multiplied with the eigenvectors from
*        the overall problem.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, the eigenvalues of the rank-1-perturbed matrix.
*         On exit, the eigenvalues of the repaired matrix.
*
*  Q      (input/output) REAL array, dimension (LDQ,N)
*         On entry, the eigenvectors of the rank-1-perturbed matrix.
*         On exit, the eigenvectors of the repaired tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input/output) INTEGER array, dimension (N)
*         On entry, the permutation which separately sorts the two
*         subproblems in D into ascending order.
*         On exit, the permutation which will reintegrate the
*         subproblems back into sorted order,
*         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
*
*  RHO    (input) REAL
*         The subdiagonal entry used to create the rank-1 modification.
*
*  CUTPNT (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) <= CUTPNT <= N/2.
*
*  WORK   (workspace) REAL array, dimension (4*N + N**2)
*
*  IWORK  (workspace) INTEGER array, dimension (4*N)
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            COLTYP, CPP1, I, IDLMDA, INDX, INDXC, INDXP,
     $                   IQ2, IS, IW, IZ, K, N1, N2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAED2, SLAED3, SLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( MIN( 1, N / 2 ).GT.CUTPNT .OR. ( N / 2 ).LT.CUTPNT ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED1', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are integer pointers which indicate
*     the portion of the workspace
*     used by a particular array in SLAED2 and SLAED3.
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      CALL SCOPY( CUTPNT, Q( CUTPNT, 1 ), LDQ, WORK( IZ ), 1 )
      CPP1 = CUTPNT + 1
      CALL SCOPY( N-CUTPNT, Q( CPP1, CPP1 ), LDQ, WORK( IZ+CUTPNT ), 1 )
*
*     Deflate eigenvalues.
*
      CALL SLAED2( K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK( IZ ),
     $             WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ),
     $             IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ),
     $             IWORK( COLTYP ), INFO )
*
      IF( INFO.NE.0 )
     $   GO TO 20
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT +
     $        ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2
         CALL SLAED3( K, N, CUTPNT, D, Q, LDQ, RHO, WORK( IDLMDA ),
     $                WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ),
     $                WORK( IW ), WORK( IS ), INFO )
         IF( INFO.NE.0 )
     $      GO TO 20
*
*     Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO 10 I = 1, N
            INDXQ( I ) = I
   10    CONTINUE
      END IF
*
   20 CONTINUE
      RETURN
*
*     End of SLAED1
*
      END
      SUBROUTINE SLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,
     $                   Q2, INDX, INDXC, INDXP, COLTYP, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQ, N, N1
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),
     $                   INDXQ( * )
      REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
     $                   W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED2 merges the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny entry in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  N1     (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) <= N1 <= N/2.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  Q      (input/output) REAL array, dimension (LDQ, N)
*         On entry, Q contains the eigenvectors of two submatrices in
*         the two square blocks with corners at (1,1), (N1,N1)
*         and (N1+1, N1+1), (N,N).
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input/output) INTEGER array, dimension (N)
*         The permutation which separately sorts the two sub-problems
*         in D into ascending order.  Note that elements in the second
*         half of this permutation must first have N1 added to their
*         values. Destroyed on exit.
*
*  RHO    (input/output) REAL
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         SLAED3.
*
*  Z      (input) REAL array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  DLAMDA (output) REAL array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         SLAED3 to form the secular equation.
*
*  W      (output) REAL array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to SLAED3.
*
*  Q2     (output) REAL array, dimension (N1**2+(N-N1)**2)
*         A copy of the first K eigenvectors which will be used by
*         SLAED3 in a matrix multiply (SGEMM) to solve for the new
*         eigenvectors.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of DLAMDA into
*         ascending order.
*
*  INDXC  (output) INTEGER array, dimension (N)
*         The permutation used to arrange the columns of the deflated
*         Q matrix into three groups:  the first group contains non-zero
*         elements only at and above N1, the second contains
*         non-zero elements only below N1, and the third is dense.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  COLTYP (workspace/output) INTEGER array, dimension (N)
*         During execution, a label which will indicate which of the
*         following types a column in the Q2 matrix is:
*         1 : non-zero in the upper half only;
*         2 : dense;
*         3 : non-zero in the lower half only;
*         4 : deflated.
*         On exit, COLTYP(i) is the number of columns of type i,
*         for i=1 to 4 only.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0E0, ZERO = 0.0E0, ONE = 1.0E0,
     $                   TWO = 2.0E0, EIGHT = 8.0E0 )
*     ..
*     .. Local Arrays ..
      INTEGER            CTOT( 4 ), PSM( 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1,
     $                   N2, NJ, PJ
      REAL               C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SLAPY2
      EXTERNAL           ISAMAX, SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLACPY, SLAMRG, SROT, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( MIN( 1, ( N / 2 ) ).GT.N1 .OR. ( N / 2 ).LT.N1 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL SSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1.  Since z is the concatenation of
*     two normalized vectors, norm2(z) = sqrt(2).
*
      T = ONE / SQRT( TWO )
      CALL SSCAL( N, T, Z, 1 )
*
*     RHO = ABS( norm(z)**2 * RHO )
*
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 10 I = N1P1, N
         INDXQ( I ) = INDXQ( I ) + N1
   10 CONTINUE
*
*     re-integrate the deflated parts from the last pass
*
      DO 20 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
   20 CONTINUE
      CALL SLAMRG( N1, N2, DLAMDA, 1, 1, INDXC )
      DO 30 I = 1, N
         INDX( I ) = INDXQ( INDXC( I ) )
   30 CONTINUE
*
*     Calculate the allowable deflation tolerance
*
      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = SLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IQ2 = 1
         DO 40 J = 1, N
            I = INDX( J )
            CALL SCOPY( N, Q( 1, I ), 1, Q2( IQ2 ), 1 )
            DLAMDA( J ) = D( I )
            IQ2 = IQ2 + N
   40    CONTINUE
         CALL SLACPY( 'A', N, N, Q2, N, Q, LDQ )
         CALL SCOPY( N, DLAMDA, 1, D, 1 )
         GO TO 190
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      DO 50 I = 1, N1
         COLTYP( I ) = 1
   50 CONTINUE
      DO 60 I = N1P1, N
         COLTYP( I ) = 3
   60 CONTINUE
*
*
      K = 0
      K2 = N + 1
      DO 70 J = 1, N
         NJ = INDX( J )
         IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            COLTYP( NJ ) = 4
            INDXP( K2 ) = NJ
            IF( J.EQ.N )
     $         GO TO 100
         ELSE
            PJ = NJ
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      NJ = INDX( J )
      IF( J.GT.N )
     $   GO TO 100
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( PJ )
         C = Z( NJ )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = SLAPY2( C, S )
         T = D( NJ ) - D( PJ )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( NJ ) = TAU
            Z( PJ ) = ZERO
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) )
     $         COLTYP( NJ ) = 2
            COLTYP( PJ ) = 4
            CALL SROT( N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S )
            T = D( PJ )*C**2 + D( NJ )*S**2
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
            D( PJ ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = PJ
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = PJ
               END IF
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
            PJ = NJ
         ELSE
            K = K + 1
            DLAMDA( K ) = D( PJ )
            W( K ) = Z( PJ )
            INDXP( K ) = PJ
            PJ = NJ
         END IF
      END IF
      GO TO 80
  100 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      DLAMDA( K ) = D( PJ )
      W( K ) = Z( PJ )
      INDXP( K ) = PJ
*
*     Count up the total number of the various types of columns, then
*     form a permutation which positions the four column types into
*     four uniform groups (although one or more of these groups may be
*     empty).
*
      DO 110 J = 1, 4
         CTOT( J ) = 0
  110 CONTINUE
      DO 120 J = 1, N
         CT = COLTYP( J )
         CTOT( CT ) = CTOT( CT ) + 1
  120 CONTINUE
*
*     PSM(*) = Position in SubMatrix (of types 1 through 4)
*
      PSM( 1 ) = 1
      PSM( 2 ) = 1 + CTOT( 1 )
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 )
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 )
      K = N - CTOT( 4 )
*
*     Fill out the INDXC array so that the permutation which it induces
*     will place all type-1 columns first, all type-2 columns next,
*     then all type-3's, and finally all type-4's.
*
      DO 130 J = 1, N
         JS = INDXP( J )
         CT = COLTYP( JS )
         INDX( PSM( CT ) ) = JS
         INDXC( PSM( CT ) ) = J
         PSM( CT ) = PSM( CT ) + 1
  130 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      I = 1
      IQ1 = 1
      IQ2 = 1 + ( CTOT( 1 )+CTOT( 2 ) )*N1
      DO 140 J = 1, CTOT( 1 )
         JS = INDX( I )
         CALL SCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ1 = IQ1 + N1
  140 CONTINUE
*
      DO 150 J = 1, CTOT( 2 )
         JS = INDX( I )
         CALL SCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
         CALL SCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ1 = IQ1 + N1
         IQ2 = IQ2 + N2
  150 CONTINUE
*
      DO 160 J = 1, CTOT( 3 )
         JS = INDX( I )
         CALL SCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ2 = IQ2 + N2
  160 CONTINUE
*
      IQ1 = IQ2
      DO 170 J = 1, CTOT( 4 )
         JS = INDX( I )
         CALL SCOPY( N, Q( 1, JS ), 1, Q2( IQ2 ), 1 )
         IQ2 = IQ2 + N
         Z( I ) = D( JS )
         I = I + 1
  170 CONTINUE
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      CALL SLACPY( 'A', N, CTOT( 4 ), Q2( IQ1 ), N, Q( 1, K+1 ), LDQ )
      CALL SCOPY( N-K, Z( K+1 ), 1, D( K+1 ), 1 )
*
*     Copy CTOT into COLTYP for referencing in SLAED3.
*
      DO 180 J = 1, 4
         COLTYP( J ) = CTOT( J )
  180 CONTINUE
*
  190 CONTINUE
      RETURN
*
*     End of SLAED2
*
      END
      SUBROUTINE SLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,
     $                   CTOT, W, S, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQ, N, N1
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            CTOT( * ), INDX( * )
      REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
     $                   S( * ), W( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED3 finds the roots of the secular equation, as defined by the
*  values in D, W, and RHO, between 1 and K.  It makes the
*  appropriate calls to SLAED4 and then updates the eigenvectors by
*  multiplying the matrix of eigenvectors of the pair of eigensystems
*  being combined by the matrix of eigenvectors of the K-by-K system
*  which is solved here.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  K       (input) INTEGER
*          The number of terms in the rational function to be solved by
*          SLAED4.  K >= 0.
*
*  N       (input) INTEGER
*          The number of rows and columns in the Q matrix.
*          N >= K (deflation may result in N>K).
*
*  N1      (input) INTEGER
*          The location of the last eigenvalue in the leading submatrix.
*          min(1,N) <= N1 <= N/2.
*
*  D       (output) REAL array, dimension (N)
*          D(I) contains the updated eigenvalues for
*          1 <= I <= K.
*
*  Q       (output) REAL array, dimension (LDQ,N)
*          Initially the first K columns are used as workspace.
*          On output the columns 1 to K contain
*          the updated eigenvectors.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  RHO     (input) REAL
*          The value of the parameter in the rank one update equation.
*          RHO >= 0 required.
*
*  DLAMDA  (input/output) REAL array, dimension (K)
*          The first K elements of this array contain the old roots
*          of the deflated updating problem.  These are the poles
*          of the secular equation. May be changed on output by
*          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
*          Cray-2, or Cray C-90, as described above.
*
*  Q2      (input) REAL array, dimension (LDQ2, N)
*          The first K columns of this matrix contain the non-deflated
*          eigenvectors for the split problem.
*
*  INDX    (input) INTEGER array, dimension (N)
*          The permutation used to arrange the columns of the deflated
*          Q matrix into three groups (see SLAED2).
*          The rows of the eigenvectors found by SLAED4 must be likewise
*          permuted before the matrix multiply can take place.
*
*  CTOT    (input) INTEGER array, dimension (4)
*          A count of the total number of the various types of columns
*          in Q, as described in INDX.  The fourth column type is any
*          column which has been deflated.
*
*  W       (input/output) REAL array, dimension (K)
*          The first K elements of this array contain the components
*          of the deflation-adjusted updating vector. Destroyed on
*          output.
*
*  S       (workspace) REAL array, dimension (N1 + 1)*K
*          Will contain the eigenvectors of the repaired matrix which
*          will be multiplied by the previously accumulated eigenvectors
*          to update the system.
*
*  LDS     (input) INTEGER
*          The leading dimension of S.  LDS >= max(1,K).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, IQ2, J, N12, N2, N23
      REAL               TEMP
*     ..
*     .. External Functions ..
      REAL               SLAMC3, SNRM2
      EXTERNAL           SLAMC3, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMM, SLACPY, SLAED4, SLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.K ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED3', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DLAMDA(I) if it is 1; this makes the subsequent
*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DLAMDA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DLAMDA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
      DO 10 I = 1, K
         DLAMDA( I ) = SLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
*
      DO 20 J = 1, K
         CALL SLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
*
*        If the zero finder fails, the computation is terminated.
*
         IF( INFO.NE.0 )
     $      GO TO 120
   20 CONTINUE
*
      IF( K.EQ.1 )
     $   GO TO 110
      IF( K.EQ.2 ) THEN
         DO 30 J = 1, K
            W( 1 ) = Q( 1, J )
            W( 2 ) = Q( 2, J )
            II = INDX( 1 )
            Q( 1, J ) = W( II )
            II = INDX( 2 )
            Q( 2, J ) = W( II )
   30    CONTINUE
         GO TO 110
      END IF
*
*     Compute updated W.
*
      CALL SCOPY( K, W, 1, S, 1 )
*
*     Initialize W(I) = Q(I,I)
*
      CALL SCOPY( K, Q, LDQ+1, W, 1 )
      DO 60 J = 1, K
         DO 40 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   40    CONTINUE
         DO 50 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
   60 CONTINUE
      DO 70 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
   70 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      DO 100 J = 1, K
         DO 80 I = 1, K
            S( I ) = W( I ) / Q( I, J )
   80    CONTINUE
         TEMP = SNRM2( K, S, 1 )
         DO 90 I = 1, K
            II = INDX( I )
            Q( I, J ) = S( II ) / TEMP
   90    CONTINUE
  100 CONTINUE
*
*     Compute the updated eigenvectors.
*
  110 CONTINUE
*
      N2 = N - N1
      N12 = CTOT( 1 ) + CTOT( 2 )
      N23 = CTOT( 2 ) + CTOT( 3 )
*
      CALL SLACPY( 'A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 )
      IQ2 = N1*N12 + 1
      IF( N23.NE.0 ) THEN
         CALL SGEMM( 'N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23,
     $               ZERO, Q( N1+1, 1 ), LDQ )
      ELSE
         CALL SLASET( 'A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ )
      END IF
*
      CALL SLACPY( 'A', N12, K, Q, LDQ, S, N12 )
      IF( N12.NE.0 ) THEN
         CALL SGEMM( 'N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q,
     $               LDQ )
      ELSE
         CALL SLASET( 'A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ )
      END IF
*
*
  120 CONTINUE
      RETURN
*
*     End of SLAED3
*
      END
      SUBROUTINE SLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     December 23, 1999
*
*     .. Scalar Arguments ..
      INTEGER            I, INFO, N
      REAL               DLAM, RHO
*     ..
*     .. Array Arguments ..
      REAL               D( * ), DELTA( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  This subroutine computes the I-th updated eigenvalue of a symmetric
*  rank-one modification to a diagonal matrix whose elements are
*  given in the array d, and that
*
*             D(i) < D(j)  for  i < j
*
*  and that RHO > 0.  This is arranged by the calling routine, and is
*  no loss in generality.  The rank-one modified system is thus
*
*             diag( D )  +  RHO *  Z * Z_transpose.
*
*  where we assume the Euclidean norm of Z is 1.
*
*  The method consists of approximating the rational functions in the
*  secular equation by simpler interpolating rational functions.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The length of all arrays.
*
*  I      (input) INTEGER
*         The index of the eigenvalue to be computed.  1 <= I <= N.
*
*  D      (input) REAL array, dimension (N)
*         The original eigenvalues.  It is assumed that they are in
*         order, D(I) < D(J)  for I < J.
*
*  Z      (input) REAL array, dimension (N)
*         The components of the updating vector.
*
*  DELTA  (output) REAL array, dimension (N)
*         If N .ne. 1, DELTA contains (D(j) - lambda_I) in its  j-th
*         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
*         contains the information necessary to construct the
*         eigenvectors.
*
*  RHO    (input) REAL
*         The scalar in the symmetric updating formula.
*
*  DLAM   (output) REAL
*         The computed lambda_I, the I-th updated eigenvalue.
*
*  INFO   (output) INTEGER
*         = 0:  successful exit
*         > 0:  if INFO = 1, the updating process failed.
*
*  Internal Parameters
*  ===================
*
*  Logical variable ORGATI (origin-at-i?) is used for distinguishing
*  whether D(i) or D(i+1) is treated as the origin.
*
*            ORGATI = .true.    origin at i
*            ORGATI = .false.   origin at i+1
*
*   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
*   if we are working with THREE poles!
*
*   MAXIT is the maximum number of iterations allowed for each
*   eigenvalue.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
      REAL               ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   THREE = 3.0E0, FOUR = 4.0E0, EIGHT = 8.0E0,
     $                   TEN = 10.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      REAL               A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW,
     $                   EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI,
     $                   RHOINV, TAU, TEMP, TEMP1, W
*     ..
*     .. Local Arrays ..
      REAL               ZZ( 3 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAED5, SLAED6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Since this routine is called in an inner loop, we do no argument
*     checking.
*
*     Quick return for N=1 and 2.
*
      INFO = 0
      IF( N.EQ.1 ) THEN
*
*         Presumably, I=1 upon entry
*
         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL SLAED5( I, D, Z, DELTA, RHO, DLAM )
         RETURN
      END IF
*
*     Compute machine epsilon
*
      EPS = SLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
*
*     The case I = N
*
      IF( I.EQ.N ) THEN
*
*        Initialize some basic variables
*
         II = N - 1
         NITER = 1
*
*        Calculate initial guess
*
         MIDPT = RHO / TWO
*
*        If ||Z||_2 is not one, then TEMP should be set to
*        RHO * ||Z||_2^2 / TWO
*
         DO 10 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
   10    CONTINUE
*
         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
   20    CONTINUE
*
         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) +
     $       Z( N )*Z( N ) / DELTA( N )
*
         IF( W.LE.ZERO ) THEN
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) +
     $             Z( N )*Z( N ) / RHO
            IF( C.LE.TEMP ) THEN
               TAU = RHO
            ELSE
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               IF( A.LT.ZERO ) THEN
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
            END IF
*
*           It can be proved that
*               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
*
            DLTLB = MIDPT
            DLTUB = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
*
*           It can be proved that
*               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
*
            DLTLB = ZERO
            DLTUB = MIDPT
         END IF
*
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
   30    CONTINUE
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $            ABS( TAU )*( DPSI+DPHI )
*
         W = RHOINV + PHI + PSI
*
*        Test for convergence
*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            DLAM = D( I ) + TAU
            GO TO 250
         END IF
*
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
*
*        Calculate the new step
*
         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W -
     $       DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( C.LT.ZERO )
     $      C = ABS( C )
         IF( C.EQ.ZERO ) THEN
*          ETA = B/A
*           ETA = RHO - TAU
            ETA = DLTUB - TAU
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
*
*        Note, eta should be positive if w is negative, and
*        eta should be negative otherwise. However,
*        if for some reason caused by roundoff, eta*w > 0,
*        we simply use one Newton step instead. This way
*        will guarantee eta*w < 0.
*
         IF( W*ETA.GT.ZERO )
     $      ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
   50    CONTINUE
*
         TAU = TAU + ETA
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $            ABS( TAU )*( DPSI+DPHI )
*
         W = RHOINV + PHI + PSI
*
*        Main loop to update the values of the array   DELTA
*
         ITER = NITER + 1
*
         DO 90 NITER = ITER, MAXIT
*
*           Test for convergence
*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               DLAM = D( I ) + TAU
               GO TO 250
            END IF
*
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
*
*           Calculate the new step
*
            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W -
     $          DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
            IF( W*ETA.GT.ZERO )
     $         ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
   70       CONTINUE
*
            TAU = TAU + ETA
*
*           Evaluate PSI and the derivative DPSI
*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )
*
*           Evaluate PHI and the derivative DPHI
*
            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $               ABS( TAU )*( DPSI+DPHI )
*
            W = RHOINV + PHI + PSI
   90    CONTINUE
*
*        Return with INFO = 1, NITER = MAXIT and not converged
*
         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250
*
*        End for the case I = N
*
      ELSE
*
*        The case for I < N
*
         NITER = 1
         IP1 = I + 1
*
*        Calculate initial guess
*
         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
         DO 100 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
  100    CONTINUE
*
         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
  110    CONTINUE
*
         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) +
     $       Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
*
         IF( W.GT.ZERO ) THEN
*
*           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
*
*           We choose d(i) as origin.
*
            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            IF( A.GT.ZERO ) THEN
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = ZERO
            DLTUB = MIDPT
         ELSE
*
*           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
*
*           We choose d(i+1) as origin.
*
            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = -MIDPT
            DLTUB = ZERO
         END IF
*
         IF( ORGATI ) THEN
            DO 130 J = 1, N
               DELTA( J ) = ( D( J )-D( I ) ) - TAU
  130       CONTINUE
         ELSE
            DO 140 J = 1, N
               DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
  140       CONTINUE
         END IF
         IF( ORGATI ) THEN
            II = I
         ELSE
            II = I + 1
         END IF
         IIM1 = II - 1
         IIP1 = II + 1
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE
*
         W = RHOINV + PHI + PSI
*
*        W is the value of the secular function with
*        its ii-th element removed.
*
         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO )
     $         SWTCH3 = .TRUE.
         ELSE
            IF( W.GT.ZERO )
     $         SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N )
     $      SWTCH3 = .FALSE.
*
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $            THREE*ABS( TEMP ) + ABS( TAU )*DW
*
*        Test for convergence
*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            IF( ORGATI ) THEN
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            END IF
            GO TO 250
         END IF
*
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
*
*        Calculate the new step
*
         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            IF( ORGATI ) THEN
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )*
     $             ( Z( I ) / DELTA( I ) )**2
            ELSE
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*
     $             ( Z( IP1 ) / DELTA( IP1 ) )**2
            END IF
            A = ( DELTA( I )+DELTA( IP1 ) )*W -
     $          DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )*
     $                   ( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )*
     $                   ( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
*
*           Interpolation using THREE most relevant poles
*
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -
     $             ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*
     $                   ( ( DPSI-TEMP1 )+DPHI )
            ELSE
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -
     $             ( D( IIP1 )-D( IIM1 ) )*TEMP1
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*
     $                   ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            CALL SLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,
     $                   INFO )
            IF( INFO.NE.0 )
     $         GO TO 250
         END IF
*
*        Note, eta should be positive if w is negative, and
*        eta should be negative otherwise. However,
*        if for some reason caused by roundoff, eta*w > 0,
*        we simply use one Newton step instead. This way
*        will guarantee eta*w < 0.
*
         IF( W*ETA.GE.ZERO )
     $      ETA = -W / DW
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
*
         PREW = W
*
  170    CONTINUE
         DO 180 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
  180    CONTINUE
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 190 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  190    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         DPHI = ZERO
         PHI = ZERO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  200    CONTINUE
*
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $            THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
*
         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN )
     $         SWTCH = .TRUE.
         ELSE
            IF( W.GT.ABS( PREW ) / TEN )
     $         SWTCH = .TRUE.
         END IF
*
         TAU = TAU + ETA
*
*        Main loop to update the values of the array   DELTA
*
         ITER = NITER + 1
*
         DO 240 NITER = ITER, MAXIT
*
*           Test for convergence
*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               IF( ORGATI ) THEN
                  DLAM = D( I ) + TAU
               ELSE
                  DLAM = D( IP1 ) + TAU
               END IF
               GO TO 250
            END IF
*
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
*
*           Calculate the new step
*
            IF( .NOT.SWTCH3 ) THEN
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DELTA( IP1 )*DW -
     $                   ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  ELSE
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*
     $                   ( Z( IP1 ) / DELTA( IP1 ) )**2
                  END IF
               ELSE
                  TEMP = Z( II ) / DELTA( II )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               END IF
               A = ( DELTA( I )+DELTA( IP1 ) )*W -
     $             DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN
                           A = Z( I )*Z( I ) + DELTA( IP1 )*
     $                         DELTA( IP1 )*( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) +
     $                         DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        END IF
                     ELSE
                        A = DELTA( I )*DELTA( I )*DPSI +
     $                      DELTA( IP1 )*DELTA( IP1 )*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            ELSE
*
*              Interpolation using THREE most relevant poles
*
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               ELSE
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -
     $                   ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*
     $                         ( ( DPSI-TEMP1 )+DPHI )
                  ELSE
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -
     $                   ( D( IIP1 )-D( IIM1 ) )*TEMP1
                     ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*
     $                         ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               CALL SLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,
     $                      INFO )
               IF( INFO.NE.0 )
     $            GO TO 250
            END IF
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
            IF( W*ETA.GE.ZERO )
     $         ETA = -W / DW
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
*
            DO 210 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
  210       CONTINUE
*
            TAU = TAU + ETA
            PREW = W
*
*           Evaluate PSI and the derivative DPSI
*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 220 J = 1, IIM1
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  220       CONTINUE
            ERRETM = ABS( ERRETM )
*
*           Evaluate PHI and the derivative DPHI
*
            DPHI = ZERO
            PHI = ZERO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  230       CONTINUE
*
            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $               THREE*ABS( TEMP ) + ABS( TAU )*DW
            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN )
     $         SWTCH = .NOT.SWTCH
*
  240    CONTINUE
*
*        Return with INFO = 1, NITER = MAXIT and not converged
*
         INFO = 1
         IF( ORGATI ) THEN
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         END IF
*
      END IF
*
  250 CONTINUE
*
      RETURN
*
*     End of SLAED4
*
      END
      SUBROUTINE SLAED5( I, D, Z, DELTA, RHO, DLAM )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            I
      REAL               DLAM, RHO
*     ..
*     .. Array Arguments ..
      REAL               D( 2 ), DELTA( 2 ), Z( 2 )
*     ..
*
*  Purpose
*  =======
*
*  This subroutine computes the I-th eigenvalue of a symmetric rank-one
*  modification of a 2-by-2 diagonal matrix
*
*             diag( D )  +  RHO *  Z * transpose(Z) .
*
*  The diagonal elements in the array D are assumed to satisfy
*
*             D(i) < D(j)  for  i < j .
*
*  We also assume RHO > 0 and that the Euclidean norm of the vector
*  Z is one.
*
*  Arguments
*  =========
*
*  I      (input) INTEGER
*         The index of the eigenvalue to be computed.  I = 1 or I = 2.
*
*  D      (input) REAL array, dimension (2)
*         The original eigenvalues.  We assume D(1) < D(2).
*
*  Z      (input) REAL array, dimension (2)
*         The components of the updating vector.
*
*  DELTA  (output) REAL array, dimension (2)
*         The vector DELTA contains the information necessary
*         to construct the eigenvectors.
*
*  RHO    (input) REAL
*         The scalar in the symmetric updating formula.
*
*  DLAM   (output) REAL
*         The computed lambda_I, the I-th updated eigenvalue.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   FOUR = 4.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               B, C, DEL, TAU, TEMP, W
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      DEL = D( 2 ) - D( 1 )
      IF( I.EQ.1 ) THEN
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
         IF( W.GT.ZERO ) THEN
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DEL
*
*           B > ZERO, always
*
            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
            DLAM = D( 1 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / TAU
            DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         ELSE
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DEL
            IF( B.GT.ZERO ) THEN
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            ELSE
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            END IF
            DLAM = D( 2 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            DELTA( 2 ) = -Z( 2 ) / TAU
         END IF
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE
*
*     Now I=2
*
         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         IF( B.GT.ZERO ) THEN
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         ELSE
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         END IF
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      END IF
      RETURN
*
*     End OF SLAED5
*
      END
      SUBROUTINE SLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      LOGICAL            ORGATI
      INTEGER            INFO, KNITER
      REAL               FINIT, RHO, TAU
*     ..
*     .. Array Arguments ..
      REAL               D( 3 ), Z( 3 )
*     ..
*
*  Purpose
*  =======
*
*  SLAED6 computes the positive or negative root (closest to the origin)
*  of
*                   z(1)        z(2)        z(3)
*  f(x) =   rho + --------- + ---------- + ---------
*                  d(1)-x      d(2)-x      d(3)-x
*
*  It is assumed that
*
*        if ORGATI = .true. the root is between d(2) and d(3);
*        otherwise it is between d(1) and d(2)
*
*  This routine will be called by SLAED4 when necessary. In most cases,
*  the root sought is the smallest in magnitude, though it might not be
*  in some extremely rare situations.
*
*  Arguments
*  =========
*
*  KNITER       (input) INTEGER
*               Refer to SLAED4 for its significance.
*
*  ORGATI       (input) LOGICAL
*               If ORGATI is true, the needed root is between d(2) and
*               d(3); otherwise it is between d(1) and d(2).  See
*               SLAED4 for further details.
*
*  RHO          (input) REAL
*               Refer to the equation f(x) above.
*
*  D            (input) REAL array, dimension (3)
*               D satisfies d(1) < d(2) < d(3).
*
*  Z            (input) REAL array, dimension (3)
*               Each of the elements in z must be positive.
*
*  FINIT        (input) REAL
*               The value of f at 0. It is more accurate than the one
*               evaluated inside this routine (if someone wants to do
*               so).
*
*  TAU          (output) REAL
*               The root of the equation f(x).
*
*  INFO         (output) INTEGER
*               = 0: successful exit
*               > 0: if INFO = 1, failure to converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 20 )
      REAL               ZERO, ONE, TWO, THREE, FOUR, EIGHT
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   THREE = 3.0E0, FOUR = 4.0E0, EIGHT = 8.0E0 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Local Arrays ..
      REAL               DSCALE( 3 ), ZSCALE( 3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, SCALE
      INTEGER            I, ITER, NITER
      REAL               A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F,
     $                   FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1,
     $                   SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4
*     ..
*     .. Save statement ..
      SAVE               FIRST, SMALL1, SMINV1, SMALL2, SMINV2, EPS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      NITER = 1
      TAU = ZERO
      IF( KNITER.EQ.2 ) THEN
         IF( ORGATI ) THEN
            TEMP = ( D( 3 )-D( 2 ) ) / TWO
            C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
            A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
            B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
         ELSE
            TEMP = ( D( 1 )-D( 2 ) ) / TWO
            C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
            A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
            B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
         END IF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            TAU = B / A
         ELSE IF( A.LE.ZERO ) THEN
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         TEMP = RHO + Z( 1 ) / ( D( 1 )-TAU ) +
     $          Z( 2 ) / ( D( 2 )-TAU ) + Z( 3 ) / ( D( 3 )-TAU )
         IF( ABS( FINIT ).LE.ABS( TEMP ) )
     $      TAU = ZERO
      END IF
*
*     On first call to routine, get machine parameters for
*     possible scaling to avoid overflow
*
      IF( FIRST ) THEN
         EPS = SLAMCH( 'Epsilon' )
         BASE = SLAMCH( 'Base' )
         SMALL1 = BASE**( INT( LOG( SLAMCH( 'SafMin' ) ) / LOG( BASE ) /
     $            THREE ) )
         SMINV1 = ONE / SMALL1
         SMALL2 = SMALL1*SMALL1
         SMINV2 = SMINV1*SMINV1
         FIRST = .FALSE.
      END IF
*
*     Determine if scaling of inputs necessary to avoid overflow
*     when computing 1/TEMP**3
*
      IF( ORGATI ) THEN
         TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
      ELSE
         TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
      END IF
      SCALE = .FALSE.
      IF( TEMP.LE.SMALL1 ) THEN
         SCALE = .TRUE.
         IF( TEMP.LE.SMALL2 ) THEN
*
*        Scale up by power of radix nearest 1/SAFMIN**(2/3)
*
            SCLFAC = SMINV2
            SCLINV = SMALL2
         ELSE
*
*        Scale up by power of radix nearest 1/SAFMIN**(1/3)
*
            SCLFAC = SMINV1
            SCLINV = SMALL1
         END IF
*
*        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
*
         DO 10 I = 1, 3
            DSCALE( I ) = D( I )*SCLFAC
            ZSCALE( I ) = Z( I )*SCLFAC
   10    CONTINUE
         TAU = TAU*SCLFAC
      ELSE
*
*        Copy D and Z to DSCALE and ZSCALE
*
         DO 20 I = 1, 3
            DSCALE( I ) = D( I )
            ZSCALE( I ) = Z( I )
   20    CONTINUE
      END IF
*
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO 30 I = 1, 3
         TEMP = ONE / ( DSCALE( I )-TAU )
         TEMP1 = ZSCALE( I )*TEMP
         TEMP2 = TEMP1*TEMP
         TEMP3 = TEMP2*TEMP
         FC = FC + TEMP1 / DSCALE( I )
         DF = DF + TEMP2
         DDF = DDF + TEMP3
   30 CONTINUE
      F = FINIT + TAU*FC
*
      IF( ABS( F ).LE.ZERO )
     $   GO TO 60
*
*        Iteration begins
*
*     It is not hard to see that
*
*           1) Iterations will go up monotonically
*              if FINIT < 0;
*
*           2) Iterations will go down monotonically
*              if FINIT > 0.
*
      ITER = NITER + 1
*
      DO 50 NITER = ITER, MAXIT
*
         IF( ORGATI ) THEN
            TEMP1 = DSCALE( 2 ) - TAU
            TEMP2 = DSCALE( 3 ) - TAU
         ELSE
            TEMP1 = DSCALE( 1 ) - TAU
            TEMP2 = DSCALE( 2 ) - TAU
         END IF
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
         B = TEMP1*TEMP2*F
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            ETA = B / A
         ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( F*ETA.GE.ZERO ) THEN
            ETA = -F / DF
         END IF
*
         TEMP = ETA + TAU
         IF( ORGATI ) THEN
            IF( ETA.GT.ZERO .AND. TEMP.GE.DSCALE( 3 ) )
     $         ETA = ( DSCALE( 3 )-TAU ) / TWO
            IF( ETA.LT.ZERO .AND. TEMP.LE.DSCALE( 2 ) )
     $         ETA = ( DSCALE( 2 )-TAU ) / TWO
         ELSE
            IF( ETA.GT.ZERO .AND. TEMP.GE.DSCALE( 2 ) )
     $         ETA = ( DSCALE( 2 )-TAU ) / TWO
            IF( ETA.LT.ZERO .AND. TEMP.LE.DSCALE( 1 ) )
     $         ETA = ( DSCALE( 1 )-TAU ) / TWO
         END IF
         TAU = TAU + ETA
*
         FC = ZERO
         ERRETM = ZERO
         DF = ZERO
         DDF = ZERO
         DO 40 I = 1, 3
            TEMP = ONE / ( DSCALE( I )-TAU )
            TEMP1 = ZSCALE( I )*TEMP
            TEMP2 = TEMP1*TEMP
            TEMP3 = TEMP2*TEMP
            TEMP4 = TEMP1 / DSCALE( I )
            FC = FC + TEMP4
            ERRETM = ERRETM + ABS( TEMP4 )
            DF = DF + TEMP2
            DDF = DDF + TEMP3
   40    CONTINUE
         F = FINIT + TAU*FC
         ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) +
     $            ABS( TAU )*DF
         IF( ABS( F ).LE.EPS*ERRETM )
     $      GO TO 60
   50 CONTINUE
      INFO = 1
   60 CONTINUE
*
*     Undo scaling
*
      IF( SCALE )
     $   TAU = TAU*SCLINV
      RETURN
*
*     End of SLAED6
*
      END
      SUBROUTINE SLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
     $                   LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,
     $                   PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,
     $                   QSIZ, TLVLS
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
     $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      REAL               D( * ), GIVNUM( 2, * ), Q( LDQ, * ),
     $                   QSTORE( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED7 computes the updated eigensystem of a diagonal
*  matrix after modification by a rank-one symmetric matrix. This
*  routine is used only for the eigenproblem which requires all
*  eigenvalues and optionally eigenvectors of a dense symmetric matrix
*  that has been reduced to tridiagonal form.  SLAED1 handles
*  the case in which all eigenvalues and eigenvectors of a symmetric
*  tridiagonal matrix are desired.
*
*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)
*
*     where Z = Q'u, u is a vector of length N with ones in the
*     CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
*
*     The eigenvectors of the original matrix are stored in Q, and the
*     eigenvalues are in D.  The algorithm consists of three stages:
*
*        The first stage consists of deflating the size of the problem
*        when there are multiple eigenvalues or if there is a zero in
*        the Z vector.  For each such occurence the dimension of the
*        secular equation problem is reduced by one.  This stage is
*        performed by the routine SLAED8.
*
*        The second stage consists of calculating the updated
*        eigenvalues. This is done by finding the roots of the secular
*        equation via the routine SLAED4 (as called by SLAED9).
*        This routine also calculates the eigenvectors of the current
*        problem.
*
*        The final stage consists of computing the updated eigenvectors
*        directly using the updated eigenvalues.  The eigenvectors for
*        the current problem are multiplied with the eigenvectors from
*        the overall problem.
*
*  Arguments
*  =========
*
*  ICOMPQ  (input) INTEGER
*          = 0:  Compute eigenvalues only.
*          = 1:  Compute eigenvectors of original dense symmetric matrix
*                also.  On entry, Q contains the orthogonal matrix used
*                to reduce the original matrix to tridiagonal form.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  QSIZ   (input) INTEGER
*         The dimension of the orthogonal matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  TLVLS  (input) INTEGER
*         The total number of merging levels in the overall divide and
*         conquer tree.
*
*  CURLVL (input) INTEGER
*         The current level in the overall merge routine,
*         0 <= CURLVL <= TLVLS.
*
*  CURPBM (input) INTEGER
*         The current problem in the current level in the overall
*         merge routine (counting from upper left to lower right).
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, the eigenvalues of the rank-1-perturbed matrix.
*         On exit, the eigenvalues of the repaired matrix.
*
*  Q      (input/output) REAL array, dimension (LDQ, N)
*         On entry, the eigenvectors of the rank-1-perturbed matrix.
*         On exit, the eigenvectors of the repaired tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (output) INTEGER array, dimension (N)
*         The permutation which will reintegrate the subproblem just
*         solved back into sorted order, i.e., D( INDXQ( I = 1, N ) )
*         will be in ascending order.
*
*  RHO    (input) REAL
*         The subdiagonal element used to create the rank-1
*         modification.
*
*  CUTPNT (input) INTEGER
*         Contains the location of the last eigenvalue in the leading
*         sub-matrix.  min(1,N) <= CUTPNT <= N.
*
*  QSTORE (input/output) REAL array, dimension (N**2+1)
*         Stores eigenvectors of submatrices encountered during
*         divide and conquer, packed together. QPTR points to
*         beginning of the submatrices.
*
*  QPTR   (input/output) INTEGER array, dimension (N+2)
*         List of indices pointing to beginning of submatrices stored
*         in QSTORE. The submatrices are numbered starting at the
*         bottom left of the divide and conquer tree, from left to
*         right and bottom to top.
*
*  PRMPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in PERM a
*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
*         indicates the size of the permutation and also the size of
*         the full, non-deflated problem.
*
*  PERM   (input) INTEGER array, dimension (N lg N)
*         Contains the permutations (from deflation and sorting) to be
*         applied to each eigenblock.
*
*  GIVPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in GIVCOL a
*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
*         indicates the number of Givens rotations.
*
*  GIVCOL (input) INTEGER array, dimension (2, N lg N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (input) REAL array, dimension (2, N lg N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  WORK   (workspace) REAL array, dimension (3*N+QSIZ*N)
*
*  IWORK  (workspace) INTEGER array, dimension (4*N)
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP,
     $                   IQ2, IS, IW, IZ, K, LDQ2, N1, N2, PTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLAED8, SLAED9, SLAEDA, SLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED7', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are for bookkeeping purposes only.  They are
*     integer pointers which indicate the portion of the workspace
*     used by a particular array in SLAED8 and SLAED9.
*
      IF( ICOMPQ.EQ.1 ) THEN
         LDQ2 = QSIZ
      ELSE
         LDQ2 = N
      END IF
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
      IS = IQ2 + N*LDQ2
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
     $             GIVCOL, GIVNUM, QSTORE, QPTR, WORK( IZ ),
     $             WORK( IZ+N ), INFO )
*
*     When solving the final problem, we no longer need the stored data,
*     so we will overwrite the data from this level onto the previously
*     used storage space.
*
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
*
*     Sort and Deflate eigenvalues.
*
      CALL SLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT,
     $             WORK( IZ ), WORK( IDLMDA ), WORK( IQ2 ), LDQ2,
     $             WORK( IW ), PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ),
     $             GIVCOL( 1, GIVPTR( CURR ) ),
     $             GIVNUM( 1, GIVPTR( CURR ) ), IWORK( INDXP ),
     $             IWORK( INDX ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         CALL SLAED9( K, 1, K, N, D, WORK( IS ), K, RHO, WORK( IDLMDA ),
     $                WORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
         IF( ICOMPQ.EQ.1 ) THEN
            CALL SGEMM( 'N', 'N', QSIZ, K, K, ONE, WORK( IQ2 ), LDQ2,
     $                  QSTORE( QPTR( CURR ) ), K, ZERO, Q, LDQ )
         END IF
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
*
*     Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
*
   30 CONTINUE
      RETURN
*
*     End of SLAED7
*
      END
      SUBROUTINE SLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,
     $                   CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, INDXP, INDX, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,
     $                   QSIZ
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      REAL               D( * ), DLAMDA( * ), GIVNUM( 2, * ),
     $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED8 merges the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny element in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  ICOMPQ  (input) INTEGER
*          = 0:  Compute eigenvalues only.
*          = 1:  Compute eigenvectors of original dense symmetric matrix
*                also.  On entry, Q contains the orthogonal matrix used
*                to reduce the original matrix to tridiagonal form.
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  QSIZ   (input) INTEGER
*         The dimension of the orthogonal matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, the eigenvalues of the two submatrices to be
*         combined.  On exit, the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  Q      (input/output) REAL array, dimension (LDQ,N)
*         If ICOMPQ = 0, Q is not referenced.  Otherwise,
*         on entry, Q contains the eigenvectors of the partially solved
*         system which has been previously updated in matrix
*         multiplies with other partially solved eigensystems.
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input) INTEGER array, dimension (N)
*         The permutation which separately sorts the two sub-problems
*         in D into ascending order.  Note that elements in the second
*         half of this permutation must first have CUTPNT added to
*         their values in order to be accurate.
*
*  RHO    (input/output) REAL
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         SLAED3.
*
*  CUTPNT (input) INTEGER
*         The location of the last eigenvalue in the leading
*         sub-matrix.  min(1,N) <= CUTPNT <= N.
*
*  Z      (input) REAL array, dimension (N)
*         On entry, Z contains the updating vector (the last row of
*         the first sub-eigenvector matrix and the first row of the
*         second sub-eigenvector matrix).
*         On exit, the contents of Z are destroyed by the updating
*         process.
*
*  DLAMDA (output) REAL array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         SLAED3 to form the secular equation.
*
*  Q2     (output) REAL array, dimension (LDQ2,N)
*         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
*         a copy of the first K eigenvectors which will be used by
*         SLAED7 in a matrix multiply (SGEMM) to update the new
*         eigenvectors.
*
*  LDQ2   (input) INTEGER
*         The leading dimension of the array Q2.  LDQ2 >= max(1,N).
*
*  W      (output) REAL array, dimension (N)
*         The first k values of the final deflation-altered z-vector and
*         will be passed to SLAED3.
*
*  PERM   (output) INTEGER array, dimension (N)
*         The permutations (from deflation and sorting) to be applied
*         to each eigenblock.
*
*  GIVPTR (output) INTEGER
*         The number of Givens rotations which took place in this
*         subproblem.
*
*  GIVCOL (output) INTEGER array, dimension (2, N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (output) REAL array, dimension (2, N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of D into ascending
*         order.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0E0, ZERO = 0.0E0, ONE = 1.0E0,
     $                   TWO = 2.0E0, EIGHT = 8.0E0 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      REAL               C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SLAPY2
      EXTERNAL           ISAMAX, SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLACPY, SLAMRG, SROT, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -10
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL SSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1
*
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL SSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL SLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
*
*     Calculate the allowable deflation tolerence
*
      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = SLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IF( ICOMPQ.EQ.0 ) THEN
            DO 50 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
   50       CONTINUE
         ELSE
            DO 60 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
               CALL SCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   60       CONTINUE
            CALL SLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ),
     $                   LDQ )
         END IF
         RETURN
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      K = 0
      GIVPTR = 0
      K2 = N + 1
      DO 70 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N )
     $         GO TO 110
         ELSE
            JLAM = J
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 100
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( JLAM )
         C = Z( J )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = SLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( J ) = TAU
            Z( JLAM ) = ZERO
*
*           Record the appropriate Givens rotation
*
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            IF( ICOMPQ.EQ.1 ) THEN
               CALL SROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                    Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            END IF
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 80
  100 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  110 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      IF( ICOMPQ.EQ.0 ) THEN
         DO 120 J = 1, N
            JP = INDXP( J )
            DLAMDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
  120    CONTINUE
      ELSE
         DO 130 J = 1, N
            JP = INDXP( J )
            DLAMDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
            CALL SCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  130    CONTINUE
      END IF
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      IF( K.LT.N ) THEN
         IF( ICOMPQ.EQ.0 ) THEN
            CALL SCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
         ELSE
            CALL SCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
            CALL SLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2,
     $                   Q( 1, K+1 ), LDQ )
         END IF
      END IF
*
      RETURN
*
*     End of SLAED8
*
      END
      SUBROUTINE SLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,
     $                   S, LDS, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N
      REAL               RHO
*     ..
*     .. Array Arguments ..
      REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),
     $                   W( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAED9 finds the roots of the secular equation, as defined by the
*  values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
*  appropriate calls to SLAED4 and then stores the new matrix of
*  eigenvectors for use in calculating the next level of Z vectors.
*
*  Arguments
*  =========
*
*  K       (input) INTEGER
*          The number of terms in the rational function to be solved by
*          SLAED4.  K >= 0.
*
*  KSTART  (input) INTEGER
*  KSTOP   (input) INTEGER
*          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
*          are to be computed.  1 <= KSTART <= KSTOP <= K.
*
*  N       (input) INTEGER
*          The number of rows and columns in the Q matrix.
*          N >= K (delation may result in N > K).
*
*  D       (output) REAL array, dimension (N)
*          D(I) contains the updated eigenvalues
*          for KSTART <= I <= KSTOP.
*
*  Q       (workspace) REAL array, dimension (LDQ,N)
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max( 1, N ).
*
*  RHO     (input) REAL
*          The value of the parameter in the rank one update equation.
*          RHO >= 0 required.
*
*  DLAMDA  (input) REAL array, dimension (K)
*          The first K elements of this array contain the old roots
*          of the deflated updating problem.  These are the poles
*          of the secular equation.
*
*  W       (input) REAL array, dimension (K)
*          The first K elements of this array contain the components
*          of the deflation-adjusted updating vector.
*
*  S       (output) REAL array, dimension (LDS, K)
*          Will contain the eigenvectors of the repaired matrix which
*          will be stored for subsequent Z vector calculation and
*          multiplied by the previously accumulated eigenvectors
*          to update the system.
*
*  LDS     (input) INTEGER
*          The leading dimension of S.  LDS >= max( 1, K ).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               TEMP
*     ..
*     .. External Functions ..
      REAL               SLAMC3, SNRM2
      EXTERNAL           SLAMC3, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAED4, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) THEN
         INFO = -2
      ELSE IF( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) )
     $          THEN
         INFO = -3
      ELSE IF( N.LT.K ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDS.LT.MAX( 1, K ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED9', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DLAMDA(I) if it is 1; this makes the subsequent
*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DLAMDA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DLAMDA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
      DO 10 I = 1, N
         DLAMDA( I ) = SLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
*
      DO 20 J = KSTART, KSTOP
         CALL SLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
*
*        If the zero finder fails, the computation is terminated.
*
         IF( INFO.NE.0 )
     $      GO TO 120
   20 CONTINUE
*
      IF( K.EQ.1 .OR. K.EQ.2 ) THEN
         DO 40 I = 1, K
            DO 30 J = 1, K
               S( J, I ) = Q( J, I )
   30       CONTINUE
   40    CONTINUE
         GO TO 120
      END IF
*
*     Compute updated W.
*
      CALL SCOPY( K, W, 1, S, 1 )
*
*     Initialize W(I) = Q(I,I)
*
      CALL SCOPY( K, Q, LDQ+1, W, 1 )
      DO 70 J = 1, K
         DO 50 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
         DO 60 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   60    CONTINUE
   70 CONTINUE
      DO 80 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
   80 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      DO 110 J = 1, K
         DO 90 I = 1, K
            Q( I, J ) = W( I ) / Q( I, J )
   90    CONTINUE
         TEMP = SNRM2( K, Q( 1, J ), 1 )
         DO 100 I = 1, K
            S( I, J ) = Q( I, J ) / TEMP
  100    CONTINUE
  110 CONTINUE
*
  120 CONTINUE
      RETURN
*
*     End of SLAED9
*
      END
      SUBROUTINE SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CURLVL, CURPBM, INFO, N, TLVLS
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),
     $                   PRMPTR( * ), QPTR( * )
      REAL               GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAEDA computes the Z vector corresponding to the merge step in the
*  CURLVLth step of the merge process with TLVLS steps for the CURPBMth
*  problem.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  TLVLS  (input) INTEGER
*         The total number of merging levels in the overall divide and
*         conquer tree.
*
*  CURLVL (input) INTEGER
*         The current level in the overall merge routine,
*         0 <= curlvl <= tlvls.
*
*  CURPBM (input) INTEGER
*         The current problem in the current level in the overall
*         merge routine (counting from upper left to lower right).
*
*  PRMPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in PERM a
*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
*         indicates the size of the permutation and incidentally the
*         size of the full, non-deflated problem.
*
*  PERM   (input) INTEGER array, dimension (N lg N)
*         Contains the permutations (from deflation and sorting) to be
*         applied to each eigenblock.
*
*  GIVPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in GIVCOL a
*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
*         indicates the number of Givens rotations.
*
*  GIVCOL (input) INTEGER array, dimension (2, N lg N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (input) REAL array, dimension (2, N lg N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  Q      (input) REAL array, dimension (N**2)
*         Contains the square eigenblocks from previous levels, the
*         starting positions for blocks are given by QPTR.
*
*  QPTR   (input) INTEGER array, dimension (N+2)
*         Contains a list of pointers which indicate where in Q an
*         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates
*         the size of the block.
*
*  Z      (output) REAL array, dimension (N)
*         On output this vector contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*
*  ZTEMP  (workspace) REAL array, dimension (N)
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            BSIZ1, BSIZ2, CURR, I, K, MID, PSIZ1, PSIZ2,
     $                   PTR, ZPTR1
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMV, SROT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -1
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAEDA', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine location of first number in second half.
*
      MID = N / 2 + 1
*
*     Gather last/first rows of appropriate eigenblocks into center of Z
*
      PTR = 1
*
*     Determine location of lowest level subproblem in the full storage
*     scheme
*
      CURR = PTR + CURPBM*2**CURLVL + 2**( CURLVL-1 ) - 1
*
*     Determine size of these matrices.  We add HALF to the value of
*     the SQRT in case the machine underestimates one of these square
*     roots.
*
      BSIZ1 = INT( HALF+SQRT( REAL( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
      BSIZ2 = INT( HALF+SQRT( REAL( QPTR( CURR+2 )-QPTR( CURR+1 ) ) ) )
      DO 10 K = 1, MID - BSIZ1 - 1
         Z( K ) = ZERO
   10 CONTINUE
      CALL SCOPY( BSIZ1, Q( QPTR( CURR )+BSIZ1-1 ), BSIZ1,
     $            Z( MID-BSIZ1 ), 1 )
      CALL SCOPY( BSIZ2, Q( QPTR( CURR+1 ) ), BSIZ2, Z( MID ), 1 )
      DO 20 K = MID + BSIZ2, N
         Z( K ) = ZERO
   20 CONTINUE
*
*     Loop thru remaining levels 1 -> CURLVL applying the Givens
*     rotations and permutation and then multiplying the center matrices
*     against the current Z.
*
      PTR = 2**TLVLS + 1
      DO 70 K = 1, CURLVL - 1
         CURR = PTR + CURPBM*2**( CURLVL-K ) + 2**( CURLVL-K-1 ) - 1
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         ZPTR1 = MID - PSIZ1
*
*       Apply Givens at CURR and CURR+1
*
         DO 30 I = GIVPTR( CURR ), GIVPTR( CURR+1 ) - 1
            CALL SROT( 1, Z( ZPTR1+GIVCOL( 1, I )-1 ), 1,
     $                 Z( ZPTR1+GIVCOL( 2, I )-1 ), 1, GIVNUM( 1, I ),
     $                 GIVNUM( 2, I ) )
   30    CONTINUE
         DO 40 I = GIVPTR( CURR+1 ), GIVPTR( CURR+2 ) - 1
            CALL SROT( 1, Z( MID-1+GIVCOL( 1, I ) ), 1,
     $                 Z( MID-1+GIVCOL( 2, I ) ), 1, GIVNUM( 1, I ),
     $                 GIVNUM( 2, I ) )
   40    CONTINUE
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         DO 50 I = 0, PSIZ1 - 1
            ZTEMP( I+1 ) = Z( ZPTR1+PERM( PRMPTR( CURR )+I )-1 )
   50    CONTINUE
         DO 60 I = 0, PSIZ2 - 1
            ZTEMP( PSIZ1+I+1 ) = Z( MID+PERM( PRMPTR( CURR+1 )+I )-1 )
   60    CONTINUE
*
*        Multiply Blocks at CURR and CURR+1
*
*        Determine size of these matrices.  We add HALF to the value of
*        the SQRT in case the machine underestimates one of these
*        square roots.
*
         BSIZ1 = INT( HALF+SQRT( REAL( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
         BSIZ2 = INT( HALF+SQRT( REAL( QPTR( CURR+2 )-QPTR( CURR+
     $           1 ) ) ) )
         IF( BSIZ1.GT.0 ) THEN
            CALL SGEMV( 'T', BSIZ1, BSIZ1, ONE, Q( QPTR( CURR ) ),
     $                  BSIZ1, ZTEMP( 1 ), 1, ZERO, Z( ZPTR1 ), 1 )
         END IF
         CALL SCOPY( PSIZ1-BSIZ1, ZTEMP( BSIZ1+1 ), 1, Z( ZPTR1+BSIZ1 ),
     $               1 )
         IF( BSIZ2.GT.0 ) THEN
            CALL SGEMV( 'T', BSIZ2, BSIZ2, ONE, Q( QPTR( CURR+1 ) ),
     $                  BSIZ2, ZTEMP( PSIZ1+1 ), 1, ZERO, Z( MID ), 1 )
         END IF
         CALL SCOPY( PSIZ2-BSIZ2, ZTEMP( PSIZ1+BSIZ2+1 ), 1,
     $               Z( MID+BSIZ2 ), 1 )
*
         PTR = PTR + 2**( TLVLS-K )
   70 CONTINUE
*
      RETURN
*
*     End of SLAEDA
*
      END
      SUBROUTINE SLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  SLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) REAL
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) REAL
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
*
*  C       (input) REAL
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) REAL
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) REAL
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) REAL
*  SN1     (output) REAL
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      REAL               AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of SLAEV2
*
      END
      SUBROUTINE SLAMRG( N1, N2, A, STRD1, STRD2, INDEX )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            N1, N2, STRD1, STRD2
*     ..
*     .. Array Arguments ..
      INTEGER            INDEX( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAMRG will create a permutation list which will merge the elements
*  of A (which is composed of two independently sorted sets) into a
*  single set which is sorted in ascending order.
*
*  Arguments
*  =========
*
*  N1     (input) INTEGER
*  N2     (input) INTEGER
*         These arguements contain the respective lengths of the two
*         sorted lists to be merged.
*
*  A      (input) REAL array, dimension (N1+N2)
*         The first N1 elements of A contain a list of numbers which
*         are sorted in either ascending or descending order.  Likewise
*         for the final N2 elements.
*
*  STRD1  (input) INTEGER
*  STRD2  (input) INTEGER
*         These are the strides to be taken through the array A.
*         Allowable strides are 1 and -1.  They indicate whether a
*         subset of A is sorted in ascending (STRDx = 1) or descending
*         (STRDx = -1) order.
*
*  INDEX  (output) INTEGER array, dimension (N1+N2)
*         On exit this array will contain a permutation such that
*         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
*         sorted in ascending order.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IND1, IND2, N1SV, N2SV
*     ..
*     .. Executable Statements ..
*
      N1SV = N1
      N2SV = N2
      IF( STRD1.GT.0 ) THEN
         IND1 = 1
      ELSE
         IND1 = N1
      END IF
      IF( STRD2.GT.0 ) THEN
         IND2 = 1 + N1
      ELSE
         IND2 = N1 + N2
      END IF
      I = 1
*     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
         IF( A( IND1 ).LE.A( IND2 ) ) THEN
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + STRD1
            N1SV = N1SV - 1
         ELSE
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + STRD2
            N2SV = N2SV - 1
         END IF
         GO TO 10
      END IF
*     end while
      IF( N1SV.EQ.0 ) THEN
         DO 20 N1SV = 1, N2SV
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + STRD2
   20    CONTINUE
      ELSE
*     N2SV .EQ. 0
         DO 30 N2SV = 1, N1SV
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + STRD1
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of SLAMRG
*
      END
      REAL             FUNCTION SLANST( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  SLANST  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real symmetric tridiagonal matrix A.
*
*  Description
*  ===========
*
*  SLANST returns the value
*
*     SLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in SLANST as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, SLANST is
*          set to zero.
*
*  D       (input) REAL array, dimension (N)
*          The diagonal elements of A.
*
*  E       (input) REAL array, dimension (N-1)
*          The (n-1) sub-diagonal or super-diagonal elements of A.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $              ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+
     $                 ABS( E( I-1 ) ) )
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL SLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL SLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      SLANST = ANORM
      RETURN
*
*     End of SLANST
*
      END
      REAL             FUNCTION SLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               X, Y
*     ..
*
*  Purpose
*  =======
*
*  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) REAL
*  Y       (input) REAL
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         SLAPY2 = W
      ELSE
         SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of SLAPY2
*
      END
      REAL             FUNCTION SLAPY3( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               X, Y, Z
*     ..
*
*  Purpose
*  =======
*
*  SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) REAL
*  Y       (input) REAL
*  Z       (input) REAL
*          X, Y and Z specify the values x, y and z.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         SLAPY3 = ZERO
      ELSE
         SLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of SLAPY3
*
      END
      SUBROUTINE SLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      REAL               CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  SLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine SROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in SBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) REAL
*          The first component of vector to be rotated.
*
*  G       (input) REAL
*          The second component of vector to be rotated.
*
*  CS      (output) REAL
*          The cosine of the rotation.
*
*  SN      (output) REAL
*          The sine of the rotation.
*
*  R       (output) REAL
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      REAL               EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = SLAMCH( 'S' )
         EPS = SLAMCH( 'E' )
         SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( SLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of SLARTG
*
      END
      SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL               CFROM, CTO
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) REAL
*  CTO     (input) REAL
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,M)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of SLASCL
*
      END
      SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SLASET initializes an m-by-n matrix A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of A is not changed.
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of A is not changed.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) REAL
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) REAL
*          The constant to which the diagonal elements are to be set.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On exit, the leading m-by-n submatrix of A is set as follows:
*
*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*
*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the strictly upper triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the strictly lower triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE
*
*        Set the leading m-by-n submatrix to ALPHA.
*
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Set the first min(M,N) diagonal elements to BETA.
*
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
*
      RETURN
*
*     End of SLASET
*
      END
      SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( * ), S( * )
*     ..
*
*  Purpose
*  =======
*
*  SLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n real matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  This version vectorises across rows of the array A when SIDE = 'L'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) REAL arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SLASR
*
      END
      SUBROUTINE SLASRT( ID, N, D, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      REAL               D( * )
*     ..
*
*  Purpose
*  =======
*
*  Sort the numbers in D in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the array to be sorted.
*          On exit, D has been sorted into increasing order
*          (D(1) <= ... <= D(N) ) or into decreasing order
*          (D(1) >= ... >= D(N) ), depending on ID.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      REAL               D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of SLASRT
*
      END
      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      REAL               X( * )
*     ..
*
*  Purpose
*  =======
*
*  SLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) REAL array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) REAL
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) REAL
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      REAL               ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of SLASSQ
*
      END
      SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  SSTEDC computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*  The eigenvectors of a full or band real symmetric matrix can also be
*  found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See SLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original dense symmetric
*                  matrix also.  On entry, Z contains the orthogonal
*                  matrix used to reduce the original matrix to
*                  tridiagonal form.
*
*  N       (input) INTEGER
*          The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) REAL array, dimension (LDZ,N)
*          On entry, if COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If  COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If eigenvectors are desired, then LDZ >= max(1,N).
*
*  WORK    (workspace/output) REAL array,
*                                         dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
*          If COMPZ = 'V' and N > 1 then LWORK must be at least
*                         ( 1 + 3*N + 2*N*lg N + 3*N**2 ),
*                         where lg( N ) = smallest integer k such
*                         that 2**k >= N.
*          If COMPZ = 'I' and N > 1 then LWORK must be at least
*                         ( 1 + 4*N + N**2 ).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
*          If COMPZ = 'V' and N > 1 then LIWORK must be at least
*                         ( 6 + 6*N + 5*N*lg N ).
*          If COMPZ = 'I' and N > 1 then LIWORK must be at least
*                         ( 3 + 5*N ).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            END, I, ICOMPZ, II, J, K, LGN, LIWMIN, LWMIN,
     $                   M, SMLSIZ, START, STOREZ, STRTRW
      REAL               EPS, ORGNRM, P, TINY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANST
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLACPY, SLAED0, SLASCL, SLASET, SLASRT,
     $                   SSTEQR, SSTERF, SSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( N.LE.1 .OR. ICOMPZ.LE.0 ) THEN
         LIWMIN = 1
         LWMIN = 1
      ELSE
         LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( ICOMPZ.EQ.1 ) THEN
            LWMIN = 1 + 3*N + 2*N*LGN + 3*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1 + 4*N + N**2
            LIWMIN = 3 + 5*N
         END IF
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -8
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN 
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN 
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 )
     $      Z( 1, 1 ) = ONE
         RETURN 
      END IF
*
      SMLSIZ = ILAENV( 9, 'SSTEDC', ' ', 0, 0, 0, 0 )
*
*     If the following conditional clause is removed, then the routine
*     will use the Divide and Conquer routine to compute only the
*     eigenvalues, which requires (3N + 3N**2) real workspace and
*     (2 + 5N + 2N lg(N)) integer workspace.
*     Since on many architectures SSTERF is much faster than any other
*     algorithm for finding eigenvalues only, it is used here
*     as the default.
*
*     If COMPZ = 'N', use SSTERF to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL SSTERF( N, D, E, INFO )
         RETURN 
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPZ.EQ.0 ) THEN
            CALL SSTERF( N, D, E, INFO )
            RETURN 
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            CALL SSTEQR( 'I', N, D, E, Z, LDZ, WORK, INFO )
            RETURN
         ELSE
            CALL SSTEQR( 'V', N, D, E, Z, LDZ, WORK, INFO )
            RETURN 
         END IF
      END IF
*
*     If COMPZ = 'V', the Z matrix must be stored elsewhere for later
*     use.
*
      IF( ICOMPZ.EQ.1 ) THEN
         STOREZ = 1 + N*N
      ELSE
         STOREZ = 1
      END IF
*
      IF( ICOMPZ.EQ.2 ) THEN
         CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
      END IF
*
*     Scale.
*
      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO )
     $   RETURN
*
      EPS = SLAMCH( 'Epsilon' )
*
      START = 1
*
*     while ( START <= N )
*
   10 CONTINUE
      IF( START.LE.N ) THEN
*
*     Let END be the position of the next subdiagonal entry such that
*     E( END ) <= TINY or END = N if no such subdiagonal exists.  The
*     matrix identified by the elements between START and END
*     constitutes an independent sub-problem.
*
         END = START
   20    CONTINUE
         IF( END.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( END ) ) )*SQRT( ABS( D( END+1 ) ) )
            IF( ABS( E( END ) ).GT.TINY ) THEN
               END = END + 1
               GO TO 20
            END IF
         END IF
*
*        (Sub) Problem determined.  Compute its size and solve it.
*
         M = END - START + 1
         IF( M.EQ.1 ) THEN
            START = END + 1
            GO TO 10
         END IF
         IF( M.GT.SMLSIZ ) THEN
            INFO = SMLSIZ
*
*           Scale.
*
            ORGNRM = SLANST( 'M', M, D( START ), E( START ) )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M,
     $                   INFO )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ),
     $                   M-1, INFO )
*
            IF( ICOMPZ.EQ.1 ) THEN
               STRTRW = 1
            ELSE
               STRTRW = START
            END IF
            CALL SLAED0( ICOMPZ, N, M, D( START ), E( START ),
     $                   Z( STRTRW, START ), LDZ, WORK( 1 ), N,
     $                   WORK( STOREZ ), IWORK, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) +
     $                MOD( INFO, ( M+1 ) ) + START - 1
               RETURN 
            END IF
*
*           Scale back.
*
            CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M,
     $                   INFO )
*
         ELSE
            IF( ICOMPZ.EQ.1 ) THEN
*
*     Since QR won't update a Z matrix which is larger than the
*     length of D, we must solve the sub-problem in a workspace and
*     then multiply back into Z.
*
               CALL SSTEQR( 'I', M, D( START ), E( START ), WORK, M,
     $                      WORK( M*M+1 ), INFO )
               CALL SLACPY( 'A', N, M, Z( 1, START ), LDZ,
     $                      WORK( STOREZ ), N )
               CALL SGEMM( 'N', 'N', N, M, M, ONE, WORK( STOREZ ), LDZ,
     $                     WORK, M, ZERO, Z( 1, START ), LDZ )
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               CALL SSTEQR( 'I', M, D( START ), E( START ),
     $                      Z( START, START ), LDZ, WORK, INFO )
            ELSE
               CALL SSTERF( M, D( START ), E( START ), INFO )
            END IF
            IF( INFO.NE.0 ) THEN
               INFO = START*( N+1 ) + END
               RETURN 
            END IF
         END IF
*
         START = END + 1
         GO TO 10
      END IF
*
*     endwhile
*
*     If the problem split any number of times, then the eigenvalues
*     will not be properly ordered.  Here we permute the eigenvalues
*     (and the associated eigenvectors) into ascending order.
*
      IF( M.NE.N ) THEN
         IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
            CALL SLASRT( 'I', N, D, INFO )
*
         ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
            DO 40 II = 2, N
               I = II - 1
               K = I
               P = D( I )
               DO 30 J = II, N
                  IF( D( J ).LT.P ) THEN
                     K = J
                     P = D( J )
                  END IF
   30          CONTINUE
               IF( K.NE.I ) THEN
                  D( K ) = D( I )
                  D( I ) = P
                  CALL SSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
               END IF
   40       CONTINUE
         END IF
      END IF
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of SSTEDC
*
      END
      SUBROUTINE SSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  SSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band symmetric matrix can also be found
*  if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to
*  tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  symmetric matrix.  On entry, Z must contain the
*                  orthogonal matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) REAL array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) REAL array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is orthogonally similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   THREE = 3.0E0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      REAL               ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANST, SLAPY2
      EXTERNAL           LSAME, SLAMCH, SLANST, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAE2, SLAEV2, SLARTG, SLASCL, SLASET, SLASR,
     $                   SLASRT, SSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = SLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = SLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL SLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL SLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL SLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL SLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL SLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL SLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL SLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL SLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL SLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
*
  190 CONTINUE
      RETURN
*
*     End of SSTEQR
*
      END
      SUBROUTINE SSTERF( N, D, E, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  SSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm failed to find all of the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   THREE = 3.0E0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
     $                   NMAXIT
      REAL               ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANST, SLAPY2
      EXTERNAL           SLAMCH, SLANST, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAE2, SLASCL, SLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'SSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = SLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*
     $       SQRT( ABS( D( M+1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = SLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use SLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL SLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = SLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $         GO TO 120
  110    CONTINUE
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use SLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL SLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = SLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
*
      END IF
*
*     Undo scaling if necessary
*
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      GO TO 180
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      CALL SLASRT( 'I', N, D, INFO )
*
  180 CONTINUE
      RETURN
*
*     End of SSTERF
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      REAL             FUNCTION SLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992 
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  SLAMCH determines single precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by SLAMCH:
*          = 'E' or 'e',   SLAMCH := eps
*          = 'S' or 's ,   SLAMCH := sfmin
*          = 'B' or 'b',   SLAMCH := base
*          = 'P' or 'p',   SLAMCH := eps*base
*          = 'N' or 'n',   SLAMCH := t
*          = 'R' or 'r',   SLAMCH := rnd
*          = 'M' or 'm',   SLAMCH := emin
*          = 'U' or 'u',   SLAMCH := rmin
*          = 'L' or 'l',   SLAMCH := emax
*          = 'O' or 'o',   SLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      SLAMCH = RMACH
      RETURN
*
*     End of SLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  SLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  SLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = SLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of SLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  SLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) REAL
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) REAL
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) REAL
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  SLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine SLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call SLAMC5 to compute EMAX and RMAX.
*
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of SLAMC2
*
      END
*
************************************************************************
*
      REAL             FUNCTION SLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B
*     ..
*
*  Purpose
*  =======
*
*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) REAL
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      SLAMC3 = A + B
*
      RETURN
*
*     End of SLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
*     ..
*
*  Purpose
*  =======
*
*  SLAMC4 is a service routine for SLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) REAL
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of SLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
*     ..
*
*  Purpose
*  =======
*
*  SLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) REAL
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of SLAMC5
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      subroutine sscal(n,sa,sx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real sa,sx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end
c***DECK CHEMV
      SUBROUTINE CHEMV (UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
C***BEGIN PROLOGUE  CHEMV
C***PURPOSE  Multiply a complex vector by a complex Hermitian matrix.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B4
C***TYPE      COMPLEX (SHEMV-S, DHEMV-D, CHEMV-C)
C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J. J., (ANL)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  CHEMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n hermitian matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the hermitian matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the hermitian matrix and the strictly
C           upper triangular part of A is not referenced.
C           Note that the imaginary parts of the diagonal elements need
C           not be set and are assumed to be zero.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - COMPLEX          array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX          array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
C                 Hanson, R. J.  An extended set of Fortran basic linear
C                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
C                 pp. 1-17, March 1988.
C***ROUTINES CALLED  LSAME, XERBLA
C***REVISION HISTORY  (YYMMDD)
C   861022  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  CHEMV
C     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
C     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C     .. Local Scalars ..
      COMPLEX            TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
C***FIRST EXECUTABLE STATEMENT  CHEMV
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + CONJG( A( I, J ) )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*REAL( A( J, J ) ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + CONJG( A( I, J ) )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*REAL( A( J, J ) ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J ) + TEMP1*REAL( A( J, J ) )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + CONJG( A( I, J ) )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY ) + TEMP1*REAL( A( J, J ) )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + CONJG( A( I, J ) )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of CHEMV .
C
      END
      subroutine  cscal(n,ca,cx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, linpack,  3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex ca,cx(*)
      integer i,incx,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
      complex function cdotc(n,cx,incx,cy,incy)
c
c     forms the dot product of two vectors, conjugating the first
c     vector.
c     jack dongarra, linpack,  3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n
c
      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotc = ctemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
   30 continue
      cdotc = ctemp
      return
      end
      subroutine caxpy(n,ca,cx,incx,cy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1),ca
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if (abs(real(ca)) + abs(aimag(ca)) .eq. 0.0 ) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
      subroutine  ccopy(n,cx,incx,cy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end
      complex function cdotu(n,cx,incx,cy,incy)
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      ctemp = (0.0,0.0)
      cdotu = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = ctemp + cx(ix)*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotu = ctemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = ctemp + cx(i)*cy(i)
   30 continue
      cdotu = ctemp
      return
      end
      real function cmach(job)
      integer job
c
c     smach computes machine parameters of floating point
c     arithmetic for use in testing only.  not required by
c     linpack proper.
c
c     if trouble with automatic computation of these quantities,
c     they can be set by direct assignment statements.
c     assume the computer has
c
c        b = base of arithmetic
c        t = number of base  b  digits
c        l = smallest possible exponent
c        u = largest possible exponent
c
c     then
c
c        eps = b**(1-t)
c        tiny = 100.0*b**(-l+t)
c        huge = 0.01*b**(u-t)
c
c     dmach same as smach except t, l, u apply to
c     double precision.
c
c     cmach same as smach except if complex division
c     is done by
c
c        1/(x+i*y) = (x-i*y)/(x**2+y**2)
c
c     then
c
c        tiny = sqrt(tiny)
c        huge = sqrt(huge)
c
c
c     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
c
c
      real eps,tiny,huge,s
c
      eps = 1.0
   10 eps = eps/2.0
      s = 1.0 + eps
      if (s .gt. 1.0) go to 10
      eps = 2.0*eps
      cmach =eps
      if( job .eq. 1) return
c
      s = 1.0
   20 tiny = s
      s = s/16.0
      if (s*1.0 .ne. 0.0) go to 20
      tiny = (tiny/eps)*100.
      s = real((1.0,0.0)/cmplx(tiny,0.0))
      if (s .ne. 1.0/tiny) tiny = sqrt(tiny)
      huge = 1.0/tiny
      if (job .eq. 1) cmach = eps
      if (job .eq. 2) cmach = tiny
      if (job .eq. 3) cmach = huge
      return
      end
      subroutine crotg(ca,cb,c,s)
      complex ca,cb,s
      real c
      real norm,scale
      complex alpha
      if (cabs(ca) .ne. 0.) go to 10
         c = 0.
         s = (1.,0.)
         ca = cb
         go to 20
   10 continue
         scale = cabs(ca) + cabs(cb)
         norm = scale * sqrt((cabs(ca/scale))**2 + (cabs(cb/scale))**2)
         alpha = ca /cabs(ca)
         c = cabs(ca) / norm
         s = alpha * conjg(cb) / norm
         ca = alpha * norm
   20 continue
      return
      end
      subroutine  csscal(n,sa,cx,incx)
c
c     scales a complex vector by a real constant.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1)
      real sa
      integer i,incx,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
   30 continue
      return
      end
      subroutine  cswap (n,cx,incx,cy,incy)
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
   30 continue
      return
      end
      integer function icamax(n,cx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      icamax = 0
      if( n .lt. 1 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end
      real function scasum(n,cx,incx)
c
c     takes the sum of the absolute values of a complex vector and
c     returns a single precision result.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1)
      real stemp
      integer i,incx,n,nincx
c
      scasum = 0.0e0
      stemp = 0.0e0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
   10 continue
      scasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
   30 continue
      scasum = stemp
      return
      end
      real function scnrm2( n, cx, incx)
      logical imag, scale
      integer          next
      real         cutlo, cuthi, hitest, sum, xmax, absx, zero, one
      complex      cx(1)
      data         zero, one /0.0e0, 1.0e0/
c
c     unitary norm of the complex n-vector stored in cx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson , 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 4.441e-16,  1.304e19 /
c
      if(n .gt. 0) go to 10
         scnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      do 210 i=1,nn,incx
         absx = abs(real(cx(i)))
         imag = .false.
         go to next,(30, 50, 70, 90, 110)
   30 if( absx .gt. cutlo) go to 85
      assign 50 to next
      scale = .false.
c
c                        phase 1.  sum is zero
c
   50 if( absx .eq. zero) go to 200
      if( absx .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 assign 110 to next
      sum = (sum / absx) / absx
  105 scale = .true.
      xmax = absx
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( absx .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( absx .le. xmax ) go to 115
         sum = one + sum * (xmax / absx)**2
         xmax = absx
         go to 200
c
  115 sum = sum + (absx/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
   85 assign 90 to next
      scale = .false.
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
      hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
   90 if(absx .ge. hitest) go to 100
         sum = sum + absx**2
  200 continue
c                  control selection of real and imaginary parts.
c
      if(imag) go to 210
         absx = abs(aimag(cx(i)))
         imag = .true.
      go to next,(  50, 70, 90, 110 )
c
  210 continue
c
c              end of main loop.
c              compute square root and adjust for scaling.
c
      scnrm2 = sqrt(sum)
      if(scale) scnrm2 = scnrm2 * xmax
  300 continue
      return
      end
      subroutine srotg(sa,sb,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      real sa,sb,c,s,roe,scale,r,z
c
      roe = sb
      if( abs(sa) .gt. abs(sb) ) roe = sa
      scale = abs(sa) + abs(sb)
      if( scale .ne. 0.0 ) go to 10
         c = 1.0
         s = 0.0
         r = 0.0
         go to 20
   10 r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
      r = sign(1.0,roe)*r
      c = sa/r
      s = sb/r
   20 z = 1.0
      if( abs(sa) .gt. abs(sb) ) z = s
      if( abs(sb) .ge. abs(sa) .and. c .ne. 0.0 ) z = 1.0/c
      sa = r
      sb = z
      return
      end
      SUBROUTINE CHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  CHER2  performs the hermitian rank 2 operation
*
*     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an n
*  by n hermitian matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX          array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX          array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the hermitian matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the hermitian matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set, they are assumed to be zero, and on exit they
*           are set to zero.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     .. Local Scalars ..
      COMPLEX            TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHER2 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in the upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*CONJG( Y( J ) )
                  TEMP2 = CONJG( ALPHA*X( J ) )
                  DO 10, I = 1, J - 1
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
                  A( J, J ) = REAL( A( J, J ) ) +
     $                        REAL( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  A( J, J ) = REAL( A( J, J ) )
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*CONJG( Y( JY ) )
                  TEMP2 = CONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J - 1
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
                  A( J, J ) = REAL( A( J, J ) ) +
     $                        REAL( X( JX )*TEMP1 + Y( JY )*TEMP2 )
               ELSE
                  A( J, J ) = REAL( A( J, J ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in the lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1     = ALPHA*CONJG( Y( J ) )
                  TEMP2     = CONJG( ALPHA*X( J ) )
                  A( J, J ) = REAL( A( J, J ) ) +
     $                        REAL( X( J )*TEMP1 + Y( J )*TEMP2 )
                  DO 50, I = J + 1, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               ELSE
                  A( J, J ) = REAL( A( J, J ) )
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1     = ALPHA*CONJG( Y( JY ) )
                  TEMP2     = CONJG( ALPHA*X( JX ) )
                  A( J, J ) = REAL( A( J, J ) ) +
     $                        REAL( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX        = JX
                  IY        = JY
                  DO 70, I = J + 1, N
                     IX        = IX        + INCX
                     IY        = IY        + INCY
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  A( J, J ) = REAL( A( J, J ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CHER2 .
*
      END
      SUBROUTINE CHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL               BETA
      COMPLEX            ALPHA
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  CHER2K  performs one of the hermitian rank 2k operations
*
*     C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
*
*  or
*
*     C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C,
*
*  where  alpha and beta  are scalars with  beta  real,  C is an  n by n
*  hermitian matrix and  A and B  are  n by k matrices in the first case
*  and  k by n  matrices in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'    C := alpha*A*conjg( B' )          +
*                                         conjg( alpha )*B*conjg( A' ) +
*                                         beta*C.
*
*              TRANS = 'C' or 'c'    C := alpha*conjg( A' )*B          +
*                                         conjg( alpha )*conjg( B' )*A +
*                                         beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
*           matrices  A and B.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - COMPLEX          array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  hermitian matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  hermitian matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set,  they are assumed to be zero,  and on exit they
*           are set to zero.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.
*     Ed Anderson, Cray Research Inc.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX            TEMP1, TEMP2
*     .. Parameters ..
      REAL               ONE
      PARAMETER        ( ONE  = 1.0E+0 )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHER2K', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.REAL( ZERO ) )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.REAL( ZERO ) )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  C( J, J ) = BETA*REAL( C( J, J ) )
                  DO 70, I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) +
*                   C.
*
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.REAL( ZERO ) )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
               ELSE
                  C( J, J ) = REAL( C( J, J ) )
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*CONJG( B( J, L ) )
                     TEMP2 = CONJG( ALPHA*A( J, L ) )
                     DO 110, I = 1, J - 1
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  110                CONTINUE
                     C( J, J ) = REAL( C( J, J ) )         +
     $                           REAL( A( J, L )*TEMP1 +
     $                                 B( J, L )*TEMP2   )
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.REAL( ZERO ) )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
               ELSE
                  C( J, J ) = REAL( C( J, J ) )
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*CONJG( B( J, L ) )
                     TEMP2 = CONJG( ALPHA*A( J, L ) )
                     DO 160, I = J + 1, N
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  160                CONTINUE
                     C( J, J ) = REAL( C( J, J ) )         +
     $                           REAL( A( J, L )*TEMP1 +
     $                                 B( J, L )*TEMP2   )
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
*
*        Form  C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A +
*                   C.
*
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + CONJG( A( L, I ) )*B( L, J )
                     TEMP2 = TEMP2 + CONJG( B( L, I ) )*A( L, J )
  190             CONTINUE
                  IF( I.EQ.J )THEN
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( J, J ) = REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     ELSE
                        C( J, J ) = BETA*REAL( C( J, J ) )         +
     $                              REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     END IF
                  ELSE
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( I, J ) = ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     ELSE
                        C( I, J ) = BETA *C( I, J ) +
     $                              ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     END IF
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + CONJG( A( L, I ) )*B( L, J )
                     TEMP2 = TEMP2 + CONJG( B( L, I ) )*A( L, J )
  220             CONTINUE
                  IF( I.EQ.J )THEN
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( J, J ) = REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     ELSE
                        C( J, J ) = BETA*REAL( C( J, J ) )         +
     $                              REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     END IF
                  ELSE
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( I, J ) = ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     ELSE
                        C( I, J ) = BETA *C( I, J ) +
     $                              ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     END IF
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CHER2K.
*
      END
      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL               TEMP
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SGEMM .
*
      END
*DECK SCOPY
      SUBROUTINE SCOPY (N, SX, INCX, SY, INCY)
C***BEGIN PROLOGUE  SCOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      SINGLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  copy of vector SX (unchanged if N .LE. 0)
C
C     Copy single precision SX to single precision SY.
C     For I = 0 to N-1, copy  SX(LX+I*INCX) to SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SCOPY
      REAL SX(*), SY(*)
C***FIRST EXECUTABLE STATEMENT  SCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I+1) = SX(I+1)
        SY(I+2) = SX(I+2)
        SY(I+3) = SX(I+3)
        SY(I+4) = SX(I+4)
        SY(I+5) = SX(I+5)
        SY(I+6) = SX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        SY(I) = SX(I)
   70 CONTINUE
      RETURN
      END
      integer function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real sx(*),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
      SUBROUTINE CGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  CGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
*
*     y := alpha*conjg( A' )*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - COMPLEX          array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX         .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX          array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
      NOCONJ = LSAME( TRANS, 'T' )
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CGEMV .
*
      END
      SUBROUTINE CGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  CGERC  performs the rank 1 operation
*
*     A := alpha*x*conjg( y' ) + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX          array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX          array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGERC ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*CONJG( Y( JY ) )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*CONJG( Y( JY ) )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of CGERC .
*
      END
      SUBROUTINE CTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX          array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
*
            IF( UPPER )THEN
               DO 120, J = 1, N
                  DO 110, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 90, K = 1, I - 1
                           TEMP = TEMP + A( K, I )*B( K, J )
   90                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*CONJG( A( I, I ) )
                        DO 100, K = 1, I - 1
                           TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  100                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  110             CONTINUE
  120          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = 1, M
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 130, K = I + 1, M
                           TEMP = TEMP + A( K, I )*B( K, J )
  130                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*CONJG( A( I, I ) )
                        DO 140, K = I + 1, M
                           TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  140                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 200, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 170, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  170             CONTINUE
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
  200          CONTINUE
            ELSE
               DO 240, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 210, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  210             CONTINUE
                  DO 230, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 220, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  220                   CONTINUE
                     END IF
  230             CONTINUE
  240          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
*
            IF( UPPER )THEN
               DO 280, K = 1, N
                  DO 260, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*CONJG( A( J, K ) )
                        END IF
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*CONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
  280          CONTINUE
            ELSE
               DO 320, K = N, 1, -1
                  DO 300, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*CONJG( A( J, K ) )
                        END IF
                        DO 290, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  290                   CONTINUE
                     END IF
  300             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*CONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 310, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  310                CONTINUE
                  END IF
  320          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of CTRMM .
*
      END
      SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  CGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX         .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX          array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX            TEMP
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.CONJA                ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.CONJB                ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
*
*           Form  C := alpha*conjg( A' )*B + beta*C.
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
*
*           Form  C := alpha*A*conjg( B' ) + beta*C.
*
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
*
*           Form  C := alpha*A*B'          + beta*C
*
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
*
*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
*
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
*
*           Form  C := alpha*conjg( A' )*B' + beta*C
*
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
*
*           Form  C := alpha*A'*conjg( B' ) + beta*C
*
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CGEMM .
*
      END
*DECK CTRMV
      SUBROUTINE CTRMV (UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
C***BEGIN PROLOGUE  CTRMV
C***PURPOSE  Multiply a complex vector by a complex triangular matrix.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B4
C***TYPE      COMPLEX (STRMV-S, DTRMV-D, CTRMV-C)
C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J. J., (ANL)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  CTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,   or   x := conjg( A')*x,
C
C  where x is an n element vector and  A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := conjg( A' )*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - COMPLEX          array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           transformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
C                 Hanson, R. J.  An extended set of Fortran basic linear
C                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
C                 pp. 1-17, March 1988.
C***ROUTINES CALLED  LSAME, XERBLA
C***REVISION HISTORY  (YYMMDD)
C   861022  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  CTRMV
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOCONJ, NOUNIT
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
C***FIRST EXECUTABLE STATEMENT  CTRMV
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOCONJ = LSAME( TRANS, 'T' )
      NOUNIT = LSAME( DIAG , 'N' )
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( LSAME( TRANS, 'N' ) )THEN
C
C        Form  x := A*x.
C
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := A'*x  or  x := conjg( A' )*x.
C
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 110, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 90, I = J - 1, 1, -1
                        TEMP = TEMP + A( I, J )*X( I )
   90                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 100, I = J - 1, 1, -1
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100                CONTINUE
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 140, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 120, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  120                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 130, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  130                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = 1, N
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 150, I = J + 1, N
                        TEMP = TEMP + A( I, J )*X( I )
  150                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 160, I = J + 1, N
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  160                CONTINUE
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               JX = KX
               DO 200, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 180, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  180                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 190, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  190                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  200          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of CTRMV .
C
      END
      REAL             FUNCTION SNRM2 ( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      REAL                              X( * )
*     ..
*
*  SNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     SNRM2 := sqrt( x'*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to SLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      REAL                  ONE         , ZERO
      PARAMETER           ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     .. Local Scalars ..
      INTEGER               IX
      REAL                  ABSXI, NORM, SCALE, SSQ
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      SNRM2 = NORM
      RETURN
*
*     End of SNRM2.
*
      END
      SUBROUTINE SROT(N,SX,INCX,SY,INCY,SC,SS)
C***BEGIN PROLOGUE  SROT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A8
C***KEYWORDS  BLAS,GIVENS ROTATION,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Apply s.p. Givens rotation
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C       SC  element of rotation matrix
C       SS  element of rotation matrix
C
C     --Output--
C       SX  rotated vector SX (unchanged if N .LE. 0)
C       SY  rotated vector SY (unchanged if N .LE. 0)
C
C     Multiply the 2 x 2 matrix  ( SC SS) times the 2 x N matrix (SX**T)
C                                (-SS SC)                        (SY**T)
C     where **T indicates transpose.  The elements of SX are in
C     SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = (-INCX)*N, and similarly for SY using LY and INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SROT
C
      REAL             SX,SY,SC,SS,ZERO,ONE,W,Z
      DIMENSION SX(1),SY(1)
      DATA ZERO,ONE/0.E0,1.E0/
C***FIRST EXECUTABLE STATEMENT  SROT
      IF(N .LE. 0 .OR. (SS .EQ. ZERO .AND. SC .EQ. ONE)) GO TO 40
      IF(.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
           NSTEPS=INCX*N
           DO 10 I=1,NSTEPS,INCX
                W=SX(I)
                Z=SY(I)
                SX(I)=SC*W+SS*Z
                SY(I)=-SS*W+SC*Z
   10           CONTINUE
           GO TO 40
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF(INCX .LT. 0) KX=1-(N-1)*INCX
           IF(INCY .LT. 0) KY=1-(N-1)*INCY
C
           DO 30 I=1,N
                W=SX(KX)
                Z=SY(KY)
                SX(KX)=SC*W+SS*Z
                SY(KY)=-SS*W+SC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SGEMV .
*
      END
      subroutine sswap (n,sx,incx,sy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
