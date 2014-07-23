        SUBROUTINE FITLINE(X,Y,I,A,B)
C
C   THIS SUBROUTINE WILL CALCULATE
C   THE BEST STRAIGHT LINE THROUGH THE DATA USING A LEAST SQUARES
C   METHOD.THE SUBROUTINE THEN CALCULATES THE ERRORS ON THE OUPUT
C   PARAMETERS A & B.
C
C   Y = AX + B   ERRORS : SA & SB
C
C
	IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 X(100),Y(100),SIGX,SIGXY,SIGX2,SIGY2,SN2,SA,SB,A,B
        REAL*8 N
C
C   TEST TO SEE IF ONLY ONE POINT, IF SO THEN SET SLOPE TO ZERO
C   AND INTERCEPT AS VALUE OF POINT
C
      IF (I.EQ.1) THEN
         A = 0.0
         B = Y(1)
         RETURN
      ENDIF
        SIGX = 0.0
        SIGY = 0.0
        SIGXY = 0.0
        SIGX2 = 0.0
        SIGY2 = 0.0
C
C    NOW CALCULATE SUMMATIONS
C
        DO 55 J=1,I
        SIGX=SIGX+X(J)
        SIGY=SIGY+Y(J)
        SIGXY=SIGXY+X(J)*Y(J)
        SIGX2=SIGX2+X(J)**2
        SIGY2=SIGY2+Y(J)**2
 55     CONTINUE
C
C    CALCULATE A & B
C
        N=FLOAT(I)
        A=(N*SIGXY-SIGX*SIGY)/(N*SIGX2-SIGX**2)
        B=(SIGX2*SIGY-SIGX*SIGXY)/(N*SIGX2-SIGX**2)
C
C    NOW THE ERRORS
C
        DELI2 = 0.0
        DO 500 J = 1,I
           DELI2 = DELI2 + ( Y(J) - A*X(J) - B )**2
 500    CONTINUE
C       SN2 = DELI2/(N-2.0)
C       SA = ( SN2*N/(N*SIGX2-SIGX**2) )**0.5
C       SB = ( SN2*SIGX2/(N*SIGX2-SIGX**2) )**0.5C
        RETURN
        END
        SUBROUTINE PARABOL(X,Y,I,A,B,C)
C
C   THIS SUBROUTINE WILL CALCULATE
C   THE BEST PARABOLA THROUGH THE DATA USING A LEAST SQUARES
C   METHOD.
C   PARAMETERS A B & C.
C
C   Y = AX**2 + BX + C   ERRORS : SA SB & SC
C
C
	IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 N,X(100),Y(100)
C
C  FIRST SUBRTACT OFF X(1) IN ORDER TO RETAIN ACCURACY
C
	XSTART = X(1)
	DO 550 J=1,I
	   X(J) = X(J) - XSTART
 550	CONTINUE
C
C  NOW NORMALISE
C
        SGX =    0.0
        SGY =    0.0
        SGX2 =   0.0
        SGX3 =   0.0
        SGX4 =   0.0
        SGXY =   0.0
        SGX2Y = 0.0
        SGX4Y2 = 0.0
        SGX3Y2 = 0.0
        SGX2Y2 = 0.0
        SGXY2  = 0.0
        SGY2   = 0.0
        DO 55 J=1,I
        SGX=SGX+X(J)
        SGY=SGY+Y(J)
        SGXY=SGXY+X(J)*Y(J)
        SGX2=SGX2+X(J)**2
        SGX2Y=SGX2Y+X(J)**2*Y(J)
        SGY2=SGY2+Y(J)**2
        SGX3=SGX3+X(J)**3
        SGX4=SGX4+X(J)**4
        SGX4Y2=SGX4Y2+(X(J)**4)*Y(J)**2
        SGX3Y2=SGX3Y2+(X(J)**3)*Y(J)**2
        SGX2Y2=SGX2Y2+(X(J)**2)*Y(J)**2
        SGXY2=SGXY2+X(J)*Y(J)**2
 55     CONTINUE
C
C    CALCULATE A B & C
C
        N=DBLE(I)
        DETERM = SGX4*(SGX2*N - SGX**2) +
     %           SGX3*(SGX2*SGX - N*SGX3) +
     %           SGX2*(SGX3*SGX - SGX2**2)
        A = SGX2Y*(N*SGX2 - SGX**2) +
     %      SGXY*(SGX*SGX2 - N*SGX3) +
     %      SGY*(SGX*SGX3 - SGX2**2)
        A = A/DETERM
        B = SGX2Y*(SGX2*SGX - N*SGX3) +
     %      SGXY*(N*SGX4 - SGX2**2) +
     %      SGY*(SGX2*SGX3 - SGX4*SGX)
        B = B/DETERM
        C = SGX2Y*(SGX3*SGX - SGX2**2) +
     %      SGXY*(SGX2*SGX3 - SGX*SGX4) +
     %      SGY*(SGX4*SGX2 - SGX3**2)
        C = C/DETERM
C
C  Alter Coeffs. to undo effect of subtracting X(1)
C
	C = C + XSTART**2 - B*XSTART
	B = B - 2.0*XSTART*A
	A = A
        RETURN
        END
      SUBROUTINE INV33M(M,MINV)
C
C   THIS ROUTINE INVERTS A 3 X 3 MATRIX USING THE TRANSPOSE COFACTOR MATRIX
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DETM, M(3,3),MINV(3,3)
      INTEGER I,J
      MINV(1,1) = M(2,2)*M(3,3) - M(2,3)*M(3,2)
      MINV(2,2) = M(1,1)*M(3,3) - M(3,1)*M(1,3)
      MINV(3,3) = M(1,1)*M(2,2) - M(2,1)*M(1,2)
      MINV(2,1) = M(3,1)*M(2,3) - M(2,1)*M(3,3)
      MINV(3,1) = M(2,1)*M(3,2) - M(2,2)*M(3,1)
      MINV(1,2) = M(3,2)*M(1,3) - M(1,2)*M(3,3)
      MINV(3,2) = M(3,1)*M(1,2) - M(1,1)*M(3,2)
      MINV(1,3) = M(1,2)*M(2,3) - M(2,2)*M(1,3)
      MINV(2,3) = M(2,1)*M(1,3) - M(1,1)*M(2,3)
      DETM = M(1,1)*( M(2,2)*M(3,3)-M(3,2)*M(2,3) ) -
     %       M(1,2)*( M(2,1)*M(3,3)-M(3,1)*M(2,3) ) +
     %       M(1,3)*( M(2,1)*M(3,2)-M(2,2)*M(3,1) )
      DO 10 I = 1,3
         DO 100 J = 1,3
            MINV(I,J) = MINV(I,J)/DETM
 100     CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE CL1(K, L, M, N, KLMD, KLM2D, NKLMD, N2D,
     * Q, KODE, TOLER, ITER, X, RES, ERROR, CU, IU, S)
C THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX
C METHOD OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION
C TO A K BY N SYSTEM OF LINEAR EQUATIONS
C             AX=B
C SUBJECT TO L LINEAR EQUALITY CONSTRAINTS
C             CX=D
C AND M LINEAR INEQUALITY CONSTRAINTS
C             EX .LE. F  .
C DESCRIPTION OF THE PARAMETERS
C K     NUMBER OF ROWS OF THE MATRIX A (K .GE. 1).
C L     NUMBER OF ROWS OF THE MATRIX C (L .GE. 0).
C M     NUMBER OF ROWS OF THE MATRIX E (M .GE. 0).
C N     NUMBER OF COLUMNS OF THE MATRICES A,C,E (N .GE. 1)
C KLMD   SET TO AT LEAST K+L+M FOR ADJUSTABLE DIMENSIONS.
C KLM2D  SET TO AT LEAST K+L+M+2 FOR ADJUSTABLE DIMENSIONS.
C NKLMD  SET TO AT LEAST N+K+L+M FOR ADJUSTABLE DIMENSONS.
C N2D    SET TO AT LEAST N+2 FOR ADJUSTABLE DIMENSIONS.
C Q      TWO DIMENSIONAL REAL*8 ARRAY WITH KLM2D ROWS AND
C        AT LEAST N2D COLUMNS.
C        ON ENTRY THE MATRICES A,C AND E, AND THE VECTORS
C        B,D AND F MUST BE STORED IN THE FIRST K+L+M ROWS
C        AND N+1 COLUMNS OF Q AS FOLLOWS
C             A B
C         Q = C D
C             E F
C        THESE VALUES ARE DESTROYED BY THE SUBROUTINE.
C KODE   A CODE USED ON ENTRY TO , AND EXIT
C        FROM, THE SUBROUTINE.
C        ON ENTRY, THIS SHOULD NORMALLY BE SET TO 0.
C        HOWEVER, IF CERTAIN NONNEGATIVITY CONSTRAINTS
C        ARE TO BE INCLUDED IMPLICITLY, RATHER THAN
C        EXPLICITLY IN THE CONSTRAINTS EX .LE. F, THEN KODE
C        SHOULD BE SET TO 1, AND THE NONNEGATIVITY
C        CONSTRAINTS INCLUDED IN THE ARRAYS X AND
C        RES (SEE BELOW).
C        ON EXIT, KODE HAS ONE OF THE
C        FOLLOWING VALUES
C             0- OPTIMAL SOLUTION FOUND,
C             1- NO FEASIBLE SOLUTION TO THE
C                CONSTRAINTS
C             2- CALCULATIONS TERMINATED
C                PREMATURELY DUE TO ROUNDING ERRORS,
C             3- MAXIMUM NUMBER OF ITERATIONS REACHED.
C TOLER  A SMALL POSITIVE TOLERANCE. EMPIRICAL
C        EVIDENCE SUGGESTS TOLER = 10**(-D*2/3),
C        WHERE D REPRESENTS THE NUMBER OF DECIMAL
C        DIGITS OF ACCURACY AVAILABLE, ESSENTIALLY,
C        THE SUBROUTINE CANNOT DISTINGUISH BETWEEN ZERO
C        AND ANY QUANTITY WHOSE MAGNITUDE DOES NOT EXCEED
C        TOLER, IN PARTICULAR, IT WILL NOT PIVOT ON ANY
C        NUMBER WHOSE MAGNITUDE DOES NOT EXCEED TOLER.
C ITER   ON ENTRY ITER MUST CONTAIN AN UPPER BOUND ON
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C        A SUGGESTED VALUE IS 10*(K+L+M). ON EXIT ITER
C        GIVES THE NUMBER OF SIMPLEX ITERATIONS.
C X      ONE DIMENSIONAL REAL*8 ARRAY OF SIZE AT LEAST N2D.
C        ON EXIT THIS ARRAY CONTAINS A
C        SOLUTION TO THE L1 PROBLEM. IF KODE=1
C        ON ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE
C        SIMPLE NONNEGATIVITY CONSTRAINTS ON THE
C        VARIABLES. THE VALUES -1, 0, OR 1
C        FOR X(J) INDICATE THAT THE J-TH VARIABLE
C        IS RESTRICTED TO BE .LE. 0, UNRESTRICTED,
C        OR .GE. 0 RESPECTIVELY.
C RES    ONE DIMENSIONAL REAL*8 ARRAY OF SIZE AT LEAST KLMD.
C        ON EXIT THIS CONTAINS THE RESIDUALS B-AX
C        IN THE FIRST K COMPONENTS, D-CX IN THE
C        NEXT L COMPONENTS (THESE WILL BE =0), AND
C        F-EX IN THE NEXT M COMPONENTS. IF KODE=1 ON
C        ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE SIMPLE
C        NONNEGATIVITY CONSTRAINTS ON THE RESIDUALS
C        B-AX. THE VALUES -1, 0, OR 1 FOR RES(I)
C        INDICATE THAT THE I-TH RESIDUAL (1 .LE. I .LE. K) IS
C        RESTRICTED TO BE .LE. 0, UNRESTRICTED, OR .GE. 0
C        RESPECTIVELY.
C ERROR  ON EXIT, THIS GIVES THE MINIMUM SUM OF
C        ABSOLUTE VALUES OF THE RESIDUALS.
C CU     A TWO DIMENSIONAL REAL*8 ARRAY WITH TWO ROWS AND
C        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.
C IU     A TWO DIMENSIONAL INTEGER ARRAY WITH TWO ROWS AND
C        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.
C S      INTEGER ARRAY OF SIZE AT LEAST KLMD, USED FOR
C        WORKSPACE.
C IF YOUR FORTRAN COMPILER PERMITS A SINGLE COLUMN OF A TWO
C DIMENSIONAL ARRAY TO BE PASSED TO A ONE DIMENSIONAL ARRAY
C THROUGH A SUBROUTINE CALL, CONSIDERABLE SAVINGS IN
C EXECUTION TIME MAY BE ACHIEVED THROUGH THE USE OF THE
C FOLLOWING SUBROUTINE, WHICH OPERATES ON COLUMN VECTORS.
C     SUBROUTINE COL(V1, V2, XMLT, NOTROW, K)
C THIS SUBROUTINE ADDS TO THE VECTOR V1 A MULTIPLE OF THE
C VECTOR V2 (ELEMENTS 1 THROUGH K EXCLUDING NOTROW).
C     DIMENSION V1(K), V2(K)
C     KEND = NOTROW - 1
C     KSTART = NOTROW + 1
C     IF (KEND .LT. 1) GO TO 20
C     DO 10 I=1, KEND
C        V1(I) = V1(I) + XMLT*V2(I)
C 10  CONTINUE
C     IF(KSTART .GT. K) GO TO 40
C 20  DO 30 I=KSTART,K
C        V1(I) = V1(I) + XMLT*V2(I)
C 30  CONTINUE
C 40  RETURN
C     END
C SEE COMMENTS FOLLOWING STATEMENT LABELLED 440 FOR
C INSTRUCTIONS ON THE IMPLEMENTATION OF THIS MODIFICATION.
C
C THIS CODE WAS OBTAINED FROM BILL JEFFERYS.  HE OR ANAND
C SIVARAMAKRISHNAN SHOULD BE PESTERED ABOUT HOW IT WORKS.
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION SUM
      DOUBLE PRECISION DBLE
      REAL*8 Q, X, Z, CU, SN, ZU, ZV, CUV, RES, XMAX, XMIN,
     * ERROR, PIVOT, TOLER, TPIVOT
      REAL*8 ABS
      INTEGER I, J, K, L, M, N, S, IA, II, IN, IU, JS, KK,
     * NK, N1, N2, JMN, JPN, KLM, NKL, NK1, N2D, IIMN,
     * IOUT, ITER, KLMD, KLM1, KLM2, KODE, NKLM, NKL1,
     * KLM2D, MAXIT, NKLMD, IPHASE, KFORCE, IINEG
      INTEGER IABS
      DIMENSION Q(KLM2D,N2D), X(N2D), RES(KLMD),
     * CU(2,NKLMD), IU(2,NKLMD), S(KLMD)
C
C INITIALIZATION.
C
      MAXIT = ITER
      N1 = N + 1
      N2 = N + 2
      NK = N + K
      NK1 = NK + 1
      NKL = NK + L
      NKL1 = NKL + 1
      KLM = K + L + M
      KLM1 = KLM + 1
      KLM2 = KLM + 2
      NKLM = N + KLM
      KFORCE = 1
      ITER = 0
      JS = 1
      IA = 0
C SET UP LABELS IN Q.
      DO 10 J=1,N
         Q(KLM2,J) = J
 10   CONTINUE
      DO 30 I=1,KLM
         Q(I,N2) = N + I
         IF (Q(I,N1) .GE. 0) GOTO 30
         DO 20 J=1,N2
            Q(I,J) = -Q(I,J)
 20      CONTINUE
 30   CONTINUE
C SET UP PHASE 1 COSTS.
      IPHASE = 2
      DO 40 J=1,NKLM
         CU(1,J) = 0.
         CU(2,J) = 0.
         IU(1,J) = 0
         IU(2,J) = 0
 40   CONTINUE
      IF (L .EQ. 0) GO TO 60
      DO 50 J=NK1,NKL
         CU(1,J) = 1.
         CU(2,J) = 1.
         IU(1,J) = 1
         IU(2,J) = 1
 50   CONTINUE
      IPHASE = 1
 60   IF (M .EQ. 0) GO TO 80
      DO 70 J=NKL1,NKLM
         CU(2,J) = 1.
         IU(2,J) = 1
         JMN = J - N
         IF (Q(JMN,N2) .LT. 0.) IPHASE = 1
 70   CONTINUE
 80   IF (KODE .EQ. 0) GO TO 150
      DO 110 J=1,N
         IF (X(J)) 90, 110, 100
 90      CU(1,J) = 1.
         IU(1,J) = 1
         GO TO 110
 100     CU(2,J) = 1.
         IU(2,J) = 1
 110  CONTINUE
      DO 140 J=1,K
         JPN = J + N
         IF (RES(J)) 120, 140, 130
 120     CU(1,JPN) = 1.
         IU(1,JPN) = 1
         IF (Q(J,N2) .GT. 0.0) IPHASE = 1
         GO TO 140
 130     CU(2,JPN) = 1.
         IU(2,JPN) = 1
         IF (Q(J,N2) .LT. 0.0) IPHASE = 1
 140  CONTINUE
 150  IF (IPHASE .EQ. 2) GO TO 500
C COMPUTE THE MARGINAL COSTS.
 160  DO 200 J=JS,N1
         SUM = 0.D0
         DO 190 I=1,KLM
            II = Q(I,N2)
            IF (II .LT. 0) GO TO 170
            Z = CU(1,II)
            GO TO 180
 170        IINEG = -II
            Z = CU(2,IINEG)
 180        SUM = SUM + DBLE(Q(I,J))*DBLE(Z)
 190     CONTINUE
         Q(KLM1,J) = SUM
 200  CONTINUE
      DO 230 J=JS,N
         II = Q(KLM2,J)
         IF (II .LT. 0) GO TO 210
         Z = CU(1,II)
         GO TO 220
 210     IINEG = -II
         Z = CU(2,IINEG)
 220     Q(KLM1,J) = Q(KLM1,J) - Z
 230  CONTINUE
C DETERMINE THE VECTOR TO ENTER THE BASIS.
 240  XMAX = 0.
      IF (JS .GT. N) GO TO 490
      DO 280 J=JS,N
         ZU = Q(KLM1,J)
         II = Q(KLM2,J)
         IF (II .GT. 0) GO TO 250
         II = -II
         ZV = ZU
         ZU = -ZU - CU(1,II) - CU(2,II)
         GO TO 260
 250     ZV = -ZU - CU(1,II) - CU(2,II)
 260     IF (KFORCE .EQ. 1 .AND. II .GT. N) GO TO 280
         IF (IU(1,II) .EQ. 1) GO TO 270
         IF (ZU .LE. XMAX) GO TO 270
         XMAX = ZU
         IN = J
 270     IF(IU(2,II) .EQ. 1) GO TO 280
         IF (ZV .LE. XMAX) GO TO 280
         XMAX = ZV
         IN = J
 280  CONTINUE
      IF (XMAX .LE. TOLER) GO TO 490
      IF (Q(KLM1,IN) .EQ. XMAX) GO TO 300
      DO 290 I=1,KLM2
         Q(I,IN) = -Q(I,IN)
 290  CONTINUE
      Q(KLM1,IN) = XMAX
C DETERMINE THE VECTOR TO LEAVE THE BASIS.
 300  IF (IPHASE .EQ. 1 .OR. IA .EQ. 0) GO TO 330
      XMAX = 0.
      DO 310 I = 1,IA
         Z = ABS(Q(I,IN))
         IF (Z .LE. XMAX) GO TO 310
         XMAX = Z
         IOUT = I
 310  CONTINUE
      IF(XMAX .LE. TOLER) GO TO 330
      DO 320 J = 1,N2
         Z = Q(IA,J)
         Q(IA,J) = Q(IOUT,J)
         Q(IOUT,J) = Z
 320  CONTINUE
      IOUT = IA
      IA = IA - 1
      PIVOT = Q(IOUT,IN)
      GO TO 420
 330  KK = 0
      DO 340 I = 1,KLM
         Z = Q(I,IN)
         IF(Z .LE. TOLER) GO TO 340
         KK = KK + 1
         RES(KK) = Q(I,N1)/Z
         S(KK) = I
 340  CONTINUE
 350  IF (KK .GT. 0) GO TO 360
      KODE = 2
      GO TO 590
 360  XMIN = RES(1)
      IOUT = S(1)
      J = 1
      IF (KK .EQ. 1) GO TO  380
      DO 370 I=2,KK
         IF (RES(I) .GE. XMIN) GO TO 370
         J = I
         XMIN = RES(I)
         IOUT = S(I)
 370  CONTINUE
      RES(J) = RES(KK)
      S(J) = S(KK)
 380  KK = KK - 1
      PIVOT = Q(IOUT,IN)
      II = Q(IOUT,N2)
      IF(IPHASE .EQ. 1) GO TO 400
      IF (II .LT. 0) GO TO 390
      IF (IU(2,II) .EQ. 1) GO TO 420
      GO TO 400
 390  IINEG = -II
      IF (IU(1,IINEG) .EQ. 1) GO TO 420
 400  II = IABS(II)
      CUV = CU(1,II) + CU(2,II)
      IF (Q(KLM1,IN)-PIVOT*CUV .LE. TOLER) GO TO 420
C BYPASS INTERMEDIATE MATRICES.
      DO 410 J=JS,N1
         Z = Q(IOUT,J)
         Q(KLM1,J) = Q(KLM1,J) - Z*CUV
         Q(IOUT,J) = -Z
 410  CONTINUE
      Q(IOUT,N2) = -Q(IOUT,N2)
      GO TO 350
C GAUSS-JORDAN ELIMINATION
 420  IF (ITER .LT. MAXIT) GO TO 430
      KODE = 3
      GO TO 590
 430  ITER = ITER + 1
      DO 440 J = JS,N1
         IF (J .NE. IN) Q(IOUT,J) = Q(IOUT,J)/PIVOT
 440  CONTINUE
C IF PERMITTED, USE SUBROUTINE COL OF THE DESCRIPTION
C SECTION AND REPLACE THE FOLLOWING SEVEN STATEMENTS DOWN
C TO AND INCLUDING STATEMENT NUMBER 460 BY..
C     DO 460 J = JS,N1
C        IF (J .EQ. IN) GO TO 460
C        Z = -Q(IOUT,J)
C        CALL COL(Q(1,J), Q(1,IN), Z, IOUT, KLM1)
C 460 CONTINUE
      DO 460 J = JS,N1
         IF (J .EQ. IN) GO TO 460
         Z = -Q(IOUT,J)
         DO 450 I = 1,KLM1
            IF (I .NE. IOUT)  Q(I,J) = Q(I,J) + Z*Q(I,IN)
 450     CONTINUE
 460  CONTINUE
      TPIVOT = -PIVOT
      DO 470 I = 1,KLM1
         IF(I .NE. IOUT) Q(I,IN) = Q(I,IN)/TPIVOT
 470  CONTINUE
      Q(IOUT,IN) = 1./PIVOT
      Z = Q(IOUT,N2)
      Q(IOUT,N2) = Q(KLM2,IN)
      Q(KLM2,IN) = Z
      II = ABS(Z)
      IF (IU(1,II) .EQ. 0 .OR. IU(2,II) .EQ. 0) GO TO 240
      DO 480 I = 1,KLM2
         Z = Q(I,IN)
         Q(I,IN) = Q(I,JS)
         Q(I,JS) = Z
 480  CONTINUE
      JS = JS + 1
      GO TO 240
C TEST FOR OPTIMALITY.
 490  IF  (KFORCE .EQ. 0) GO TO 580
      IF (IPHASE .EQ. 1 .AND. Q(KLM1,N1) .LE. TOLER) GO TO 500
      KFORCE = 0
      GO TO 240
C SET UP PHASE 2 COSTS
 500  IPHASE = 2
      DO 510 J = 1,NKLM
         CU(1,J)= 0.
         CU(2,J) = 0.
 510  CONTINUE
      DO 520 J = N1,NK
         CU(1,J) = 1.
         CU(2,J) = 1.
 520  CONTINUE
      DO 560 I = 1,KLM
         II = Q(I,N2)
         IF (II .GT. 0) GO TO 530
         II = -II
         IF (IU(2,II) .EQ. 0) GO TO 560
         CU(2,II) = 0.
         GO TO 540
 530     IF (IU(1,II) .EQ. 0) GO TO 560
         CU(1,II) = 0.
 540     IA = IA + 1
         DO 550 J = 1,N2
            Z = Q(IA,J)
            Q(IA,J) = Q(I,J)
            Q(I,J) = Z
 550     CONTINUE
 560  CONTINUE
      GO TO 160
 570  IF (Q(KLM1,N1) .LE. TOLER) GO TO 500
      KODE = 1
      GO TO 590
 580  IF (IPHASE .EQ. 1) GO TO 570
C PREPARE OUTPUT.
      KODE = 0
 590  SUM = 0.D0
      DO 600 J = 1,N
         X(J) = 0.
 600  CONTINUE
      DO 610 I = 1,KLM
         RES(I) = 0.
 610  CONTINUE
      DO 640 I = 1,KLM
         II = Q(I,N2)
         SN = 1.
         IF(II .GT. 0) GO TO 620
         II = -II
         SN = -1.
 620     IF (II .GT. N) GO TO 630
         X(II) = SN * Q(I,N1)
         GO TO 640
 630     IIMN = II - N
         RES(IIMN) = SN * Q(I,N1)
         IF (II .GE. N1 .AND. II .LE. NK) SUM = SUM +
     *    DBLE(Q(I,N1))
 640  CONTINUE
      ERROR = SUM
      RETURN
      END
      SUBROUTINE CROSS(RFSPEC,SPECTRM,NPTS,SEARCH,ISHIFT)
C   THIS PERFORMS A CROSS-CORRELATION BETWEEN RFSPEC AND SPECTRM
C   STOLEN FROM JEFF BROWN WHO STOLE IT FROM FRITZ BENEDIT ...
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RFSPEC(10000),SPECTRM(10000)
      INTEGER RFNPTS,SEARCH
      DIMENSION LAG(401), RXY(401)
      RFNPTS = NPTS
C     IF(NPTS .EQ. RFNPTS) GO TO 20
C     MINPT = MAX0(101 - ISHIFT, 150)
C     MAXPT = 100 * ((MIN0(NPTS,RFNPTS) / 100) - 1)
C     GO TO 10
20    MINPT = SEARCH + 1 + ISHIFT
      IF (ISHIFT-SEARCH .LT. 0) MINPT = SEARCH-ISHIFT+1
      MINPT = MAX(MINPT,1)
      MAXPT = NPTS - MINPT
10    IF(MINPT .GT. MAXPT) GO TO 7734
C
      DO 100 I=1,2*SEARCH+1
         XYSUM = 0.
         XSUM = 0.
         YSUM = 0.
         X2SUM = 0.
         Y2SUM = 0.
         LAG(I) = ISHIFT - (SEARCH+1) + I
         DO 30 J = MINPT, MAXPT
            Y = SPECTRM(J + LAG(I))
            X = RFSPEC(J)
            XSUM = XSUM + X
            YSUM = YSUM + Y
            X2SUM = X2SUM + X * X
            Y2SUM = Y2SUM + Y * Y
            XYSUM = XYSUM + X * Y
30       CONTINUE
         RANGE = FLOAT(MAXPT - MINPT + 1)
         COVARXY = RANGE * XYSUM - XSUM * YSUM
         VARX = SQRT(RANGE * X2SUM - XSUM * XSUM)
         VARY = SQRT(RANGE * Y2SUM - YSUM * YSUM)
         IF (COVARXY.NE.0.0  .AND. VARX.NE.0.0  .AND. VARY.NE.0.0) THEN
            RXY(I) = COVARXY / VARX / VARY
         ELSE
            RXY(I) = 0.0
         ENDIF
100   CONTINUE
      WRITE(8,120) (LAG(I), RXY(I), I=1,2*SEARCH+1)
120   FORMAT(1X,I4,1X,F6.3)
      IMAX = 1
      DO 150 I=1,2*SEARCH+1
         IF(RXY(I) .LE. RXY(IMAX)) GO TO 150
         IMAX = I
150   CONTINUE
      ISHIFT = LAG(IMAX)
      WRITE(1,180)ISHIFT,RXY(IMAX)
180   FORMAT(' THE SPECTRUM IS SHIFTED BY ',I4,' DIODES, WITH A CROSS '
     %      ,'PRODUCT OF ',1PE15.7)
      RETURN
7734  WRITE(1,7740)
7740  FORMAT('0 DEATH IN CROSSC: MAX + MIN OF POINTS TO BE',
     1 ' CROSS-CORRELATED ARE BUNGED.')
      WRITE(1,7750) MAXPT, MINPT
7750  FORMAT('  MAX=',I5,'   MIN =',I5)
      STOP
      END

      SUBROUTINE INTERP(TEMP1,TEMP2,TEMP3,FILE1,FILE2,FILE3)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 NE1,NE2,NE3,MU1,MU2,MU3,KROSS1,KROSS2,KROSS3
      CHARACTER*(*) FILE1,FILE2,FILE3
      CALL OPNFLE(FILE1,FILE2,FILE3)
      DO 5 I=1,1000
      READ(7,*,END=99)TAU1,T1,PGAS1,NE1,KROSS1,MU1
      READ(8,*,END=99)TAU2,T2,PGAS2,NE2,KROSS2,MU2
      IF ( (TEMP1.NE.TEMP2).AND.(TEMP1.NE.TEMP3)
     %     .AND.(TEMP2.NE.TEMP3) )GOTO 200
         TAU3 = TAU1
         T3 = T1
         PGAS3 = PGAS1
         NE3 = NE1
         KROSS3 = KROSS1
         MU3 = MU1
         GOTO 300
 200  CONTINUE
      FRACTN=(TEMP1-TEMP3)/(TEMP1-TEMP2)
      TAU3=TAU1 - (TAU1-TAU2)*FRACTN
      T3=T1 - (T1-T2)*FRACTN
      PGAS3=PGAS1 - (PGAS1-PGAS2)*FRACTN
      NE3=NE1 - (NE1-NE2)*FRACTN
      KROSS3 = KROSS1 - (KROSS1-KROSS2)*FRACTN
      MU3 = MU1 - (MU1-MU2)*FRACTN
 300  CONTINUE
      WRITE(9,100)TAU3,T3,PGAS3,NE3,KROSS3,MU3
 100  FORMAT(6G13.6)
 5    CONTINUE
 99   REWIND 7
      REWIND 8
      REWIND 9
      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      RETURN
      END

      SUBROUTINE KINTERP(TEMP1,TEMP2,TEMP3,FILE1,FILE2,FILE3)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 NE1,NE2,NE3
      CHARACTER*(*) FILE1,FILE2,FILE3
      CALL OPNFLE(FILE1,FILE2,FILE3)
      DO 5 I=1,1000
      READ(7,*,END=99)TAU1,T1,PGAS1,NE1
      READ(8,*,END=99)TAU2,T2,PGAS2,NE2
      IF ( (TEMP1.NE.TEMP2).AND.(TEMP1.NE.TEMP3)
     %     .AND.(TEMP2.NE.TEMP3) )GOTO 200
         TAU3 = TAU1
         T3 = T1
         PGAS3 = PGAS1
         NE3 = NE1
         GOTO 300
 200  CONTINUE
      FRACTN=(TEMP1-TEMP3)/(TEMP1-TEMP2)
      TAU3=TAU1 - (TAU1-TAU2)*FRACTN
      T3=T1 - (T1-T2)*FRACTN
      PGAS3=PGAS1 - (PGAS1-PGAS2)*FRACTN
      NE3=NE1 - (NE1-NE2)*FRACTN
 300  CONTINUE
      WRITE(9,100)TAU3,T3,PGAS3,NE3
 100  FORMAT(G15.8,F9.1,1PE10.3,1PE10.3)
 5    CONTINUE
 99   REWIND 7
      REWIND 8
      REWIND 9
      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      RETURN
      END


      SUBROUTINE OPNFLE(FILE1,FILE2,FILE3)
      CHARACTER*(*) FILE1,FILE2,FILE3
      OPEN(FILE=FILE1,UNIT=7,BLANK='NULL')
      REWIND 7
      OPEN(FILE=FILE2,UNIT=8,BLANK='NULL')
      REWIND 8
      OPEN(UNIT=9,FILE=FILE3)
      RETURN
      END





      REAL*8 FUNCTION HELENA(X,F,A,N)                                  
C  *********************************************************************
C  RUDOLF LOESER, 8 SEP 64                                             
C  A(I) = INTEGRAL FROM X(I-1) TO X(I) OF F(X), FOR I FROM 1 TO N.     
C  UPON RETURN, THE VALUE OF THE FUNCTION WILL BE SUM OF ALL A'S.      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION X(5000),F(5000),A(5000)
      EQUIVALENCE (ISS,SS),(ISP,SP),(ISM,SM)                           
C  SET VALUE OF FIRST COMPONENT INTEGRAL.                              
      A(1) = 0.                                                        
C  COMPUTE SECOND COMPONENT INTEGRAL,USING SPECIALLY FITTED PARABOLIC  
C  SEGMENT                                                             
      TL = X(2) - X(1)                                                 
      TH = X(3) - X(1)                                                 
      A(2) = TL*TL*((4.*F(2)+2.*F(1))/TL-(F(3)-F(1))/TH)/6.            
C  INITIALIZE HELENA                                                   
      ANS = A(2)                                                       
C  LOOP OVER THE REMAINING INTEGRALS                                   
      DO 110 I=3,N                                                     
C  IF THIS IS THE LAST INTEGRAL, USE LOWER PARABOLIC SEGMENT.          
      IF(I .EQ. N) GO TO 106                                           
C  DECIDE WHETHER TO USE THE UPPER OR LOWER PARABOLIC SEGMENT FOR      
C  THE CURRENT INTEGRAL                                                
      SP = (F(I+1) - F(I))/(X(I+1) - X(I))                             
      SS = (F(I) - F(I-1))/(X(I) - X(I-1))                             
      SM = (F(I-1) - F(I-2))/(X(I-1) - X(I-2))                         
      IF(ISS - ISM) 101,103,103                                        
101   IF(ISS - ISP) 105,102,102                                        
102   IF(SS) 106,105,105                                               
103   IF(ISS - ISP) 104,104,106                                        
104   IF(SS) 105,106,106                                               
C  SET INDEX TO THE APPROPRIATE PARABOLA                               
105   J = I                                                            
      GO TO 107                                                        
106   J = I - 1                                                        
107   TL = X(J) - X(J-1)                                               
      TH = X(J+1) - X(J-1)                                             
      TU = X(J+1) - X(J)                                               
      FL = F(J) - F(J-1)                                               
C  SELECT FORMULA                                                      
      IF(J .NE. I) GO TO 109                                           
      FX = F(I+1) - F(I-1)                                             
      A(I)=TL*(F(I-1)+(3.*FL*TH*TH-2.*FL*TH*TL-FX*TL*TL)/(6.*TH*TU))   
      GO TO 110                                                        
C  COMPUTE WITH LOWER PARABOLA                                         
109   FX = F(I) - F(I-1)                                               
      A(I)=TU*(F(I-1)+(FL*TU*TU+2.*FX*TL*TU+3.*FX*TL*TL)/(6.*TH*TL))   
C  UPDATE HELENA                                                       
110   ANS = ANS + A(I)                                                 
      HELENA = ANS                                                     
      RETURN                                                           
      END                                                              





	subroutine polyfit(x,y,sig,ma,a,chisq,ndata)
c
c  A program to perform a non-linear least squares fit to data to a polynomial!
c  ma is order of polynomial (ie number of coefficients to be found)
c  only laziness keeps me from removing undecessary code
c
        implicit real*8 (a-h,o-z)
	dimension lista(50)
	real*8 x(1000),y(1000),sig(1000),a(50),covar(50,50),
     %            alpha(50,50)
	data tol/0.05/
        data nca/1000/
	call minmax(x,ndata,amin,amax)
c       call fixsig(sig,ndata)             *** to be used when weights are equal
	call fixlst(lista,ma,mfit)
	call guessa(a,ma)
	alamda = -1.0
	ochisq = 0.0
	do 100 i = 1,100
	   call mrqmin(x,y,sig,ndata,a,ma,lista,mfit,covar,alpha,nca,
     %                 chisq,alamda)
	   delchi = abs(ochisq-chisq)
	   if ( (delchi.le.tol).and.(i.ne.1) ) then
	      alamda = 0.0
              call mrqmin(x,y,sig,ndata,a,ma,lista,mfit,covar,alpha,
     %                    nca,chisq,alamda)
	      return
	   endif
	   ochisq = chisq
 100	continue
	return
	end


	subroutine minmax(x,ndata,amin,amax)
        implicit real*8 (a-h,o-z)
	real*8 x(1000)
	amax = x(1)
	amin = x(1)
	do 100 i = 2,ndata
	   if ( x(i).gt.amax ) amax = x(i)
	   if ( x(i).lt.amin ) amin = x(i)
 100	continue
	return
	end


	subroutine fixsig(sig,ndata)
	implicit real*8 (a-h,o-z)
	real*8 sig(1000)
	do 100 i = 1,ndata
 100	   sig(i) = 100.0
	return
	end

	subroutine fixlst(lista,ma,mfit)
        implicit real*8 (a-h,o-z)
c
c  Routine to set up the sorting list lista.  Presently the simple case
c  where all parameters are to be iterated on is set.
c
	dimension lista(50)
	do 100 i = 1,ma
	   lista(i) = i
 100	continue
	mfit = ma
	return
	end

	subroutine guessa(a,ma)
	implicit real*8 (a-h,o-z)
	real*8 a(50)
	a(1) = 1.0
        do 100 i = 2,ma
           a(i) = 0.0 
 100	continue
	return
	end

	subroutine wrdat(a,ma,covar,chisq,alpha,amin,amax)
	implicit real*8 (a-h,o-z)
c
c  Output subroutine
c
        real*8 a(50),covar(50,50),alpha(50,50),chisq,amin,amax,range
	write(6,'(17h A values follow )')
	write(6,*)(a(i),i=1,ma)
	write(6,'(10h chi**2 = ,g12.5)')chisq
	write(6,'(27h covarience matrix follows )')
	write(6,'(4g12.5)')((covar(j,i),i=1,ma),j=1,ma)
	amin = aint( 10*amin )/10.0
	amax = aint( 10*amax )/10.0 + 0.1
	range = amax - amin
        do 100 i = 1,101
       	   xi = amin + dble(i-1)*range/100.0
       	   yi = a(1) + a(2)*xi + a(3)*xi**2 + a(4)*xi**3 
       	   write(6,*)yi,xi
 100    continue
	return
	end

	subroutine mrqmin(x,y,sig,ndata,a,ma,lista,mfit,
     %                    covar,alpha,nca,chisq,alamda)
	implicit real*8 (a-h,o-z)
c
c  Levenberg-Marquardt method, attempting to reduce chi**2 of a fit
c  between ndata points x(i),Y(i) with individual std deviations sig(i),
c  and nonlinear function dependent on ma coefficients a.  Array lista
c  numbers the parameters a such that the first mfit elements
c  correspond to values actually being adjusted.  The remaining ma-mfit
c  parameters are kept at the input value.  The routine returns current
c  best fit values for the ma fit parameters a, and chi**2 chisq.
c  The arrays covar(nca,nca0, alpha(nca,nca) with physical dimensions
c  nca (>= mfit) are used as working space during most iterations.
c
c  A user subroutine is required, FUNCS(x,a,yfit,dyda,ma) that evaluates
c  the fitting function yfit, and it's derivatives dyda with respect to 
c  the fitting parameters a at x.  On the first call provide an initial
c  guess for the parameters a, and set alamda < 0 (which then sets
c  alamda=0.001).  If a step succeeds chisq becomes smaller and alamda
c  decreases by a factor of 10.  If the step fails alamda is increased by
c  a factor of 10.
c
c  To use the subroutine one must make repeated calls until convergence is
c  achieved, then a final call with alamda=0 so that covar(i,j) is returned
c  with the covarience matrix and alpha(i,j) with the curvature matrix.
c
	parameter (mmax=50)
	dimension x(1000),y(1000),sig(1000),a(50),lista(50),
     %  covar(1000,1000),alpha(1000,1000),atry(50),beta(50),da(50)
	if (alamda.lt.0.0) then
	   kk = mfit + 1
	   do 12 j=1,ma
	      ihit=0
	      do 11 k=1,mfit
	         if (lista(k).eq.j)ihit=ihit+1
 11              continue
	      if (ihit.eq.0) then
	         lista(kk) = j
                 kk = kk + 1
	      elseif (ihit.gt.1) then
	         pause 'Improper permutation in LISTA'
	      endif
 12	      continue
	   if (kk.ne.(ma+1)) pause 'Improper permutation in LISTA'
	   alamda = 100.0
           call mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,nca,
     %                 chisq)
	   ochisq = chisq
	   do 13 j = 1,ma
	      atry(j) = a(j)
 13           continue
	endif
	do 15 j = 1,mfit
	   do 14 k = 1,mfit
	      covar(j,k) = alpha(j,k)
 14           continue
           covar(j,j) = alpha(j,j)*(1.0+alamda)
           da(j) = beta(j)
 15        continue
c
c  Matrix solution
c
	call gaussj(covar,mfit,nca,da,1,1)
	if (alamda.eq.0.0) then
c  diagnostic for covar
	   call covsrt(covar,nca,ma,lista,mfit)
	   return
        endif
c
c  See if the trial succeeded
c
	do 16 j = 1,mfit
           atry(lista(j)) = a(lista(j)) + da(j)
 16	   continue
           call mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,
     %                 nca,chisq)
	if (chisq.lt.ochisq) then
	   alamda = 0.1*alamda
	   ochisq = chisq
	   do 18 j = 1,mfit
	      do 17 k = 1,mfit
	         alpha(j,k) = covar(j,k)
 17	         continue
 	      beta(j) = da(j)
	      a(lista(j)) = atry(lista(j))
 18	      continue
	else
	   alamda = 10.0*alamda
	   chisq = ochisq
	endif
	return
	end

	subroutine gaussj(a,n,np,b,m,mp)
	implicit real*8 (a-h,o-z)
c
c  Linear equation solving by Gauss-Jordan elimination. a is the input
c  matrix of n by n elements, stored in a physical array of np by np
c  elements.  b is a n by m input matrix containing m right hand side
c  vectors, stored in an array of physical dimension np by mp.  On output
c  a is replaced by its matrix inverse, and b is replaced by the 
c  corresponding set of solution vectors.
c  
	parameter (nmax=50)
        dimension a(1000,1000),b(1000,1),ipiv(50),indxr(50),indxc(50)
c
c  ipiv, indxr and indxc are used for pivoting bookkeeping.  nmax
c  should be as large as the largest anticipated value of n.
c
	do 11 j = 1,n
	   ipiv(j) = 0
 11	continue
	do 22 i = 1,n
	   big = 0.0
	   do 13 j = 1,n
	      if (ipiv(j).ne.1) then
	         do 12 k = 1,n
	            if (ipiv(k).eq.0) then
	               if (abs(a(j,k)).ge.big) then
	                  big = abs(a(j,k))
	                  irow = j
                          icol = k
                       endif
                    elseif (ipiv(k).gt.1) then
                       pause ' Singular matrix'
                    endif
 12	         continue
              endif
 13	   continue
	   ipiv(icol) = ipiv(icol) + 1
	   if (irow.ne.icol) then
	      do 14 l = 1,n
	         dum = a(irow,l)
	         a(irow,l) = a(icol,l)
	         a(icol,l) = dum
 14	      continue
	      do 15 l = 1,m
	         dum = b(irow,l)
	         b(irow,l) = b(icol,l)
	         b(icol,l) = dum
 15	      continue
	   endif
	   indxr(i) = irow
	   indxc(i) = icol
	   if (a(icol,icol).eq.0) pause 'Singular matrix'
	   pivinv = 1.0/a(icol,icol)
	   a(icol,icol) = 1.0
	   do 16 l = 1,n
	      a(icol,l) = a(icol,l)*pivinv
 16	   continue
	   do 17 l = 1,m
	      b(icol,l) = b(icol,l)*pivinv
 17	   continue
	   do 21 ll = 1,n
	      if (ll.ne.icol) then
	         dum = a(ll,icol)
	         a(ll,icol) = 0.0
	         do 18 l = 1,n
	            a(ll,l) = a(ll,l) - a(icol,l)*dum
 18	         continue
	         do 19 l = 1,m
	            b(ll,l) = b(ll,l) - b(icol,l)*dum
 19	         continue
	      endif
 21	   continue
 22	continue
	do 24 l = n,1,-1
	   if (indxr(l).ne.indxc(l)) then 
	      do 23 k = 1,n
	         dum = a(k,indxr(l))
	         a(k,indxr(l)) = a(k,indxc(l))
	         a(k,indxc(l)) = dum
 23	      continue
	   endif
 24	continue
	return
	end

	subroutine mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,
     %                    nalp,chisq)
	implicit real*8(a-h,o-z)
c
c  This routine is used by mrqmin to evaluate the linearized fitting matrix
c  alpha, and vector beta.
c
	parameter (mmax=50)
        dimension x(1000),y(1000),sig(1000),alpha(1000,1000),
     %            beta(50),dyda(50),lista(50),a(50)
c
c  Initialize
c
	do 12 j = 1,mfit
	   do 11 k = 1,j
	      alpha(j,k) = 0.0
 11	      continue
 	   beta(j) = 0.0
 12	   continue
	chisq = 0.0
	do 15 i = 1,ndata
	   call fcos(x(i),a,ymod,dyda,ma)
	   sig2i = 1.0/(sig(i)*sig(i))
	   dy = y(i) - ymod
	   do 14 j = 1,mfit
	      wt = dyda(lista(j))*sig2i
	      do 13 k = 1,j
	         alpha(j,k) = alpha(j,k) + wt*dyda(lista(k))
 13	         continue
	      beta(j) = beta(j) + dy*wt
 14	      continue
 	   chisq = chisq + dy*dy*sig2i
 15	   continue
	do 17 j = 2,mfit
	   do 16 k = 1,j-1
	      alpha(k,j) = alpha(j,k)
 16	      continue
 17	   continue
	return
	end

	subroutine covsrt(covar,ncvm,ma,lista,mfit)
	implicit real*8 (a-h,o-z)
c
c   Given the covarience matrix covar of a fit mfit of ma parameters,
c   and their ordering lista(i), repack the covarience matrix to the true order
c   of the parameters.  Elements associated with fixed parameters are zero.  
c   ncvm is the physical dimension of covar.
	dimension covar(1000,1000),lista(50)
	do 12 j = 1,ma-1
	   do 11 i = j+1,ma
	      covar(i,j) = 0.0
 11	      continue
 12	   continue
	do 14 i = 1,mfit-1
	   do 13 j = i+1,mfit
	      if (lista(j).gt.lista(i)) then
	         covar(lista(j),lista(i)) = covar(i,j)
	      else
	         covar(lista(i),lista(j)) = covar(i,j)
	      endif
 13	      continue
 14	   continue
	swap = covar(1,1)
	do 15 j = 1,ma
	   covar(1,j) = covar(j,j)
	   covar(j,j) = 0.0
 15	   continue
	covar(lista(1),lista(1)) = swap
 	do 16 j = 2,mfit
	   covar(lista(j),lista(j)) = covar(1,j)
 16	   continue
	do 18 j = 2,ma
	   do 17 i = 1,j-1
	      covar(i,j) = covar(j,i)
 17	      continue
 18	   continue
	return
	end 
	subroutine fcos(xi,a,yi,dyda,npar)
	implicit real*8 (a-h,o-z)
c
c  A subroutine to return the y value and derivative of a polynomial order npar
c   y = a1 + a2.x + a3.x**2 + a4.x**3 + a5.x**4 + .... + anpar.x**(npar-1)
c
	real*8 a(50),dyda(50),xi,yi
        yi = 0.0
        do 100 i = 1,npar
           yi = yi + a(i)*xi**(i-1)
           dyda(i) = xi**(i-1)
 100	continue
	return
	end

	subroutine polylin(x,y,sig,ma,a,covar,chisq,ndata)
c
c  Given a set of ndata points x(i), y(i) with individual deviations sig(i)
c  use Chi**2 minimization to determine mfit of ma coefficients a, of a 
c  function that depends linearly only on a, y=Sigma a.func(x).  The array lista
c  renumbers the parameters so that the first mfit elements correspond to the
c  parameters actually being determined; the remaining ma - mfit elements are
c  held fixed at their input value.  The program returns values for the ma fit
c  parameters a, Chi**2, CHI2, and the covariance matrix COVAR(I,J).  NCVM is 
c  the physical dimension of COVAR(NVCM,NVCM) in the calling routine.  The
c  user supplies a subroutine FUNCS(X,AFUNC,MA) that returns the MA basis
c  functions evaluated at x=A in the array AFUNC.
c
	implicit real*8(a-h,o-z)
	dimension x(1000),y(1000),sig(1000),a(50),lista(50),
     %            covar(50,50),beta(50),afunc(50)
        ncvm = 50
        mmax = 50
c
c  make all list variable
c
	do 100 i = 1,ma
	   lista(i) = i
 100	continue
c
c  don't fix any coefficients, ie mfit = ma
c
        mfit = ma
	kk = mfit+1
	do 12 j=1,ma
	   ihit = 0
	   do 11 k = 1,mfit
	      if (lista(k).eq.j)ihit=ihit+1
 11	   continue
	   if (ihit.eq.0) then
	      lista(kk)=j
	      kk = kk+1
	   elseif (ihit.gt.1) then
	      pause 'Improper set in lista'
           endif
 12	continue
	if (kk.ne.(ma+1)) pause 'Improper set in lista'
	do 14 j = 1,mfit
	   do 13 k=1,mfit
	      covar(j,k)=0.0
 13	   continue
	   beta(j) = 0.0
 14	continue
	do 18 i=1,ndata
	   call funcs(x(i),afunc,ma)
	   ym = y(i)
	   if (mfit.lt.ma) then
	      do 15 j=mfit+1,ma
	         ym = ym-a(lista(j))*afunc(lista(j))
 15	      continue
	   endif
	   sig2i = 1.0/sig(i)**2
	   do 17 j = 1,mfit
	      wt = afunc(lista(j))*sig2i
	      do 16 k=1,j
	         covar(j,k) = covar(j,k)+wt*afunc(lista(k))
 16	      continue
	      beta(j) = beta(j)+ym*wt
 17	   continue
 18	continue
	if (mfit.gt.1) then
	   do 21 j=2,mfit
	      do 19 k=1,j-1
	         covar(k,j) = covar(j,k)
 19	      continue
 21	   continue
	endif
	call gaussj2(covar,mfit,ncvm,beta,1,1)
	do 22 j=1,mfit
	   a(lista(j)) = beta(j)
 22	continue
	chisq = 0.0
	do 24 i = 1,ndata
	   call funcs(x(i),afunc,ma)
	   sum = 0.0
	   do 23 j=1,ma
	      sum = sum + a(j)*afunc(j)
 23	   continue
	   chisq = chisq + ((y(i)-sum)/sig(i))**2
 24	continue
	call covsrt2(covar,ncvm,ma,lista,mfit)
	return
	end
	subroutine gaussj2(a,n,np,b,m,mp)
c
c  Linear equation solving by Gauss-Jordan elimination. a is the input
c  matrix of n by n elements, stored in a physical array of np by np
c  elements.  b is a n by m input matrix containing m right hand side
c  vectors, stored in an array of physical dimension np by mp.  On output
c  a is replaced by its matrix inverse, and b is replaced by the 
c  corresponding set of solution vectors.
c  
	implicit real*8 (a-h,o-z)
	parameter (nmax=50)
        dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
c
c  ipiv, indxr and indxc are used for pivoting bookkeeping.  nmax
c  should be as large as the largest anticipated value of n.
c
	do 11 j = 1,n
	   ipiv(j) = 0
 11	continue
	do 22 i = 1,n
	   big = 0.0
	   do 13 j = 1,n
	      if (ipiv(j).ne.1) then
	         do 12 k = 1,n
	            if (ipiv(k).eq.0) then
	               if (abs(a(j,k)).ge.big) then
	                  big = abs(a(j,k))
	                  irow = j
                          icol = k
                       endif
                    elseif (ipiv(k).gt.1) then
                       pause ' Singular matrix'
                    endif
 12	         continue
              endif
 13	   continue
	   ipiv(icol) = ipiv(icol) + 1
	   if (irow.ne.icol) then
	      do 14 l = 1,n
	         dum = a(irow,l)
	         a(irow,l) = a(icol,l)
	         a(icol,l) = dum
 14	      continue
	      do 15 l = 1,m
	         dum = b(irow,l)
	         b(irow,l) = b(icol,l)
	         b(icol,l) = dum
 15	      continue
	   endif
	   indxr(i) = irow
	   indxc(i) = icol
	   if (a(icol,icol).eq.0) pause 'Singular matrix'
	   pivinv = 1.0/a(icol,icol)
	   a(icol,icol) = 1.0
	   do 16 l = 1,n
	      a(icol,l) = a(icol,l)*pivinv
 16	   continue
	   do 17 l = 1,m
	      b(icol,l) = b(icol,l)*pivinv
 17	   continue
	   do 21 ll = 1,n
	      if (ll.ne.icol) then
	         dum = a(ll,icol)
	         a(ll,icol) = 0.0
	         do 18 l = 1,n
	            a(ll,l) = a(ll,l) - a(icol,l)*dum
 18	         continue
	         do 19 l = 1,m
	            b(ll,l) = b(ll,l) - b(icol,l)*dum
 19	         continue
	      endif
 21	   continue
 22	continue
	do 24 l = n,1,-1
	   if (indxr(l).ne.indxc(l)) then 
	      do 23 k = 1,n
	         dum = a(k,indxr(l))
	         a(k,indxr(l)) = a(k,indxc(l))
	         a(k,indxc(l)) = dum
 23	      continue
	   endif
 24	continue
	return
	end
	subroutine covsrt2(covar,ncvm,ma,lista,mfit)
c
c   Given the covarience matrix covar of a fit mfit of ma parameters,
c   and their ordering lista(i), repack the covarience matrix to the true order
c   of the parameters.  Elements associated with fixed parameters are zero.  
c   ncvm is the physical dimension of covar.
	implicit real*8 (a-h,o-z)
	dimension covar(ncvm,ncvm),lista(mfit)
	do 12 j = 1,ma-1
	   do 11 i = j+1,ma
	      covar(i,j) = 0.0
 11	      continue
 12	   continue
	do 14 i = 1,mfit-1
	   do 13 j = i+1,mfit
	      if (lista(j).gt.lista(i)) then
	         covar(lista(j),lista(i)) = covar(i,j)
	      else
	         covar(lista(i),lista(j)) = covar(i,j)
	      endif
 13	      continue
 14	   continue
	swap = covar(1,1)
	do 15 j = 1,ma
	   covar(1,j) = covar(j,j)
	   covar(j,j) = 0.0
 15	   continue
	covar(lista(1),lista(1)) = swap
 	do 16 j = 2,mfit
	   covar(lista(j),lista(j)) = covar(1,j)
 16	   continue
	do 18 j = 2,ma
	   do 17 i = 1,j-1
	      covar(i,j) = covar(j,i)
 17	      continue
 18	   continue
	return
	end 
	subroutine funcs(x,p,np)
c
c  fitting routine for polynomial of degree np-1 with np coefficients
c
	implicit real*8 (a-h,o-z)
	dimension p(np)
	p(1) = 1.0
	do 11 j=2,np
	   p(j) = p(j-1)*x
 11	continue
	return
	end
