	program makekurucz3
c
c a routine to generate an interpolated Kurucz model
c for Sneden's MOOG
c
c A.McWilliam

        real*8 teff,logg,feh,micro
        character*10 mtype

	write(6,'(a)')'Enter Teff, logg, [Fe/H] and micro.'
        read(5,*)teff,logg,feh,micro

c
c interpolate model of specified type
c
        write(6,*)'Enter Kurucz model type: KUROLD, ODFNEW, AODFNEW, NOVER'
        read(5,'(a)')mtype
        if (mtype.ne.'KUROLD' .and.  mtype.ne.'ODFNEW' .and. 
     %     mtype.ne.'AODFNEW' .and. mtype.ne.'NOVER' )  then
           write(6,*)'Invalid type: must be KUROLD, ODFNEW, AODFNEW, or NOVER'
           stop
        endif
c

        call make_kurucz(teff,logg,feh,micro,mtype)
	stop
	end





c
c This routine produces interpolated models in the format expected by
c Chris Sneden's 1998 MOOG.
c
c When changing systems you need to change assignments to variable DIR 
c in subroutine PRINT_MODEL
c
c A.McWilliam
c

      SUBROUTINE MAKE_KURUCZ(TEFF,LOGG,FE,MICRO,MTYPEIN)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 TEFF,LOGG,FE,MICRO
      LOGICAL ERROR
      CHARACTER*80 TITLE,MTYPEIN*10
      COMMON/MODTYPE/ MTYPE
      CHARACTER*10 MTYPE

      TITLE = '          '
      ERROR = .FALSE.
      MTYPE = MTYPEIN

c write parameters
      WRITE(TITLE,100)TEFF,LOGG,FE,MICRO
 100  FORMAT(F7.0,1h/,F7.2,1h/,F7.2,5x,7h mic = ,F6.4)


      IF (MTYPE.EQ.'KUROLD') THEN
         CALL GET_MODEL(TEFF,LOGG,FE,ERROR)
      ELSEIF (MTYPE.EQ.'ODFNEW') THEN
         CALL GET_MODEL_ODFNEW(TEFF,LOGG,FE,ERROR)
      ELSEIF (MTYPE.EQ.'AODFNEW') THEN
         CALL GET_MODEL_AODFNEW(TEFF,LOGG,FE,ERROR)
      ELSEIF (MTYPE.EQ.'NOVER') THEN
         CALL GET_MODEL_NOVER(TEFF,LOGG,FE,ERROR)
      ENDIF
c
      IF (ERROR) THEN
          WRITE(6,'(a)')'ERROR IN MAKE_KURUCZ SUBROUTINE '
          STOP 
      ENDIF
      CALL SET_UP_LINES_INPUT_FILE(MICRO,FE,TITLE)
      CLOSE(UNIT=7)
      RETURN
      END

      



      SUBROUTINE GET_MODEL(TEFF,LOGG,FE,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 LOGG,TEFF,FE,METAL(10),M1,M2
      LOGICAL ERROR
      CHARACTER*8 M1FIL,M2FIL,MODEL
      DATA METAL/+1.0D0,+0.5D0,0.0D0,-0.5D0,-1.0D0,-1.5D0,-2.0D0,-2.5D0,
     %           -3.0D0,-3.5D0/
      DATA M1FIL,M2FIL,MODEL/'M1     ','M2     ','model  '/
      
      IF (FE .GT. +1.0D0) THEN   
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL RICH THAN +1.0: +1.0 MODEL USED'
         FE = 1.0D0   
      ENDIF
      IF (FE .LT. -3.5D0)   THEN
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL POOR THAN -3.5: -3.5 MODEL USED'
         FE = -3.5D0  
      ENDIF
      IF (FE .EQ. +1.0D0 .OR. FE .EQ. +0.5D0 .OR. FE .EQ.  0.0D0 .OR. 
     %    FE .EQ. -0.5D0 .OR. FE .EQ. -1.0D0 .OR. FE .EQ. -1.5D0 .OR.
     %    FE .EQ. -2.0D0 .OR. FE .EQ. -2.5D0 .OR. FE. EQ. -3.0D0 .OR.
     %    FE .EQ. -3.5D0) THEN

         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,FE,MODEL,ERROR)   
      ELSE
         DO I = 2, 10   
            IF (FE .GT. METAL(I)) THEN
               M1 = METAL(I-1)
               M2 = METAL(I)
               GOTO 800
            ENDIF
         ENDDO
 800     CONTINUE
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M1,M1FIL,ERROR)
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M2,M2FIL,ERROR)
         IF (ERROR)  RETURN
         CALL KINTERP(M1,M2,FE,M1FIL,M2FIL,MODEL)
      ENDIF
      RETURN
      END




      SUBROUTINE GET_MODEL_ODFNEW(TEFF,LOGG,FE,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 LOGG,TEFF,FE,METAL(7),M1,M2
      LOGICAL ERROR
      CHARACTER*8 M1FIL,M2FIL,MODEL
      DATA METAL/+0.5D0,0.0D0,-0.5D0,-1.0D0,-1.5D0,-2.0D0,-2.5D0/

      DATA M1FIL,M2FIL,MODEL/'M1     ','M2     ','model  '/
      
      IF (FE .GT. +0.5D0) THEN   
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL RICH THAN +0.5: +0.5 MODEL USED'
         FE = 1.0D0   
      ENDIF
      IF (FE .LT. -2.5D0)   THEN
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL POOR THAN -2.5: -2.5 MODEL USED'
         FE = -2.5D0  
      ENDIF
      IF (FE .EQ. +0.5D0 .OR. FE .EQ.  0.0D0 .OR. 
     %    FE .EQ. -0.5D0 .OR. FE .EQ. -1.0D0 .OR. FE .EQ. -1.5D0 .OR.
     %    FE .EQ. -2.0D0 .OR. FE .EQ. -2.5D0) THEN

         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,FE,MODEL,ERROR)   
      ELSE
         DO I = 2, 7   
            IF (FE .GT. METAL(I)) THEN
               M1 = METAL(I-1)
               M2 = METAL(I)
               GOTO 800
            ENDIF
         ENDDO
 800     CONTINUE
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M1,M1FIL,ERROR)
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M2,M2FIL,ERROR)
         IF (ERROR)  RETURN
         CALL KINTERP(M1,M2,FE,M1FIL,M2FIL,MODEL)
      ENDIF
      RETURN
      END





      SUBROUTINE GET_MODEL_AODFNEW(TEFF,LOGG,FE,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 LOGG,TEFF,FE,METAL(8),M1,M2
      LOGICAL ERROR
      CHARACTER*8 M1FIL,M2FIL,MODEL
      DATA METAL/+0.5D0,0.0D0,-0.5D0,-1.0D0,-1.5D0,-2.0D0,-2.5D0,-4.0D0/

      DATA M1FIL,M2FIL,MODEL/'M1     ','M2     ','model  '/
      
      IF (FE .GT. +0.5D0) THEN   
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL RICH THAN +0.5: +0.5 MODEL USED'
         FE = 1.0D0   
      ENDIF
      IF (FE .LT. -4.0D0)   THEN
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL POOR THAN -4.0: -4.0 MODEL USED'
         FE = -4.0D0  
      ENDIF
      IF (FE .EQ. +0.5D0 .OR. FE .EQ.  0.0D0 .OR. 
     %    FE .EQ. -0.5D0 .OR. FE .EQ. -1.0D0 .OR. FE .EQ. -1.5D0 .OR.
     %    FE .EQ. -2.0D0 .OR. FE .EQ. -2.5D0 .OR. FE .EQ. -4.0D0) THEN

         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,FE,MODEL,ERROR)   
      ELSE
         DO I = 2, 8
            IF (FE .GT. METAL(I)) THEN
               M1 = METAL(I-1)
               M2 = METAL(I)
               GOTO 800
            ENDIF
         ENDDO
 800     CONTINUE
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M1,M1FIL,ERROR)
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M2,M2FIL,ERROR)
         IF (ERROR)  RETURN
         CALL KINTERP(M1,M2,FE,M1FIL,M2FIL,MODEL)
      ENDIF
      RETURN
      END



      SUBROUTINE GET_MODEL_NOVER(TEFF,LOGG,FE,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 LOGG,TEFF,FE,METAL(7),M1,M2
      LOGICAL ERROR
      CHARACTER*8 M1FIL,M2FIL,MODEL
      DATA METAL/+0.5D0,0.0D0,-0.5D0,-1.0D0,-1.5D0,-2.0D0,-2.5D0/

      DATA M1FIL,M2FIL,MODEL/'M1     ','M2     ','model  '/
      
      IF (FE .GT. +0.5D0) THEN   
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL RICH THAN +0.5: +0.5 MODEL USED'
         FE = 1.0D0   
      ENDIF
      IF (FE .LT. -4.0D0)   THEN
         WRITE(6,'(a)')
     %      'STAR IS MORE METAL POOR THAN -2.5: -2.5 MODEL USED'
         FE = -4.0D0  
      ENDIF
      IF (FE .EQ. +0.5D0 .OR. FE .EQ.  0.0D0 .OR. 
     %    FE .EQ. -0.5D0 .OR. FE .EQ. -1.0D0 .OR. FE .EQ. -1.5D0 .OR.
     %    FE .EQ. -2.0D0 .OR. FE .EQ. -2.5D0 ) THEN

         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,FE,MODEL,ERROR)   
      ELSE
         DO I = 2, 7
            IF (FE .GT. METAL(I)) THEN
               M1 = METAL(I-1)
               M2 = METAL(I)
               GOTO 800
            ENDIF
         ENDDO
 800     CONTINUE
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M1,M1FIL,ERROR)
         CALL COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,M2,M2FIL,ERROR)
         IF (ERROR)  RETURN
         CALL KINTERP(M1,M2,FE,M1FIL,M2FIL,MODEL)
      ENDIF
      RETURN
      END





      SUBROUTINE COMPUTE_MODEL_ATMOSPHERE(TEFF,LOGG,FE,OUTFILE,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 OUTFILE
      REAL*8 TEFF,TEFFS(2),LOGG,GRAVS(2),FE
      LOGICAL OKG,OKT,ERROR

      CALL CHECK_PARAMETERS(OKG,OKT,TEFF,LOGG,TEFFS,GRAVS,ERROR)
      IF (ERROR) RETURN
      IF (OKG .AND. OKT) THEN
         CALL PRINT_MODEL(TEFF,LOGG,OUTFILE,FE,ERROR)
         IF (ERROR) RETURN
      ELSEIF (OKG) THEN
         CALL PRINT_MODEL(TEFFS(1),LOGG,'MOD1    ',FE,ERROR)
         CALL PRINT_MODEL(TEFFS(2),LOGG,'MOD2    ',FE,ERROR)
         IF (ERROR) RETURN
         CALL KINTERP(TEFFS(1),TEFFS(2),TEFF,'MOD1    ','MOD2    ',OUTFILE)
      ELSEIF (OKT)   THEN
         CALL PRINT_MODEL(TEFF,GRAVS(1),'MOD1    ',FE,ERROR)
         CALL PRINT_MODEL(TEFF,GRAVS(2),'MOD2    ',FE,ERROR)
         IF (ERROR) RETURN  
         CALL KINTERP(GRAVS(1),GRAVS(2),LOGG,'MOD1    ','MOD2    ',OUTFILE)
      ELSE   
         CALL PRINT_MODEL(TEFFS(1),GRAVS(1),'MOD1    ',FE,ERROR)
         CALL PRINT_MODEL(TEFFS(2),GRAVS(1),'MOD2    ',FE,ERROR)
         CALL PRINT_MODEL(TEFFS(1),GRAVS(2),'MOD3    ',FE,ERROR)
         CALL PRINT_MODEL(TEFFS(2),GRAVS(2),'MOD4    ',FE,ERROR)
         IF (ERROR) RETURN
         CALL KINTERP(GRAVS(1),GRAVS(2),LOGG,'MOD1    ','MOD3    ','MOD5    ')
         CALL KINTERP(GRAVS(1),GRAVS(2),LOGG,'MOD2    ','MOD4    ','MOD6    ')
         CALL KINTERP(TEFFS(1),TEFFS(2),TEFF,'MOD5    ','MOD6    ',OUTFILE)   
      ENDIF
      RETURN
      END


      SUBROUTINE CHECK_PARAMETERS(OKG,OKT,TEFF,GRAV,TEFFS,GRAVS,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKG,OKT,ERROR
      INTEGER I
      REAL*8 TEFF,TEFFS(2),GRAV,GRAVS(2)

      IF (TEFF.LT.3500.0D0  .OR. TEFF.GT.10000.0D0) THEN
         WRITE(6,'(a)')'ERROR: TEMPERATURE IS OUT OF BOUNDS'
         STOP   
      ENDIF

c  check to see if Teff is a multiple of 250K  

      I = INT(TEFF/250.0D0)
      IF (DABS(TEFF - DBLE(I)*250.0D0) .LT. 0.5D0) THEN
         TEFFS(1) = DBLE(I)*250.0D0
         TEFFS(2) = DBLE(I)*250.0D0
         TEFF = DBLE(NINT(TEFF))
         OKT = .TRUE.   
      ELSE   
         TEFFS(1) = DBLE(I)*250.0D0
         TEFFS(2) = DBLE(I+1)*250.0D0
         OKT = .FALSE.   
      ENDIF
      
c  Now logg  

      IF (GRAV.LT.0.0D0  .OR.  GRAV.GT.5.0D0) THEN
         WRITE(6,'(a)')'ERROR: GRAVITY IS OUT OF BOUNDS'
         STOP   
      ENDIF

c  check to see if logg is a multiple of 0.5  

      I = INT(GRAV/0.5D0)
      IF (DABS(GRAV - DBLE(I)*0.5D0) .LT. 0.005D0) THEN
         GRAVS(1) = DBLE(I)*0.5D0
         GRAVS(2) = DBLE(I)*0.5D0
         GRAV = DBLE(NINT(GRAV*2.0))*0.5
         OKG = .TRUE.   
      ELSE   
         GRAVS(1) = DBLE(I)*0.5D0
         GRAVS(2) = DBLE(I+1)*0.5D0
         OKG = .FALSE.
      ENDIF
      RETURN
      END
      


      SUBROUTINE PRINT_MODEL(TEFF,LOGG,MODEL,FE,ERROR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MODTYPE/ MTYPE
      CHARACTER*10 MTYPE
      LOGICAL ERROR
      CHARACTER*(4) MODEL*(8),FILES(11)
      CHARACTER*(80) DIR,LINE
      CHARACTER*(80) FNAME
      REAL*8 TEMP,GRAV,LOGG,TEFF,FE
      DATA DIR/'                                              '/
      DATA FILES /'ap10','ap05','ap00','am05','am10','am15','am20',
     %            'am25', 'am30','am35', 'am40'/

      ERROR = .TRUE.
c
c Set the directory according to the type of Kurucz model being used
c
      IF (MTYPE.EQ.'KUROLD') THEN
      DIR  = 'kurold/'
      ELSEIF (MTYPE.EQ.'ODFNEW') THEN
      DIR  = 'odfnew/'
      ELSEIF (MTYPE.EQ.'AODFNEW') THEN
      DIR  = 'aodfnew/'
      ELSEIF (MTYPE.EQ.'NOVER') THEN
      DIR  = 'nover/'
      ENDIF

c
      IPOSN = INDEX(DIR,' ')
      IF (FE .EQ. 1.0D0) THEN  
         IFILE = 1
      ELSEIF (FE .EQ. 0.5D0) THEN
          IFILE = 2
      ELSEIF (FE .EQ.  0.0D0)  THEN
         IFILE = 3  
      ELSEIF (FE .EQ. -0.5D0)  THEN  
         IFILE = 4  
      ELSEIF (FE .EQ. -1.0D0)  THEN  
         IFILE = 5  
      ELSEIF (FE .EQ. -1.5D0)  THEN  
         IFILE = 6  
      ELSEIF (FE .EQ. -2.0D0)  THEN  
         IFILE = 7  
      ELSEIF (FE .EQ. -2.5D0)  THEN  
         IFILE = 8  
      ELSEIF (FE .EQ. -3.0D0)  THEN  
         IFILE = 9  
      ELSEIF (FE .EQ. -3.5D0)  THEN  
         IFILE = 10  
      ELSEIF (FE .EQ. -4.0D0)  THEN  
         IFILE = 11
      ELSE  
         ERROR = .TRUE.
         WRITE(6,'(a)') 'FE/H NOT +1.0, +0.5, 0.0, -0.5, -1.0, -1.5,',
     %' -2.0, -2.5, -3.0, -3.5, OR -4.0'
         RETURN   
      ENDIF
      FNAME = DIR
      FNAME(IPOSN:) = FILES(IFILE)
      OPEN(FILE=FNAME,UNIT=11,STATUS='OLD')
      OPEN(FILE=MODEL,UNIT=4)
      REWIND 11
      REWIND 4
      DO I = 1,1000000
         READ(11,'(A80)',END=900)LINE
         IF (LINE(:5).EQ.'MODEL') THEN
            READ(11,*)TEMP,GRAV
            IF (LOGG.EQ.GRAV  .AND. TEMP.EQ.TEFF) THEN
               ERROR = .FALSE.
               DO J = 1,1000000
                  READ(11,'(A80)',END=900)LINE
                  IF (LINE(:5).EQ.'MODEL') GOTO 900  
                  WRITE(4,'(A80)')LINE
               ENDDO
            ENDIF
         ENDIF
      ENDDO
 900  REWIND 4
      REWIND 11
      CLOSE(UNIT=11)
      CLOSE(UNIT=4)
      IF (ERROR) THEN
         WRITE(6,100)TEFF,LOGG,FNAME
         STOP
      ENDIF
 100  FORMAT('COULD NOT FIND FILE T = ',f6.0,'  LOGG = ',f6.2,' IN ',A80)
      RETURN
      END
      
      

      SUBROUTINE SET_UP_LINES_INPUT_FILE(MICRO,FE,TITLE)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MODTYPE/ MTYPE
      CHARACTER*10 MTYPE
      CHARACTER*80 TITLE,LINE
      REAL*8 MICRO

c I must insert code for LINES to scale the SAD as required for RHOX models "

      OPEN(UNIT=4,FILE='model')
      OPEN(UNIT=7,FILE='FINALMODEL')
      REWIND 4
      WRITE(7,'(A)')'KURTYPE'
      WRITE(7,'(A80)')TITLE

      IF (MTYPE.EQ.'KUROLD') THEN
         WRITE(7,'(13X,2H64)')
      ELSE
         WRITE(7,'(13X,2H72)')
      ENDIF

      WRITE(7,'(6H5000.0)')
      DO I = 1,1000000
         READ(4,'(A44)',END=900)LINE(:44)
         WRITE(7,'(A44)')LINE(:44)
      ENDDO
 900  CONTINUE
      WRITE(7,'(F13.2)')MICRO
      WRITE(7,'(6HNATOMS,8X,1H0,10X,F5.2)')FE
      WRITE(7,'(4HNMOL,10X,1H0)')
      REWIND 4
      CLOSE(UNIT=4,STATUS='DELETE')
      RETURN
      END


