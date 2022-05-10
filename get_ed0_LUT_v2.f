!	This routine reads returns Ed0- LUT which has 4 dimensions for
!		1. Wavelength = 290 : 700 : 5
!		2. ThetaS = 0 : 90 : 5
!		3. Ozone = 200 : 550 : 50
!		4. Cloud optical Thickness = 0 to 64 = c(0,1,2,4,8,16,32,64)
!		5. Surface Albedo = 0.05 : 0.95 : 0.15
!
        INTEGER FUNCTION get_ed0(LUT,choice)
        USE iso_c_binding
!       Define parameters
        PARAMETER (nwl=83, nthetas=19, nozone=10, ntaucl=8, nalb=7)
        LOGICAL THERE
        LOGICAL OK
!       Declare Look-up-table
        CHARACTER PATH*50, fn*100
        INTEGER i,j,k,l,m
        REAL(c_float) LUT(nwl, nthetas, nozone, ntaucl, nalb)
        INTEGER choice

        CALL getenv('OCDATAROOT', PATH )
        PATH=TRIM(PATH)//'/common'
        IF (choice == 1) THEN
          fn=TRIM(PATH)//'/Ed0plus_LUT_5nm_v2.dat'
        ELSE
          fn=TRIM(PATH)//'/Ed0moins_LUT_5nm_v2.dat'
        ENDIF

        INQUIRE( FILE=TRIM(fn), EXIST=THERE )
        IF ( THERE ) THEN
            OPEN(10, FILE=TRIM(fn),status="old",form="formatted",access="sequential")
            INQUIRE( UNIT=10, OPENED=OK )
            IF ( OK ) THEN
                DO i=1,nthetas
                    DO j=1,nozone
                        DO k=1,ntaucl
                           DO l=1,nalb
                                READ(10, *) (LUT(m,i,j,k,l), m=1,nwl)
                           ENDDO
                        ENDDO
                    ENDDO
                ENDDO
                CLOSE(10)
                get_ed0=0
            ELSE
                get_ed0=2
                PRINT *,"Cannot open file to read Ed LUT: ",fn
            ENDIF
        ELSE
            get_ed0=1
            PRINT *,"Cannot find file to read Ed LUT: ",fn
        ENDIF
        RETURN
        END
