!#########################################################################
! GLOBAL THREE DIMENSIONAL DISPERSION MODEL VERSION 3.8
!-------------------------------------------------------------------------
! Finite difference solution of the advection-diffusion equation on a global
! latitude-longitude grid. The model can run on any 1 or 2.5 deg grid with a
! regular grid spacing. The concentration grid is automatically configured
! to match the meteorological grid.  The grid system indicies increase from
! south to north and from west to east. The longitude coordinate system can
! go from 000 to 360 or from -180 to 180. The first and last longitude grid
! points are adjacent to each other but do not overlap (cyclic boundary 
! condition). The south to north coordinate system would run from -90 to +90
! degrees latitude.  The grid's vertical coordinate system is on pressure-
! sigma surfaces. Pollutants are dispersed and advected in terms of their 
! mass mixing ratio which is converted to concentration at STP for output.
! see https://www.hysplit.com/html/kr85_global.html for more information.
!-------------------------------------------------------------------------
! Use the following to compile for a single processor with gfortran:
!     gfortran -ogblmdl38 -O2 -fconvert=big-endian -frecord-marker=4 
!     -ffree-form gblmdl38.f
! Use the following to compile for multiple processors with gfortran:
!     gfortran -ogblmdl38 -O2 -fopenmp -fconvert=big-endian -frecord-marker=4 
!    -ffree-form gblmdl38.f
!-------------------------------------------------------------------------
! Author: roland.draxler@meteozone.com
! This code is distributed under the GNU General Public license 
! https://www.gnu.org/licenses/gpl-3.0.en.html
!-------------------------------------------------------------------------
! REVISION	Version	MODIFICATION ABSTRACT
! 04 Aug 1989 - 1.0	Initial version   
! 19 Feb 1992 - 1.x	380 to 190 km conformal projection data    
! 19 Jun 2003 - 1.y	Hysplit compatible output files 
! 03 May 2004 - 2.0	Non-homogeneous initialization
! 10 Aug 2005 - 3.0	Free format and lat-lon projection meteorology
! 12 Sep 2005 - 3.1	Enhanced convection module
! 15 Sep 2005 - 3.2	Step terrain equations
! 16 Dec 2005 - 3.3	Pressure-sigma coordinate system
! 22 Jun 2006 - 3.4     Monthly emission factors defined by year
! 18 Sep 2006 - 3.5	Delta-height limit test & monthly emission factors
! 05 Dec 2007 - 3.6	Fix to dew point init, output freq, w&rh init
! 28 Mar 2008 - 3.7	Adjoint version
! 18 Apr 2025 - 3.8	Update most input files to standard CSV formats
!-------------------------------------------------------------------------

MODULE constant

! Module with all the constant values such as input/output file
! unit numbers, standard atmosphere values, and mapping conversions

                                             ! FILE UNIT NUMBERS
  INTEGER*4, PARAMETER :: KUMET = 10         ! meteorological data input
  INTEGER*4, PARAMETER :: KSINP = 15         ! simulation parameters input
  INTEGER*4, PARAMETER :: KCOUT = 20         ! packed concentration output
  INTEGER*4, PARAMETER :: KC3DD = 25         ! daily 3D concentration dump
  INTEGER*4, PARAMETER :: KSTNS = 30         ! station file for initialization
  INTEGER*4, PARAMETER :: KQVAL = 35         ! source location and rate data
  INTEGER*4, PARAMETER :: KDIAG = 50         ! diagnositic message output

                                             ! PHYSICAL CONSTANTS
  REAL*8, PARAMETER :: ROW     = 1.17        ! STP air density (g/kg)
  REAL*8, PARAMETER :: REARTH  = 6378140.0   ! radius of the earth (m)
  REAL*8, PARAMETER :: PI      = 3.141592654 ! circumference over diameter
  REAL*8, PARAMETER :: RADPDEG = PI/180.0    ! radians per degree
  REAL*8, PARAMETER :: GRAV    = 9.80665     ! acceleration of gravity (m/s2)
  REAL*8, PARAMETER :: VONK    = 0.4         ! von Karman's constant
  REAL*8, PARAMETER :: P2JM    = 100.0       ! convert mb to J/m3 
  REAL*8, PARAMETER :: RDRY    = 287.04      ! dry air constant (J/Kg-K)

  REAL*8, PARAMETER :: PSFC    = 1013.0      ! default surface pressure
  REAL*8, PARAMETER :: PTOP    =   10.0      ! top of the model pressure
  REAL*8, PARAMETER :: VWLIM   =    0.02     ! vertical velocity limits (m/s)
  REAL*8, PARAMETER :: TWGHT   =    0.40     ! central point terrain smoothing
  REAL*8, PARAMETER :: DZMIN   =  100.00     ! minimum height difference

! standard atmosphere heights at 25 hPa intervals from 1000 hpa to 500 hPa
!                             at 50 hPa intervals from  500 hPa to 100 hPa
!                             at 10 hPa intervals from  100 hPa to  10 hPa

  REAL*8, PARAMETER :: STDATM(38)  = (/                                  & 
                       111.,323.,540.,762.,988.,1220.,1457.,1700.,       &
                       1949.,2204.,2466.,2735.,3012.,3297.,3591.,3894.,  &
                       4206.,4530.,4865.,5213.,5574.,6344.,7185.,8117.,  &
                       9164.,10363.,11784.,13608.,16180.,16848.,17595.,  &
                       18442.,19419.,20576.,22000.,23849.,26481.,31055. /)

  REAL*8, PARAMETER :: STDPPP(38)  = (/                                  & 
                       1000.,975.,950.,925.,900.,875.,850.,825.,800.,    &
                       775.,750.,725.,700.,675.,650.,625.,600.,575.,     &
                       550.,525.,500.,450.,400.,350.,300.,250.,200.,     &
                       150.,100.,90.,80.,70.,60.,50.,40.,30.,20.,10.    /)

END MODULE constant

!###############################################################################

MODULE metarray

! METARRAY - Module that defines the meteorological data arrays
! that are dependent upon the dimensions of the input data

! PPP - PRESSURE                mb
! UUU - U WIND COMPONENTS       m/s
! VVV - V WIND COMPONENTS       m/s
! TTT - TEMPERATURE             Kelvin
! MMM - MOISTURE (RH)           percent
! HHH - HEIGHT MSL OF FIELD     meters
! WWW - VERTICAL VELOCITY       mb/s to m/s
! KKK - VERTICAL MIXING         m2/s
! RRR - LOCAL AIR DENSITY       Kg/m3
! QQQ - MIXING RATIO            g/Kg
! GSX - WE GRID DIMENSIONS      meters
! GSY - SN GRID DIMENSIONS      meters
! GXY - AREA OF GRID CELL       m2
! SFC - SURFACE PRESSURE        mb
! NGP - HORIZONTAL POINTS       number

! meteorological file variables
  INTEGER*4, ALLOCATABLE :: NVAR(:)            ! number variables per level
  INTEGER*4, ALLOCATABLE :: KNDX(:)            ! standard atm index for input

! one-dimensional vertical grid variables
  REAL*8,    ALLOCATABLE :: SIG(:)             ! pressure-sigma coordinate
  REAL*8,    ALLOCATABLE :: PPP(:),DPP(:)      ! pressure and delta-pressure
  REAL*8,    ALLOCATABLE :: TCF(:),POT(:)      ! potential temp conversion and temp

! two-dimensional horizontal grid variables
  REAL*8,    ALLOCATABLE :: SFC(:,:), AVG(:,:) ! surface pressure, average pressure
  REAL*8,    ALLOCATABLE :: GSX(:,:), GSY(:,:) ! grid spacing
  REAL*8,    ALLOCATABLE :: GXY(:,:)           ! grid cell area

! other variables
  INTEGER*4, ALLOCATABLE :: NGP(:)             ! number of X grid points by latitude
  REAL*8,    ALLOCATABLE :: VB4(:)             ! temporary profile holding variable
  REAL*8,    ALLOCATABLE :: VB8(:)             ! temporary profile holding variable

! three-dimensional meteorological variables
  REAL*8,    ALLOCATABLE :: VAL(:,:)           ! temporary variables
  REAL*8,    ALLOCATABLE :: UUU(:,:,:), VVV(:,:,:), TTT(:,:,:)
  REAL*8,    ALLOCATABLE :: MMM(:,:,:), HHH(:,:,:), WWW(:,:,:)
  REAL*8,    ALLOCATABLE :: KKK(:,:,:), RRR(:,:,:)

! three-dimensional pollutant mass
  REAL*8,    ALLOCATABLE :: QQQ(:,:,:)         ! final mixing ratio
  REAL*8,    ALLOCATABLE :: CCC(:,:,:)         ! temporary mixing ratio array

! two-dimensional output concentration
  REAL*4,    ALLOCATABLE :: CONC(:,:)          ! mixing ratio converted to concentration

  SAVE

END MODULE metarray

!###############################################################################

PROGRAM gblmdlv3

  USE constant
  USE metarray

  IMPLICIT NONE

  LOGICAL       :: BACK = .FALSE.
  CHARACTER(4)  :: TRACER        ! pollutant identification
  CHARACTER(80) :: FNAME         ! initial meteo dir/file

  REAL*8        :: QTOT          ! total system masss
  REAL*8        :: DECAY = 0.0   ! decay rate per time step
  REAL*8        :: HMIX          ! horizontal mixing coefficient

  REAL*8        :: CFACT         ! output units (1 Bq/SCM = 27 pCi/SCM)
  REAL*8        :: HALFL         ! pollutant half life (days)
  REAL*8        :: CMIN,CMAX     ! valid concentration range 
  REAL*8        :: VMIX          ! vertical mixing method
  REAL*8        :: WFACT         ! vertical velocity scaling
  REAL*8        :: QBASE = 1.0   ! base source term rate adjustment

  INTEGER*4     :: KDIVG         ! vertical velocity computation method
  INTEGER*4     :: KOUNT = 0     ! integration counter
  INTEGER*4     :: KDSK          ! concentration initialization method
  INTEGER*4     :: KREC          ! meteorology file record number
  INTEGER*4     :: KMET          ! number of records per time period
  INTEGER*4     :: DELTA         ! integration time step (min)
  INTEGER*4     :: LEVEL         ! output level index or hPa height (output:meters)
  INTEGER*4     :: KOUT          ! output level index value
  INTEGER*4     :: LFREQ         ! conc output frequency (hrs)
  INTEGER*4     :: IFREQ         ! initial hour for counting output frequency (0-23)
  INTEGER*4     :: KHRS          ! temporal meteo data frequency (hrs)
  INTEGER*4     :: KMASS         ! mass conservation (0=skip 1=show 2=conserve)
  INTEGER*4     :: QFREQ         ! source term frequency (0:none 1:once >1:step)

  INTEGER       :: DATE_TIME(8)
  CHARACTER(12) :: REAL_CLOCK(3)

  INTEGER*4 :: yr1,mo1,da1       ! beginning simulation date
  INTEGER*4 :: yr2,mo2,da2       ! ending simulation date
  INTEGER*4 :: myr,mmo,mda,mhr   ! last input meteorology date
  INTEGER*4 :: iyr, imo, ida     ! current integration date
  INTEGER*4 :: ihr = 0           ! current integration hour always starts at zero
  INTEGER*4 :: imn = 0           ! current integration minute always starts at zero

! emission rate conversion: 0=none; 1=pBq/y->Ci/min; 2=MCi/y->Bq/min; 3=pBq/y->Bq/min
  INTEGER*4 :: kfact = 3         ! units conversion factor source input

  INTEGER*4 :: i,j,k,l
  INTEGER*4 :: nx,ny,nz
  COMMON /SIZEGRID/ nx,ny,nz     ! meteorology and concentration grid dimensions 

!----------------------------------------------------------
! subroutine interfaces required with dynamic allocation

  INTERFACE
    SUBROUTINE stblty (vmix,delta)
    REAL*8,    INTENT(IN)  :: vmix
    INTEGER*4, INTENT(OUT) :: delta
    END SUBROUTINE stblty

    SUBROUTINE diverg (wfact,kdivg)
    REAL*8,    INTENT(IN) :: wfact
    INTEGER*4, INTENT(IN) :: kdivg
    END SUBROUTINE diverg

    SUBROUTINE dskset (tracer,iyr,imo,ida,ihr,level)
    CHARACTER(4), INTENT(IN) :: tracer
    INTEGER*4,    INTENT(IN) :: iyr,imo,ida,ihr,level
    END SUBROUTINE dskset
 
    SUBROUTINE dskout (tracer,level,kout,iy,im,id,ih,cfact)
    CHARACTER(4), INTENT(IN) :: tracer
    INTEGER*4,    INTENT(IN) :: level,kout,iy,im,id,ih
    REAL*8,       INTENT(IN) :: cfact
    END SUBROUTINE dskout

    SUBROUTINE dayout (IYR,IMO,IDA,IHR,IMN,CFACT)
    INTEGER*4, INTENT(IN) :: iyr,imo,ida,ihr,imn
    REAL*8,    INTENT(IN) :: cfact
    END SUBROUTINE dayout

    SUBROUTINE calend (iyr,imo,ida,ihr,imn)
    INTEGER*4, INTENT(INOUT) :: iyr,imo,ida,ihr,imn
    END SUBROUTINE calend

    SUBROUTINE metinp (fname,myr,mmo,mda,mhr,krec)
    CHARACTER(80), INTENT(INOUT) :: FNAME  
    INTEGER*4,     INTENT(OUT)   :: myr,mmo,mda,mhr
    INTEGER*4,     INTENT(INOUT) :: krec
    END SUBROUTINE metinp

    SUBROUTINE posmet (iyr,imo,ida,ihr,krec)
    INTEGER*4,  INTENT(IN)  :: iyr,imo,ida,ihr
    INTEGER*4,  INTENT(OUT) :: krec
    END SUBROUTINE posmet

    SUBROUTINE chkmet (fname,khrs,kmet)
    CHARACTER(80), INTENT(IN)  :: fname
    INTEGER*4,     INTENT(OUT) :: khrs
    INTEGER*4,     INTENT(OUT) :: kmet
    END SUBROUTINE chkmet

    SUBROUTINE conset (kdsk,iyr,imo,ida,cfact,cmin,cmax)
    INTEGER*4, INTENT(IN) :: kdsk,iyr,imo,ida
    REAL*8,    INTENT(IN) :: cfact,cmin,cmax
    END SUBROUTINE conset

    SUBROUTINE eqntns (decay,hmix,delta,qtot,kmass)
    REAL*8,    INTENT(IN)  :: decay,hmix
    INTEGER*4, INTENT(IN)  :: delta,kmass
    REAL*8,    INTENT(OUT) :: qtot
    END SUBROUTINE eqntns

    SUBROUTINE source (iyr,imo,ida,ihr,imn,delta,kfact,qbase)
    INTEGER*4,  INTENT(IN) :: iyr,imo,ida,ihr,imn
    INTEGER*4,  INTENT(IN) :: delta
    INTEGER*4,  INTENT(IN) :: kfact
    REAL*8,     INTENT(IN) :: qbase
    END SUBROUTINE source

    SUBROUTINE xcross (cfact)
    REAL*8, INTENT(IN) :: cfact
    END SUBROUTINE xcross

  END INTERFACE

!---------------------------------------------------------
! configure simulation from input file

  OPEN  (KDIAG, FILE='message.txt')

  OPEN  (KSINP, FILE='default.dat')

  READ  (KSINP,*)  YR1, MO1, DA1       ! START DATE (inclusive)
  READ  (KSINP,*)  YR2, MO2, DA2       ! stop date  (exclusive)

  READ  (KSINP,*)  KDSK                ! INITIALIZATION
                                       ! =0 initial values all zero
                                       ! =1 avg latitude bands: startup.txt
                                       ! =2 interpolation from: startup.txt
                                       ! =3 from 3D dump file: gbl3dim.bin
                                       ! >3 the value given by KDSK in tenths
  READ  (KSINP,*)  CMIN,CMAX           ! valid concentration range (kdsk = 1 or 2)

  READ  (KSINP,*)  HMIX                ! HORIZONTAL MIXING COEFFICIENT
  READ  (KSINP,*)  WFACT               ! vertical velocity scaling (0.0 -> 1.0)
  READ  (KSINP,*)  KDIVG               ! vertical velocity method (0=data 1=divergence)
  READ  (KSINP,*)  VMIX                ! VERTICAL mixing scaling value (~50)

  READ  (KSINP,*)  KMASS               ! Mass (0=no 1=show 2=conserve)
  READ  (KSINP,*)  QFREQ               ! Source term frequency (0:none 1:once >=2:step)
  READ  (KSINP,*)  QBASE               ! Baseline source rate adjustment (default = 1.0)

  READ  (KSINP,*)  LFREQ               ! CONCENTRATION OUTPUT FREQUENCY (hours <=24)
  READ  (KSINP,*)  IFREQ               ! initial hour for output (0Z-23Z)
  ifreq = MAX(0, MIN(23,ifreq))
  IF (lfreq.GT.24) lfreq=MOD(lfreq-1,24)+1

  READ  (KSINP,*)  LEVEL               ! output values at index level or hPa pressure
  READ  (KSINP,*)  CFACT               ! units conversion from 'X' to output
  READ  (KSINP,*)  HALFL               ! pollutant half life (days)
  READ  (KSINP,*)  TRACER              ! pollutant 4-character identification

  READ  (KSINP,*)  FNAME               ! INITIAL METEOROLOGY DATA DIR/FILE (quoted)
                                       ! new line required for each file
 
  iyr=yr1                              ! Set the internal integration date from input
  imo=mo1                              ! starting date. Integration proceeds until
  ida=da1                              ! internal date equals the end date.

! very simple test may need to be disabled at times
! IF(yr1*366+mo1*31+da1.GT.yr2*366+mo2*31+da2) back=.TRUE.                 

!-----------------------------------------------------
! open meteorology file and create data arrays

  CALL chkmet (fname,khrs,kmet)        ! determine if meteo file can be used
  CALL setgrd                          ! compute horizontal grid spacing 
  CALL posmet (iyr,imo,ida,ihr,krec)   ! position meteorology to first record

! read all the data for the first period
  CALL metinp (fname,myr,mmo,mda,mhr,krec) 
  IF(iyr.NE.myr.OR.imo.NE.mmo.OR.ida.NE.mda.OR.ihr.NE.mhr)THEN
     WRITE(kdiag,'(A,4I3)')'ERROR: meteorology time not correct - ',myr,mmo,mda,mhr
     CLOSE(kdiag)
     STOP 900
  ELSE
     WRITE(kdiag,'(A,4I3)')'Meteorology positioned:',iyr,imo,ida,ihr
     WRITE(kdiag,'(A,I6)') '..... to record number:',krec
  END IF

  CALL metdat
  CALL diverg (wfact,kdivg)            ! convert vertical velocity from p to z
  CALL stblty (vmix,delta)             ! compute vertical mixing and time step
  IF(back) delta=-delta

! compute the output index level and convert pressure to height ...
! the height is computed as the top-of-the-box and the coordinate system is such
! that the meteorology points always represent the center of the grid box
  IF(level.LT.nz)THEN
     kout=level
  ELSE
     kout=1
     DO WHILE (ppp(kout).GE.level.AND.kout.LT.nz)
        kout=kout+1
     END DO
  END IF
! level=0.5*(SUM(hhh(:,:,kout))+SUM(hhh(:,:,MIN(kout+1,nz))))/nx/ny
  level=SUM(hhh(:,:,kout))/nx/ny
  CALL dskset (tracer,iyr,imo,ida,ihr,level)  ! open the concentration output files

! daily model dump of 3D concentration fields
  OPEN(kc3dd,FILE='gbl3dim.bin',FORM='UNFORMATTED')

! initialization of the concentration field
  CALL conset (kdsk,iyr,imo,ida,cfact,cmin,cmax)
  CALL xcross (cfact)

! write initial surface concentration distribution to file
  CALL dskout (tracer,level,kout,iyr,imo,ida,ihr,cfact)

! decay (1/s) ... delta (min) ... half-life (days)
  IF(halfl.GT.0.0) decay=(1.0-DEXP(abs(DBLE(delta))*DLOG(DBLE(0.5)) &
                         /abs(halfl)/DBLE(1440.0)))                 &
                         /abs(DBLE(delta))/DBLE(60.0)

! current processor clock time
  WRITE(KDIAG,'(A)')'Simulation start ...'
  CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),DATE_TIME)
  WRITE(kdiag,'(2A)')'Start Date (CCYYMMDD): ',REAL_CLOCK(1)
  WRITE(kdiag,'(2A)')'Start Time (HHMMSS.S): ',REAL_CLOCK(2)

!-------------------------------------------------------------
! Main time integration loop

  DO WHILE (.NOT.(IYR.EQ.YR2.AND.IMO.EQ.MO2.AND.IDA.EQ.DA2))
     KOUNT=KOUNT+1

!    meteorological data input and processing
     IF(IMN.EQ.0.AND.MOD(IHR+(KHRS/2),KHRS).EQ.0)THEN
        IF(back) krec=krec-2*kmet
        CALL metinp (fname,myr,mmo,mda,mhr,krec)
        CALL metdat
        CALL diverg (wfact,kdivg)
        CALL stblty (vmix,delta)
        IF(back) delta=-delta
        IF(halfl.GT.0.0) decay=(1.0-DEXP(abs(DBLE(delta))*DLOG(DBLE(0.5)) &
                               /abs(halfl)/DBLE(1440.0)))                 &
                               /abs(DBLE(delta))/DBLE(60.0)
     END IF

!    finite difference equations
     IF(qfreq.GE.2.OR.(qfreq.EQ.1.AND.kount.EQ.1))   &
        CALL source (iyr,imo,ida,ihr,imn,delta,kfact,qbase)
     CALL eqntns (decay,hmix,delta,qtot,kmass)

!    increment time
     imn=imn+delta
     CALL calend (iyr,imo,ida,ihr,imn)

!    disk output frequency tests and also open a new file each year
     IF(imo.eq.1.AND.ida.EQ.1.AND.ihr.EQ.0.AND.imn.EQ.0)THEN
        CLOSE(kcout)
        CALL dskset (tracer,iyr,imo,ida,ihr,level)
     END IF
     IF(imn.EQ.0.AND.lfreq.GT.0.AND.MOD(ihr-ifreq,lfreq).EQ.0) &
        CALL dskout (tracer,level,kout,iyr,imo,ida,ihr,cfact)

!    once per month diagnostic dump of the concentration latitudinal gradients
     IF(ida.EQ.1.AND.ihr.EQ.0.AND.imn.EQ.0) CALL xcross (cfact)
    
!    end-of-day processing
     IF(ihr.EQ.0.AND.imn.EQ.0)THEN
!       daily mass inventory
        WRITE(kdiag,'(A,5I2,A,4I2,A,I3,A,E15.8 )') 'Model: ',iyr,imo,ida,ihr,imn,   &
        '  Meteo: ',myr,mmo,mda,mhr,'    Step: ',delta,'         Mass: ',qtot
        CALL dayout (IYR,IMO,IDA,IHR,IMN,CFACT)
     END IF

!    once per hour message to standard output
     IF(imn.EQ.0) WRITE(*,'(A,5I2,E20.8)')'Finished: ',iyr,imo,ida,ihr,imn,qtot
  END DO

  WRITE(KDIAG,'(A)')'Simulation complete'
  CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),DATE_TIME)
  WRITE(kdiag,'(2A)')'Date (CCYYMMDD): ',REAL_CLOCK(1)
  WRITE(kdiag,'(2A)')'Time (HHMMSS.S): ',REAL_CLOCK(2)

  DEALLOCATE (nvar,ppp,sig,pot,vb4,vb8,gsx,gsy,gxy,sfc)
  DEALLOCATE (uuu,vvv,ttt,mmm,hhh,www,kkk,rrr,qqq,ccc,conc)
  CLOSE(kdiag)

END PROGRAM gblmdlv3


!###############################################################################
! SOURCE - Sets the emission rate each time step at appropriate locations. The
! emission rate data array is initialized from a text file containg the lat,
! lon, height, and emission rate in mass units per minute for each location. 
! If the emission rate file is not found, the model will be run with no 
! emissions.
!-------------------------------------------------------------------------------
! LAST REVISED: 19 Feb 1992 - initial version
!               16 Aug 2005 - fortran90 update
!               07 Sep 2005 - multi year emission files
!               26 Dec 2005 - monthly correction (annual total unchanged)
!               18 Sep 2006 - monthly emissions factors to input file
!               28 Mar 2008 - absolute value of time step for backward
!               26 Jun 2008 - species conversion MCi - PBq
!               08 Sep 2008 - half grid point correction lat,lon to i,j
!               13 Aug 2013 - revised test of monthly emission factors
!               27 Mar 2025 - convert to reading annual data as CSV 
!-------------------------------------------------------------------------------

SUBROUTINE source (iyr,imo,ida,ihr,imn,delta,kfact,qbase)

  USE constant
  USE metarray

  IMPLICIT NONE

  INTEGER*4,  INTENT(IN) :: iyr,imo,ida,ihr,imn
  INTEGER*4,  INTENT(IN) :: delta                 ! integration time step (min)
  INTEGER*4,  INTENT(IN) :: kfact                 ! emission rate conversion
                                                  ! 0 = none
                                                  ! 1 = pBq/y -> Ci/min
                                                  ! 2 = MCi/y -> Bq/min
                                                  ! 3 = pBq/y -> Bq/min 
  REAL*8,     INTENT(IN) :: qbase

  INTEGER*4, ALLOCATABLE :: iq(:),jq(:),kq(:)     ! indices of emission  point
  INTEGER*4, ALLOCATABLE :: monq(:,:)             ! monthly q adjustment index
  REAL*8,    ALLOCATABLE :: rate(:,:)             ! emission rate (/min) by year
  REAL*8,    ALLOCATABLE :: qadj(:,:)             ! monthly emission factor profiles

  CHARACTER*512 :: line

  CHARACTER*4 :: site
  REAL*8      :: depth, volume
  INTEGER*4   :: i,j,k,l,numq,nump,kret
  INTEGER*4   :: year,year1,year2

  LOGICAL     :: mfile  = .false. ! monthly emission adjustment factors
  LOGICAL     :: qfile  = .false. ! annual emission totals
  LOGICAL     :: poslon = .true.  ! all positive longitudes

  INTEGER*4   :: clat,clon,nx,ny,nz
  REAL*4      :: clat1,clon1,dlat,dlon

  COMMON /MAINGRID/ clat1,clon1,dlat,dlon
  COMMON /SIZEGRID/ nx,ny,nz

  SAVE qfile,numq,iq,jq,kq,rate,monq,qadj

! compute the base emission rate at the start of the simulation
  year=1900+iyr
  if(iyr.LT.40)year=2000+iyr

!-------------------------------------------------------------------
! one time initialization to read yearly emission values by location

  IF(.NOT.qfile)THEN

!    monthly emission factors

     nump=0
     INQUIRE(FILE='monthly.csv',EXIST=mfile)
     IF(.NOT.mfile)THEN
        WRITE(kdiag,'(A)')'Monthly emission factors file not found ... using default!'
        ALLOCATE (qadj(12,1))
        qadj=1.0

     ELSE
        OPEN (kqval,FILE='monthly.csv')
        kret=0
        DO WHILE (kret.EQ.0)
           READ(kqval,*,IOSTAT=KRET)
           IF(kret.EQ.0) nump=nump+1
        END DO
        REWIND(kqval)
        WRITE(kdiag,'(A,I5)')'Number of monthly emission profiles defined: ', nump
        ALLOCATE (qadj(12,nump))

        DO L=1,nump
           READ(kqval,*) (qadj(K,L),K=1,12)
           IF(SUM(qadj(:,L)).NE.12.0) WRITE(kdiag,'(A,I5)')'Emission factor sum<>12 on line: ',L   
        END DO
        CLOSE(kqval)
     END IF

!    yearly emission totals (format depends upon availability of monthly)

     INQUIRE(FILE='emission.csv',EXIST=qfile)
     IF(.NOT.qfile)THEN
        WRITE(kdiag,'(A)')'ERROR: emission.csv rate file not found!'
        CLOSE(kdiag)
        STOP 905
     END IF

!    determine the number of columns (i.e. years)
     OPEN (kqval,FILE='emission.csv')
     READ(kqval,'(A)') line

     year1=0      ! start year
     year2=0      ! end year
     k=999        ! comma character position
     l=0          ! number of data columns found
     i=1          ! start character number
     DO WHILE (k.GT.0)
        k=INDEX(line,',')
        IF(k.GT.0)THEN
           l=l+1
           IF(l.GE.4) READ(line,'(I4)')year
           IF(l.EQ.4) year1=year
           i=k+1 
           line=line(i:)  
        END IF
     END DO
     READ(line,'(I4)')year2

     kret=0
     numq=0
     DO WHILE (kret.EQ.0)
        READ(kqval,*,IOSTAT=KRET) site
        IF(kret.EQ.0) numq=numq+1
     END DO
     REWIND(kqval)

     ALLOCATE (iq(numq),jq(numq),kq(numq))
     ALLOCATE (rate(numq,year1:year2))
     ALLOCATE (monq(numq,year1:year2))
     monq = 1  

!    set the longitude system flag
     IF(clon1.LT.0.0)poslon=.false.

     READ(kqval,*)     ! skip first record
     DO L=1,numq
        IF(mfile)THEN
           READ(kqval,*,IOSTAT=KRET) site, clat, clon, (monq(L,K),rate(L,K),K=year1,year2)
           IF(KRET.NE.0)THEN
              WRITE(kdiag,'(A)')'ERROR: reading emissions ... monthly factors required!'
              CLOSE(kdiag)
              STOP
           END IF

        ELSE
           READ(kqval,*,IOSTAT=KRET) site, clat, clon, (rate(L,K),K=year1,year2)
           IF(KRET.NE.0)THEN
              WRITE(kdiag,'(A)')'ERROR: reading emissions ... monthly factors not required!'
              CLOSE(kdiag)
              STOP
           END IF
        END IF
        
!       adjust the longitude system to reflect meteo data
        IF(poslon)THEN
           IF(clon.LT.0)clon=clon+360
        ELSE
           IF(clon.GT.180)clon=360-clon
        END IF

!       save grid indicies of each source location
        jq(L)=1+(clat-clat1+dlat/2.0)/dlat
        iq(L)=1+(clon-clon1+dlon/2.0)/dlon

!       emissions forced from the lowest grid box
        kq(L)=1

     END DO
     CLOSE(kqval)

!    test if there is a monthly emission profile for each specified in emission.csv
     IF(nump.GT.MAXVAL(monq))THEN
        WRITE(*,*)'Insufficient number of monthly emission profiles defined!'
        WRITE(*,*)'Number of definitions defined in monthly factors file:',nump
        WRITE(*,*)'Maximum monthly index defined in yearly emission file:',MAXVAL(monq)
        CLOSE(kdiag)
        STOP
     END IF  

!    emission rate conversion: 0=none; 1=pBq/y->Ci/min; 2=MCi/y->Bq/min; pBq/y->Bq/min
!    1 Ci  = 37 x 10^9 Bq ....... 1 Bq = 27 x 10^-12 Ci
!    1 MCi = 37 pBq ............. 1 Bq = 27 pCi
         
     IF(kfact.EQ.1)THEN
        rate = (1.0E+06)*rate/37.0  ! pBq/y -> Ci/year
     ELSEIF(kfact.EQ.2)THEN
        rate = (1.0E+15)*rate*37.0  ! MCi/y -> Bq/year
     ELSEIF(kfact.EQ.3)THEN
        rate = (1.0E+15)*rate       ! pBq/y -> Bq/year
     ELSE
        CONTINUE
     END IF
     rate = qbase * rate / DBLE(525960.0)        ! units/year -> units/minute

     WRITE(kdiag,'(A,E10.3)')'Emission file total (mass/min): ',SUM(rate(:,year))
     WRITE(kdiag,'(A,I5)')'Total number of sources: ', numq
     IF(SUM(rate).EQ.0.0) numq=0

  END IF

!--------------------------------------------------
! add emissions to each grid point

  IF(numq.EQ.0)RETURN

  DO L=1,numq
     i=iq(L)
     j=jq(L)
     k=kq(L)

     IF(kq(L).LT.nz)THEN
        depth=HHH(i,j,k+1)-HHH(i,j,k)
     ELSE
        depth=HHH(i,j,k)-HHH(i,j,k-1)
     END IF

!    grid cell volume
     volume=GXY(i,j)*depth

!    add emissions to each source grid point
     QQQ(i,j,k)=QQQ(i,j,k)+rate(L,year)*abs(DBLE(DELTA))*qadj(imo,monq(L,year))/(volume*RRR(i,j,k))
   END DO

!  once a month diagnostic dump of the total emissions 
   IF(ida.EQ.1.AND.ihr.EQ.0.AND.imn.EQ.0)THEN
      WRITE(kdiag,'(A,E10.3)')'Emission file total (mass/min): ',SUM(rate(:,year))
      WRITE(kdiag,'(A,50F5.2)')'Q-factors: ',(qadj(imo,monq(L,year)),L=1,numq)
   END IF
      
END SUBROUTINE source


!###############################################################################
! EQNTNS - Main grid finite difference equations gives the solution for all
! grid points in x except the first (south pole) and last (north pole) y points.
! The south and north pole regions are solved in subroutine bounds.
!-------------------------------------------------------------------------------
! PRIMARY METEOROLOGICAL VARIABLES
! UUU - U WIND COMPONENTS    M/S
! VVV - V WIND COMPONENTS    M/S
! TTT - TEMPERATURE          KELVIN
! MMM - MOISTURE (RH)        PERCENT
! HHH - HEIGHT MSL OF FIELD  METERS
! WWW - INPUT AS DIVERGENCE  1/SEC
! KKK - VERTICAL MIXING      M2/SEC
! RRR - LOCAL AIR DENSITY    KG/M3
! QQQ - MIXING RATIO         g/KG
!-------------------------------------------------------------------------------
! LAST REVISED: 12 Apr 1992 - initial version
!               16 Aug 2005 - fortran90 update
!               18 Sep 2006 - delta height tests
!               28 Mar 2008 - absolute value of time step for backward
!                           - integration direction test for gradients
!-------------------------------------------------------------------------------

SUBROUTINE eqntns (DECAY,HMIX,DELTA,QTOT,KMASS)

  USE constant
  USE metarray

! this allows use of omp library functions and subroutines
  use omp_lib

  IMPLICIT NONE

  REAL*8,    INTENT(IN)  :: decay,hmix
  INTEGER*4, INTENT(IN)  :: delta,kmass
  REAL*8,    INTENT(OUT) :: qtot

  REAL*8    :: HM                     ! adjusted horizontal mixing
  REAL*8    :: XDYD                   ! ratio of (E-W)/(N-S) grid distances  
  REAL*8    :: QSUM                   ! mass summation
  REAL*8    :: DQDT                   ! rate of mass change
  REAL*8    :: DIST,DZI,DZT,DZB       ! various grid distances
  REAL*8    :: TRX,TRY,TRZ            ! transport (advection) fluxes
  REAL*8    :: DFX,DFY,DFT,DFB        ! diffusion fluxes

  INTEGER*4 :: numb                   ! flux averaging grid point counter
  INTEGER*4 :: i,j,k,m                ! basic loop indicies for 3D grid
  INTEGER*4 :: ip,ii                  ! loop indices for sub-grid summation
  INTEGER*4 :: im1,ip1                ! plus and minus one indicies
  INTEGER*4 :: nx,ny,nz               ! 3D grid dimensions

  COMMON /SIZEGRID/ nx,ny,nz

! initial system mass prior to computation
  qtot = 0.0
  qsum = 0.0   

!------------------------------
! in parallelization section below, to force 8 threads... 
! PARALLEL NUM_THREADS(8) 
!------------------------------
!$OMP PARALLEL &
!$OMP SHARED(CCC,UUU,QQQ,VVV,WWW,KKK,HHH) &
!$OMP SHARED(GSX,GSY,NGP) &
!$OMP SHARED(HMIX,NX,NY,NZ,DELTA) &
!$OMP PRIVATE(K,J,I,ii,ip,IM1,IP1) &
!$OMP PRIVATE(TRX,DFX,TRY,DFY,TRZ,DFT,DFB,DZT,DZI,DZB) &
!$OMP PRIVATE(XDYD,HM,DIST,numb,DQDT) 
!$OMP DO REDUCTION(+:qsum,qtot)
!------------------------------
       
  DO K=1,nz
  DO J=1,ny

! skip to the first horizontal grid point of each multi-meteo concentration cell
! process the first and last point for east-west fluxes

  DO I=1,nx,ngp(j)

!    ratio is just COS(latitude) - horizontal mixing reduces with latitude
     XDYD=GSX(i,j)/GSY(i,j)
     HM=HMIX*XDYD

!    tropical convection enhancement (0.87 = 30 deg lat)
     if(xdyd.ge.0.87)hm=hm*2.0

!    Except near the poles, the number of concentration grid points (ngp) per
!    meteorological grid point is one.

     TRX=0.0
     DFX=0.0

     IF(j.GT.1.AND.j.LT.ny)THEN
!       Horizontal advection use upwind gradient and for multi-meteorology 
!       concentration grid cells, advection is computed from the W-E ends.
!       The polar cell is excluded for W-E advection.

        DIST=GSX(I,J)*DBLE(NGP(J))     ! adjusted grid distance

        IM1=I-1                        ! index minus one
        IP1=I+1                        ! index plus one
        IF(I.EQ.1)IM1=nx               ! cyclic boundary adjustments
        IF(I.EQ.nx)IP1=1

!       west to east advection
        IF(UUU(I,J,K)*DELTA.GT.0.0)THEN
           TRX=UUU(I,J,K)*(QQQ(I,J,K)-QQQ(IM1,J,K))/DIST
        ELSEIF(UUU(I,J,K)*DELTA.LT.0.0)THEN
           TRX=UUU(I,J,K)*(QQQ(IP1,J,K)-QQQ(I,J,K))/DIST
        END IF

!       west to east horizontal diffusion
        DFX=HM*(QQQ(IP1,J,K)-2.0*QQQ(I,J,K)+QQQ(IM1,J,K))/DIST/DIST

!       set the number of meteo cells to average per conc cell
        numb=ngp(j)
     ELSE
!       at the poles process only one cell
        numb=1
     END IF

     TRY=0.0
     DFY=0.0
     TRZ=0.0
     DFT=0.0
     DFB=0.0

     ii=0
     DO WHILE(ii.LT.numb)

!       Near the poles where ngp>1 the fluxes are averaged for each
!       concentration grid cell.

        ii=ii+1
        ip=ii+i-1

!       south to north advection and diffusion

        IF(j.GT.1.AND.j.LT.ny)THEN
!          interior of the grid       
           IF(VVV(IP,J,K)*DELTA.GT.0.0)THEN
              TRY=TRY+VVV(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP,J-1,K))/GSY(IP,J)
           ELSEIF(VVV(I,J,K)*DELTA.LT.0.0)THEN
              TRY=TRY+VVV(IP,J,K)*(QQQ(IP,J+1,K)-QQQ(IP,J,K))/GSY(IP,J)
           END IF
           DFY=DFY+HM*(QQQ(IP,J+1,K)-2.0*QQQ(IP,J,K)+QQQ(IP,J-1,K))/GSY(IP,J)/GSY(IP,J)

!       At the poles advection and diffusion is computed only for one grid cell. At other
!       grid cells rate terms are computed as the average of all meteo cells that are in
!       that concentration grid cell. The gradients are computed across the pole such that
!       a point at (I,J+1) or (I,J-1) is the same as a point at (I+NX/2,J) 

        ELSEIF(j.EQ.1.AND.IP.EQ.1)THEN
!          south pole
           IF(VVV(IP,J,K)*DELTA.GT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP+NX/2,J+1,K))/GSY(IP,J)
           ELSEIF(VVV(I,J,K)*DELTA.LT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP,J+1,K)-QQQ(IP,J,K))/GSY(IP,J)
           END IF
           DFY=HM*(QQQ(IP,J+1,K)-2.0*QQQ(IP,J,K)+QQQ(IP+NX/2,J+1,K))/GSY(IP,J)/GSY(IP,J)

        ELSEIF(J.EQ.NY.AND.IP.EQ.1)THEN
!         north pole
           IF(VVV(IP,J,K)*DELTA.GT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP,J-1,K))/GSY(IP,J)
           ELSEIF(VVV(I,J,K)*DELTA.LT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP+NX/2,J-1,K)-QQQ(IP,J,K))/GSY(IP,J)
           END IF
           DFY=HM*(QQQ(IP,J-1,K)-2.0*QQQ(IP,J,K)+QQQ(IP+NX/2,J-1,K))/GSY(IP,J)/GSY(IP,J)
        END IF
 
!       vertical diffusion and advection

        IF (k.GT.1.AND.k.LT.nz) THEN
!          interior of the grid
           DZT=MAX(DZMIN,HHH(IP,J,K+1)-HHH(IP,J,K))
           DZI=MAX(DZMIN,0.5*(HHH(IP,J,K+1)-HHH(IP,J,K-1)))
           DZB=MAX(DZMIN,HHH(IP,J,K)-HHH(IP,J,K-1))

           IF(WWW(IP,J,K)*DELTA.GT.0.0)THEN
              TRZ=TRZ+WWW(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP,J,K-1))/DZB
           ELSE
              TRZ=TRZ+WWW(IP,J,K)*(QQQ(IP,J,K+1)-QQQ(IP,J,K))/DZT
           END IF
           DFT=DFT+KKK(IP,J,K)*(QQQ(IP,J,K+1)-QQQ(IP,J,K))/DZT/DZI
           DFB=DFB+KKK(IP,J,K-1)*(QQQ(IP,J,K)-QQQ(IP,J,K-1))/DZB/DZI

        ELSEIF (k.EQ.1) THEN
!          bottom boundary
           DZI=MAX(DZMIN,HHH(IP,J,2)-HHH(IP,J,1))
           TRZ=TRZ+WWW(IP,J,2)*(QQQ(IP,J,2)-QQQ(IP,J,1))/DZI
           DFT=DFT+KKK(IP,J,2)*(QQQ(IP,J,2)-QQQ(IP,J,1))/DZI/DZI

        ELSE
!          top boundary
           DZI=MAX(DZMIN,HHH(IP,J,nz)-HHH(IP,J,nz-1))
           TRZ=TRZ+WWW(IP,J,nz)*(QQQ(IP,J,nz)-QQQ(IP,J,nz-1))/DZI
           DFB=DFB+KKK(IP,J,nz-1)*(QQQ(IP,J,nz)-QQQ(IP,J,nz-1))/DZI/DZI
        END IF

     END DO

!    add all final derivative terms
     IF(DELTA.GT.0.0)THEN
        DQDT=((DFX-TRX)+(DFY-TRY-TRZ+DFT-DFB)/DBLE(numb))*DBLE(DELTA)*DBLE(60.0)
     ELSE
        DQDT=((DFX+TRX)+(DFY+TRY+TRZ+DFT-DFB)/DBLE(numb))*abs(DBLE(DELTA))*DBLE(60.0)
     END IF

!    and update the concentration working array
     CCC(I,J,K)=QQQ(I,J,K)+DQDT

!    any #instability# go back to the orignal value
     IF(CCC(I,J,K).LT.0.0) CCC(I,J,K)=QQQ(I,J,K)

!    for computational purposes multi-meteo-cell concentration grid cells 
!    near the poles all have the same concentration values
     IF(ngp(j).GT.1) CCC(I+1:I+NGP(J)-1,J,K) = CCC(i,j,k)

!    diagnostic mass summations
     IF(kmass.GT.0)THEN
        DO IP=I,I+NGP(J)-1
           IF(k.LT.nz)THEN
              DZI=HHH(ip,j,k+1)-HHH(ip,j,k)
           ELSE
              DZI=HHH(ip,j,k)-HHH(ip,j,k-1)
           END IF
           qtot=qtot+qqq(ip,j,k)*RRR(ip,j,k)*GXY(ip,j)*DZI   ! before differencing     
           qsum=qsum+ccc(ip,j,k)*RRR(ip,j,k)*GXY(ip,j)*DZI   ! after  differencing 

!          error check (#effectively disabled with previous instability test#)
           IF(ccc(ip,j,k).LT.0.0)THEN
              WRITE(*,*)'Concentration <0: ',ccc(ip,j,k)
              WRITE(*,*)'At grid location: ',ip,j,k
              WRITE(*,*)'Surface pressure: ',sfc(ip,j)
              WRITE(*,'(7A10)')'sig','ppp','hhh','ttt','www','mmm','rrr'
              DO m=nz,1,-1
                 WRITE(*,'(7E10.3)') sig(m),ppp(m),hhh(ip,j,m),ttt(ip,j,m),www(ip,j,m), &
                                                   mmm(ip,j,m),rrr(ip,j,m)
              END DO 
              READ(*,*)
           END IF
        END DO
     END IF

  END DO
  END DO
  END DO

!----------------------
!$OMP END DO
!$OMP END PARALLEL
!----------------------

! computation mass adjustment factor to maintain system mass
  IF(kmass.EQ.2.AND.qsum.NE.0.0) ccc = qtot*ccc/qsum

! copy back to master array with decay adjustment  
  qqq  = ccc  - decay*abs(DBLE(delta))*DBLE(60.0)*ccc
      
END SUBROUTINE eqntns


!###############################################################################
! CONSET - Set the intial concentration values according to several methods.
! =0 no initialization all zero
! =1 by latitude bands: station.txt
! =2 pos interpolation: station.txt
! =3 from 3D dump file: gbl3dim.bin
! >3 from single value: float(kdsk)
!-------------------------------------------------------------------------------
! LAST REVISED: 05 Mar 2004 - initial version
!               16 Aug 2005 - fortran90 upgrade
!               21 Feb 2006 - 3D initialization patch
!               22 Sep 2006 - array correction
!               27 Nov 2007 - bands and dew pnt init correction
!               29 Mar 2025 - constant initialization input as milli units
!-------------------------------------------------------------------------------

SUBROUTINE conset (KDSK,IYR,IMO,IDA,CFACT,CMIN,CMAX)

  USE constant
  USE metarray

  IMPLICIT NONE

  INTEGER*4, INTENT(IN)  :: KDSK             ! initialization method
  INTEGER*4, INTENT(IN)  :: IYR,IMO,IDA      ! initialization time
  REAL*8,    INTENT(IN)  :: CFACT            ! concentration units conversion
  REAL*8,    INTENT(IN)  :: CMIN,CMAX        ! concentration range limits

  REAL*8,    ALLOCATABLE :: xcon(:)          ! concentration sampling data
  REAL*8,    ALLOCATABLE :: xvar(:), yvar(:) ! working array for regression
  INTEGER*4, ALLOCATABLE :: clat(:), clon(:) ! sampler positions

  LOGICAL   :: poslon = .true.               ! all positive longitudes
  LOGICAL   :: blend  = .false.              ! set the data blending flag

  REAL*8    :: height,conval                 ! binary input file values
  REAL*8    :: es,ea,dp                      ! dewpoint computation variables
  REAL*8    :: yintc_nh,slope_nh             ! linear regression
  REAL*8    :: yintc_sh,slope_sh
  REAL*8    :: cini,cini_nh,cini_sh
  REAL*8    :: xlat,xlon                     ! temporary position

  INTEGER*4 :: i,j,k,l,ii,jj,kk              ! grid indicies
  INTEGER*4 :: kret,numb                     ! return code and number of input records
  INTEGER*4 :: kyr,kmo,kda,khr,kmn           ! input data date

  INTEGER*4 :: nxg,nyg,nzg                   ! dimensions of 3D input grid
  REAL*8    :: clatg,clong,glat,glon         ! grid corners and spacing
  REAL*8,   ALLOCATABLE :: ggg(:)            ! vertical levels for input data
  REAL*8,   ALLOCATABLE :: gqc(:,:,:)        ! temp 3D mixing ratio array

  INTEGER*4 :: nx,ny,nz                      ! dimensions of current simulation grid
  REAL*4    :: clat1,clon1,dlat,dlon

  REAL*8            :: pres
  REAL*8, PARAMETER :: cbot = 0.87           ! lower stratosphere mass correction (0.87)
  REAL*8, PARAMETER :: ctop = 0.54           ! upper stratosphere mass correction (0.54)

  COMMON /MAINGRID/ clat1,clon1,dlat,dlon
  COMMON /SIZEGRID/ nx,ny,nz

!---------------------------------------------------

  INTERFACE
    SUBROUTINE linreg (xvar,yvar,slope,yintc)
    REAL*8, INTENT(IN)  :: xvar(:), yvar(:)
    REAL*8, INTENT(OUT) :: slope, yintc
    END SUBROUTINE linreg
  END INTERFACE

!---------------------------------------------------
! read station file for initialization

  IF(kdsk.EQ.1.OR.kdsk.EQ.2)THEN

     OPEN (kstns,FILE='startup.txt')
     kret=0
     numb=0
     DO WHILE (kret.EQ.0)
        READ(kstns,*,IOSTAT=KRET) glat, glon
        IF(kret.EQ.0) numb=numb+1
     END DO
     REWIND(kstns)
     ALLOCATE (xcon(numb),clat(numb),clon(numb))

!    set the longitude system flag
     IF(clon1.LT.0.0)poslon=.false.

     jj=0
     kk=0
!    concentration always the dependent variable
     DO L=1,numb
        READ(kstns,*)CLAT(L),CLON(L),XCON(L)
        IF(clat(L).LE.0)jj=jj+1 
        IF(clat(L).GT.0)kk=kk+1

!       adjust the longitude system to reflect meteo data
        IF(poslon)THEN
           IF(clon(L).LT.0)clon(L)=clon(L)+360
        ELSE
           IF(clon(L).GT.180)clon(L)=360-clon(L)
        END IF

     END DO
     CLOSE(kstns)

     WRITE(kdiag,'(A,I5)')'Station initialization from file (#): ',numb
     IF(jj.LT.2.OR.kk.LT.2)THEN
        WRITE(*,*)'Insufficient number in {startup.txt} for initialization!'
        WRITE(*,*)'Northern hemisphere: ',kk
        WRITE(*,*)'Southern hemisphere: ',jj
        CLOSE(kdiag)
        STOP
     END IF
  END IF

!---------------------------------------------------
! no initialization

  IF(kdsk.LE.0)THEN
     qqq = 0.0
     WRITE(kdiag,'(A)')'All zero initialization for concentration'

!---------------------------------------------
! initialization from station file by latitude

  ELSEIF(kdsk.EQ.1)THEN

     ALLOCATE (xvar(kk),yvar(kk))
!    northern hemisphere regression
     kk=0
     DO L=1,numb
        IF(clat(L).GT.0)THEN
           kk=kk+1
           xvar(kk)=clat(L)
           yvar(kk)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_nh,yintc_nh)
     DEALLOCATE (xvar,yvar)

     ALLOCATE (xvar(jj),yvar(jj))
!    southern hemisphere regression
     jj=0
     DO L=1,numb
        IF(clat(L).LE.0)THEN
           jj=jj+1
           xvar(jj)=clat(L)
           yvar(jj)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_sh,yintc_sh)
     DEALLOCATE (xvar,yvar)

     DO j=1,ny
        xlat=(j-1)*dlat+clat1-dlat/2.0
        IF(xlat.GT.0.0)THEN
           CINI=MIN(CMAX,MAX(CMIN,yintc_nh+slope_nh*xlat))/CFACT/ROW
        ELSE
           CINI=MIN(CMAX,MAX(CMIN,yintc_sh+slope_sh*xlat))/CFACT/ROW
        END IF 
        DO i=1,nx
           DO k=1,nz
              PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP
              IF(pres.GT.150.0)THEN
                 qqq(i,j,k)=CINI
              ELSEIF(pres.GT.70.0)THEN
                 qqq(i,j,k)=CBOT*CINI
              ELSE
                 qqq(i,j,k)=CTOP*CINI
              END IF
           END DO
        END DO
     END DO
     DEALLOCATE (clat,clon,xcon)
     WRITE(kdiag,'(A)')'Latitudinal initialization'

!----------------------------------------------
! dewpoint interpolation from station file

  ELSEIF(kdsk.EQ.2)THEN

     ALLOCATE (xvar(kk),yvar(kk))
!    northern hemisphere regression
     kk=0
     DO L=1,numb
        IF(clat(L).GT.0)THEN
           kk=kk+1
           j=1+(clat(L)-clat1+dlat/2.0)/dlat
           i=1+(clon(L)-clon1+dlon/2.0)/dlon
           ES=10.0**(9.4051-2353.0/TTT(I,J,1))
           EA=MAX(DBLE(1.0),MMM(I,J,1))*ES/100.0
           xvar(kk)=2353.0/(9.4051-DLOG10(EA))
           yvar(kk)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_nh,yintc_nh)
     DEALLOCATE (xvar,yvar)

     ALLOCATE (xvar(jj),yvar(jj))
!    southern hemisphere regression
     jj=0
     DO L=1,numb
        IF(clat(L).LE.0)THEN
           jj=jj+1
           j=1+(clat(L)-clat1+dlat/2.0)/dlat
           i=1+(clon(L)-clon1+dlon/2.0)/dlon
           ES=10.0**(9.4051-2353.0/TTT(I,J,1))
           EA=MAX(DBLE(1.0),MMM(I,J,1))*ES/100.0
           xvar(jj)=2353.0/(9.4051-DLOG10(EA))
           yvar(jj)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_sh,yintc_sh)
     DEALLOCATE (xvar,yvar)

     DO J=1,ny
     DO I=1,nx

        ES=10.0**(9.4051-2353.0/TTT(I,J,1))
        EA=MAX(DBLE(1.0),MMM(I,J,1))*ES/100.0
        DP=2353.0/(9.4051-DLOG10(EA))

        xlat=(j-1)*dlat+clat1-dlat/2.0
        IF(xlat.GE.5.0)THEN
           CINI=MIN(CMAX,MAX(CMIN,yintc_nh+slope_nh*dp))/CFACT/ROW
        ELSEIF(xlat.LT.5.0.AND.xlat.GT.-5.0)THEN
           CINI_NH=MIN(CMAX,MAX(CMIN,yintc_nh+slope_nh*dp))/CFACT/ROW
           CINI_SH=MIN(CMAX,MAX(CMIN,yintc_sh+slope_sh*dp))/CFACT/ROW
           CINI=0.5*(CINI_NH+CINI_SH)
        ELSE
           CINI=MIN(CMAX,MAX(CMIN,yintc_sh+slope_sh*dp))/CFACT/ROW
        END IF 

        DO k=1,nz
           PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP
           IF(pres.GT.150.0)THEN
              qqq(i,j,k)=CINI
           ELSEIF(pres.GT.70.0)THEN
              qqq(i,j,k)=CBOT*CINI
           ELSE
              qqq(i,j,k)=CTOP*CINI
           END IF
        END DO
		
     END DO
     END DO
     DEALLOCATE (clat,clon,xcon)
     WRITE(kdiag,'(A)')'Dew point initialization'

!----------------------------------------------
! initialization from previous simulation

  ELSEIF(kdsk.EQ.3)THEN

     READ(kc3dd) KYR,KMO,KDA,KHR,KMN
     READ(kc3dd) NXG,NYG,NZG
     READ(kc3dd) CLATG,CLONG,GLAT,GLON
     ALLOCATE (GGG(nzg))
     READ(kc3dd) GGG

     IF(IYR.NE.KYR.OR.IMO.NE.KMO.OR.IDA.NE.KDA)THEN
        WRITE(kdiag,'(A)')    'GBL3DM.BIN header inconsistent with simulation start date'
        WRITE(kdiag,'(A,3I3)')'FILE date: ',KYR,KMO,KDA
        WRITE(kdiag,'(A,3I3)')'PROG date: ',IYR,IMO,IDA
        CLOSE(kdiag)
        STOP 910
     END IF

!    if grids are not the same, set the data blending flag
     IF(nxg.NE.nx.OR.nyg.NE.ny.OR.nzg.NE.nz) THEN
        BLEND=.true.
        WRITE(kdiag,'(A,3I3)')   'Input data grid dimensions: ',nxg,nyg,nzg
        WRITE(kdiag,'(A,4F10.1)')'  Horizontal domain limits: ',clatg,clong,glat,glon
        WRITE(kdiag,'(A,6F10.1)')'  Vertical levels (1 to 6): ',ggg(1:6)
     END IF

     IF(blend)THEN
        ALLOCATE (gqc(nxg,nyg,nzg))
        READ(kc3dd) GQC

!       loop through the current grid and find position on input grid
        DO K=1,nz

!          find the corresponding height on the input grid (ggg)
           kk=nzg
           DO WHILE (kk.GT.1.AND.0.5*(ggg(kk)+ggg(kk-1)).LT.sig(k))
              kk=kk-1
           END DO
           kk=MAX(1,kk)
       
        DO J=1,ny
        DO I=1,nx

!          location on the current grid
           xlat=(j-1)*dlat+clat1-dlat/2.0
           xlon=(i-1)*dlon+clon1-dlon/2.0

!          index on the old grid
           jj=1+(xlat-clatg+glat/2.0)/glat
           ii=1+(xlon-clong+glon/2.0)/glon
           jj=MAX(1,MIN(jj,nyg))
           ii=MAX(1,MIN(ii,nxg))

!          fill in concentration value
           QQQ(i,j,k)=GQC(ii,jj,kk)/CFACT/ROW

        END DO
        END DO
        END DO
        WRITE(kdiag,'(A)')'Blended initialization from file: GBL3DM.BIN'
        DEALLOCATE (gqc)

     ELSE
!       input grid is the same as the current simulation grid
        READ(kc3dd) CCC
        QQQ = CCC/CFACT/ROW
        WRITE(kdiag,'(A)')'Direct initialization from file: GBL3DM.BIN'
     END IF

     REWIND(kc3dd)
     DEALLOCATE(ggg)

!---------------------------------------------
! single value initialization (kdsk>3 in tenths)
! e.g. an input of 15 would be evaluated as 1.5

  ELSE
     CINI = DBLE(REAL(kdsk)/10.0)/CFACT/ROW/1000.0
     DO j=1,ny	 
        DO i=1,nx
           DO k=1,nz
              PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP
              IF(pres.GT.150.0)THEN
                 qqq(i,j,k)=CINI
              ELSEIF(pres.GT.70.0)THEN
                 qqq(i,j,k)=CBOT*CINI
              ELSE
                 qqq(i,j,k)=CTOP*CINI
              END IF
           END DO
        END DO
     END DO
     WRITE(kdiag,'(A,I3)')'Constant value initialization: ',kdsk

  END IF

END SUBROUTINE conset


!###############################################################################
! CHKMET - Analysis of the meteorological data file: Opens file, determines 
! grid dimensions, and checks for appropriate data contents.
!-------------------------------------------------------------------------------
! LAST REVISED: 11 AUG 2005 - initial version from program chk_rec
!               24 Feb 2006 - standard atm testing
!               28 Mar 2008 - number of records per time period
!-------------------------------------------------------------------------------
 
SUBROUTINE chkmet (fname,khrs,kmet)

  USE constant
  USE metarray

  IMPLICIT NONE

  CHARACTER(80), INTENT(IN)  :: FNAME
  INTEGER*4,     INTENT(OUT) :: KHRS
  INTEGER*4,     INTENT(OUT) :: KMET

  CHARACTER(50)   :: LABEL
  CHARACTER(4)    :: KVAR, MODEL
  CHARACTER(3072) :: HEADER

  LOGICAL    :: FTEST
  INTEGER*4  :: K,L,KOL,NREC,KRET
  INTEGER*4  :: IYR,IMO,IDA,IHR,IFH
  INTEGER*4  :: IFHX,MINUTES,KFLAG,LENH

  REAL*8     :: CLAT2, CLON2
  REAL*8     :: SIZE, ORIENT, CONE, SYNCXP, SYNCYP, DUMMY
 
  INTEGER*4 :: nx,ny,nz
  REAL*4    :: clat1,clon1,dlat,dlon

  COMMON /MAINGRID/ clat1,clon1,dlat,dlon
  COMMON /SIZEGRID/ nx,ny,nz

! test for meteo file existence
  INQUIRE(FILE=FNAME,EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     WRITE(KDIAG,'(A)')'ERROR: Unable to find meteorological data file'
     WRITE(KDIAG,'(A)') FNAME
     CLOSE(kdiag)
     STOP 915
  END IF

! open file to decode the standard label (50) plus the fixed portion (108) of the header
  OPEN(KUMET,FILE=FNAME,RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED')

! decode the standard portion of the index record
  READ(KUMET,REC=1)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,4X,A4)',IOSTAT=kret)IYR,IMO,IDA,IHR,IFH,KVAR
  WRITE(KDIAG,'(A,4I3)')'Opened file: ',IYR,IMO,IDA,IHR

  IF(KVAR.NE.'INDX'.OR.KRET.NE.0)THEN
     WRITE(KDIAG,'(A)')'ERROR: Meteorological data file has an unsupported format'
     CLOSE(kdiag)
     STOP 920
  END IF

! decode extended portion of the header
  READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)',IOSTAT=kret)            &
     MODEL,IFHX,MINUTES,CLAT2,CLON2,DLAT,DLON,SIZE,ORIENT,CONE,            &
     SYNCXP,SYNCYP,CLAT1,CLON1,DUMMY,NX,NY,NZ,KFLAG,LENH
  NZ=NZ-1  ! first level are surface fields not used in 3D arrays

  IF(KRET.NE.0)THEN
     WRITE(KDIAG,'(A)')'ERROR: Decoding meteorological data header record'
     WRITE(KDIAG,'(A)') HEADER
     CLOSE(kdiag)
     STOP 925
  END IF

! close file 
  CLOSE (KUMET)

! final checks to determine if the data are of the proper type
  IF(SIZE.NE.0.0.OR.ORIENT.NE.0.0.OR.CONE.NE.0.0)THEN
     WRITE(KDIAG,'(A)')'ERROR: meteorological data grid is not lat-lon'
     CLOSE(kdiag)
     STOP 930

  ELSEIF(KFLAG.NE.2)THEN
     WRITE(KDIAG,'(A,I3)')'ERROR: meteorological coordinate not pressure - ',kflag
     CLOSE(kdiag)
     STOP 935

  ELSEIF(DLON.NE.1.0.AND.DLON.NE.2.5)THEN
     WRITE(KDIAG,'(A)')'ERROR: only 1.0 deg or 2.5 deg lat-lon grids are supported!'
     WRITE(KDIAG,'(A,2F10.2)')'Grid delta  (lat,lon): ',DLAT,DLON
     CLOSE(kdiag)
     STOP 938

  ELSEIF((CLON2+DLON-CLON1.EQ.360.0).OR.(CLON2+DLON-CLON1.EQ.0.0).AND.  &
          CLAT2-CLAT1.EQ.180.0)THEN
!    the meteorology grid is a global lat-lon grid  
     WRITE(KDIAG,'(A,2F10.2)')'Lower left  (lat,lon): ',CLAT1,CLON1
     WRITE(KDIAG,'(A,2F10.2)')'Upper right (lat,lon): ',CLAT2,CLON2
     WRITE(KDIAG,'(A,2F10.2)')'Grid delta  (lat,lon): ',DLAT,DLON

  ELSE
     WRITE(KDIAG,'(A)')'ERROR: not a global lat-lon meteorological grid'
     CLOSE(kdiag)
     STOP 940
  END IF

! allocate all the meteorological data array space
  CALL setmet
      
! open meteorological file with correct length
  OPEN(KUMET,FILE=FNAME,RECL=(50+NX*NY),ACCESS='DIRECT',FORM='UNFORMATTED')

! read the entire extended header record now that the record length is known
  READ(KUMET,REC=1)LABEL,HEADER(1:LENH)

! loop through and decode the remainder of the index string
  KOL=109
  NREC=1   ! number of records per time period

! note that level=0 represents surface fields not used in the calculation
  DO L=0,NZ
     READ(HEADER(KOL:KOL+7),'(F6.2,I2)') PPP(L), NVAR(L)
     KOL=KOL+8
     DO K=1,NVAR(L)
        READ(HEADER(KOL:KOL+7),'(A4,I3)') KVAR, KRET
        KOL=KOL+8
        NREC=NREC+1
     END DO
  END DO
  KMET=NREC
  WRITE(KDIAG,'(A,I3)')'Numb records per time: ',KMET

  DO K=1,nz
!    define sigma coordinate
     SIG(K)=(PPP(K)-PTOP)/(PSFC-PTOP)

!    set the standard atm index value for each level
     KNDX(K)=0
     stdlvl : DO L=1,38
        IF(PPP(K).EQ.STDPPP(L))THEN
           KNDX(K)=L
           EXIT stdlvl
        END IF
     END DO stdlvl
     IF(KNDX(K).EQ.0)THEN
        WRITE(KDIAG,'(A,F10.1)')'Input level not in std atmosphere: ',PPP(K)
        CLOSE(kdiag)
        STOP 945
     END IF
  END DO

! determine the temporal interval
  NREC=NREC+1
  READ(KUMET,REC=NREC)LABEL
  READ(LABEL,'(6X,I2)')IFH
  KHRS=IFH-IHR
  WRITE(KDIAG,'(A,I3)')'Grid delta time (hrs): ',KHRS

END SUBROUTINE chkmet


!###############################################################################
! SETMET - Allocate the meteorological data arrays
!-------------------------------------------------------------------------------
! UUU - U WIND COMPONENTS       m/s   
! VVV - V WIND COMPONENTS       m/s
! TTT - TEMPERATURE             Kelvin
! MMM - MOISTURE (RH)           percent
! HHH - HEIGHT MSL OF FIELD     meters
! WWW - DIVERGENCE & VELOCITY   1/s to m/s
! KKK - VERTICAL MIXING         m2/s
! RRR - LOCAL AIR DENSITY       Kg/m3
! QQQ - MIXING RATIO            g/Kg
! GSX - WE GRID DIMENSIONS      meters
! GSY - SN GRID DIMENSIONS      meters
! SFC - SURFACE PRESSURE        mb
!-------------------------------------------------------------------------------
! LAST REVISED: 11 AUG 2005 - initial version
!               06 SEP 2005 - added surface pressure
!-------------------------------------------------------------------------------
 
SUBROUTINE setmet

  USE constant
  USE metarray

  IMPLICIT NONE

  INTEGER*4 :: nx,ny,nz,kret,ktot
  COMMON /SIZEGRID/ nx,ny,nz

  ktot=0

  ALLOCATE (NVAR(0:nz), NGP(ny), KNDX(nz), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #1 - ',kret
  ktot=ktot+kret

  ALLOCATE (PPP(0:nz), SIG(nz), VAL(nz,7), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #2 - ',kret
  ktot=ktot+kret

  ALLOCATE (POT(nz), VB4(nz), VB8(nz), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #3 - ',kret
  ktot=ktot+kret

  ALLOCATE (UUU(nx,ny,nz),VVV(nx,ny,nz),WWW(nx,ny,nz), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #4 - ',kret
  ktot=ktot+kret

  ALLOCATE (TTT(nx,ny,nz),MMM(nx,ny,nz),HHH(nx,ny,nz), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #5 - ',kret
  ktot=ktot+kret

  ALLOCATE (KKK(nx,ny,nz),RRR(nx,ny,nz), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #6 - ',kret
  ktot=ktot+kret

  ALLOCATE (GSX(nx,ny),GSY(nx,ny),GXY(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #7 - ',kret
  ktot=ktot+kret

  ALLOCATE (SFC(nx,ny), AVG(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #8 - ',kret
  ktot=ktot+kret

  ALLOCATE (QQQ(nx,ny,nz), CCC(nx,ny,nz), CONC(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(kdiag,'(A,I3)')'ERROR: memory allocation #9 - ',kret
  ktot=ktot+kret

  IF(ktot.gt.0) THEN
     CLOSE(kdiag)
     STOP 945
  END IF

! initialize diagnostic variable arrays (may not always be present at all levels)
  www = 0.0
  mmm = 0.0

END SUBROUTINE setmet


!###############################################################################
! SETGRD - Creates the latitude/longitude dependent grid dimensions for
! computing the finite difference fluxes and converting cell mass to air 
! concentration. In the current global lat-lon grid configuration the 1,1 grid
! point is at the southwest corner and the maximum points are at the northwest
! corner.  Lat-lon grids can start at either 000 or 180E longitude.
!-------------------------------------------------------------------------------
! GRID COORDINATE SYSTEM
! 
!     *      *      *        where * are the meteorlogical grid points
!                            of dimensions nx,ny     
!         +-----+
!         |     |            and the box represents the concentration
!     *   |  *  |   *        grid cell whose dimensions are computed 
!         |     |            between the + verticies
!         +-----+
!
!     *      *      *
! ------------------------------------------------------------------------------
! LAST REVISED: 16 Oct 1989 - initial version
!               12 Aug 2005 - converted to lat/lon grid system
!               26 Sep 2008 - grid size computation needed clat1 added
!-------------------------------------------------------------------------------

SUBROUTINE setgrd

  USE constant
  USE metarray

  IMPLICIT NONE

! establish the default computational horzontal grid spacing for latitudes 70-90
! (near the poles) in terms of the number of meteorological grid points comprise
! each concentration grid cell for both the 1-deg and 2.5 deg grids

  INTEGER*4 :: nh10(181),nh25(73)
  DATA nh10 /360,2*180,4*90,2*45,5*15,6*3,141*1,6*3,5*15,2*45,4*90,2*180,360/
  DATA nh25 /144,72,36,18,6,6,2,2,57*1,2,2,6,6,18,36,72,144/                          
                          
  INTEGER*4 :: j
  REAL*8    :: clat,base,pole,abot,atop

  INTEGER*4 :: nx,ny,nz
  REAL*4    :: clat1,clon1,dlat,dlon

  COMMON /MAINGRID/ clat1,clon1,dlat,dlon
  COMMON /SIZEGRID/ nx,ny,nz

! establish the horizontal grid aggregation factors to limit the number
! of computational points near the poles
  
  IF(ny.EQ.181)THEN
     ngp = nh10
  ELSEIF(ny.EQ.73)THEN
     ngp = nh25
  ELSE
     ngp = 1
     WRITE(kdiag,'(A)')'ERROR: grid spacing aggregation not defined - ',ny
     CLOSE(kdiag)
     STOP 948
  END IF

! Distance on the earth's surface is given by the earth's circumference
! converted to distance per degree times the cell spacing. The grid spacing
! is used in the finite difference equations and therefore should apply
! at the meteorological grid points.

  base = 2.0*pi*rearth/360.0 ! base distance on earth in meters/deg-latitude

  gsy  = dlat*base           ! spacing (meters) same all latitudes

! Distance in the horizontal uses a similar approach by computing the 
! circumference at each latitude. Horizontal distance = 0 at the poles.

  gsx = 0.0

  DO j=2,(ny-1)   
     clat = (j-1)*dlat+clat1
     gsx(:,j) = dlon*base*cos(clat*radpdeg)  ! the same for all longitudes
  END DO

! The area of the cell surrounding each meteorological grid point is computed
! from the difference in the sector areas as measured from the pole to the 
! the base and top of each cell. The sector area is defined as (pi P**2),
! where P is the linear distance from the pole to the latitude. Northern
! hemisphere distances are the same as southern hemisphere distances.

  DO j=2,ny/2
     clat = (j-1)*dlat-(0.5*dlat)+clat1   ! bot of box lat
     base=rearth*cos(clat*radpdeg)
     pole=rearth*(1.0-sin(clat*radpdeg))
     abot=pi*(base*base+pole*pole)          

     clat = j*dlat-(0.5*dlat)+clat1       ! top of box lat
     base=rearth*cos(clat*radpdeg)
     pole=rearth*(1.0-sin(clat*radpdeg))
     atop=pi*(base*base+pole*pole)          

     gxy(:,j)      = abs(atop-abot)/nx    ! area (sq-meters) per grid cell
     gxy(:,ny-j+1) = abs(atop-abot)/nx
  END DO

! special case at the equator

  gxy(:,1+ny/2) = gsx(:,1+ny/2) * gsy(:,1+ny/2) 

! The area of the polar grid cell (there is only one, but all are the same in
! the array) is computed from the radius using the above equations.

  clat=90.0-0.5*dlat
  base=rearth*cos(clat*radpdeg)
  pole=rearth*(1.0-sin(clat*radpdeg))
  gxy(:,1)  = pi*(base*base+pole*pole)          
  gxy(:,ny) = pi*(base*base+pole*pole)          

END SUBROUTINE setgrd


!###############################################################################
! POSMET - Position the meteorology file to the index record number for the 
! designated time.  Output is the record number of the index record which is
! always the first record of a time period group of records.
!-------------------------------------------------------------------------------
! LAST REVISED: 19 Feb 1992 - initial version
!               12 Aug 2005 - fortran90 update
!-------------------------------------------------------------------------------

SUBROUTINE posmet (iyr,imo,ida,ihr,krec)

   USE constant

   IMPLICIT NONE

   INTEGER*4,  INTENT(IN)  :: iyr,imo,ida,ihr
   INTEGER*4,  INTENT(OUT) :: krec

   CHARACTER(50) :: LABEL
   INTEGER*4     :: KYR,KMO,KDA,KHR,KRET

   KREC=1
   KRET=0

   DO WHILE (kret.EQ.0)
      READ(KUMET,REC=KREC,IOSTAT=kret)LABEL

      IF(kret.NE.0)THEN
         WRITE(KDIAG,'(A)')'ERROR: positioning data file to start hour'
         WRITE(KDIAG,'(A,3I3)')'Requested start time: ',iyr,imo,ida,ihr
         WRITE(KDIAG,'(A,3I3)')'Last valid file time: ',kyr,kmo,kda,khr
         WRITE(KDIAG,'(A, I3)')'Current record  numb: ',krec
         CLOSE(kdiag)
         STOP 950
      END IF

      READ(LABEL,'(4I2)')KYR,KMO,KDA,KHR
      IF(IYR.EQ.KYR.AND.IMO.EQ.KMO.AND.IDA.EQ.KDA.AND.IHR.EQ.KHR)RETURN
      KREC=KREC+1
   END DO

END SUBROUTINE posmet


!###############################################################################
! METINP - Meteorology Input for one time period 
!-------------------------------------------------------------------------------
! LAST REVISED: 19 Feb 1992 - initial version
!               04 May 2004 - data file format upgrade
!               15 Aug 2005 - global lat/lon grids
!               06 Sep 2005 - added surface pressure
!               16 Feb 2006 - auto file name construct or input line
!               27 Nov 2007 - zero out w and rh before input
!-------------------------------------------------------------------------------

SUBROUTINE metinp (fname,myr,mmo,mda,mhr,krec)

  USE constant
  USE metarray

  IMPLICIT NONE

  CHARACTER(80), INTENT(INOUT) :: FNAME  
  INTEGER*4,     INTENT(OUT)   :: myr,mmo,mda,mhr
  INTEGER*4,     INTENT(INOUT) :: krec

  LOGICAL       :: AUTO = .FALSE. ! auto file name construct

  LOGICAL       :: FTEST
  CHARACTER(80) :: FNAME2
  CHARACTER(50) :: LABEL
  CHARACTER(4)  :: VARB
  CHARACTER(3)  :: MONTH
  CHARACTER(3)  :: MON(12) = (/ 'jan','feb','mar','apr','may','jun',  &
                                'jul','aug','sep','oct','nov','dec'  /)

  INTEGER*4     :: K,L,NX,NY,NZ,NXY,MFH,LEV,NEXP,KRET,KEND,YEAR,WEEK
  REAL*8        :: VAR1,PREC

  CHARACTER(1), ALLOCATABLE :: BUF(:)

  SAVE buf, auto, mon

  INTERFACE
    SUBROUTINE unpack (RVAR,CVAR,NX,NY,NXY,NEXP,VAR1)
    REAL*8,       INTENT(OUT) :: RVAR(:,:)
    CHARACTER(1), INTENT(IN)  :: CVAR(:)
    INTEGER*4,    INTENT(IN)  :: NX,NY,NXY,NEXP
    REAL*8,       INTENT(IN)  :: VAR1
    END SUBROUTINE unpack
  END INTERFACE

  COMMON /SIZEGRID/ nx,ny,nz

  nxy=nx*ny
  IF(.NOT.ALLOCATED(buf)) ALLOCATE (BUF(nxy))

!------------------------------------------------------------
! read the index record for testing (krec always points to index)

  READ(KUMET,REC=KREC,IOSTAT=kend)LABEL

  IF(kend.NE.0)THEN
!    closing old meteo input file
     CLOSE(kumet)

!    meteorology directory and file from input
!    multiple input lines turn off file name auto construct 
     IF(.NOT.auto)THEN
        READ(ksinp,*,IOSTAT=kret) FNAME2
        IF(kret.NE.0)THEN
           AUTO=.TRUE.
        ELSE
!          name of next meteo file to process
           FNAME=FNAME2
        END IF
     END IF

!    auto file name construct section
     IF(auto)THEN
        IF(ny.EQ.181)THEN
!          one-degree resolution files
           kend=INDEX(fname,' ')-1
           READ(fname(kend-7:),'(a3,i2,2x,i1)')month,year,week
           DO K=1,12
              IF(month.EQ.mon(k))mmo=k
           END DO
           week=week+1
           IF(week.GT.5.OR.(week.EQ.5.AND.mmo.EQ.2.AND.MOD(year,4).NE.0))THEN
              week=1
              mmo=mmo+1
              if(mmo.GT.12)THEN
                 year=MOD(year+1,100)
                 mmo=1
              END IF
           END IF
           WRITE(fname(kend-7:),'(a3,i2.2,a2,i1)')mon(mmo),year,'.w',week
   
        ELSE
!          two and half degree resolution files
           kend=INDEX(fname,' ')-1
           READ(fname(kend-9:),'(i4,i2)')year,mmo
           mmo=mmo+1
           if(mmo.GT.12)THEN
              year=year+1
              mmo=1
           END IF
           WRITE(fname(kend-9:),'(I4,I2.2,A)')year,mmo,'.gbl'
        END IF
     END IF

!    check for data file
     INQUIRE(FILE=fname,EXIST=ftest)
     IF(.NOT.ftest)THEN    
        WRITE(kdiag,'(2A)')'Seeking file: ',FNAME(1:70)
        WRITE(kdiag,'(A)')'No more meteorological files to open!'
        CLOSE(kdiag)
        STOP 955
     END IF

!    open new meteo data file
     OPEN(kumet,FILE=FNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=(50+nxy))
     krec=1
     READ(KUMET,REC=KREC)LABEL
     WRITE(kdiag,'(2A)')'Opened file: ',FNAME(1:70)
  END IF

! insure that the first record is the index record
  READ(LABEL,'(14X,A4)')VARB
  IF(VARB.NE.'INDX')THEN
     WRITE(KDIAG,'(A,I5)')'ERROR: index record not found at rec number - ',KREC
     WRITE(KDIAG,'(A)') LABEL
     CLOSE(kdiag)
     STOP 958
  END IF

  krec=1+krec ! skip past the index record

! zero out variables that may not be complete in the vertical
  WWW=0.0
  MMM=0.0
    
  DO L=0,nz
  DO K=1,nvar(L)

!    read one data record (variable at a level)
     READ(KUMET,REC=KREC)LABEL,BUF

!    decode index record for packing information
     READ(LABEL,'(6I2,2X,A4,I4,2E14.7)',IOSTAT=kret)       &
          MYR,MMO,MDA,MHR,MFH,LEV,VARB,NEXP,PREC,VAR1

     IF(kret.NE.0)THEN
!       error decoding label field
        WRITE(KDIAG,'(A)')'ERROR: decoding meteo label field'
        WRITE(KDIAG,'(A)') LABEL
        CLOSE(kdiag)
        STOP 960
     END IF

     IF(mfh.LT.0.OR.varb.EQ.'NULL')THEN
!       missing data not permitted
        WRITE(KDIAG,'(A)')'ERROR: missing meteorological data'
        WRITE(KDIAG,'(A)') LABEL
        CLOSE(kdiag)
        STOP 965
     END IF

     IF(lev.NE.L)THEN
!       level mismatch
        WRITE(KDIAG,'(A)')'ERROR: input data level mismatch'
        WRITE(KDIAG,'(3(A,I3))')'Loop: ',L,'  File: ',lev,'  Record: ',krec 
        CLOSE(kdiag)
        STOP 968
     END IF

     IF(L.GT.0)THEN
!       upper level variables
        IF(VARB.EQ.'HGTS') CALL UNPACK(HHH(:,:,L),BUF,NX,NY,NXY,NEXP,VAR1)
        IF(VARB.EQ.'TEMP') CALL UNPACK(TTT(:,:,L),BUF,NX,NY,NXY,NEXP,VAR1)
        IF(VARB.EQ.'UWND') CALL UNPACK(UUU(:,:,L),BUF,NX,NY,NXY,NEXP,VAR1)
        IF(VARB.EQ.'VWND') CALL UNPACK(VVV(:,:,L),BUF,NX,NY,NXY,NEXP,VAR1)
        IF(VARB.EQ.'WWND') CALL UNPACK(WWW(:,:,L),BUF,NX,NY,NXY,NEXP,VAR1)
        IF(VARB.EQ.'RELH') CALL UNPACK(MMM(:,:,L),BUF,NX,NY,NXY,NEXP,VAR1)

     ELSE
!       surface variables
        IF(VARB.EQ.'PRSS') CALL UNPACK(SFC(:,:),BUF,NX,NY,NXY,NEXP,VAR1)
     END IF
     krec=krec+1

  END DO
  END DO
! WRITE(KDIAG,'(A,4I2)')' Loaded meteorology data: ',MYR,MMO,MDA,MHR

END SUBROUTINE metinp


!###############################################################################
! UNPACK - Meteorology data record unpacker convert 1-byte packed to REAL*8
!-------------------------------------------------------------------------------
! LAST REVISED: 16 Oct 1989 - initial version
!               15 Aug 2005 - fortran90 upgrade
!-------------------------------------------------------------------------------

SUBROUTINE unpack (RVAR,CVAR,NX,NY,NXY,NEXP,VAR1)

  IMPLICIT NONE

  REAL*8,       INTENT(OUT) :: RVAR(:,:)       ! unpacked data array
  CHARACTER(1), INTENT(IN)  :: CVAR(:)         ! packed data 
  INTEGER*4,    INTENT(IN)  :: NX,NY,NXY,NEXP  ! grid dimensions and exponent
  REAL*8,       INTENT(IN)  :: VAR1            ! data value at 1,1

  REAL*8    :: sfact,qold
  INTEGER*4 :: i,j,indx,itemp

  CHARACTER(1)  :: mychr
  INTEGER       :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)

  SFACT=2.0**(7-NEXP)
  QOLD=DBLE(VAR1)
  INDX=0
  DO J=1,NY
     DO I=1,NX
        INDX=INDX+1
        ITEMP=JCHAR(CVAR(INDX))
        RVAR(I,J)=((ITEMP-127.)/SFACT)+QOLD
        QOLD=RVAR(I,J)
     END DO
     QOLD=RVAR(1,J)
   END DO

END SUBROUTINE unpack


!###############################################################################
! DSKSET - Initialize packed concentration output file in standard HYSPLIT
! binary format. Horizontal grid is the same as the meteorological grid with the
! exception that the cells are centered over the meteo grid points and the pole 
! grid cell is circular about the pole. Also open the daily 3D binary dump file 
! which can be used for initialization.
!-------------------------------------------------------------------------------
! LAST REVISED: 20 Jun 2003 - initial version
!               12 Aug 2005 - fortran90 update
!-------------------------------------------------------------------------------

SUBROUTINE dskset (tracer,iyr,imo,ida,ihr,level)

  USE constant

  IMPLICIT NONE

  CHARACTER(4), INTENT(IN) :: tracer
  INTEGER*4,    INTENT(IN) :: iyr,imo,ida,ihr,level

  CHARACTER(11) :: fname
  INTEGER*4     :: nx,ny,nz, year
  REAL*4        :: clat1,clon1,dlat,dlon

  COMMON /MAINGRID/ clat1,clon1,dlat,dlon
  COMMON /SIZEGRID/ nx,ny,nz

! compute the four digit year
  year=1900+iyr
  if(iyr.LT.40)year=2000+iyr

! construct file name
  WRITE(fname,'(a3,i4,a4)')'sfc',year,'.bin'

! two-dimensional concentration output file
  OPEN(kcout,FILE=fname,FORM='UNFORMATTED')

! Binary header record
! Global Tracer Model,    Start date, Forecast, Sources, Packing
  WRITE(kcout)'GTM3',IYR,IMO,IDA,IHR,        0,       1,       0

! calculation start by location ....... lat   long  level
  WRITE(kcout)IYR,IMO,IDA,IHR,          0.0,  0.0,  REAL(level)

! horizontal grid 
  WRITE(kcout)NY,NX,DLAT,DLON,CLAT1,CLON1

! vertical grid index record (only one level is output)
! ............ number  level
  WRITE(kcout)   1,    LEVEL

! pollutant identification record (only one pollutant is output)
! ............ number   type
  WRITE(kcout)   1,    TRACER

END SUBROUTINE dskset


!###############################################################################
! CALEND - Simple calendar routine
!-------------------------------------------------------------------------------
! LAST REVISED: 15 Aug 2005 - initial version
!               28 Mar 2008 - added backward calendar
!-------------------------------------------------------------------------------

SUBROUTINE calend (IYR,IMO,IDA,IHR,IMN)

  IMPLICIT NONE

  INTEGER*4, INTENT(INOUT) :: IYR,IMO,IDA,IHR,IMN

! number of days per month
  INTEGER :: NDM(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31/)
  SAVE NDM

  NDM(2)=28
  IF(MOD(IYR,4).EQ.0)NDM(2)=29

! forward

  IF(IMN.GE.60)THEN
!    new hour
     IMN=0
     IHR=IHR+1

  IF(IHR.GE.24)THEN
!    new day
     IHR=0
     IDA=IDA+1

  IF(IDA.GT.NDM(IMO)) THEN
!    new month
     IMO=IMO+1
     IDA=1

  IF(IMO.GT.12) THEN
!    new year
     IYR=MOD(IYR+1,100)
     IMO=1

  END IF
  END IF
  END IF
  END IF

! backward

  IF(IMN.LT.0)THEN
!    new hour
     IMN=60+IMN
     IHR=IHR-1

  IF(IHR.LT.0)THEN
!    new day
     IHR=24+IHR
     IDA=IDA-1

  IF(IDA.LT.1) THEN
!    new month
     IMO=IMO-1
     IDA=NDM(IMO)

  IF(IMO.LT.1) THEN
!    new year
     IYR=MOD(IYR-1,100)
     IMO=12

  END IF
  END IF
  END IF
  END IF

END SUBROUTINE calend


!###############################################################################
! DSKOUT - Write concentrations to disk in standard HYSPLIT format
!-------------------------------------------------------------------------------
! LAST REVISED: 23 Jun 2003 - initial version
!               15 Aug 2005 - fortran90 upgrade
!-------------------------------------------------------------------------------

SUBROUTINE dskout (TRACER,LEVEL,KOUT,IY,IM,ID,IH,CFACT)

  USE constant
  USE metarray

  IMPLICIT NONE

  CHARACTER(4), INTENT(IN) :: tracer
  INTEGER*4,    INTENT(IN) :: level,kout,iy,im,id,ih
  REAL*8,       INTENT(IN) :: cfact

  INTEGER*4  :: i,j,k,nx,ny,nz

  COMMON /SIZEGRID/ nx,ny,nz

! snapshot concentrations valid at output time
  WRITE(kcout)IY,IM,ID,IH,0,0
  WRITE(kcout)IY,IM,ID,IH,0,0

! convert to output units at STP
  DO J=1,NY
  DO I=1,NX
     CONC(I,J) = QQQ(i,j,kout)*CFACT*ROW
  END DO
  END DO

! write all grid points according to the old format ...
! packing not efficient when entire grid has non-zero values
  WRITE(kcout)tracer,level,CONC

END SUBROUTINE dskout


!###############################################################################
! DAYOUT - Write concentrations to disk in binary format for initialization
!-------------------------------------------------------------------------------
! LAST REVISED: 24 Aug 2005 - initial version
!               14 May 2008 - save annual file
!-------------------------------------------------------------------------------

SUBROUTINE dayout (IYR,IMO,IDA,IHR,IMN,CFACT)

  USE constant
  USE metarray

  IMPLICIT NONE

  INTEGER*4, INTENT(IN) :: iyr,imo,ida,ihr,imn
  REAL*8,    INTENT(IN) :: cfact

  CHARACTER(11) :: fname
  REAL*4        :: clat1,clon1,dlat,dlon
  INTEGER*4     :: i,j,k,nx,ny,nz,year

  COMMON /SIZEGRID/ nx,ny,nz
  COMMON /MAINGRID/ clat1,clon1,dlat,dlon

! special case for annual initialization
  IF(imo.EQ.1.AND.ida.EQ.1)THEN
     CLOSE(kc3dd)
     year=1900+iyr
     if(iyr.LT.40)year=2000+iyr
     WRITE(fname,'(a3,i4,a4)')'gbl',year,'.bin'
     OPEN(kc3dd,FILE=fname,FORM='UNFORMATTED')
  END IF

! dump all array contents to file
  WRITE(kc3dd) iyr,imo,ida,ihr,imn
  WRITE(kc3dd) nx,ny,nz
  WRITE(kc3dd) DBLE(clat1),DBLE(clon1),DBLE(dlat),DBLE(dlon)
  WRITE(kc3dd) sig
  CCC = QQQ*CFACT*ROW
  WRITE(kc3dd) CCC

  IF(imo.EQ.1.AND.ida.EQ.1)THEN
!    close special annual initialization file and reopen daily
     CLOSE(kc3dd)
     OPEN(kc3dd,FILE='gbl3dim.bin',FORM='UNFORMATTED')
  ELSE
!    over-write daily file
     REWIND(kc3dd)
  END IF

END SUBROUTINE dayout


!###############################################################################
! DIVERG - Convert vertical velocity from mb/s to m/s
!-------------------------------------------------------------------------------
! LAST REVISED: 19 Feb 1992 - initial version
!               15 Aug 2005 - fortran90 upgrade
!               15 Sep 2005 - added divegence or data option
!-------------------------------------------------------------------------------

SUBROUTINE diverg (WFACT,KDIVG)

  USE constant
  USE metarray

  IMPLICIT NONE

  REAL*8,    INTENT(IN) :: wfact   ! fraction of vertical velocity in each cell
  INTEGER*4, INTENT(IN) :: kdivg   ! when 0 use data, when 1 compute from divergence

  REAL*8    :: dzi,delu,delv,divg,wvel
  INTEGER*4 :: i,j,k,nx,ny,nz,im1,ip1

  COMMON /SIZEGRID/ nx,ny,nz

  DO I=1,nx

     IM1=I-1                        ! index minus one
     IP1=I+1                        ! index plus one
     IF(I.EQ.1)IM1=nx               ! cyclic boundary adjustments
     IF(I.EQ.nx)IP1=1

     DO J=2,ny-1

        DIVG=0.0
        DO K=1,nz

!          set vertical box sizes between interfaces (DZI)
           IF(K.EQ.1) THEN
              DZI=HHH(I,J,K+1)-HHH(I,J,K)
           ELSEIF (K.EQ.NZ) THEN
              DZI=HHH(I,J,K)-HHH(I,J,K-1)
           ELSE
              DZI=0.5*(HHH(I,J,K+1)-HHH(I,J,K-1))
           END IF

           IF(kdivg.EQ.1)THEN
!             divergence in sec-1
              DELU=0.5*(UUU(IP1,J,K)-UUU(IM1,J,K))/GSX(I,J)
              DELV=0.5*(VVV(I,J+1,K)-VVV(I,J-1,K))/GSY(I,J)
              DIVG=DIVG+(DELU+DELV)*DZI

!             positive divergence implies downward motion or a negative vertical velocity
              WVEL=-DIVG*WFACT

!             vertical velocity limits
              IF(wvel.LT.-vwlim) wvel=-vwlim
              IF(wvel.GT.+vwlim) wvel=+vwlim

!             only apply divergence to those levels that originally had a defined velocity
              IF(www(i,j,k).NE.0.0) www(i,j,k)=wvel

           ELSE
!             convert to m/s decreasing pressure with height is a positive vertical velocity
              WVEL=-WFACT*P2JM*WWW(I,J,K)/RRR(I,J,K)/GRAV
              IF(wvel.LT.-vwlim) wvel=-vwlim
              IF(wvel.GT.+vwlim) wvel=+vwlim
              WWW(I,J,K)=wvel
           END IF

        END DO
     END DO

!    pole values the same as adjacent "j" grid point
     WWW(I,1,:)=WWW(I,2,:)
     WWW(I,NY,:)=WWW(I,NY-1,:)

  END DO

END SUBROUTINE diverg


!###############################################################################
! STBLTY - Compute the maximum permitted time step to maintain computational
! stability according to the maximum wind speed. The routine also computes the 
! vertical mixing coefficient array. Note that the vertical coordinate system is
! defined such that the meteorological data represent values at the box centers.
! The heights of the bottom and top cells are the same as their adjacent cells.
!-------------------------------------------------------------------------------
! LAST REVISED: 12 Apr 1991 - initial version
!               16 Aug 2005 - fortran90 upgrade
!               07 Sep 2005 - corrected Louis and Media equations
!               20 Sep 2006 - enhanced limits tests
!-------------------------------------------------------------------------------

SUBROUTINE stblty (VMIX,DELTA)

  USE constant
  USE metarray

  IMPLICIT NONE

  REAL*8,    INTENT(IN)  :: VMIX   ! vertical mixing method
  INTEGER*4, INTENT(OUT) :: DELTA  ! time step in minutes

  REAL*8  :: HSTEP, VSTEP   ! temporary internal time step
  REAL*8  :: FRACZ = 0.50   ! stability limit for vertical mixing
  REAL*8  :: FRACH = 0.75   ! stability limit for horizontal advection  
  REAL*8  :: UMAX  = 150.0  ! maximum permitted horizontal wind speed  

  INTEGER*4 :: i,j,k,nx,ny,nz
  REAL*8    :: ri,fri,tbar,delu,delt,eddy,edsq,pres
  REAL*8    :: htzi,con1,con2,con3,cons,dzc,zmix

  COMMON /SIZEGRID/ nx,ny,nz

  HSTEP = 10800.0  !  default intitial time step to 3h (in sec)
  VSTEP = 10800.0

  DO J=1,ny
  DO I=1,nx
  DO K=1,nz

!    delta box sizes between centers
     IF(K.EQ.1)THEN
        DZC=MAX(DZMIN,HHH(I,J,K+1)-HHH(I,J,K))
     ELSE
        DZC=MAX(DZMIN,HHH(I,J,K)-HHH(I,J,K-1))
     END IF

!    potential temperature and scalar wind speed profile
     PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP        
     POT(K)=TTT(I,J,K)*(1000.0/PRES)**0.286
     VB4(K)=SQRT(UUU(I,J,K)*UUU(I,J,K)+VVV(I,J,K)*VVV(I,J,K))

!    limit maximum scalar wind speed but keep vector orientation
     IF(VB4(K).GT.UMAX)THEN
        UUU(I,J,K)=UUU(I,J,K)*UMAX/VB4(K)
        VVV(I,J,K)=VVV(I,J,K)*UMAX/VB4(K)
     END IF

!    time step due to advection
     IF(j.GT.1.AND.j.LT.ny) &
        HSTEP=MIN(HSTEP,FRACH*REAL(GSX(I,J))*NGP(J)/MIN(MAX(VB4(K),1.0),UMAX))

!    vertical mixing
     IF(K.GT.1)THEN

!       mixing coefficient (ri method) at top of each box
        TBAR=2.0/(POT(K)+POT(K-1))
        DELT=(POT(K)-POT(K-1))/DZC
        DELU=(VB4(K)-VB4(K-1))/DZC
        RI=DSIGN(998.0D00,DELT)
        IF (DELU.NE.0.0) RI=GRAV*DELT*TBAR/(DELU*DELU)
        RI=MAX(-1.0,MIN(RI,998.0))

!       mixing via nmc-ngm boundary layer method
        ZMIX=VMIX/(2.0+RI)

!       time step due to vertical mixing
        VSTEP=MIN(VSTEP,FRACZ*DZC*DZC/ZMIX)

!       set mixing limits
        IF(PRES.GE.250.0)THEN
           KKK(i,j,k-1)=MAX(0.03, MIN(100.0, ZMIX))
        ELSE
           KKK(i,j,k-1)=MAX(0.03, MIN( 10.0, ZMIX))
        END IF
     END IF

  END DO
  KKK(i,j,nz)=0.0

  END DO
  END DO

! WRITE(kdiag,'(A,2F10.1)')'Minimum time steps (s): ',hstep,vstep
  HSTEP=MIN(HSTEP,VSTEP)

  IF(HSTEP.GE.3600.0)THEN
!    time step in even hours evenly divisible into 24
     DELTA=MAX(1.0,HSTEP/3600.0)
     DO WHILE (MOD(24,DELTA).NE.0)
        DELTA=DELTA-1
     END DO
     DELTA=DELTA*60
  ELSE
!    time step in even minutes evenly divisible into 60
     DELTA=MAX(1.0,HSTEP/60.0)
     DO WHILE (MOD(60,DELTA).NE.0)
        DELTA=DELTA-1
     END DO
     DELTA=MAX(1,DELTA)    
  END IF

END SUBROUTINE stblty


!###############################################################################
! LINREG - Linear regression module used in concentration initialization
!-------------------------------------------------------------------------------
! LAST REVISED: 15 AUG 2005 - initial version
!-------------------------------------------------------------------------------

SUBROUTINE linreg (xvar,yvar,slope,yintc)

  USE constant

  IMPLICIT NONE

  REAL*8, INTENT(IN)  :: xvar(:), yvar(:)
  REAL*8, INTENT(OUT) :: slope, yintc

  INTEGER*4 :: numb
  REAL*8    :: sumx,sumy,sumxy,sumx2

! summations
  NUMB=SIZE(xvar,1)
  SUMX=SUM(xvar)
  SUMY=SUM(yvar)
  SUMXY=SUM(xvar*yvar)
  SUMX2=SUM(xvar*xvar)

! final coefficients
  IF(numb.GE.2)THEN
     SLOPE=(NUMB*SUMXY-SUMX*SUMY)/(NUMB*SUMX2-SUMX*SUMX)
     YINTC=(SUMY-SLOPE*SUMX)/NUMB
  ELSEIF(numb.EQ.1)THEN
     SLOPE=0.0
     YINTC=SUMY
  ELSE
     SLOPE=0.0
     YINTC=0.0
  END IF

  WRITE(kdiag,'(A)')'Regression for concentration initialziation'
  WRITE(kdiag,'(A,2E10.3,I5)')'Slope, Intercept, Number: ',SLOPE,YINTC,NUMB

END SUBROUTINE linreg


!###############################################################################
! XCROSS - Cross-section diagnostic output
!-------------------------------------------------------------------------------
! LAST REVISED: 09 Sep 2005 - initial version
!-------------------------------------------------------------------------------

SUBROUTINE xcross (cfact)

  USE constant
  USE metarray

  IMPLICIT NONE

  REAL*8, INTENT(IN)  :: cfact
  INTEGER*4           :: j,k,nx,ny,nz
  REAL*4              :: clat1,clon1,dlat,dlon

  COMMON /MAINGRID/ clat1,clon1,dlat,dlon
  COMMON /SIZEGRID/ nx,ny,nz

! latitude band labels
  WRITE(kdiag,'(A)')' '
  WRITE(kdiag,'(7X,50F5.0)') (((j-1)*dlat+clat1-dlat/2.0),j=ny,(ny-2),-1),   &
                             (((j-1)*dlat+clat1-dlat/2.0),j=(ny-8),3,-5),    &
                             (((j-1)*dlat+clat1-dlat/2.0),j=2,1,-1)

! print latitude average concentrations by level from the top
  DO K=NZ,1,-1
     DO J=1,NY
        CONC(1,J)=SUM(QQQ(:,j,k))*CFACT*ROW/NX
     END DO
     WRITE(kdiag,'(F6.1,1X,50F5.1)')PPP(k),   &
          (CONC(1,j),j=ny,(ny-2),-1),         &
          (CONC(1,j),j=(ny-8),3,-5),          & 
          (CONC(1,j),j=2,1,-1)
  END DO
  WRITE(kdiag,'(A)')' '

END SUBROUTINE xcross


!###############################################################################
! METDAT - METeorological DATa processing
!-------------------------------------------------------------------------------
! LAST REVISED: 26 Sep 2005 - initial version
!               22 Dec 2005 - sigma surface interpolation
!               16 Feb 2006 - additional terrain smoothing
!               14 Jan 2009 - moved terrain correction after k-loop
!-------------------------------------------------------------------------------

SUBROUTINE metdat

  USE constant
  USE metarray

  IMPLICIT NONE

  REAL*8    :: dzi,dsp
  REAL*8    :: pres,fact,zsfc
  INTEGER*4 :: i,j,k,kb,kt,nx,ny,nz

  COMMON /SIZEGRID/ nx,ny,nz

!--------------------------------------------------------------
! average terrain by averaging surface pressure field
!--------------------------------------------------------------

! main grid 
  DO I=2,(nx-1)
  DO J=2,(ny-1)
     avg(i,j)=0.25*(1.0-twght)*(sfc(i-1,j)+sfc(i+1,j)+sfc(i,j-1)+sfc(i,j+1))+twght*sfc(i,j)
  END DO
  END DO

! east side of edge
  I=1
  DO J=2,(ny-1)
     avg(i,j)=0.25*(1.0-twght)*(sfc(nx,j)+sfc(i+1,j)+sfc(i,j-1)+sfc(i,j+1))+twght*sfc(i,j)
  END DO

! west side of edge
  I=nx
  DO J=2,(ny-1)
     avg(i,j)=0.25*(1.0-twght)*(sfc(i-1,j)+sfc(1,j)+sfc(i,j-1)+sfc(i,j+1))+twght*sfc(i,j)
  END DO

! poles
  avg(:,1) =SUM(sfc(:,1))/nx
  avg(:,ny)=SUM(sfc(:,ny))/nx

! replace surface pressure with averaged field
  sfc=avg

!--------------------------------------------------------------
! interpolation to sigma surfaces
!--------------------------------------------------------------

  DO J=1,ny
  DO I=1,nx

!---------------------------------
!    compute the terrain height

     IF(sfc(i,j).LT.ppp(1))THEN
!       interpolate height from existing points
        zloop : DO k=2,nz
           IF(sfc(i,j).GE.ppp(k))THEN
              IF(sfc(i,j).EQ.ppp(k))THEN
                 zsfc=hhh(i,j,k)
              ELSE
                 fact=(ppp(k-1)-sfc(i,j))/(ppp(k-1)-ppp(k))
                 zsfc=fact*(hhh(i,j,k)-hhh(i,j,k-1))+hhh(i,j,k-1)
              END IF
              EXIT zloop
           END IF
        END DO zloop

     ELSEIF(sfc(i,j).GT.ppp(1))THEN
!       compute height to the ground
        zsfc=HHH(I,J,1)-0.5*(HHH(I,J,2)-HHH(I,J,1))

     ELSE
!       coincidently perfect match
        zsfc=hhh(i,j,1)
     END IF

!----------------------------------
!    compute the density profile

!##  DZI=0.0        
!##  DSP=1013.0

     DO K=1,nz

!##     thermodynamic method to recompute height fields
!##     RRR(I,J,K)=P2JM*PPP(K)/(TTT(I,J,K)*RDRY)
!##     DZI=DZI+P2JM*(DSP-PPP(K))/(RRR(I,J,K)*GRAV)
!##     HHH(I,J,K)=DZI
!##     DSP=PPP(K)

!       check height field for consistency (needed for error check on external data)
        IF((K.GT.1).AND.(HHH(I,J,K).LT.HHH(I,J,K-1)))THEN
           WRITE(KDIAG,'(A,3I5,F10.1)')'Height inversion at pressure level: ',I,J,K,ppp(K)
           WRITE(KDIAG,'(F8.1,A,F8.1)')HHH(I,J,K),' should be closer to',STDATM(KNDX(K))
           HHH(I,J,K)=STDATM(KNDX(K))
        END IF

!       use vertical box sizes between interfaces (DZI) to compute density
        IF(K.EQ.1) THEN
           DZI=HHH(I,J,K+1)-HHH(I,J,K)
           DSP=PPP(K)-PPP(K+1)
        ELSEIF (K.EQ.NZ) THEN
           DZI=HHH(I,J,K)-HHH(I,J,K-1)
           DSP=PPP(K-1)-PPP(K)
        ELSE
           DZI=0.5*(HHH(I,J,K+1)-HHH(I,J,K-1))
           DSP=0.5*(PPP(K-1)-PPP(K+1))
        END IF

!       average density (g/kg) within each vertical box
        RRR(I,J,K)=P2JM*DSP/DZI/GRAV

     END DO

!    optional adjustment for  heights relative to terrain
!    does not affect subsequent computations based upon height differences
     HHH(i,j,:)=HHH(i,j,:)-zsfc

!-----------------------------------
!    interpolate to sigma surfaces

     DO K=1,nz
!       pressure at the interpolation sigma level
        PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP        

        kt=1
        DO WHILE (kt.LT.nz.AND.PPP(kt).GT.PRES)
!          find first index above sigma level
           kt=kt+1
        END DO
        kb=MAX(1,kt-1)

!       vertical interpolation factor
        IF(kt.GT.kb)THEN
           FACT=(PPP(kb)-PRES)/(PPP(kb)-PPP(kt))
        ELSE
           FACT=0.0
        END IF
  
        VAL(k,1)=HHH(i,j,kb)-FACT*(HHH(i,j,kb)-HHH(i,j,kt))
        VAL(k,2)=TTT(i,j,kb)-FACT*(TTT(i,j,kb)-TTT(i,j,kt)) 
        VAL(k,3)=UUU(i,j,kb)-FACT*(UUU(i,j,kb)-UUU(i,j,kt)) 
        VAL(k,4)=VVV(i,j,kb)-FACT*(VVV(i,j,kb)-VVV(i,j,kt)) 
        VAL(k,5)=WWW(i,j,kb)-FACT*(WWW(i,j,kb)-WWW(i,j,kt)) 
        VAL(k,6)=MMM(i,j,kb)-FACT*(MMM(i,j,kb)-MMM(i,j,kt)) 
        VAL(k,7)=RRR(i,j,kb)-FACT*(RRR(i,j,kb)-RRR(i,j,kt)) 
     END DO

!    put data back into 3D array
     HHH(i,j,:)=VAL(:,1)
     TTT(i,j,:)=VAL(:,2)
     UUU(i,j,:)=VAL(:,3)
     VVV(i,j,:)=VAL(:,4)
     WWW(i,j,:)=VAL(:,5)
     MMM(i,j,:)=VAL(:,6) 
     RRR(i,j,:)=VAL(:,7) 

  END DO
  END DO

END SUBROUTINE metdat
