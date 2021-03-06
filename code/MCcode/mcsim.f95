!---------------------------------------------------------------*
      
!     This subroutine performs a Monte Carlo simulation on the 
!     polymer chain.
      
      SUBROUTINE MCsim(R,U,NT,N,NP,NSTEP,BROWN,INTON,IDUM,PARA, &
	  MCAMP,SUCCESS,SUCCESS_TOTAL,MOVEON,WINDOW,RING,TWIST,Lk,LT,LP,L)

      use mt19937, only : grnd, sgrnd, rnorm, mt, mti
      IMPLICIT NONE

      
      REAL, PARAMETER :: PI=3.141592654 ! Value of pi
      
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION RP(NT,3)  ! Bead positions
      DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
      INTEGER N,NP,NT           ! Number of beads
      INTEGER NSTEP             ! Number of MC steps
      INTEGER BROWN            ! Turn on fluctuations
      INTEGER INTON             ! Include polymer interactions
      INTEGER RING              ! Is polymer a ring?
      INTEGER TWIST             ! Include twist?
      INTEGER Lk                ! Linking number
      DOUBLE PRECISION L        ! Contour length

      DOUBLE PRECISION ETWIST
!     Variables for the simulation
      
      INTEGER ISTEP             ! Current MC step index
      DOUBLE PRECISION PROB     ! Calculated test prob
      DOUBLE PRECISION TEST     ! Random test variable
      INTEGER IB                ! Test bead 
      INTEGER IP                ! Test polymer 
      INTEGER IB1               ! Test bead position 1
      INTEGER IT1               ! Index of test bead 1
      INTEGER IB2               ! Test bead position 2
      INTEGER IT2               ! Index of test bead 2
      INTEGER DIB               ! Number of segments in portion of chain moved
      INTEGER TEMP
      REAL ran1                 ! Random number generator
      INTEGER IDUM              ! Seed for the generator
      INTEGER NOW(3)            ! Time now (hr,min,sec)
      INTEGER I,J
      DOUBLE PRECISION R0(3)
      
!     Energy variables
      
      DOUBLE PRECISION DEELAS   ! Change in bending energy
      DOUBLE PRECISION DESELF   ! Change in self energy
      DOUBLE PRECISION ESELF    ! Self energy of current structure
      DOUBLE PRECISION ESELFP   ! Self energy of test structure
      DOUBLE PRECISION DEEX     ! Change in external energy
      DOUBLE PRECISION ENERGY
           
!     MC adaptation variables
      
      DOUBLE PRECISION MCAMP(6) ! Amplitude of random change      
      INTEGER MCTYPE            ! Type of MC move
      INTEGER NADAPT(6)         ! Num steps btwn adapt
      DOUBLE PRECISION PHIT     ! % hits per total steps
      DOUBLE PRECISION PDESIRE(6) ! Desired hit rate
      INTEGER SUCCESS(6),SUCCESS_TOTAL(6)        ! Number of successes
      DOUBLE PRECISION MINAMP(6) ! Minimum amp to stop
      DOUBLE PRECISION MAXAMP(6) ! Minimum amp to stop
      INTEGER MOVEON(6)			! Is the move active
      INTEGER WINDOW(6)			! Size of window for bead selection
      INTEGER WINDOW_MAX(6)             ! Maximum size of window for bead selection
!     Alexander Polynomial Variables
      DOUBLE PRECISION, ALLOCATABLE :: CROSS(:,:)   !Matrix of information for crossings in a 2-D projection of the polymer
      DOUBLE PRECISION, ALLOCATABLE :: CROSSP(:,:)  !Matrix of crossings for the trial configuration
      INTEGER NCROSS
      INTEGER NCROSSP
      INTEGER CrossSize
      INTEGER DELTA             !Alexander polynomial evaluated at t=-1; used for knot checking
      INTEGER DELTAP            !Alexandper polynomial of trial configuration

!     Variables in the simulation
      
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LT       ! Twist persistence length
      DOUBLE PRECISION LP       ! Bend persistence length
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION WR       ! Writhe of current structure
      DOUBLE PRECISION WRP      ! Writhe of test structure
      
!     Load the input parameters
      
      EB=PARA(1)
      EPAR=PARA(2)
      EPERP=PARA(3)
      GAM=PARA(4)
      ETA=PARA(5)
      XIR=PARA(6)
      XIU=PARA(7)
      LBOX=PARA(8)
      LHC=PARA(9)
      VHC=PARA(10)

!     Set-up MC Adaptation variables

      MINAMP(1)=0.07*2*PI
      MINAMP(2)=(1.0E-1)*L/N
      MINAMP(3)=0.0*PI
      MINAMP(4)=0.0*PI
      MINAMP(5)=0.1*PI
      MINAMP(6)=0.01
	  
      MAXAMP(1)=2.*PI
      MAXAMP(2)=LBOX
      MAXAMP(3)=2.*PI
      MAXAMP(4)=2.*PI
      MAXAMP(5)=2.*PI
      MAXAMP(6)=LBOX

      WINDOW_MAX=N/2
	  
      NADAPT(1)=1000
      NADAPT(2)=1000
      NADAPT(3)=1000
      NADAPT(4)=1000
      NADAPT(5)=1000
      NADAPT(6)=1000


      if (NSTEP.LE.NADAPT(1)) then
         NADAPT(1)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(2)) then
         NADAPT(2)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(3)) then
         NADAPT(3)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(4)) then
         NADAPT(4)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(5)) then
         NADAPT(5)=NSTEP
      endif
	  
      PDESIRE(1)=0.5
      PDESIRE(2)=0.5
      PDESIRE(3)=0.5
      PDESIRE(4)=0.5
      PDESIRE(5)=0.5
      PDESIRE(6)=0.5

      SUCCESS(1)=0
      SUCCESS(2)=0
      SUCCESS(3)=0
      SUCCESS(4)=0
      SUCCESS(5)=0
      SUCCESS(6)=0
      SUCCESS_TOTAL=0

      DEELAS=0.
      DESELF=0.
      DEEX=0.

!     Get initial Writhe
      call WRITHE(R,N,Wr)

!     Initialize the Cross matrix

      CrossSize=N**2
      ALLOCATE(Cross(CrossSize,6))
      ALLOCATE(CrossP(CrossSize,6))

!     Get initial value of Alexander polynomial and Cross matrix
      NCross=0
      CALL ALEXANDERP(R,N,DELTA,Cross,CrossSize,NCross)
!     Begin Monte Carlo simulation

          
      ISTEP=1
           
      DO WHILE (ISTEP.LE.NSTEP)
         
         DO 10 MCTYPE=1,6

            if (MOVEON(MCTYPE).EQ.0) then
               goto 60
            endif
		             
            CALL MC_move_new(R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,IDUM,MCTYPE,MCAMP,WINDOW,RING,DIB)
     
!     If chain is a ring, calculate the alexander polynomial evaluated at t=-1. Reject the move if the chain becomes knotted (DELTA.NE.1)
!     This is currently only set-up to handle one chain. Need to make modifications to handle multiple chain systems
            
            IF (RING.EQ.1) THEN
               CrossP=Cross
               NCrossP=NCross
               IF (MCTYPE.EQ.1) THEN
                  CALL alexanderp_crank(RP,N,DELTA,CrossP,CrossSize,NCrossP,IT1,IT2,DIB)
               ELSEIF (MCTYPE.EQ.2) THEN
                  IF (DIB.NE.N) THEN
                     CALL alexanderp_slide(RP,N,DELTA,CrossP,CrossSize,NCrossP,IT1,IT2,DIB)
                  ENDIF
               ELSE
                  CALL ALEXANDERP(RP,N,DELTA,CrossP,CrossSize,NCrossP)
               ENDIF
               IF (DELTA.NE.1) THEN
                  PROB=-1.
                  GOTO 80   !Skip check for move acceptance (reject move)
               ENDIF
            ENDIF
!     Calculate the change in compression and bending energy
            IF (MCTYPE.NE.5.AND.MCTYPE.NE.6) THEN
             CALL  MC_eelas(DEELAS,R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,EB,EPAR,EPERP,GAM,ETA,RING,TWIST,Lk,lt,LP,L,MCTYPE,WR,WRP)
            ENDIF
            
!     Calculate the change in the self-interaction energy
            
            IF (INTON.EQ.1) THEN
               IF (MCTYPE.EQ.1) THEN
                  CALL DE_SELF_CRANK(DESELF,R,RP,NT,N,NP,PARA,RING,IB1,IB2)
            
               ELSEIF (MCTYPE.EQ.2) THEN
                  CALL ENERGY_SELF_SLIDE(ESELF,R,NT,N,NP,PARA,RING,IB1,IB2)
                  CALL ENERGY_SELF_SLIDE(ESELFP,RP,NT,N,NP,PARA,RING,IB1,IB2)

            
                  DESELF=ESELFP-ESELF
               ELSE
                  DESELF=0.
               ENDIF
            ENDIF
           
          
            
!     Change the position if appropriate
            
            ENERGY=DEELAS+DESELF
           
            PROB= exp(-ENERGY)

            if (BROWN.EQ.1) then
               TEST=grnd()
            else
               TEST=1.

            endif
            if (TEST.LE.PROB) then
              
               DO 20 I=N*(IP-1)+1,N*IP
                  R(I,1)=RP(I,1)
                  R(I,2)=RP(I,2)
                  R(I,3)=RP(I,3)
                  U(I,1)=UP(I,1)
                  U(I,2)=UP(I,2)
                  U(I,3)=UP(I,3)
 20            CONTINUE 
                WR=WRP
                NCross=NCrossP
                Cross=CrossP
               
               SUCCESS(MCTYPE)=SUCCESS(MCTYPE)+1
               SUCCESS_TOTAL(MCTYPE)=SUCCESS_TOTAL(MCTYPE)+1
               
            endif

80         CONTINUE
            
!     Adapt the amplitude of step every NADAPT steps
            
            if (mod(ISTEP,NADAPT(MCTYPE)).EQ.0) then
               PHIT=real(SUCCESS(MCTYPE))/real(NADAPT(MCTYPE))
               if (PHIT.GT.PDESIRE(MCTYPE)) then
                  MCAMP(MCTYPE)=MCAMP(MCTYPE)*1.05
               else
                  MCAMP(MCTYPE)=MCAMP(MCTYPE)*0.95
               endif
               if (MCAMP(MCTYPE).GT.MAXAMP(MCTYPE)) then
                  MCAMP(MCTYPE)=MAXAMP(MCTYPE)
                  WINDOW(MCTYPE)=WINDOW(MCTYPE)+1
               endif
               if (MCAMP(MCTYPE).LT.MINAMP(MCTYPE)) then
                  MCAMP(MCTYPE)=MINAMP(MCTYPE)
                  WINDOW(MCTYPE)=WINDOW(MCTYPE)-1
               endif

               IF (WINDOW(MCTYPE).LT.1) THEN
                  WINDOW(MCTYPE)=1
               ENDIF

               IF (WINDOW(MCTYPE).GT.WINDOW_MAX(MCTYPE)) THEN
                  WINDOW(MCTYPE)=WINDOW_MAX(MCTYPE)
               ENDIF
			   
               SUCCESS(MCTYPE)=0
               
               IB=1
               DO 40 I=1,NP
                  R0(1)=nint(R(IB,1)/LBOX-0.5)*LBOX
                  R0(2)=nint(R(IB,2)/LBOX-0.5)*LBOX
                  R0(3)=nint(R(IB,3)/LBOX-0.5)*LBOX
                  DO 50 J=1,N
                     R(IB,1)=R(IB,1)-R0(1)
                     R(IB,2)=R(IB,2)-R0(2)
                     R(IB,3)=R(IB,3)-R0(3)
                     IB=IB+1
 50              CONTINUE
 40           CONTINUE
               
            endif

 60         CONTINUE


 10      CONTINUE
         
         ISTEP=ISTEP+1   
         
      ENDDO

         
         
      RETURN      
      END
      
!---------------------------------------------------------------*
