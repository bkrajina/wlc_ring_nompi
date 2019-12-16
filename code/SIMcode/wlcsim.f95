!---------------------------------------------------------------*
      
      PROGRAM wlcsim
      
!     This simulation models the equilibrium conformation statistics of DNA molecules modeled
!     as a discrete,stretchable, shearable worm-like chain (see Koslover, 2013) using a Metropolis
!     Monte Carlo procedure. Brownian dynamics simulations are also included.

!     The program involves performing a series of monte carlo simulations with a specified number of monte
!     carlo steps. After each monte carlo simulation, the configuration, energy, and structural parameters
!     of the resulting polymer configuration are saved for future analysis.

!     This version of the simulation allows for a parallel tempering scheme in which multiple "replicas"
!     with different linking numbers are simulated in parallel, periodically testing for exchange between
!     replicas. Note that this parallel tempering scheme applies only to closed-circular DNA.

!     DNA molecules may be linear or circular. Topology of circular DNA molecules are specified
!     by their linking number. This allows the properties of supercoiled DNA molecules to be studied.

!     Brad Krajina, Andrew J. Spakowitz.
!     Last modified: 10/2/2014     


      use mt19937, only : grnd, sgrnd, rnorm, mt, mti
      IMPLICIT NONE

      REAL, PARAMETER :: PI=3.141592654 ! Value of pi

      !Variables used in simulation 

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R	 ! Conformation of polymer chains
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U	 ! Conformation of polymer chains
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R0	 ! Conformation of polymer chains
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U0	 ! Conformation of polymer chains
      
   
      DOUBLE PRECISION L0       ! Equilibrium segment length
      DOUBLE PRECISION DEL      ! Segment contour length
      DOUBLE PRECISION ENERGY   ! Total energy
      DOUBLE PRECISION TIME     ! Current time
      DOUBLE PRECISION TSAVE     ! Time of save point
      DOUBLE PRECISION T0,TF    ! Initial/final times
      DOUBLE PRECISION DT       ! Time step size
      INTEGER I,J,IB            ! Index
      INTEGER INDMAX            ! Maximum index in series
      INTEGER IND               ! Ind in series
      INTEGER TENS              ! Decimal of index
      character*5 fileind       ! Index of output
      character*16 snapnm       ! File for output
   
     
!     Simulation input variables

      INTEGER NT                 ! Number of beads in simulation
      INTEGER N                 ! Number of beads in simulation
      INTEGER NP                ! Number of polymers in simulation
      DOUBLE PRECISION L        ! Total contour length
      DOUBLE PRECISION LP       ! Bending persistencce length
      DOUBLE PRECISION LT       ! Twist persistence length
      INTEGER LK                ! Linking number
      INTEGER FRMFILE           ! Initial condition
      INTEGER BROWN             ! Include Brownian forces
      INTEGER INTON             ! Include polymer interactions
      INTEGER LOGTIME           ! Is data recorded in log time?
      INTEGER RING              ! Is polymer a ring?
      INTEGER TWIST             ! Include twist?
      DOUBLE PRECISION DT0      ! Initial time step size
      INTEGER NSTEP,NINIT

!     Monte Carlo variables

      DOUBLE PRECISION MCAMP(6) ! Amplitude of random change
      INTEGER MOVEON(6)			! Is the move active
      INTEGER WINDOW(6)			! Size of window for bead selection
      INTEGER SUCCESS(6),SUCCESS_TOTAL(6)        ! Number of successes

!     Parallel tempering variables

      INTEGER ParTempOn         !Use parallel tempering?
      INTEGER LkMin             !Minimum LK for parallel tempering
      INTEGER LkMax             !Maximum LK for parallel tempering
      INTEGER LkStep            !Step size between LK values
      INTEGER LkMinus           !Lk for replica one Lk step below current
      INTEGER LkPlus            !Lk for replica one Lk step above current
      INTEGER LkSwap            !LK value to test for replica exchange with current simulation
      INTEGER LkPlusInd         !Simulation index of LK replica one LK step above current replica
      INTEGER LkMinusInd          !Simulation index of LK replica one LK step below current replica
      CHARACTER*4 LkPlusStr     !String variable for Lk replica one Lk step above current replica
      CHARACTER*4 LkMinusStr      !String variable for Lk replica one LK step below current replica
      CHARACTER*4 LkSwapStr     !String variable for Lk replica to test for exchange
      INTEGER TensMinus        ! Number of characters requried to specify Lk in string format for LkMinus
      INTEGER TensPlus         !Number of characters required to specify Lk in string format for LkPlus
      INTEGER TensSwap
      DOUBLE PRECISION EReplica(5)  !Energy vector of replica to test for exchange
      DOUBLE PRECISION ECurrent(5)   !Current energy vector for current replica
      DOUBLE PRECISION EReplicaTotOld  !Total energy of replica configuration to test for exchange when evaluated under its hamiltonian
      DOUBLE PRECISION EReplicaTotNew  !Total energy of replica configuration to test for exchange when evaluated under this hamiltonian
      DOUBLE PRECISION ETotNew         !Total energy of this replica's configuration when evaluated in the replica's Hamiltonian
      DOUBLE PRECISION WrReplica !Writhe of replica to test for exchnage
      DOUBLE PRECISION TwReplica !Twist of replica to test for exchange when evaluated with the current linking number
      DOUBLE PRECISION TwNew     !Twist of the current conformation when evaluated with the replica's linking number
      DOUBLE PRECISION DEnergy  !Difference in energy between current replica and replica to test for exchange
      DOUBLE PRECISION Prob     !Probability of swapping current replica with test replica 
      DOUBLE PRECISION Test     !Test number for determining whether to swap (random number between 0. and 1.)
      INTEGER NSwapPlus         !Number of times replica with Lk one step above is swapped
      INTEGER NSwapMinus        !Number of times replica with Lk one step below is swapped
      INTEGER NTrialPlus        !Number of trials for exchanging with the replica on Lk step above 
      INTEGER NTrialMinus       !Number of trials for exchanging with the replica on Lk step below
      INTEGER IOStatus          !IO status variable for read operations
      INTEGER LkMinusComp       !Number of MC simulations completed (including replica exchange) for replica one Lk step below
      INTEGER LkPlusComp        !Number of MC simulations completed (including replica exchange) for replica one Lk step above
      INTEGER Accept            !Accept conformation of replica (1 or 0)?
      INTEGER ReplicaAccept     !Does replica to test for exchange accept conformation of this replica (1 or 0)?
      INTEGER UpOrDown          !Swap with replica with Lk above or below the current replica?
      INTEGER ReplicaIndex     !Index of this replica in relation to others (replicas are indexed starting at 0 for Lkmin)
      INTEGER ReplicaExchangeIndex !Number of exchanges tested by the current trial replica. This is used for synchronization
      INTEGER ExchangeMinusIndex   !Number of exchanges tested by the replica with one Lk step above the current
      INTEGER ExchangePlusIndex    !Number of exchanges tested by the replica with one Lk step below the current
      INTEGER ExchangeIndex
      INTEGER MC_Completed
      INTEGER Restart              !Is simulation restarting from a previously interrupted run (1 or 0)?
      DOUBLE PRECISION, ALLOCATABLE ::  Swap(:)     !Vector of swaps with other replicas (+1,-1,or 0)
!     Energy variables
      
      DOUBLE PRECISION EELAS(4) ! Elastic energy
      DOUBLE PRECISION EPONP    ! Poly-poly energy
      DOUBLE PRECISION ETOT     ! Total chain energy
      DOUBLE PRECISION, ALLOCATABLE :: EELASALL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EPONPALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: ETotAll(:)
      DOUBLE PRECISION  ETotAvg
      DOUBLE PRECISION  ETotStdev
      DOUBLE PRECISION  EtotStderr
      
!     Structure analysis
      
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: RCOM ! Center of mass
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  RCOMSQ  ! Center of mass squared
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  RGYRSQ  ! Radius of gyration squared
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  RGYRSQALL ! Radius of gyration squared matrix (all polymers and MC steps)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RGYR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: RGYRALL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: R2  !Second moment of end-to-end displacement 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: R4  !Fourth moment of end-to-end displacement 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: R6  !Sixth moment of end-to-end displacement 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DR  !End-to-end displacement 


      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R2ALL !Second moment matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R4ALL !Fourth moment matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R6ALL !Sixth moment matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRALL !Displacement matrix


      DOUBLE PRECISION RGYRSQ_AVG
      DOUBLE PRECISION RGYRSQ_STDEV
      DOUBLE PRECISION RGYRSQ_STDER
      DOUBLE PRECISION RGYR_AVG
      DOUBLE PRECISION RGYR_STDEV
      DOUBLE PRECISION RGYR_STDER
      DOUBLE PRECISION R2_AVG
      DOUBLE PRECISION R2_STDEV
      DOUBLE PRECISION R2_STDER
      DOUBLE PRECISION R4_AVG
      DOUBLE PRECISION R4_STDEV
      DOUBLE PRECISION R4_STDER
      DOUBLE PRECISION R6_AVG
      DOUBLE PRECISION R6_STDEV
      DOUBLE PRECISION R6_STDER
      DOUBLE PRECISION DR_AVG
      DOUBLE PRECISION DR_STDEV
      DOUBLE PRECISION DR_STDER

      DOUBLE PRECISION WR_AVG
      DOUBLE PRECISION WR_STDEV
      DOUBLE PRECISION WR_STDER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WrAll !Writhe of all iterations
      DOUBLE PRECISION Wr       ! Writhe
      DOUBLE PRECISION Tw       ! Twist
      DOUBLE PRECISION DELR(3)  ! Mag of gyration tensor
      DOUBLE PRECISION RCOM0(3) ! Init val RCOM
      DOUBLE PRECISION DELR0(3) ! Init val DELR
      DOUBLE PRECISION DRCOM    ! Change in RCOM
      DOUBLE PRECISION SIG(3,3)
      DOUBLE PRECISION COR
     
!     Algorithm analysis variables
      DOUBLE PRECISION, ALLOCATABLE ::  MCAMP1ALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: MCAMP2ALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: WINDOW1ALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: WINDOW2ALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: MCSTEPCUM(:)
      DOUBLE PRECISION, ALLOCATABLE :: SUCCESSALL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PHITALL(:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  RGYRSQ_AUTO(:)  !auto-correlation of radius of gyration
      DOUBLE PRECISION, ALLOCATABLE :: RSQ_AUTO(:)     !auto-correlation of RSQ
      DOUBLE PRECISION, ALLOCATABLE :: Wr_AUTO(:)                      !auto-correlation of Wr
      DOUBLE PRECISION, ALLOCATABLE :: Energy_AUTO(:)                      !auto-correlation of Wr

      INTEGER N_auto !maximum difference in indices between elements on which auto-correlation is computed

!     Variables in the simulation
      
      DOUBLE PRECISION PARA(10)

!     Variables for the random number generators

      INTEGER IDUM              ! Seed for the generator
      DOUBLE PRECISION MOM(6)


!     Load in the parameters for the simulation


      open (unit=5, file='input/input')
      read (unit=5, fmt='(32(/))')
      read (unit=5, fmt=*) N
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) TF
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) INDMAX
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) DT
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) FRMFILE
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) BROWN
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) INTON
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) RING
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) TWIST
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LOGTIME
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NINIT
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) NSTEP
      close(5)

!     Read from the parallel tempering input file

      open (unit=1,file='input/partemp')
      read (unit=1, fmt='(4(/))')
      read (unit=1,fmt=*) ParTempOn
      read (unit=1,fmt='(2(/))')
      read (unit=1,fmt=*) LkMin
      read (unit=1,fmt='(2(/))')
      read (unit=1,fmt=*) LkMax
      read (unit=1,fmt='(2(/))')
      read (unit=1,fmt=*) LkStep
      close(1)

!     Get elastic parameters and allocate vectors

      call getpara(PARA,DT,DEL,L,LP,LT,Lk,RING)
      DT0=DT
	  
      NT=N*NP
      ALLOCATE(R(NT,3))
      ALLOCATE(U(NT,3))
      ALLOCATE(R0(NT,3))
      ALLOCATE(U0(NT,3))
      ALLOCATE(RCOM(NP,3))
      ALLOCATE(RCOMSQ(NP))
      ALLOCATE(RGYRSQ(NP))
      ALLOCATE(RGYR(NP))
      ALLOCATE(RGYRSQALL(INDMAX+1,NP))
      ALLOCATE(RGYRALL(INDMAX+1,NP))
      ALLOCATE(R2(NP))
      ALLOCATE(R4(NP))
      ALLOCATE(R6(NP))
      ALLOCATE(DR(NP))
    
      ALLOCATE(R2ALL(INDMAX+1,NP))
      ALLOCATE(R4ALL(INDMAX+1,NP))
      ALLOCATE(R6ALL(INDMAX+1,NP))
      ALLOCATE(DRALL(INDMAX+1,NP))

      ALLOCATE(WrAll(INDMAX+1))
      ALLOCATE(MCAMP1ALL(INDMAX+1))
      ALLOCATE(MCAMP2ALL(INDMAX+1))
      ALLOCATE(WINDOW1ALL(INDMAX+1))
      ALLOCATE(WINDOW2ALL(INDMAX+1))
      ALLOCATE(MCSTEPCUM(INDMAX+1))
      ALLOCATE(SUCCESSALL(INDMAX+1,6))
      ALLOCATE(PHITALL(INDMAX+1,6))
      ALLOCATE(EELASALL(INDMAX,4))
      ALLOCATE(EPONPALL(INDMAX))
      ALLOCATE(ETotAll(INDMAX))

      ALLOCATE(Swap(INDMAX))
      
      N_auto=INDMAX/10

      ALLOCATE(RGYRSQ_AUTO(N_auto))
      ALLOCATE(RSQ_AUTO(N_auto))
      ALLOCATE(Wr_AUTO(N_auto))
      ALLOCATE(Energy_AUTO(N_auto))
      
!     Set which MC moves to use
      MOVEON(1)=1
      MOVEON(2)=1
      MOVEON(3)=1
      MOVEON(4)=1
      
      if (RING.EQ.1) then
         MOVEON(3)=0
      endif
   
      if (INTON.EQ.1.AND.NP.GT.1) then
         MOVEON(5)=1
         MOVEON(6)=1
      else
         MOVEON(5)=0
         MOVEON(6)=0
      endif
   

!     Setup the initial condition and seed the randomm number generator

      call initcond(R,U,NT,N,NP,IDUM,FRMFILE,PARA,RING)


!     Parallel tempering: Determine Lk values for replica with one Lk step above and below
!     The current Lk. Also write them to character variables
      
      IF (ParTempOn.EQ.1) THEN
         IF (Lk.EQ.LkMax) THEN
            LkPlus=Lk
            LkMinus=Lk-LkStep
         ELSEIF (Lk.EQ.LkMin) THEN
            LkPlus=Lk+LkStep
            LkMinus=Lk
         ELSE
            LkPlus=Lk+LkStep
            LkMinus=Lk-Lkstep
         ENDIF
         write(LkMinusStr,fmt='((I4))') LkMinus
         write(LkPlusStr,fmt='((I4))') LkPlus
         TensMinus=nint(log10(1.*abs(LkMinus))-0.49999)+1
         TensPlus=nint(log10(1.*abs(LkPlus))-0.49999)+1
         IF (LkMinus.LT.0) THEN
            TensMinus=TensMinus+1
         ENDIF
         IF (LkPlus.LT.0) THEN
            TensPlus=TensPlus+1
         ENDIF

         !Determine the index of this replica in relation to the others
         !Initialize up or down variable depending on whether the index is even or odd.
         !This ensures that replicas are paired with one another for exchange

         ReplicaIndex=(Lk-LkMin)/LkStep
         IF (mod(ReplicaIndex,2).EQ.0) THEN
            UpOrDown=1
         ELSE
            UpOrDown=-1
         ENDIF

         !Restart code: Check if MC_IND has been written to. If so, read the conformation as the current
         !conformation and resume the MC simulation

         OPEN(unit=1,file='data/MC_IND',IOSTAT=IOStatus,status='old')
         IF (IOStatus.EQ.0) THEN
            Restart=1
            READ(unit=1,fmt=*) IND
            CLOSE(unit=1)
            
        
            !Read the exchange index and the MC_complete index to determine where to resume.
            OPEN(unit=1,file='data/ExchangeIndex',IOSTAT=IOStatus,status='old')
            READ(unit=1,fmt=*) ExchangeIndex
            CLOSE(unit=1)

            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF

            OPEN(unit=1,file='data/MC_Completed',IOSTAT=IOStatus,status='old')
            READ(unit=1,fmt=*) MC_Completed
            CLOSE(unit=1)

            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF

            OPEN(unit=1,file='data/UpOrDown',IOSTAT=IOStatus,status='old')
            READ(unit=1,fmt=*) UpOrDown
            CLOSE(unit=1)
            
            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF

            !Read the polymer energy,writhe, and configuration
            
            OPEN(unit=1,file='data/E_Current',IOSTAT=IOStatus,status='old')
            READ(unit=1,fmt=*) ECurrent
            CLOSE(unit=1)

            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF


            ETOT=SUM(ECurrent)

            OPEN(unit=1,file='data/Wr_Current',IOSTAT=IOStatus,status='old')
            READ(unit=1,fmt=*) Wr
            CLOSE(unit=1)
            
            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF

            
            TENS=nint(log10(1.*IND)-0.499999)+1
            write (fileind,'(I5)'), IND
            snapnm= 'data/r'//fileind((5-TENS+1):5)
            OPEN(unit=1,file=snapnm,IOSTAT=IOStatus,status='old')
           
            DO I=1,NT
               READ(unit=1,fmt=*) R(I,:)
            ENDDO
            CLOSE(unit=1)
            
            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF


            snapnm= 'data/u'//fileind((5-TENS+1):5)
           
            OPEN(unit=1,file=snapnm,IOSTAT=IOStatus,status='old')
            DO I=1,NT
               READ(unit=1,fmt=*) U(I,:)
            ENDDO
            CLOSE(unit=1)
            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF

            OPEN (unit=1,file='data/WINDOW',IOSTAT=IOStatus,status='old')
            DO I=1,6
               READ(1,*) WINDOW(I)
            ENDDO
            CLOSE(1)

            IF (IOStatus.NE.0) THEN
                 Restart=0
            ENDIF

            OPEN (unit=1,file='data/MCAMP',IOSTAT=IOStatus,status='old')
            DO I=1,6
               READ(1,*) MCAMP(I)
            ENDDO
            CLOSE(1)

            IF (IOStatus.NE.0) THEN
               Restart=0
            ENDIF

            !Resume the simulation 
            IF (Restart.EQ.1) THEN
               GOTO 10
            ENDIF
         ELSE
            Restart=0
            CLOSE(unit=1)
         ENDIF
         !Write the MCIND to file
         
         OPEN (unit=1,file='data/MC_IND',status='REPLACE')
         WRITE(unit=1,fmt=*) 0
         CLOSE (unit=1)

      ENDIF
      
	  
!     Perform an initialization MC simulation

      MCAMP(1)=1.
      MCAMP(2)=1.
      MCAMP(3)=1.
      MCAMP(4)=1.
      MCAMP(5)=1.
      MCAMP(6)=1.

      WINDOW=N
     
   
      CALL    MCsim(R,U,NT,N,NP,NINIT,BROWN,INTON,IDUM,PARA, MCAMP,SUCCESS,SUCCESS_TOTAL,MOVEON,WINDOW,RING,TWIST,Lk,LT,LP,L)
     

!     Save the conformation and PSI angles 
      
      OPEN (UNIT = 1, FILE = 'data/r0', STATUS = 'REPLACE')      
      IB=1
      DO  I=1,NP
         DO  J=1,N
            R0(IB,1)=R(IB,1)
            R0(IB,2)=R(IB,2)
            R0(IB,3)=R(IB,3)
            U0(IB,1)=U(IB,1)
            U0(IB,2)=U(IB,2)
            U0(IB,3)=U(IB,3)
            WRITE(1,*) R(IB,1),R(IB,2),R(IB,3)
            IB=IB+1
         ENDDO
      ENDDO
      CLOSE(1)
      
      OPEN (UNIT = 1, FILE = 'data/u0', STATUS = 'REPLACE')      
      IB=1
      DO  I=1,NP
         DO  J=1,N
            WRITE(1,*) U(IB,1),U(IB,2),U(IB,3)
            IB=IB+1
         ENDDO
      ENDDO
      CLOSE(1)

!     Structural Analysis
 
      CALL getdim(N,NP,NT,R,RCOM,RCOMSQ,RGYRSQ,R2,R4,R6,DR)
      
      call WRITHE(R,N,Wr)
      WrALL(1)=Wr
      RGYRSQALL(1,:)=RGYRSQ
      RGYRALL(1,:)=RGYRSQ**0.5
      R2ALL(1,:)=R2
      R4ALL(1,:)=R4
      R6ALL(1,:)=R6
      DRALL(1,:)=DR


!     Algorithm analysis: write to MCAMP and MCSTEP vectors
      MCAMP1ALL(1)=MCAMP(1)
      MCAMP2ALL(1)=MCAMP(2)
      WINDOW1ALL(1)=WINDOW(1)
      WINDOW2ALL(1)=WINDOW(2)
      MCSTEPCUM(1)=NINIT
      SUCCESSALL(1,:)=SUCCESS_TOTAL
      PHITALL(1,:)=SUCCESS_TOTAL/NINIT
     
     
!     Begin simulation
      
      IND=1
      TIME=0.
      NSwapPlus=0
      NSwapMinus=0
      LkPlusInd=0
      LkMinusInd=0
      NTrialMinus=0
      NTrialPlus=0
      LkMinusComp=0
      LkPlusComp=0
      ExchangePlusIndex=0
      ExchangeMinusIndex=0

      10 CONTINUE

      DO WHILE (IND.LE.INDMAX)

!     Restart code:If simulation is restarting, go to restart continue point

         IF (Restart.EQ.1) THEN
            GOTO 20
         ENDIF

!     Perform a MC simulation

         CALL  MCsim(R,U,NT,N,NP,NSTEP,BROWN,INTON,IDUM,PARA, MCAMP,SUCCESS,SUCCESS_TOTAL,MOVEON,WINDOW,RING,TWIST,Lk,LT,LP,L)
         
         
!     Perform a time step
         
         if (LOGTIME.EQ.0) then
            TSAVE = TF*IND/INDMAX
         else
            TSAVE = DT0*exp((IND-1.)/(INDMAX-1.)*log(TF/DT0))		 
         endif
         if (NSTEP.EQ.0) then
            call BDsim(R,U,NT,N,NP,TIME,TSAVE,DT,BROWN,INTON,IDUM,PARA)
         endif
         
!     Save the conformation and the metrics
         
         TENS=nint(log10(1.*IND)-0.499999)+1
        
         write (fileind,'(I5)'), IND
         snapnm= 'data/r'//fileind((5-TENS+1):5)
         OPEN (UNIT = 1, FILE = snapnm, STATUS = 'REPLACE')
         IB=1
         DO  I=1,NP
            DO  J=1,N
               WRITE(1,*) R(IB,1),R(IB,2),R(IB,3)
               IB=IB+1
            ENDDO
         ENDDO
         CLOSE(1)
         
         snapnm= 'data/u'//fileind((5-TENS+1):5)
         OPEN (UNIT = 1, FILE = snapnm, STATUS = 'REPLACE')
         IB=1
         DO  I=1,NP
            DO  J=1,N
               WRITE(1,*) U(IB,1),U(IB,2),U(IB,3)
               IB=IB+1
            ENDDO
         ENDDO
         CLOSE(1)

         OPEN (unit=1,file='data/WINDOW',status='REPLACE')
         DO I=1,6
            WRITE(1,*) WINDOW(I)
         ENDDO
         CLOSE(1)

         OPEN (unit=1,file='data/MCAMP',status='REPLACE')
         DO I=1,6
            WRITE(1,*) MCAMP(I)
         ENDDO
         CLOSE(1)

         
         !Get dimensions of polymer. Write writhe to file for other replicas to see
         CALL getdim(N,NP,NT,R,RCOM,RCOMSQ,RGYRSQ,R2,R4,R6,DR)
        
         call WRITHE(R,N, Wr)
         
         
         RGYRSQALL(IND+1,:)=RGYRSQ
         RGYRALL(IND+1,:)=RGYRSQ**0.5
         R2ALL(IND+1,:)=R2
         R4ALL(IND+1,:)=R4
         R6ALL(IND+1,:)=R6
         DRALL(IND+1,:)=DR
         WrAll(IND+1)=Wr 

         OPEN(unit=1,file='data/Wr_Current',Status='replace')
         WRITE(unit=1,fmt=*) Wr
         CLOSE(unit=1)
       
        
         ! Get elastic energy and interaction energy and write to file 

         CALL energy_elas(EELAS,R,U,NT,N,NP,PARA,RING,TWIST,Lk,lt,LP,L)
                        
         EPONP=0.
         if (INTON.EQ.1) then
            call  ENERGY_SELF_CHAIN(EPONP,R,NT,N,NP,PARA,RING)
         endif

         ETOT=EPONP +SUM(EELAS)
         ECurrent(1:4)=EELAS
         ECurrent(5)=EPONP

         OPEN(unit=1,file='data/E_Current',Status='replace')
         WRITE(unit=1,fmt=*) ECurrent
         CLOSE(unit=1)
         
         !Construct energy vectors

         EELASALL(IND,:)=EELAS
         EPONPALL(IND)=EPONP
         ETotAll(IND)=ETOT

         ! Construct MCAMP and MCSTEP CUM vectors
         MCAMP1ALL(IND+1)=MCAMP(1)
         MCAMP2ALL(IND+1)=MCAMP(2)
         WINDOW1ALL(IND+1)=WINDOW(1)
         WINDOW2ALL(IND+1)=WINDOW(2)
         MCSTEPCUM(IND+1)=(IND)*NSTEP+NINIT
         SUCCESSALL(IND+1,:)=SUCCESS_TOTAL
         PHITALL(IND+1,:)=real(SUCCESS_TOTAL)/real(NSTEP)
        
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Parallel Tempering Algorithm
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         PRINT*, '________________________________________'
         PRINT*, 'Time point ',IND, ' out of', INDMAX


         !For simulation restart, resume here
         20 CONTINUE

         IF (ParTempOn.EQ.1) THEN         
            !Write the Current MC Index to a file for other replicas to see

            !Restart code: If MC_Completed=IND, then skip the remainder of the exchange algorithm (it was already performed)
            IF (Restart.EQ.1.AND.MC_Completed.EQ.IND) THEN
               GOTO 50
            ENDIF

            
            OPEN (unit=1,file='data/MC_IND',status='REPLACE')
            WRITE(unit=1,fmt=*) IND
            CLOSE (unit=1)


            !Check if LK replicas one step above and one step below have finished this MC step
            !If not, loop until they have finished

            print *, 'Starting synchronization'
            DO WHILE (LkMinusInd.NE.IND.OR.LkPlusInd.NE.IND) 
               OPEN(unit=1,file='../LK_'//LkMinusStr((4-TensMinus+1):4)//'/data/MC_IND',IOSTAT=IOStatus,status='old')
               READ(unit=1,fmt=*,IOSTAT=IOStatus) LkMinusInd
               CLOSE(unit=1,IOSTAT=IOStatus)

               OPEN(unit=1,file='../LK_'//LkPlusStr((4-TensPlus+1):4)//'/data/MC_IND',IOSTAT=IOStatus,status='old')
               READ(unit=1,fmt=*,IOSTAT=IOStatus) LkPlusInd
               CLOSE(unit=1,IOSTAT=IOStatus)
               
               !Wait for 1 second before checking again to avoid excessive I/O
               CALL SLEEP(1)

            ENDDO
            print *, 'Synchronization complete'
            
         
            !Once replicas are synchronized, determine whether to swap with the replica above or below
            !based on the UpOrDown varaible
            
            IF (UpOrDown.EQ.1) THEN
               LkSwap=LkPlus
               LkSwapStr=LkPlusStr
               TensSwap=TensPlus
               NTrialPlus=NTrialPlus+1
               Swap(IND)=1
            ELSE
               LkSwap=LkMinus
               LkSwapStr=LkMinusStr
               TensSwap=TensMinus
               NTrialMinus=NTrialMinus+1
               Swap(IND)=-1
            ENDIF

            !By Convention, the replica with UpOrdown=1 will decide whether the exchange is made
            !If if the current replica is not paired with any replicas during this iteration, skip
            !the test for exchange

            print *, 'LkSwap', LkSwap
            IF (Lk.EQ.LkSwap) THEN
               OPEN(unit=1,file='data/ExchangeIndex',status='replace')
               WRITE(unit=1,fmt=*) IND
               CLOSE(unit=1)

               Accept=1
               OPEN(unit=1,file='data/Accept',status='replace')
               WRITE(unit=1,fmt=*) Accept
               CLOSE(unit=1)

               GOTO 30
            ENDIF


            !Restart Code: If IND=ExchangeIndex, then the acceptance status has already
            !been determined or is ready to be read from the replica. Read it from the 
            !file and resume the simulation from the point after the simulation exchange 
            !status has been determined

            IF (Restart.EQ.1.AND.IND.EQ.ExchangeIndex) THEN
               IF (UpOrDown.EQ.1) THEN

                  OPEN(unit=1,file='data/Accept',status='old')
                  READ(unit=1,fmt=*) Accept
                  CLOSE(unit=1)
               ENDIF
               GOTO 40
            ENDIF

            !If UpOrDown==1, this replica decides whether the exchange is made

            IF (UpOrDown.EQ.1) THEN

               !Read the energy from Lkswap
               
               OPEN(unit=1,file='../LK_'//LkSwapStr((4-TensSwap+1):4)//'/data/E_Current')
               READ(unit=1,fmt=*) EReplica
               CLOSE(unit=1)


               !Get the Writhe for LkSwap

               OPEN(unit=1,file='../LK_'//LkSwapStr((4-TensSwap+1):4)//'/data/Wr_Current')
               READ(unit=1,fmt=*) WrReplica
               CLOSE(unit=1)

               !Calculate the twist energy and total energy for the replica conformation (using the current Lk)
               TwReplica=Lk-WrReplica
             
               EReplicaTotOld=sum(EReplica)
               EReplicaTotNew=EReplicaTotOld+(((2*PI*TwReplica)**2)*LT/(2*L))-EReplica(4)

               !Calculate the total energy of the current conformation under the replica's linking number
               TwNew=LkSwap-Wr
               ETotNew=ETot+(((2*PI*TwNew)**2)*LT/(2*L))-EELAS(4)
             
               !Determine whether to accept the conformation of the replica based on a Metropolois-Hastings weighting
        
               print *,'LkSwap',LkSwap
               print *, 'WrRep',WrReplica
               print *, 'WrCurrent',Wr
               Prob=exp(-((EReplicaTotNew+ETotNew)-(EReplicaTotOld+ETot)))
               print *,'Prob',Prob

               Test=grnd()
               IF (Test.LE.Prob) THEN


                  !Accept the conformation of the replica and write this to the accept file
                  !for the other replica it is paired with to see. 
                  
                  Accept=1
                  OPEN(unit=1,file='data/Accept',status='replace')
                  WRITE(unit=1,fmt=*) Accept
                  CLOSE(unit=1)

               ELSE

                  !Reject the conformation of the replica and write this to the accept file
                  !for the other replica it is paired with to see.
                  Accept=0
                  Swap(IND)=0
                  OPEN(unit=1,file='data/Accept',status='replace')
                  WRITE(unit=1,fmt=*) Accept
                  CLOSE(unit=1)

               ENDIF
            ENDIF
            !Write the number of exchange tests that have been performed by this replica for the other
            !replica to see to ensure synchronization

            OPEN(unit=1,file='data/ExchangeIndex',status='replace')
            WRITE(unit=1,fmt=*) IND
            CLOSE(unit=1)

40          CONTINUE
            
            !If UpOrDown=-1, check if the other replica accepts the exchange

            IF (UpOrDown.EQ.-1) THEN
               !Wait for the other replica to determine its acceptance status 
               !to ensure that you are not reading from an old acceptance status

               DO WHILE (ReplicaExchangeIndex.NE.IND) 
                  OPEN(unit=1,file='../LK_'//LkSwapStr((4-TensSwap+1):4)//'/data/ExchangeIndex',IOSTAT=IOStatus,status='old')
                  READ(unit=1,fmt=*,IOSTAT=IOStatus) ReplicaExchangeIndex
                  CLOSE(unit=1,IOSTAT=IOStatus)
        
                  !Wait 1 s to avoid excessive I/O
                  CALL SLEEP(1)
        
               ENDDO
       
               !Read the acceptance status of the replica

               OPEN(unit=1,file='../LK_'//LkSwapStr((4-TensSwap+1):4)//'/data/Accept',IOSTAT=IOStatus,status='old')
               READ(unit=1,fmt=*,IOSTAT=IOStatus) Accept
               CLOSE(unit=1,IOSTAT=IOStatus)
             
            ENDIF
             !If the replica exchange is accepted, update the conformation of this replica
                 !with the trial replica conformation
            IF (Accept.EQ.1) THEN

               print *, 'accept'
               OPEN(unit=1,file='../LK_'//LkSwapStr((4-TensSwap+1):4)//'/data/r'//fileind((5-TENS+1):5),status='old')
               DO I=1,NT
                  READ(unit=1,fmt=*) R(I,:)
               ENDDO
               CLOSE(unit=1)
               
               OPEN(unit=1,file='../LK_'//LkSwapStr((4-TensSwap+1):4)//'/data/u'//fileind((5-TENS+1):5),status='old')
               DO I=1,NT
                  READ(unit=1,fmt=*) U(I,:)
               ENDDO
               CLOSE(unit=1)
               
               !Sum the number of accepts
               
               IF (UpOrDown.EQ.-1) THEN
                  NSwapMinus=NSwapMinus+1
               ELSE
                  NSwapPlus=NSwapPlus+1
               ENDIF

            ELSE
               Swap(IND)=0
            ENDIF
         

30          CONTINUE
            
      

            !Alternate the UpOrDown variable so that in the next step the pairing of replicas is alternated
            UpOrDown=UpOrDown*(-1)

            OPEN (unit=1,file='data/UpOrDown',status='replace')
            WRITE(unit=1,fmt=*) UpOrDown
            CLOSE(unit=1)
 

50          CONTINUE
            !Write the number of MC steps completed for other replicas to see
            OPEN (unit=1,file='data/MC_Completed',status='replace')
            WRITE(unit=1,fmt=*) IND
            CLOSE(unit=1)

            !Wait for the other replicas to complete their simulation before moving on 

            DO WHILE (LkMinusComp.NE.IND.OR.LkPlusComp.NE.IND) 
               OPEN(unit=1,file='../LK_'//LkMinusStr((4-TensMinus+1):4)//'/data/MC_Completed',IOSTAT=IOStatus,status='old')
               READ(unit=1,fmt=*,IOSTAT=IOStatus) LkMinusComp
               CLOSE(unit=1,IOSTAT=IOStatus)

               OPEN(unit=1,file='../LK_'//LkPlusStr((4-TensPlus+1):4)//'/data/MC_Completed',IOSTAT=IOStatus,status='old')
               READ(unit=1,fmt=*,IOSTAT=IOStatus) LkPlusComp
               CLOSE(unit=1,IOSTAT=IOStatus)
               
               !Wait 1s to avoid excessive I/O
               CALL SLEEP(1)

            ENDDO
              
      
         ENDIF


         PRINT*, 'Current time ',TIME
         PRINT*, 'Bending energy ', EELAS(1)
         PRINT*, 'Par compression energy ', EELAS(2)
         PRINT*, 'Perp compression energy ', EELAS(3)
         PRINT*, 'Twist Energy', EELAS(4)
         PRINT*, 'Writhe', Wr
         PRINT*, 'Polymer-polymer energy ', EPONP
         PRINT*, 'Current number of beads ', N
         PRINT*, 'Time step ', DT
         PRINT *, 'MCAMP(1)', MCAMP(1)
         PRINT *, 'MCAMP(2)', MCAMP(2)
         PRINT *, 'WINDOW(1)', WINDOW(1)
         PRINT *, 'WINDOW(2)', WINDOW(2)
         !     Decimate the chain if time exceeds cutoff
         
!         if (TIME.GT.(500.*DT).AND.N.GT.5.AND.NSTEP.EQ.0) then
!            call decim(R0,U0,NT,N,NP,PARA,DT)
!            N=2*N-1
!            call decim(R,U,NT,N,NP,PARA,DT)
!         endif
         
         IND=IND+1
         Restart=0
      ENDDO

     ! Output the RGYRSQ and moments to file
      OPEN(UNIT = 1, FILE = 'data/RGYRSQ',STATUS='REPLACE') 
      DO IND=1,INDMAX+1
         WRITE(1,fmt=*) RGYRSQALL(IND,:)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/RGYR',STATUS='REPLACE') 
      DO IND=1,INDMAX+1
         WRITE(1,fmt=*) RGYRALL(IND,:)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/RSQ',STATUS='REPLACE') 
      DO IND=1,INDMAX+1
         WRITE(1,fmt=*) R2ALL(IND,:)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/R4th',STATUS='REPLACE') 
      DO IND=1,INDMAX+1
         WRITE(1,fmt=*) R4ALL(IND,:)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/R6th',STATUS='REPLACE') 
      DO IND=1,INDMAX+1
         WRITE(1,fmt=*) R6ALL(IND,:)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/DR',STATUS='REPLACE') 
      DO IND=1,INDMAX+1
         WRITE(1,fmt=*) DRALL(IND,:)
      ENDDO
      CLOSE(1)


    ! Output structure statistics to file (average, standard error)

      RGYRSQ_AVG=SUM(RGYRSQALL(2:INDMAX+1,:))/REAL(INDMAX)
      RGYRSQ_STDEV=SQRT(SUM((RGYRSQALL(2:INDMAX+1,:)-RGYRSQ_AVG)**2))/SQRT(REAL(INDMAX))
      RGYRSQ_STDER=SQRT(SUM((RGYRSQALL(2:INDMAX+1,:)-RGYRSQ_AVG)**2))/REAL(INDMAX)

      RGYR_AVG=SUM(RGYRALL(2:INDMAX+1,:))/REAL(INDMAX)
      RGYR_STDEV=SQRT(SUM((RGYRALL(2:INDMAX+1,:)-RGYR_AVG)**2))/SQRT(REAL(INDMAX))
      RGYR_STDER=SQRT(SUM((RGYRALL(2:INDMAX+1,:)-RGYR_AVG)**2))/REAL(INDMAX)


      R2_AVG=SUM(R2ALL(2:INDMAX+1,:))/REAL(INDMAX)
      R2_STDEV=SQRT(SUM((R2ALL(2:INDMAX+1,:)-R2_AVG)**2))/SQRT(REAL(INDMAX))
      R2_STDER=SQRT(SUM((R2ALL(2:INDMAX+1,:)-R2_AVG)**2))/REAL(INDMAX)

      R4_AVG=SUM(R4ALL(2:INDMAX+1,:))/REAL(INDMAX)
      R4_STDEV=SQRT(SUM((R4ALL(2:INDMAX+1,:)-R4_AVG)**2))/SQRT(REAL(INDMAX))
      R4_STDER=SQRT(SUM((R4ALL(2:INDMAX+1,:)-R4_AVG)**2))/REAL(INDMAX)

      R6_AVG=SUM(R6ALL(2:INDMAX+1,:))/REAL(INDMAX)
      R6_STDEV=SQRT(SUM((R6ALL(2:INDMAX+1,:)-R6_AVG)**2))/SQRT(REAL(INDMAX))
      R6_STDER=SQRT(SUM((R6ALL(2:INDMAX+1,:)-R6_AVG)**2))/REAL(INDMAX)

      DR_AVG=SUM(DRALL(2:INDMAX+1,:))/REAL(INDMAX)
      DR_STDEV=SQRT(SUM((DRALL(2:INDMAX+1,:)-DR_AVG)**2))/SQRT(REAL(INDMAX))
      DR_STDER=SQRT(SUM((DRALL(2:INDMAX+1,:)-DR_AVG)**2))/REAL(INDMAX)


      Wr_AVG=SUM(WrALL(2:INDMAX+1))/REAL(INDMAX)
      Wr_STDEV=SQRT(SUM((WrAll-Wr_AVG)**2))/SQRT(REAL(INDMAX))
      Wr_STDER=SQRT(SUM((WrAll-Wr_AVG)**2))/REAL(INDMAX)

      ETotAvg=SUM(ETotAll)/REAL(INDMAX)
      ETotStdev=SQRT(SUM((ETotAll-ETotAvg)**2))/SQRT(REAL(INDMAX))
      ETotStderr=SQRT(SUM((EtotAll-ETotAvg)**2))/REAL(INDMAX)
    
      OPEN(UNIT = 1, FILE = 'data/RGYRSQ_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) RGYRSQ_AVG, RGYRSQ_STDEV, RGYRSQ_STDER
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/RGYR_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) RGYR_AVG, RGYR_STDEV, RGYR_STDER
      CLOSE(1)

      
      OPEN(UNIT = 1, FILE = 'data/R2_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) R2_AVG, R2_STDEV, R2_STDER
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/R4_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) R4_AVG, R4_STDEV, R4_STDER
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/R6_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) R6_AVG, R6_STDEV, R6_STDER
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/DR_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) DR_AVG, DR_STDEV, DR_STDER
      CLOSE(1)


      OPEN(UNIT = 1, FILE = 'data/Wr_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) Wr_AVG,Wr_STDEV, Wr_STDER
      CLOSE(1)

  
      OPEN (UNIT = 1, FILE = 'data/Wr', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) WrAll(IND)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 1, FILE = 'data/Energy_STAT',STATUS='REPLACE') 
      WRITE(1,fmt=*) ETotAvg, ETotStdev, ETotStderr
      CLOSE(1)


      OPEN (UNIT = 1, FILE = 'data/MCAMP1', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), MCAMP1ALL(IND)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/MCAMP2', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), MCAMP2ALL(IND)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/WrPlot', STATUS='REPLACE')
      DO IND=1,INDMAX
         WRITE(1, fmt=*) MCSTEPCUM(IND), WrALL(IND)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/window1', STATUS='REPLACE')
      DO IND=1,INDMAX
         WRITE(1, fmt=*) MCSTEPCUM(IND), WINDOW1ALL(IND)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/window2', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), WINDOW2ALL(IND)
      ENDDO
      CLOSE(1)

   
      OPEN (UNIT = 1, FILE = 'data/EELAS', STATUS='REPLACE')
      DO IND=1,INDMAX
         WRITE(1, fmt=*) EELASALL(IND,:)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/ETot', STATUS='REPLACE')
      DO IND=1,INDMAX
         WRITE(1, fmt=*) ETotAll(IND)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/EPONP', STATUS='REPLACE')
      DO IND=1,INDMAX
         WRITE(1, fmt=*) MCSTEPCUM(IND), EPONPALL(IND)
      ENDDO
      CLOSE(1)


      OPEN (UNIT = 1, FILE = 'data/PHIT1', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), PHITALL(IND,1)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/PHIT2', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), PHITALL(IND,2)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/PHIT3', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), PHITALL(IND,3)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/PHIT4', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), PHITALL(IND,4)
      ENDDO
      CLOSE(1)

      OPEN (UNIT = 1, FILE = 'data/PHIT5', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), PHITALL(IND,5)
      ENDDO
      CLOSE(1)
      
      OPEN (UNIT = 1, FILE = 'data/PHIT6', STATUS='REPLACE')
      DO IND=1,INDMAX+1
         WRITE(1, fmt=*) MCSTEPCUM(IND), PHITALL(IND,6)
      ENDDO
      CLOSE(1)


     ! Algorithm Analysis: Get auto_correlations and write to file

      CALL auto_correlation_vector(R2ALL(2:INDMAX+1,1),INDMAX,1,N_auto,RSQ_AUTO)
      OPEN (UNIT =1, FILE='data/RSQ_auto',STATUS='REPLACE')

      DO IND=1,N_auto
         WRITE(1,fmt=*) RSQ_AUTO(IND)
      ENDDO


      CALL auto_correlation_vector(RGYRSQALL(2:INDMAX+1,1),INDMAX,1,N_auto,RGYRSQ_AUTO)
      OPEN (UNIT =1, FILE='data/RGYRSQ_auto',STATUS='REPLACE')

      DO IND=1,N_auto
         WRITE(1,fmt=*) RGYRSQ_AUTO(IND)
      ENDDO


      CALL auto_correlation_vector(WrALL(2:INDMAX+1),INDMAX,1,N_auto,Wr_AUTO)
      OPEN (UNIT =1, FILE='data/Wr_auto',STATUS='REPLACE')

      DO IND=1,N_auto
         WRITE(1,fmt=*) Wr_AUTO(IND)
      ENDDO

      CALL auto_correlation_vector(ETotAll,INDMAX,1,N_auto,Energy_AUTO)
      OPEN (UNIT =1, FILE='data/Energy_auto',STATUS='REPLACE')
           DO IND=1,N_auto
         WRITE(1,fmt=*) Energy_AUTO(IND)
      ENDDO


      !Parallel tempering analysis: Save Fraction of Exchanges to file

      OPEN (UNIT=1,file='data/PSwapPlus')
      WRITE(unit=1,fmt=*) real(NSwapPlus)/real(NTrialPlus)
      CLOSE(UNIT=1)

      OPEN (UNIT=1,file='data/PSwapMinus')
      WRITE(unit=1,fmt=*) real(NSwapMinus)/real(NtrialMinus)
      CLOSE(unit=1)
        
      OPEN (UNIT=1,file='data/Swap')
      DO IND=1,INDMAX

         WRITE(unit=1,fmt=*) Swap(IND) 
      ENDDO
      CLOSE(unit=1)

      ENDPROGRAM
      
    

     

      
!---------------------------------------------------------------*
      
