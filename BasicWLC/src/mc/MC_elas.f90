!---------------------------------------------------------------*

!     subroutine MC_eelas
!
!     Calculate the change in the polymer elastic energy
!     due to the displacement from a MC move
SUBROUTINE MC_eelas(DEELAS,R,U,RP,UP,&
                    NT,NB,IB1,IB2,&
                    IT1,IT2,EB,EPAR,EPERP,GAM,ETA, &
                    RING,TWIST,Lk,lt,LP,L,MCTYPE,WR,WRP)

use setPrecision
IMPLICIT NONE
DOUBLE PRECISION, intent(in) :: R(NT,3)  ! Bead positions
DOUBLE PRECISION, intent(in) :: U(NT,3)  ! Tangent vectors
DOUBLE PRECISION, intent(in) :: RP(NT,3)  ! Bead positions
DOUBLE PRECISION, intent(in) :: UP(NT,3)  ! Tangent vectors
INTEGER, intent(in) :: NB                ! Number of beads in a polymer
INTEGER, intent(in) :: NT                ! Total number of beads
INTEGER, intent(in) :: IB1               ! Test bead position 1
INTEGER, intent(in) :: IT1               ! Index of test bead 1
INTEGER, intent(in) :: IB2               ! Test bead position 2
INTEGER, intent(in) :: IT2               ! Index of test bead 2

DOUBLE PRECISION, intent(out) :: DEELAS(4)   ! Change in ECOM

!     Polymer properties

double precision, intent(in) :: EB
double precision, intent(in) :: EPAR
double precision, intent(in) :: EPERP
double precision, intent(in) :: GAM
double precision, intent(in) :: ETA
INTEGER  Lk               ! Linking number
DOUBLE PRECISION L        ! Contour length
DOUBLE PRECISION lt       ! Twist persistence length
DOUBLE PRECISION LP       ! Bend persistence length
DOUBLE PRECISION Tw       ! Twist
DOUBLE PRECISION TWP      ! Twist of test structure
DOUBLE PRECISION WR       ! Writhe
DOUBLE PRECISION WRP      ! Writhe of test structure
DOUBLE PRECISION DWR      ! Change in Writhe
DOUBLE PRECISION WRM,WRMP ! Component of writhe affected by move
INTEGER RING              ! Is polymer a ring?
INTEGER TWIST             ! Include twist?
INTEGER IT1M1
INTEGER IT2P1
INTEGER MCTYPE            ! MC move type

DOUBLE PRECISION DR(3),DRPAR,DRPERP(3)
DOUBLE PRECISION GI(3)

! Setup parameters

      DEELAS(1) = 0.0_dp ! bending energy
      DEELAS(2) = 0.0_dp ! parallel stretch energy
      DEELAS(3) = 0.0_dp ! perpendicular stretch energy
      DEELAS(4) = 0.0_dp ! twist energy

!     Calculate the change in the energy

      if ((IB1.NE.1).or.(RING.EQ.1)) then
          ! so if RING.EQ.1, i.e.
          if (IB1.EQ.1) then
              ! then the bead to the left of IB1 is actually the end bead due to the ring inducing periodic boundaries on the array
              IT1M1 = NB
          else
              IT1M1 = IT1 - 1
          endif

         ! MC move only affects energy if it's interior to the polymer, since there are only crankshaft moves, and end
         ! crankshafts don't affect polymer
         if (SIMTYPE.EQ.1.AND.(IB1.NE.NB.OR.RING.EQ.1)) then
             if (IB1.EQ.NB) then
                 IT1P1 = 1
             else
                 IT1P1 = IT1 + 1
             endif

            U(IT1M1,1)=R(IT1,1)-R(IT1M1,1)
            U(IT1M1,2)=R(IT1,2)-R(IT1M1,2)
            U(IT1M1,3)=R(IT1,3)-R(IT1M1,3)
            UNORM=sqrt(U(IT1M1,1)**2.+U(IT1M1,2)**2.+U(IT1M1,3)**2.)
            U(IT1M1,1)=U(IT1M1,1)/UNORM
            U(IT1M1,2)=U(IT1M1,2)/UNORM
            U(IT1M1,3)=U(IT1M1,3)/UNORM

            U(IT1,1)=R(IT1+1,1)-R(IT1,1)
            U(IT1,2)=R(IT1+1,2)-R(IT1,2)
            U(IT1,3)=R(IT1+1,3)-R(IT1,3)
            UNORM=sqrt(U(IT1,1)**2.+U(IT1,2)**2.+U(IT1,3)**2.)
            U(IT1,1)=U(IT1,1)/UNORM
            U(IT1,2)=U(IT1,2)/UNORM
            U(IT1,3)=U(IT1,3)/UNORM

            UP(IT1,1)=RP(IT1+1,1)-RP(IT1,1)
            UP(IT1,2)=RP(IT1+1,2)-RP(IT1,2)
            UP(IT1,3)=RP(IT1+1,3)-RP(IT1,3)
            UNORM=sqrt(UP(IT1,1)**2.+UP(IT1,2)**2.+UP(IT1,3)**2.)
            UP(IT1,1)=UP(IT1,1)/UNORM
            UP(IT1,2)=UP(IT1,2)/UNORM
            UP(IT1,3)=UP(IT1,3)/UNORM

            U1U2=U(IT1M1,1)*U(IT1,1)+U(IT1M1,2)*U(IT1,2)+U(IT1M1,3)*U(IT1,3)
            ! only bending energy in WLC, others are init to zero, so sum in total energy calculation later will be no-op
            DEELAS(1)=DEELAS(1)+EB*U1U2
            U1U2=U(IT1M1,1)*UP(IT1,1)+U(IT1M1,2)*UP(IT1,2)+U(IT1M1,3)*UP(IT1,3)
            DEELAS(1)=DEELAS(1)-EB*U1U2

         elseif (SIMTYPE.EQ.2) then

            DR(1)=R(IT1,1)-R(IT1M1,1)
            DR(2)=R(IT1,2)-R(IT1M1,2)
            DR(3)=R(IT1,3)-R(IT1M1,3)
            DRPAR=DR(1)*U(IT1M1,1)+DR(2)*U(IT1M1,2)+DR(3)*U(IT1M1,3)

            DRPERP(1)=DR(1)-DRPAR*U(IT1M1,1)
            DRPERP(2)=DR(2)-DRPAR*U(IT1M1,2)
            DRPERP(3)=DR(3)-DRPAR*U(IT1M1,3)
            !U1U2=U(IT1M1,1)*U(IT1,1)+U(IT1M1,2)*U(IT1,2)+U(IT1M1,3)*U(IT1,3)

            GI(1)=(U(IT1,1)-U(IT1M1,1)-ETA*DRPERP(1))
            GI(2)=(U(IT1,2)-U(IT1M1,2)-ETA*DRPERP(2))
            GI(3)=(U(IT1,3)-U(IT1M1,3)-ETA*DRPERP(3))

            DEELAS(1)=DEELAS(1)-0.5_dp*EB*(GI(1)**2+GI(2)**2+GI(3)**2)
            DEELAS(2)=DEELAS(2)-0.5_dp*EPAR*(DRPAR-GAM)**2
            DEELAS(3)=DEELAS(3)-0.5_dp*EPERP*(DRPERP(1)**2+DRPERP(2)**2.+DRPERP(3)**2)

            DR(1)=RP(IT1,1)-R(IT1M1,1)
            DR(2)=RP(IT1,2)-R(IT1M1,2)
            DR(3)=RP(IT1,3)-R(IT1M1,3)
            DRPAR=DR(1)*U(IT1M1,1)+DR(2)*U(IT1M1,2)+DR(3)*U(IT1M1,3)

            DRPERP(1)=DR(1)-DRPAR*U(IT1M1,1)
            DRPERP(2)=DR(2)-DRPAR*U(IT1M1,2)
            DRPERP(3)=DR(3)-DRPAR*U(IT1M1,3)
            !U1U2=U(IT1M1,1)*UP(IT1,1)+U(IT1M1,2)*UP(IT1,2)+U(IT1M1,3)*UP(IT1,3)

            GI(1)=(UP(IT1,1)-U(IT1M1,1)-ETA*DRPERP(1))
            GI(2)=(UP(IT1,2)-U(IT1M1,2)-ETA*DRPERP(2))
            GI(3)=(UP(IT1,3)-U(IT1M1,3)-ETA*DRPERP(3))

            ! in ssWLC, we have, bending, parallel extension, and perpendicular extension energy, which we must sum later to get
            ! total energy
            DEELAS(1)=DEELAS(1)+0.5_dp*EB*(GI(1)**2+GI(2)**2+GI(3)**2)
            DEELAS(2)=DEELAS(2)+0.5_dp*EPAR*(DRPAR-GAM)**2
            DEELAS(3)=DEELaS(3)+0.5_dp*EPERP*(DRPERP(1)**2+DRPERP(2)**2+DRPERP(3)**2)

         elseif (SIMTYPE.EQ.3) then

            DR(1)=R(IT1,1)-R(IT1M1,1)
            DR(2)=R(IT1,2)-R(IT1M1,2)
            DR(3)=R(IT1,3)-R(IT1M1,3)
            ! in gaussian chain, there's only parallel stretching energy. DEELAS init'd to zeros, so sum(DEELAS) == DEELAS(2) later
            DEELAS(2)=DEELAS(2)-0.5*EPAR*(DR(1)**2.+DR(2)**2.+DR(3)**2.)
            DR(1)=RP(IT1,1)-R(IT1M1,1)
            DR(2)=RP(IT1,2)-R(IT1M1,2)
            DR(3)=RP(IT1,3)-R(IT1M1,3)
            DEELAS(2)=DEELAS(2)+0.5*EPAR*(DR(1)**2.+DR(2)**2.+DR(3)**2.)

         endif
      endif

      if ((IB2.NE.NB).or.(RING.EQ.1)) then
         if (IB2.EQ.NB) then
             IT2P1 = 1
         else
             IT2P1 = IT2 + 1
         endif

         ! if we're talking about a WLC, if we crankshaft a single bead, that's a no-op, since the u's are directly
         ! determined by the r's. Thus we're not worried about double counting the energy change here since the energy change
         ! should be zero by definition if IB1==IB2.
         if (SIMTYPE.EQ.1.AND.((IB2.NE.1).OR.(RING.EQ.1))) then
            if (IB2.EQ.1) then
                IT2M1 = NB
            else
                IT2M1 = IT2 - 1
            endif

            U(IT2-1,1)=R(IT2,1)-R(IT2-1,1)
            U(IT2-1,2)=R(IT2,2)-R(IT2-1,2)
            U(IT2-1,3)=R(IT2,3)-R(IT2-1,3)
            UNORM=sqrt(U(IT2-1,1)**2.+U(IT2-1,2)**2.+U(IT2-1,3)**2.)
            U(IT2-1,1)=U(IT2-1,1)/UNORM
            U(IT2-1,2)=U(IT2-1,2)/UNORM
            U(IT2-1,3)=U(IT2-1,3)/UNORM

            U(IT2,1)=R(IT2+1,1)-R(IT2,1)
            U(IT2,2)=R(IT2+1,2)-R(IT2,2)
            U(IT2,3)=R(IT2+1,3)-R(IT2,3)
            UNORM=sqrt(U(IT2,1)**2.+U(IT2,2)**2.+U(IT2,3)**2.)
            U(IT2,1)=U(IT2,1)/UNORM
            U(IT2,2)=U(IT2,2)/UNORM
            U(IT2,3)=U(IT2,3)/UNORM

            UP(IT2-1,1)=RP(IT2,1)-RP(IT2-1,1)
            UP(IT2-1,2)=RP(IT2,2)-RP(IT2-1,2)
            UP(IT2-1,3)=RP(IT2,3)-RP(IT2-1,3)
            UNORM=sqrt(UP(IT2-1,1)**2.+UP(IT2-1,2)**2.+UP(IT2-1,3)**2.)
            UP(IT2-1,1)=UP(IT2-1,1)/UNORM
            UP(IT2-1,2)=UP(IT2-1,2)/UNORM
            UP(IT2-1,3)=UP(IT2-1,3)/UNORM

            U1U2=U(IT2-1,1)*U(IT2,1)+U(IT2-1,2)*U(IT2,2)+U(IT2-1,3)*U(IT2,3)
            DEELAS(1)=DEELAS(1)+EB*U1U2
            U1U2=UP(IT2-1,1)*U(IT2,1)+UP(IT2-1,2)*U(IT2,2)+UP(IT2-1,3)*U(IT2,3)
            DEELAS(1)=DEELAS(1)-EB*U1U2

         elseif (SIMTYPE.EQ.2) then

            DR(1)=R(IT2+1,1)-R(IT2,1)
            DR(2)=R(IT2+1,2)-R(IT2,2)
            DR(3)=R(IT2+1,3)-R(IT2,3)
            DRPAR=DR(1)*U(IT2,1)+DR(2)*U(IT2,2)+DR(3)*U(IT2,3)

            DRPERP(1)=DR(1)-DRPAR*U(IT2,1)
            DRPERP(2)=DR(2)-DRPAR*U(IT2,2)
            DRPERP(3)=DR(3)-DRPAR*U(IT2,3)
            !U1U2=U(IT2,1)*U(IT2+1,1)+U(IT2,2)*U(IT2+1,2)+U(IT2,3)*U(IT2+1,3)

            GI(1)=(U(IT2+1,1)-U(IT2,1)-ETA*DRPERP(1))
            GI(2)=(U(IT2+1,2)-U(IT2,2)-ETA*DRPERP(2))
            GI(3)=(U(IT2+1,3)-U(IT2,3)-ETA*DRPERP(3))

            DEELAS(1)=DEELAS(1)-0.5_dp*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.)
            DEELAS(2)=DEELAS(2)-0.5_dp*EPAR*(DRPAR-GAM)**2.
            DEELAS(3)=DEELAS(3)-0.5_dp*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

            DR(1)=R(IT2+1,1)-RP(IT2,1)
            DR(2)=R(IT2+1,2)-RP(IT2,2)
            DR(3)=R(IT2+1,3)-RP(IT2,3)
            DRPAR=DR(1)*UP(IT2,1)+DR(2)*UP(IT2,2)+DR(3)*UP(IT2,3)

            DRPERP(1)=DR(1)-DRPAR*UP(IT2,1)
            DRPERP(2)=DR(2)-DRPAR*UP(IT2,2)
            DRPERP(3)=DR(3)-DRPAR*UP(IT2,3)
            !U1U2=UP(IT2,1)*U(IT2+1,1)+UP(IT2,2)*U(IT2+1,2)+UP(IT2,3)*U(IT2+1,3)

            GI(1)=(U(IT2+1,1)-UP(IT2,1)-ETA*DRPERP(1))
            GI(2)=(U(IT2+1,2)-UP(IT2,2)-ETA*DRPERP(2))
            GI(3)=(U(IT2+1,3)-UP(IT2,3)-ETA*DRPERP(3))

            DEELAS(1)=DEELAS(1)+0.5_dp*EB*(GI(1)**2.+GI(2)**2+GI(3)**2)
            DEELAS(2)=DEELAS(2)+0.5_dp*EPAR*(DRPAR-GAM)**2.
            DEELAS(3)=DEELAS(3)+0.5_dp*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

         elseif (SIMTYPE.EQ.3) then

            DR(1)=R(IT2+1,1)-R(IT2,1)
            DR(2)=R(IT2+1,2)-R(IT2,2)
            DR(3)=R(IT2+1,3)-R(IT2,3)
            DEELAS(2)=DEELAS(2)-0.5*EPAR*(DR(1)**2.+DR(2)**2.+DR(3)**2.)
            DR(1)=R(IT2+1,1)-RP(IT2,1)
            DR(2)=R(IT2+1,2)-RP(IT2,2)
            DR(3)=R(IT2+1,3)-RP(IT2,3)
            DEELAS(2)=DEELAS(2)+0.5*EPAR*(DR(1)**2.+DR(2)**2.+DR(3)**2.)

         endif

      endif

      IF (RING.EQ.1.AND.TWIST.EQ.1) THEN
          IF (MCTYPE.EQ.1) THEN
              CALL WRITHECRANK(R,IT1,IT2,N,WRM)
              CALL WRITHECRANK(RP,IT1,IT2,N,WRMP)
              DWR=WRMP-WRM
          ELSEIF (MCTYPE.EQ.2) THEN
              CALL WRITHESLIDE(R,IT1,IT2,N,WRM)
              CALL WRITHESLIDE(RP,IT1,IT2,N,WRMP)
              DWR=WRMP-WRM
          ELSE
              DWR=0.
          ENDIF
          WRP=WR+DWR
          TW=REAL(LK)-WR
          TWP=REAL(LK)-WRP
          DEELAS(4)=DEELAS(4)+(((2.*pi*TWP)**2.)*LT/(2.*L))-(((2.*pi*TW)**2.)*LT/(2.*L))
      ENDIF

      RETURN
      END

!---------------------------------------------------------------*
