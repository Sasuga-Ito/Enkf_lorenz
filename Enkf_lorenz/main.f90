PROGRAM LORENZ
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP     = KIND(0.0D0)
 INTEGER, PARAMETER :: NSKIND = DP

 INTEGER                       :: I, M, L1, L2, L3

 REAL(NSKIND)                  :: DT, A, B
 REAL(NSKIND), DIMENSION(3)    :: K1, K2, K3, K4
 REAL(NSKIND), DIMENSION(3)    :: R0, W0, X0
 REAL(NSKIND), DIMENSION(3,3000) :: XF, DXF
 REAL(NSKIND), DIMENSION(3,3000) :: DXO
 REAL(NSKIND), DIMENSION(3,3000) :: W

 REAL(NSKIND), DIMENSION(3)    :: E
 REAL(NSKIND), DIMENSION(3)    :: P0, XAH, DXAH
 REAL(NSKIND), DIMENSION(3,3)  :: PX, PDX, K
 REAL(NSKIND), DIMENSION(3,6)  :: PI
 REAL(NSKIND), DIMENSION(3,100):: XA, DXA

  DT = 0.01
  R0 = (/0.5D0, 0.2D0, 0.1D0/)
  W0 = (/2.0D0, 4.0D0, 3.0D0/)

  XF(1,1) = 5.0D0
  XF(2,1) = 3.0D0
  XF(3,1) = 3.0D0
 
  DXF(:,1) = XF(:,1)

  CALL RANDOM_NUMBER(W(:,1))
  W(:,1) = W0(:)*W(:,1)
  
  DO I = 2, 3000
     K1(:) = CALDT(XF(:,I-1))
     K2(:) = CALDT(XF(:,I-1)+DT*K1*0.5D0)
     K3(:) = CALDT(XF(:,I-1)+DT*K2*0.5D0)
     K4(:) = CALDT(XF(:,I-1)+DT*K3)

     CALL RANDOM_NUMBER(W(:,I))
     W(:,I) = W0(:)*W(:,I)

     XF(:,I) = XF(:,I-1) + DT*(K1(:)+K2(:)+K3(:)+K4(:))/6.0D0 + W(:,I-1)*DT
     DXF(:,I) = XF(:,I)
  ENDDO

  DO I = 1, 3000
     CALL RANDOM_NUMBER(DXO(:,I))
     DXO(:,I) = R0(:)*DXO(:,I)
     DXO(:,I) = DXO(:,I) + DXF(:,I)
  ENDDO

  OPEN(1,FILE='observation.dat',FORM='FORMATTED')
  DO I = 1,3000
  WRITE(1,'(6f14.4)') XF(1,I), XF(2,I), XF(3,I),DXO(1,I),DXO(2,I),DXO(3,I)
  ENDDO
  CLOSE(1)
  
  X0 = (/2.0D0,10.0D0,1.0D0/)
  P0 = (/1.0D0,3.0D0,3.0D0/)


  DO M = 1, 100
     CALL RANDOM_NUMBER(XA(:,M))
     XA(:,M) = W0(:)*XA(:,M)
     XA(:,M) = X0(:) + XA(:,M)
  ENDDO

  DO M = 1, 100
     CALL RANDOM_NUMBER(DXA(:,M))
     DXA(:,M) = R0(:)*DXA(:,M)
     DXA(:,M) = XA(:,M) + DXA(:,M)
  ENDDO

    
  !OPEN(3,FILE='forecast.dat',FORM='FORMATTED')
  OPEN(3,FILE='assimulation.dat',FORM='FORMATTED')
  DO I = 1, 3000

  XAH = 0.0D0; DXAH = 0.0D0
  DO M = 1, 100
     XAH(:)  = XAH(:) + XA(:,M)
     DXAH(:)  = DXAH(:) + DXA(:,M)
  ENDDO
  XAH = XAH/100; DXAH = DXAH/100
  
  WRITE(3,'(6f14.4)') XAH(1), XAH(2), XAH(3),DXAH(1),DXAH(2),DXAH(3)

  PX = 0.0D0; PDX = 0.0D0
  DO M = 1, 100
  DO L1 = 1,3; DO L2 = 1,3
     PX(L1,L2)  = PX(L1,L2)  + (XA(L1,M)-XAH(L1))*(DXA(L2,M)-DXAH(L2))/DBLE(100-1)
     PDX(L1,L2) = PDX(L1,L2) + (DXA(L1,M)-DXAH(L1))*(DXA(L2,M)-DXAH(L2))/DBLE(100-1)
  ENDDO; ENDDO
  ENDDO

  
  PI = 0.0D0
  PI(1:3,1:3) = PDX(1:3,1:3)
  PI(1,4) = 1.0D0; PI(2,5) = 1.0D0; PI(3,6) = 1.0D0


  DO L1 = 1, 3
     A = PI(L1,L1)
     PI(L1,:) = PI(L1,:)/A
     DO L2 = 1, 3
        IF(L2 /= L1)THEN
           B = PI(L2,L1)
           PI(L2,:) = PI(L2,:) - B*PI(L1,:)
        ENDIF
     ENDDO
  ENDDO

  PDX(1:3,1:3) = PI(1:3,4:6)

  K = 0.0D0
  DO L1 = 1,3; DO L2 = 1,3; DO L3 = 1,3
     K(L1,L2) = K(L1,L2) + PX(L1,L3)*PDX(L3,L2)
  ENDDO; ENDDO; ENDDO

  DO M = 1, 100
  DO L1 = 1,3; DO L2 = 1,3
     XA(L1,M) = XA(L1,M) + K(L1,L2)*(DXO(L2,I)-DXA(L2,M))
  ENDDO; ENDDO
     CALL RANDOM_NUMBER(E)
     XA(:,M) = XA(:,M) + CALDT(XA(:,M))*DT+W0(:)*E(:)
     CALL RANDOM_NUMBER(E)
     DXA(:,M) = XA(:,M) + R0(:)*E(:)
  ENDDO
  ENDDO

  CLOSE(3)

  

 CONTAINS
 
  FUNCTION CALDT(X) RESULT(DX)
  IMPLICIT NONE
  REAL(NSKIND)               :: S = 10.0D0
  REAL(NSKIND)               :: R = 28.0D0
  REAL(NSKIND)               :: B = 8.0D0/3.0D0
  REAL(NSKIND), DIMENSION(3) :: X, DX

   DX(1) = -S*(X(1)-X(2))
   DX(2) = R*X(1)-X(2)-X(1)*X(3)
   DX(3) = X(1)*X(2)-B*X(3)

  END FUNCTION CALDT
END PROGRAM
