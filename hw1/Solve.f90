! SYSTEM FORTRAN-90
! LAST MODIFICATION BY Matt McGarry, 13 jan2012 - Changed to F90 format
!------------------------------------------------------------
! B = LHS MATRIX , DIMENSIONED (NDIM,MDIM) IN CALLING PROGRAM
! R = RIGHT-HAND-SIDE, DIMENSIONED (NDIM) IN CALLING PROGRAM 
! NEQ= # OF EQUATIONS                                        
! IHALFB= HALF-BANDWIDTH; 2*IHALFB+1 = BANDWIDTH             
! SOLUTION RETURNS IN R                                      
! KKK = 1 PERFORMS LU DECOMPOSITION DESTRUCTIVELY            
! KKK = 2 PERFORMS BACK SUBSTITUTION
! KKK = 3 PERFORMS OPTIONS 1 AND 2
!
      SUBROUTINE SOLVE(KKK,B,R,NEQ,IHALFB,NDIM,MDIM)                           
!                                                                               
!   ASYMMETRIC BAND MATRIX EQUATION SOLVER (SIGNLE PRECISION)                                      
!   DOCTORED TO IGNORE ZEROS IN LU DECOMP. STEP                                 
!                                                                               
        INTEGER KKK,NEQ,IHALFB,NDIM,MDIM
        REAL*4 B(NDIM,MDIM),R(NDIM)                                            
        NRS=NEQ-1                                                                 
        IHBP=IHALFB+1                                                             
        IF (KKK.EQ.2) GO TO 30                                                    
!                                                                               
!  TRIANGULARIZE MATRIX A USING DOOLITTLE METHOD                                
!                                                                               
      DO 10 K=1,NRS                                                             
      PIVOT=B(K,IHBP)                                                           
      KK=K+1                                                                    
      KC=IHBP                                                                   
      DO 21 I=KK,NEQ                                                            
      KC=KC-1                                                                   
      IF(KC.LE.0) GO TO 10                                                      
      C=-B(I,KC)/PIVOT                                                          
      IF (C.EQ.0.0) GO TO 21                                                    
      B(I,KC)=C                                                                 
      KI=KC+1                                                                   
      LIM=KC+IHALFB                                                             
      DO 20 J=KI,LIM                                                            
      JC=IHBP+J-KC                                                              
   20 B(I,J)=B(I,J)+C*B(K,JC)                                                   
   21 CONTINUE                                                                  
   10 CONTINUE                                                                  
      IF(KKK.EQ.1) GO TO 100
!                                                                               
!  MODIFY LOAD VECTOR R                                                         
!                                                                               
   30 NN=NEQ+1                                                                  
      IBAND=2*IHALFB+1                                                          
      DO 40 I=2,NEQ                                                             
      JC=IHBP-I+1                                                               
      JI=1                                                                      
      IF (JC.LE.0) GO TO 50
      GO TO 60
   50 JC=1
      JI=I-IHBP+1
   60 SUM=0.0
      DO 70 J=JC,IHALFB
      SUM=SUM+B(I,J)*R(JI)
   70 JI=JI+1
   40 R(I)=R(I)+SUM
!                                                                               
!   BACK SOLUTION                                                               
!                                                                               
      R(NEQ)=R(NEQ)/B(NEQ,IHBP)
      DO 80 IBACK=2,NEQ
      I=NN-IBACK
      JP=I
      KR=IHBP+1
      MR=MIN0(IBAND,IHALFB+IBACK)
      SUM=0.0
      DO 90 J=KR,MR
      JP=JP+1
   90 SUM=SUM+B(I,J)*R(JP)
   80 R(I)=(R(I)-SUM)/B(I,IHBP)
  100 RETURN
      END SUBROUTINE