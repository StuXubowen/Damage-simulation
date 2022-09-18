*deck,usermat      USERDISTRIB  parallel                       gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,statev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c******************************************************************

#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), statev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c


      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,EMOD,ENU,EA,EB,EC
	DOUBLE PRECISION Etf,Etu,Ecf,Ecu
	double precision  Damage

	INTEGER      i,j,K1,K2
c
c************************ User defined part*****************************
c

c-----------------------Calculated element stiffness coefficient---------------------------

	PARAMETER (ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0)
c-------------EMOD(Elastic modulus)£¬ENU(Poisson's ratio)--------------------------
c-------------Etf(Tensile damage threshold strain)£¬Etu(Tensile damage limit strain)-----------
c-------------Ecf(Compression damage threshold strain)£¬Ecu(Compression damage limit strain)-----------
	EMOD=prop(1)
	ENU=prop(2)
	Etf=prop(3)
	Etu=prop(4)
	Ecf=prop(5)
	Ecu=prop(6)
c     
c     Call getDamage to calculate the damage factor
c
      call getDamage(Strain,ncomp,Etf,Etu,Ecf,Ecu,Damage)
	EA=EMOD*ENU
	EA=EA/(ONE+ENU)/(ONE-TWO*ENU)
	EA=EA+(EMOD*Damage)/(THREE*(ONE+ENU))
	EB=EMOD*(ONE-ENU)/(ONE+ENU)/(ONE-TWO*ENU)
	EB=EB-(EMOD*TWO*Damage)/(THREE*(ONE+ENU))
      EC=EMOD/(ONE+ENU)/TWO*(ONE-Damage)
c
C     Uniform tangent operator matrix
C
	DO K1=1, nDirect
	  DO K2=1, nDirect
	    dsdePl(K2, K1)=EA
	  END DO
	  dsdePl(K1, K1)=EB
	END DO
	DO K1=nDirect+1, ncomp
	  dsdePl(K1 ,K1)=EC
	END DO

C------------------------------------------------------------------
C--------------------------Update stress--------------------------------
C------------------------------------------------------------------
	DO K1=1, ncomp
	DO K2=1, ncomp
	stress(K2)=stress(K2)+dsdePl(K2,K1)*dStrain(K1)
	END DO
	END DO

C------------------------------------------------------------------
C------------------------Update state variable damage factor----------------------
C------------------------------------------------------------------
      statev(1)=Damage
      RETURN
      END

c------------------------------------------------------------------
c------------------Get damage factor D subprogram getDamage--------------------
c------------------------------------------------------------------ 
      SUBROUTINE getDamage(Strain,ncomp,Etf,Etu,Ecf,Ecu,Damage)
#include "impcom.inc"
	INTEGER ncomp
	double precision  Strain(ncomp),Etf,Etu,Ecf,Ecu
	double precision  equalStrain
      double precision  pStrain(3)
      double precision  eqTstrain
	double precision  eqCStrain
	double precision  zeroStrain(3)
	double precision  tt,cc
	integer I,J
	external prinst
	double precision tDamage,cDamage,Damage

c------------Get principal strain-------------  
      INTEGER          MSVAR
      parameter       (MSVAR = 12)
      double precision  svar(MSVAR)
	do i=1, MSVAR
         svar(i) =0.0
      end do
      do i=1,ncomp
         svar(i) = Strain(i)
      end do
      call prinst(svar(1)) 
        pStrain(1)=svar(7)
        pStrain(2)=svar(8)
        pStrain(3)=svar(9)
c-------------equalstrain-------------  
	equalStrain=pStrain(1)**2
	equalStrain=equalStrain+pStrain(2)**2
	equalStrain=equalStrain+pStrain(3)**2
c------------Denominator cannot be zero-------------
	if (equalStrain .EQ.0 ) equalStrain=1.0D-20
      equalStrain=sqrt(equalStrain)
c---------eqTstrain---------
      do i=1,3
	    IF(pStrain(i) .GT.0)  then 
		zeroStrain(i) =pStrain(i)
	      else 
		  zeroStrain(i) =0
	    End if
	end do
	eqTstrain=zeroStrain(1)**2
	eqTstrain=eqTstrain+zeroStrain(2)**2
	eqTstrain=eqTstrain+zeroStrain(3)**2
	if (eqTstrain .EQ.0 ) eqTstrain=1.0D-20
	eqTstrain=sqrt(eqTstrain)
      do i=1,3
          zeroStrain(i) =0
	end do
c---------eqCstrain---------
      do i=1,3
	    IF(pStrain(i) .LT.0) then 
		zeroStrain(i) =pStrain(i)
	      else 
		  zeroStrain(i) =0
	    End if
	end do
      eqCstrain=zeroStrain(1)**2
	eqCstrain=eqCstrain+zeroStrain(2)**2
	eqCstrain=eqCstrain+zeroStrain(3)**2
       if (eqCstrain .EQ.0 ) eqCstrain=1.0D-20
	eqCstrain=sqrt(eqCstrain)
c----------tDamage--------
      IF ((eqTstrain.GT.0).AND.(eqTstrain.LT.Etf)) then 
	       tDamage=0
      END IF
      IF ((eqTstrain.GT.Etf).AND.(eqTstrain.LT.Etu)) then 
	       tDamage=Etu*(eqTstrain-Etf)
		   tDamage=tDamage/eqTstrain
		   tDamage=tDamage/(Etu-Etf)
      END IF

c----------cDamage--------
      IF ((eqCstrain.GT.0).AND.(eqCstrain.LT.Ecf)) then 
             cDamage=0
      END IF
      IF ((eqCstrain.GT.Ecf).AND.(eqCstrain.LT.Ecu)) then 
             cDamage=Ecu*(eqCstrain-Ecf)
		   cDamage=cDamage/eqCstrain
		   cDamage=cDamage/(Ecu-Ecf)
      END IF
c---------------total Damage------------
      cc=eqCstrain/equalStrain
	cc=cc**2
      tt=eqTstrain/equalStrain
	tt=tt**2
      Damage=cDamage*cc+tDamage*tt
	RETURN
	END




