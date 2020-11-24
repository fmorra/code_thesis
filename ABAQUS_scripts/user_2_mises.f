       SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
C      INCLUDE 'ABA_PARAM.INC'
C       INCLUDE 'param_umat.txt'
       COMMON /SHARE/ Bdiff,Bdisl,QMISES
C
       INTEGER ii,e
       DIMENSION TIME(2)
       REAL*8 Bdiff(218705),Bdisl(218705),QMISES(218705)
       CHARACTER(256) filename

C
       if (LOP.EQ.0) then

       filename='/home/fabri/Earth_model_abaqus_SLE0/results_run_22/e.dat'
           open(16,file=filename)
	     do ii=1,218705
	       read(16,*) e,Bdiff(ii),Bdisl(ii),QMISES(ii)
	     enddo
           close(16)
       endif
       return
       end 

	 
      SUBROUTINE CREEP(DECRA,DESWA,STATEV,SERD,EC,ESW,P,QTILD,
     1 TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,CMNAME,LEXIMP,LEND,
     2 COORDS,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C     INCLUDE 'ABA_PARAM.INC'
C       INCLUDE 'param_umat.txt'
      COMMON /SHARE/ Bdiff,Bdisl,QMISES
C
      CHARACTER*80 CMNAME
	INTEGER els(218705),index,ii
      REAL*8 A,ALIN
      REAL*8 Bdiff(218705),Bdisl(218705),QMISES(218705)
      DIMENSION DECRA(5),DESWA(5),STATEV(*),PREDEF(*),DPRED(*),
     1 TIME(2),COORDS(*),EC(2),ESW(2)
C
      AN = 3.5
      ANM1 = 2.5

	  ALIN = Bdiff(NOEL)
	  A = Bdisl(NOEL)
	  QTILD = QMISES(NOEL)
	  
	DECRA(1) = ALIN*QTILD*DTIME + A*QTILD**AN*DTIME 
	IF(LEXIMP.EQ.1) THEN
	  DECRA(5) = ALIN*DTIME + AN*A*QTILD**(ANM1)*DTIME
	END IF
      RETURN
      END
