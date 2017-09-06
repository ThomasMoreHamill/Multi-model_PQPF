      SUBROUTINE DOY(JDATE,JY,JM,JD,JH,MDOY) 
C 
C        NOVEMBER 1996   GLAHN   TDL   MOS-2000
C
C        PURPOSE 
C            TO PARSE THE DATE/TIME IN YYYYMMDDHH FORMAT AND
C            TO COMPUTE THE DAY OF THE YEAR.
C 
C        DATA SET USE 
C            NONE. 
C 
C        VARIABLES 
C 
C            INPUT
C              JDATE = BASIC DATE OF FORMAT YYYYMMDDHH
C                      (E.G., 1993120108).
C
C            OUTPUT               
C                  JY = YEAR (E.G., 1993). 
C                  JM = MONTH (E.G., 12). 
C                  JD = DAY (E.G., 01). 
C                  JH = HOUR (E.G., 08). 
C                MDOY = DAY OF THE YEAR.
C
C            INTERNAL
C            MTEST(J) = NUMBER OF DAYS IN ALL MONTHS UP
C                       TO, BUT NOT INCLUDING, THE CURRENT
C                       ONE J (J=1,12). 
C 
C        NONSYSTEM SUBROUTINES CALLED 
C            NONE. 
C 
      DIMENSION MTEST(12)
      DATA MTEST/0,31,59,90,120,151,181,212,243,273,304,334/ 
C 
C        BREAK DATE INTO COMPONENT PARTS. 
C 
      JY=JDATE/1000000
      JM=JDATE/10000-JY*100
      JD=JDATE/100-JY*10000-JM*100
      JH=MOD(JDATE,100)
C
C        COMPUTE DOY.
C
      MDOY=MTEST(JM)+JD
C 
C        HANDLE LEAP YEAR, WHICH INCLUDES YEAR 2000
C        BUT NOT OTHER CENTENIAL YEARS.  ADD A DAY
C        IN MONTHS FOLLOWING FEBRUARY.
C
      IF(MOD(JY,4).EQ.0.AND.JM.GT.2)MDOY=MDOY+1 
      RETURN 
      END 
