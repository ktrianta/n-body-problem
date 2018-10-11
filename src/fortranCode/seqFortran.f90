PROGRAM nBodySimulation_Sequential
IMPLICIT NONE 

INTEGER, PARAMETER :: N = 5
DOUBLE PRECISION, PARAMETER :: G = 1.0D0, EPSILON = 0.001
DOUBLE PRECISION, PARAMETER :: Time=10.0D0, DT = 0.00001
INTEGER :: numTimeSteps
DOUBLE PRECISION, DIMENSION(0:N-1) :: x=0.0D0,y=0.0D0,ux=0.0D0,uy=0.0D0,m=1.0D0,ax=0.0D0,ay=0.0D0
INTEGER :: i,k,t,j

numTimeSteps = INT(Time / DT + 1)
! Main program
CALL InitialPositions()

OPEN(UNIT=1,FILE="data.txt")

DO i=0,N-1
	WRITE(1,*) x(i),achar(9),y(i),achar(9),ux(i),achar(9),uy(i)
END DO

DO k=0,N-1
	CALL computeF(k,x,y,ax,N)
	CALL computeF(k,y,x,ay,N)
END DO


timeLoop:DO t=0,numTimeSteps-1
			
			DO j=0,N-1
				ux(j) = ux(j) + 0.5*ax(j)*DT
				uy(j) = uy(j) + 0.5*ay(j)*DT

				x(j) = x(j) + ux(j)*DT
				y(j) = y(j) + uy(j)*DT
			END DO

			DO j=0,N-1
				ax(j) = 0.0D0
				ay(j) = 0.0D0
			END DO

			DO j=0,N-1
				CALL computeF(j,x,y,ax,N)
				CALL computeF(j,y,x,ay,N)

				ux(j) = ux(j) + 0.5*ax(j)*DT
				uy(j) = uy(j) + 0.5*ay(j)*DT

				IF (MOD(t,200).EQ.0) THEN
					WRITE(1,*) x(j),achar(9),y(j),achar(9),ux(j),achar(9),uy(j)
				END IF

			END DO

		END DO timeLoop

CLOSE(UNIT=1)

CONTAINS 

SUBROUTINE InitialPositions()

	INTEGER i;

	DO i=0,N-1
		x(i) = RAND(i)
		y(i) = RAND(i)
	END DO

END SUBROUTINE InitialPositions

SUBROUTINE computeF(j,x,y,ax,N)

	INTEGER,INTENT(IN) :: j,N
	DOUBLE PRECISION,DIMENSION(0:N-1) :: x,y,ax
	INTEGER :: i
	DOUBLE PRECISION :: d,denom,force

	DO i=j+1,N-1
		IF (i.NE.j) THEN
			d = sqrt( ( x(i)-x(j) )**2 + ( y(i)-y(j) )**2 )
			denom = sqrt( (d**2+EPSILON**2)**3 )
			force = -G*m(j)*m(i)*( x(j)-x(i) )/denom;
			ax(j) = ax(j) + force / m(j)
			ax(i) = ax(i) - force / m(i)
		END IF
	END DO

END SUBROUTINE computeF

END PROGRAM 