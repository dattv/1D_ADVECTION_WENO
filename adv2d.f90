MODULE ADV2D_MODULE
    
    USE IR_Precision
    USE LIMITER_MODULE
    
    IMPLICIT NONE
    
    CONTAINS
    
    SUBROUTINE ADV2D()
    
    IMPLICIT NONE
    
    INTEGER(I4P)    :: X_CELL, Y_CELL
    REAL(R8P)       :: X_MIN, X_MAX, Y_MIN, Y_MAX
    REAL(R8P)       :: DX, DY
    REAL(R8P), DIMENSION(:), ALLOCATABLE    :: X_COORD, Y_COORD
    REAL(R8P), DIMENSION(:), ALLOCATABLE    :: VAR
    INTEGER(I4P), DIMENSION(:,:), ALLOCATABLE  :: ELECONNECT
    INTEGER(I4P)    :: I, J
    REAL(R8P)       :: MEANX, MEANY
    
    ! BODY
    
    !> CREATE 2D CATESIAN GRID
    
    X_CELL = 64
    Y_CELL = 64
    
    X_MIN = 0._R8P
    Y_MIN = X_MIN
    
    X_MAX = 1._R8P
    Y_MAX = X_MAX
    
    DX = (X_MAX-X_MIN)/X_CELL
    DY = (Y_MAX-Y_MIN)/Y_CELL
    
    ALLOCATE(X_COORD((X_CELL+1)*(Y_CELL+1)))
    ALLOCATE(Y_COORD((X_CELL+1)*(Y_CELL+1)))
    
    DO I = 1, X_CELL+1
        DO J = 1, Y_CELL+1
            X_COORD((J-1)*(X_CELL+1) + I) = (I-1)*DX
            Y_COORD((J-1)*(X_CELL+1) + I) = (J-1)*DY
        END DO
    END DO
    
    ALLOCATE(ELECONNECT(4, (X_CELL)*(Y_CELL)))
    
    DO I = 1, X_CELL
        DO J = 1, Y_CELL
            ELECONNECT(1,(J-1)*X_CELL+I) = (J-1)*(X_CELL+1)+I
            ELECONNECT(2,(J-1)*X_CELL+I) = (J-1)*(X_CELL+1)+I+1
            ELECONNECT(3,(J-1)*X_CELL+I) = J*(X_CELL+1)+I+1
            ELECONNECT(4,(J-1)*X_CELL+I) = J*(X_CELL+1)+I
        END DO
    END DO
    
    !> INITIAL CONDITION
    ALLOCATE(VAR(X_CELL*Y_CELL))
    DO I = 1, X_CELL
        DO J = 1, Y_CELL
            MEANX = ( X_COORD((J-1)*(X_CELL+1) + I) +  X_COORD((J-1)*(X_CELL+1) + I+1) +  X_COORD((J)*(X_CELL+1) + I) +  X_COORD((J)*(X_CELL+1) + I+1))*0.25_R8P
            MEANY = ( Y_COORD((J-1)*(X_CELL+1) + I) +  Y_COORD((J-1)*(X_CELL+1) + I+1) +  Y_COORD((J)*(X_CELL+1) + I) +  Y_COORD((J)*(X_CELL+1) + I+1))*0.25_R8P
            
            IF (MEANX >= 0.35 .AND. MEANX <= 0.65 .AND. MEANY >= 0.6 .AND. MEANY <= 0.9) THEN 
                VAR((J-1)*X_CELL + I) = 1._R8P
            ELSE
                VAR((J-1)*X_CELL + I) = 0._R8P
            END IF
            
        END DO
    END DO
    
    
    CALL OUTPUT2D("INITIAL2D.TEC", X_CELL, Y_CELL, X_COORD, Y_COORD, ELECONNECT, VAR)
    END SUBROUTINE ADV2D
    
!
!--------------------------------------------------------------------
!
    SUBROUTINE OUTPUT2D(FILENAME, X_CELL, Y_CELL, X_COORD, Y_COORD, ELECONNECT, VAR)
    IMPLICIT NONE
    
    INTEGER(I4P), INTENT(IN)    :: X_CELL, Y_CELL
    REAL(R8P), DIMENSION((X_CELL+1)*(Y_CELL+1)), INTENT(INOUT)    :: X_COORD, Y_COORD
    REAL(R8P), DIMENSION(X_CELL*Y_CELL), INTENT(IN) :: VAR
    CHARACTER(*), INTENT(IN)    :: FILENAME
    INTEGER(I4P), DIMENSION(4,X_CELL*Y_CELL), INTENT(IN)    :: ELECONNECT
    INTEGER(I4P)    :: I, J
    
    ! BODY 
    OPEN(UNIT = 100, FILE = TRIM(FILENAME), ACTION = 'WRITE')
    WRITE(100, *) 'TITLE="FLOW_YZ"'
    WRITE(100, *) 'VARIABLES = "X","Y"'

    
    WRITE(100, *) 'ZONE N = ', (X_CELL+1)*(Y_CELL+1), ',E = ', X_CELL*Y_CELL
    WRITE(100, *) 'ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK '
    DO I = 1, X_CELL+1
        DO J = 1, Y_CELL+1
            WRITE(100, *) X_COORD((J-1)*(X_CELL+1)+I)
        END DO
    END DO      
    DO I = 1, X_CELL+1
        DO J = 1, Y_CELL+1
            WRITE(100, *) Y_COORD((J-1)*(X_CELL+1)+I)
        END DO
    END DO  
    DO I = 1, X_CELL+1
        DO J = 1, Y_CELL+1
            IF (I == 1 .OR. J == 1) THEN 
                WRITE(100, *) VAR((J-1)*(X_CELL)+I)
            ELSE IF (I <= X_CELL .AND. J <= Y_CELL) THEN               
                WRITE(100, *) (VAR((J-1)*(X_CELL)+I) + VAR((J-1)*(X_CELL)+I+1) + VAR((J)*(X_CELL)+I) + VAR((J)*(X_CELL)+I+1))*0.25_R8P
            ELSE
                WRITE(100, *) VAR((J-1)*(X_CELL)+I)
            END IF
        END DO
    END DO
    !END DO      
    DO I = 1, X_CELL
        DO J = 1, Y_CELL
            WRITE(100, *) ELECONNECT(1,(J-1)*X_CELL+I), ELECONNECT(2,(J-1)*X_CELL+I), ELECONNECT(3,(J-1)*X_CELL+I), ELECONNECT(4,(J-1)*X_CELL+I)
        END DO
    END DO
    
    !WRITE(100, *) (VAR(I),I = 1, X_CELL*Y_CELL)
    CLOSE(UNIT = 100)
    END SUBROUTINE OUTPUT2D
END MODULE ADV2D_MODULE    