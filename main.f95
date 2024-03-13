MODULE global

IMPLICIT NONE

    !definition des tailles des reels(RKind) et des entiers(IKind)
    INTEGER, PARAMETER :: RKind = SELECTED_REAL_KIND(10,200)
    INTEGER, PARAMETER :: IKind = SELECTED_INT_KIND(5)
    
    !Constante de taille pour le nom des fichiers
    INTEGER, PARAMETER :: StrLen = 40
    
    
    
    !-----------------------------------------------------
    
    
    
    !Variables d'entree
    
    !Taille du maillage spatial
    INTEGER(KIND = IKind) :: n_x, n_y
    
    !nombre d'itérations temporelles imposé
    INTEGER(KIND = IKind) :: n_t
    
    !dimension spatiale du maillage [m]
    REAL(KIND = RKind) :: l_x, l_y
    
    !viscosité [m2/s]
    REAL(KIND = RKind) :: viscosity

    !CFL nombre de courant adimensionnel en input, en x et y
    REAL(KIND = RKind) :: cfl

    !Fo nombre de Fourier en x et y
    REAL(KIND = RKind) :: fo

    !Conditions initiales pour la solution
    REAL(KIND = RKind) :: x_max, x_min, y_max, y_min, u_in, u_out, v_in, v_out

    !conditions limites a gauche et a droite [m/s] (même unité que u)
    REAL(KIND = RKind) :: boundary_condition_left, boundary_condition_right
    REAL(KIND = RKind) :: boundary_condition_up, boundary_condition_down
    
    !definition des composantes de la vitesse [m/s] comme vars globales, tableaux de scalaire pour la partie 3, on ne stocke pas les itérations
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: u, v

    !tableau de stockage intermediaire
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: u_temp, v_temp
    
    !structure qui contient les coordonnées en x et y d'un point du maillage
    TYPE COORDS
        REAL(KIND = RKind) :: x, y
    END TYPE

    !definition du maillage spatial, exprime les coordonnées en [m] de la case
    TYPE(COORDS), DIMENSION(:, :), ALLOCATABLE :: space_grid
    
    !Variable temporelle
    REAL(KIND = RKind) :: t
    
    !variable de fréquence d'affichage
    INTEGER(KIND = IKind) :: frame
    
    !definition du pas spatial [m] et du pas de temps [s]
    REAL(KIND = RKind) :: dx, dy, dt

    
    
    
CONTAINS
    


    !maj à l'Etape 3, 2D, retiré les variables plus utilisées
    SUBROUTINE read_input_file(name)
    IMPLICIT NONE
        
        CHARACTER(LEN = StrLen), INTENT(IN) :: name
        
        !Ouverture du fichier
        OPEN(10, FILE = name)
        
        
        !Saut des deux premieres lignes
        READ(10, *)
        READ(10, *)
        
        !Taille du maillage spatial
        READ(10, *)
        READ(10, *) n_x, n_y
        
        READ(10, *)

        !Nombre d'itérations imposées
        READ(10, *)
        READ(10, *) n_t
        
        READ(10, *)
        
        !longueur espace
        READ(10, *)
        READ(10, *) l_x, l_y
        
        READ(10, *)
        
        !fréquence d'affichage
        READ(10, *)
        READ(10, *) frame
        
        READ(10, *)
        
        !Viscosité
        READ(10, *)
        READ(10, *) viscosity
        !CFL
        READ(10, *)
        READ(10, *) cfl
        !Nombre de fourrier
        READ(10, *)
        READ(10, *) fo
        
        READ(10, *)
        
        !Variables d'initialisations
        READ(10, *)
        READ(10, *) x_min, x_max
        READ(10, *)
        READ(10, *) y_min, y_max
        READ(10, *)
        READ(10, *) u_in, u_out, v_in, v_out
        
        READ(10, *)
        
        !Conditions aux limites spatiales
        READ(10, *)
        READ(10, *) boundary_condition_left
        READ(10, *)
        READ(10, *) boundary_condition_right
        READ(10, *)
        READ(10, *) boundary_condition_up
        READ(10, *)
        READ(10, *) boundary_condition_down
        
        
        
        !Fermeture du fichier
        CLOSE(10)
        
    END SUBROUTINE read_input_file
    
    
    
    !maj à l'Etape 2, 2D
    !subroutine pour l'ecriture des données dans un fichier
    SUBROUTINE write_output_file(iteration)
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), INTENT(IN) :: iteration
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = RKind) :: i, j
        
        WRITE(name, '(I0)') iteration
        name = 'output/resTECPLOT_' // TRIM(name) // '.dat'
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        
        WRITE(11, *) 'TITLE = "ETAPE2"'
        WRITE(11, *) 'VARIABLES = "X", "Y", "U", "V"'
        WRITE(11, '(1X, A, ES20.13, A, I4, A, I4, A)') 'ZONE T="', t, &
            '   seconds", I=', n_x, ', J=', n_y, ', DATAPACKING=POINT'
        
        !Ecriture pour chaque position
        DO i = 1, n_x
            DO j = 1, n_y
                WRITE(11, '(4(ES20.13, 1X))') space_grid(i, j)%x, space_grid(i, j)%y, u(i, j), v(i, j)
            END DO
        END DO
        
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE write_output_file
    
    
    !maj Etape 2
    !Permet de rapidement tester la valeur de certaines variables
    SUBROUTINE Debug(iteration)
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), INTENT(IN) :: iteration
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = RKind) :: i, j
        
        WRITE(name, '(I0)') iteration
        name = 'debug/var' // TRIM(name) // '.dat'
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        
        WRITE(11, *) 'case finale', space_grid(n_x, n_y)
        WRITE(11, *) 'dx', dx, 'dy', dy
        WRITE(11, *) 'n_x', n_x, 'n_y', n_y
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE Debug
    
        
    
    
    !maj à l'Etape 2, 2D
    !subroutine de creation du maillage spatial, contient les coordonnées exactes de chaque pt, en 2D
    SUBROUTINE create_space_grid()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j
    
        
        ALLOCATE(space_grid(n_x, n_y))
        
        !calcul du pas spatial en x
        dx = l_x / REAL(n_x - 1)

        !calcul du pas spatial en y
        dy = l_y / REAL(n_y - 1)

        !assignation des coordonnées
        DO j = 1, n_y
            DO i = 1, n_x
                space_grid(i,j)%x = dx * REAL(i-1)
                space_grid(i,j)%y = dy * REAL(j-1) 
            END DO
        END DO
    END SUBROUTINE create_space_grid
    
    
    
    !maj à l'Etape 3, 2D
    !Initialisation aux conditions initiales, elles s'appliquent aux deux composantes u et v
    SUBROUTINE init_solution()
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        ALLOCATE(u(n_x, n_y))
        ALLOCATE(v(n_x, n_y))
        
        !Assignation pour chaque position
        DO i= 1, n_x
            DO j = 1, n_y
                IF ((space_grid(i,j)%x >= x_min) .AND. (space_grid(i,j)%x <= x_max) .AND. &
                    (space_grid(i,j)%y <= y_max) .AND. (space_grid(i,j)%y >= y_min)) THEN
                    u(i, j) = u_in
                    v(i, j) = v_in
                ELSE
                    u(i, j) = u_out
                    v(i, j) = v_out
                END IF
            END DO
        END DO
        
    END SUBROUTINE init_solution
    
    
    
    !maj à l'étape 3, 2D
    !Réalise toute les initialisations
    SUBROUTINE initialisation()
    
    IMPLICIT NONE
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = IKind) :: i
        
        !Nettoie le dossier output
        CALL EXECUTE_COMMAND_LINE("rm ./output/*.dat")
        !Nettoie le dossier debug
        CALL EXECUTE_COMMAND_LINE("rm ./debug/*.dat")
        
        
        t = 0
        
        !récupération des données du problème
        name = 'input.dat'
        CALL read_input_file(name)
        
        !création du maillage spatial
        CALL create_space_grid()
    
        !initalisation de la solution grâce aux conditions initiales
        CALL init_solution()
        
        i = 0
        !Ecriture de la solution initiale
        CALL write_output_file(i)
        
    END SUBROUTINE initialisation



    !créé à l'étape 3
    !calcul du pas de temps pour une itération, cfl imposé constant, mais dépend de u ou de v 
    SUBROUTINE compute_time_step()

    IMPLICIT NONE

        INTEGER(KIND = IKIND) :: i, j
        REAL(KIND = RKIND) :: dt_min, dt_temp

        !On cherche le dt le plus contraignant, càd le plus petit

        !initialisation du min en utilisant la condition du fourrier
        dt_min = fo/viscosity*dx**2

        !recherche des u_max, v_max
        DO i = 1, n_x
            DO j = 1, n_y
                dt_temp = cfl*dx/ABS(u(i,j))     !calcul du minimum potentiel sur u

                IF (dt_temp < dt_min) THEN
                    dt_min = dt_temp
                END IF

                dt_temp = cfl*dy/ABS(v(i,j))     !calcul du minimum potentiel sur v

                IF (dt_temp < dt_min) THEN
                    dt_min = dt_temp
                END IF
            END DO
        END DO
        
        dt = dt_min
        
    END SUBROUTINE compute_time_step
    
    
    
    !maj à l'Etape 2, 2D
    !calcul du profil de vitesse pour une itération temporelle
    SUBROUTINE burgers_2D()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j !compteur
        INTEGER(KIND = IKIND) :: upwind_x, upwind_y     !permet de déterminer la direction du upwind (1 si backward, 0 sinon)
        
        !Application des conditions limites sur u
        u(1,:) = boundary_condition_left
        u(n_x,:) = boundary_condition_right
        u(:,1) = boundary_condition_down
        u(:,n_y) = boundary_condition_up
        
        !Application des conditions limites sur v
        v(1,:) = boundary_condition_left
        v(n_x,:) = boundary_condition_right
        v(:,1) = boundary_condition_down
        v(:,n_y) = boundary_condition_up
        
        
        !calcul du n+1
        DO j = 2, n_y-1
            DO i = 2, n_x-1
                
                !Ajustement du sens upwind pour u
                IF (u(i, j) > 0) THEN
                    upwind_x = 1
                ELSE
                    upwind_x = 0
                END IF
                
                !Ajustement du sens upwind pour v
                IF (v(i, j) > 0) THEN
                    upwind_y = 1
                ELSE
                    upwind_y = 0
                END IF
                
                !Calcul de u
                u_temp(i,j) = u(i,j)*(1.0 - 2.0*viscosity*dt*(1.0/dx**2 + 1.0/dy**2) &
                - u(i,j)*dt/dx*2.0*(REAL(upwind_x) - 0.5) - v(i,j)*dt/dy*2.0*(REAL(upwind_y) - 0.5)) &
                + u(i-1,j)*dt*(viscosity/dx**2 + u(i,j)/dx*REAL(upwind_x)) &
                + u(i,j-1)*dt*(viscosity/dy**2 + v(i,j)/dy*REAL(upwind_y)) &
                + u(i+1,j)*dt*(viscosity/dx**2 - u(i,j)/dx*(1 - REAL(upwind_x))) &
                + u(i,j+1)*dt*(viscosity/dy**2 - v(i,j)/dy*(1 - REAL(upwind_y)))
                
                !Calcul de v
                v_temp(i,j) = v(i,j)*(1.0 - 2.0*viscosity*dt*(1.0/dx**2 + 1.0/dy**2) &
                - u(i,j)*dt/dx*2.0*(REAL(upwind_x) - 0.5) - v(i,j)*dt/dy*2.0*(REAL(upwind_y) - 0.5)) &
                + v(i-1,j)*dt*(viscosity/dx**2 + u(i,j)/dx*REAL(upwind_x)) &
                + v(i,j-1)*dt*(viscosity/dy**2 + v(i,j)/dy*REAL(upwind_y)) &
                + v(i+1,j)*dt*(viscosity/dx**2 - u(i,j)/dx*(1 - REAL(upwind_x))) &
                + v(i,j+1)*dt*(viscosity/dy**2 - v(i,j)/dy*(1 - REAL(upwind_y)))
            END DO
        END DO
        
        !assignation du n+1
        DO j = 2, n_y-1
            DO i = 2, n_x-1
                u(i,j) = u_temp(i,j)
                v(i,j) = v_temp(i,j)
            END DO
        END DO
        
        

    END SUBROUTINE burgers_2D
    
    
    
    !maj à l'étape 2, 2D
    !Réalise la boucle temporelle
    SUBROUTINE resolution_loop()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i
        
        i = 0
        
        !Allocation du tableau permettant de stocker les valeurs intermediaires de calcul
        ALLOCATE(u_temp(n_x, n_y))
        ALLOCATE(v_temp(n_x, n_y))
        
        !Boucle temporelle du calcul
        DO WHILE (i < n_t)
            
            CALL compute_time_step()
            
            !appelle le calcul pour cette itération
            CALL burgers_2D()
            i = i + 1
            t = t + dt
            
            !écrit dans un fichier toute les frames
            IF (MOD(i, frame) == 0) THEN
                CALL write_output_file(i)
            END IF 
            
            !IF (i>=200 .AND. i <= 210) THEN
            !    CALL debug(i)
            !END IF
            
            !PRINT*, i
        END DO
        
        DEALLOCATE(u_temp)
        DEALLOCATE(v_temp)
        
    END SUBROUTINE resolution_loop
    
    
    
    

END MODULE global






PROGRAM main

USE global

IMPLICIT NONE
    
    CALL initialisation()
    
    CALL resolution_loop()
    
    DEALLOCATE(space_grid)
    DEALLOCATE(u)
    DEALLOCATE(v)
    
    PRINT*, n_t
    
END PROGRAM main