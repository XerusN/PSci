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
    
    !nombre d'itérations temporelles
    INTEGER(KIND = IKind) :: n_t
    
    !dimension spatiale du maillage [m]
    REAL(KIND = RKind) :: l_x, l_y

    !vitesse de convection [m/s]
    REAL(KIND = RKind) :: c
    
    !viscosité [m2/s]
    REAL(KIND = RKind) :: viscosite

    !temps total d'observation du phénomène [s]
    REAL(KIND = RKind) :: t_f

    !CFL nombre de courant adimensionnel en input, en x et y
    REAL(KIND = RKind) :: cfl, cfl_x, cfl_y

    !Fo nombre de Fourier en x et y
    REAL(KIND = RKind) :: fo_x, fo_y

    !Conditions initiales
    REAL(KIND = RKind) :: x_max, x_min, y_max, y_min, u_in, u_out

    !conditions limites a gauche et a droite [m/s] (même unité que u)
    REAL(KIND = RKind) :: boundary_condition_left, boundary_condition_right
    REAL(KIND = RKind) :: boundary_condition_up, boundary_condition_down
    
    !definition de la vitesse [m/s] comme var globale, tableau de scalaire pour la partie 1, on ne stocke pas les itérations
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: u
    !tableau de stockage intermediaire
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: u_temp
    
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
    


    !maj à l'Etape 2, 2D
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
        
        !longueur espace
        READ(10, *)
        READ(10, *) l_x, l_y
        
        READ(10, *)
        
        !Temps de fin
        READ(10, *)
        READ(10, *) t_f
        !fréquence d'affichage
        READ(10, *)
        READ(10, *) frame
        
        READ(10, *)
        
        !Vitesse de convection
        READ(10, *)
        READ(10, *) c
        !Viscosité
        READ(10, *)
        READ(10, *) viscosite
        !CFL
        READ(10, *)
        READ(10, *) cfl
        
        READ(10, *)
        
        !Variables d'initialisations
        READ(10, *)
        READ(10, *) x_min, x_max
        READ(10, *)
        READ(10, *) y_min, y_max
        READ(10, *)
        READ(10, *) u_in, u_out
        
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
        
        CHARACTER(LEN = StrLen) :: name, format_iteration
        INTEGER(KIND = RKind) :: i, j
        
        WRITE(name, '(I0)') iteration
        name = 'output/resTECPLOT_' // TRIM(name) // '.dat'
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        
        WRITE(11, *) 'TITLE = "ETAPE2"'
        WRITE(11, *) 'VARIABLES = "X", "Y", "U"'
        WRITE(11, '(A, ES20.13, A, I4, A, I4, A)') 'ZONE T="', t, &
            '   seconds", I=', n_x, ', J=', n_y, ', DATAPACKING=POINT'
        
        !Ecriture pour chaque position
        DO i = 1, n_x
            DO j = 1, n_y
                WRITE(11, '(3(ES20.13, 1X))') space_grid(i, j)%x, space_grid(i, j)%y, u(i, j)
            END DO
        END DO
        
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE write_output_file
        
    
    
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
    
    
    
    !maj à l'Etape 2, 2D
    !Initialisation aux conditions initiales
    SUBROUTINE init_solution()
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        ALLOCATE(u(n_x, n_y))
        
        
        
        !Calcul pour chaque position
        DO i= 1, n_x
            DO j = 1, n_y
                IF ((space_grid(i,j)%x >= x_min) .AND. (space_grid(i,j)%x <= x_max) .AND. &
                    (space_grid(i,j)%y <= y_max) .AND. (space_grid(i,j)%y >= y_min)) THEN
                    u(i, j) = u_in
                ELSE
                    u(i, j) = u_out
                END IF
            END DO
        END DO
        
    END SUBROUTINE init_solution
    
    
    
    !maj à l'étape 2, 2D
    !Réalise toute les initialisations
    SUBROUTINE initialisation()
    
    IMPLICIT NONE
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = IKind) :: i
        
        t = 0
        
        !récupération des données du problème
        name = 'input.dat'
        CALL read_input_file(name)
        
        !création du maillage spatial
        CALL create_space_grid()
        
        !calcul du pas de temps à partir du CFL
        dt = ABS(cfl * dx / c)
    
        !initalisation de la solution grâce aux conditions initiales
        CALL init_solution()
        
        i = 0
        !Ecriture de la solution initiale
        CALL write_output_file(i)
        
    END SUBROUTINE initialisation
    
    
    
    !maj à l'Etape 2, 2D
    !calcul du profil de vitesse pour une itération temporelle
    SUBROUTINE convection_diffusion_2D()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j !compteur
        INTEGER :: upwind   !permet de maintenir la configuration upwind selon c
        
        !determine le sens upwind pour du/dt + c du/dx = 0
        IF (c > 0) THEN
            upwind = -1
        ELSE
            upwind = 1
        END IF
        
        !definition des conditions limites a gauche et a droite, inchangées
        u(1,:) = boundary_condition_left
        u(n_x,:) = boundary_condition_right
        u(:,1) = boundary_condition_down
        u(:,n_y) = boundary_condition_up

        !Calcul des CFL, pour cette étape ils sont égaux
        cfl_x = cfl
        cfl_y = cfl

        !calcul des nombres de Fourier
        fo_x = viscosite*dt/(dx**2)
        fo_y = viscosite*dt/(dy**2)
        
        !calcul du n+1
        DO j = 2, n_y-1
            DO i = 2, n_x-1
                u_temp(i,j) = u(i,j)*(1 - dt*((dx**2*(2*viscosite + c*dy) + dy**2*(2*viscosite + c*dx))/(dx**2 * dy**2))) &
                + u(i-1,j)*dt*((c*dx + viscosite)/(dx**2)) + u(i,j-1)*dt*((c*dy + viscosite)/(dy**2)) &
                + u(i+1,j)*((viscosite*dt)/(dx**2)) + u(i,j+1)*((viscosite*dt)/(dy**2))
            END DO
        END DO
        
        !assignation du n+1
        DO j = 2, n_y-1
            DO i = 2, n_x-1
                u(i,j) = u_temp(i,j)
            END DO
        END DO

    END SUBROUTINE convection_diffusion_2D
    
    
    
    !maj à l'étape 2, 2D
    !Réalise la boucle temporelle
    SUBROUTINE boucle_resolution()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i
        LOGICAL :: last_iteration
        
        i = 0
        
        !Allocation du tableau permettant de stocker les valeurs intermediaires de calcul
        ALLOCATE(u_temp(n_x, n_y))
        
        !permet de savoir si on a atteint la dernière iteration
        last_iteration = .FALSE.
        
        !Boucle temporelle du calcul
        DO WHILE (t < t_f)
            !verifie si on a atteint la dernière itération
            IF (t + dt > t_f) THEN
                dt = t_f - t
                last_iteration = .TRUE.
            END IF
            
            !appelle le calcul pour cette itération
            CALL convection_diffusion_2D()
            i = i + 1
            t = t + dt
            
            !écrit dans un fichier toute les frames
            IF (MOD(i, frame) == 0) THEN
                CALL write_output_file(i)
            END IF 
            
            IF (last_iteration .EQV. .TRUE.) THEN
                EXIT
            END IF
        END DO
        
        DEALLOCATE(u_temp)
        
    END SUBROUTINE boucle_resolution
    
    
    
    

END MODULE global






PROGRAM main

USE global

IMPLICIT NONE
    
    CALL initialisation()
    
    CALL boucle_resolution()
    
    DEALLOCATE(space_grid)
    DEALLOCATE(u)

END PROGRAM main