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
    !masse volumique
    REAL(KIND = RKind) :: density

    !Nombre de Reynolds
    REAL(KIND = RKind) :: re

    !CFL nombre de courant adimensionnel en input, en x et y
    REAL(KIND = RKind) :: cfl

    !Fo nombre de Fourier en x et y
    REAL(KIND = RKind) :: fo

    ! !conditions limites a gauche et a droite [m/s] (même unité que u)
    ! REAL(KIND = RKind) :: boundary_condition_left, boundary_condition_right
    ! REAL(KIND = RKind) :: boundary_condition_up, boundary_condition_down
    
    !definition des composantes de la vitesse [m/s] comme vars globales
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: u, v
    !champ de pression
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: p

    !tableau de stockage intermediaire
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: u_temp, v_temp
    
    !matrice A de l'equation de poisson
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: a, a_loc
    !vecteur b de l'équation de poisson
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: b
    !Stockage de p sous forme vectorielle et residu de la méthode de jacobi
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: p_vec, p_vec_temp, jacobi_r
    
    !structure qui contient les coordonnées en x et y d'un point du maillage
    TYPE COORDS
        REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: x, y
    END TYPE

    !definition du maillage spatial, exprime les coordonnées en [m] de la case
    TYPE(COORDS) :: space_grid
    
    !Variable temporelle et temps final
    REAL(KIND = RKind) :: t, t_f
    
    !variable de fréquence d'affichage
    INTEGER(KIND = IKind) :: frame
    
    !definition du pas spatial [m] et du pas de temps [s]
    REAL(KIND = RKind) :: dx, dy, dt

    !variable pour le choix de discretisation des termes convectifs
    CHARACTER(LEN = 3) :: scheme

    
    
    
CONTAINS
    

    !maj à l'Etape 7, 2D, retiré les variables plus utilisées et ajout du choix de scheme pour les termes convectifs
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

        !Temps final
        READ(10, *)
        READ(10, *) t_f
        
        READ(10, *)
        
        !longueur espace
        READ(10, *)
        READ(10, *) l_x, l_y
        
        READ(10, *)
        
        !fréquence d'affichage
        READ(10, *)
        READ(10, *) frame
        
        READ(10, *)
        
        !masse volumique
        READ(10, *)
        READ(10, *) density
        
        !CFL
        READ(10, *)
        READ(10, *) cfl
        !Nombre de fourrier
        READ(10, *)
        READ(10, *) fo
        
        READ(10, *)
        
        !Nombre de Reynolds
        READ(10, *)
        READ(10, *) re
        viscosity = l_x/re
        READ(10, *)
        
        !Scheme choisi pour termes convectifs
        READ(10, *)
        READ(10, *) scheme

        
        ! READ(10, *)
        
        ! !Conditions limites
        ! READ(10, *)
        ! READ(10, *) boundary_condition_up
        ! READ(10, *)
        ! READ(10, *) boundary_condition_down
        ! READ(10, *)
        ! READ(10, *) boundary_condition_left
        ! READ(10, *)
        ! READ(10, *) boundary_condition_right
        
        
        
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
        WRITE(11, *) 'VARIABLES = "X", "Y", "U", "V", "MAG", "P"'
        WRITE(11, '(1X, A, ES20.13, A, I4, A, I4, A)') 'ZONE T="', t, &
            '   seconds", I=', n_x, ', J=', n_y, ', DATAPACKING=POINT'
        
        !Ecriture pour chaque position
        DO i = 1, n_x
            DO j = 1, n_y
                WRITE(11, '(6(ES20.13, 1X))') space_grid%x(i), space_grid%y(j), u(i, j), v(i, j), &
                SQRT(u(i, j)**2.0_RKIND + v(i, j)**2.0_RKIND), p(i, j)
            END DO
        END DO
        
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE write_output_file
    
    
    !maj Etape 2
    !Permet de rapidement tester la valeur de certaines variables
    SUBROUTINE debug(iteration)
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), INTENT(IN) :: iteration
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = RKind) :: i, j
        
        WRITE(name, '(I0)') iteration
        name = 'debug/b_' // TRIM(name) // '.dat'
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        DO i = 1, n_x*n_y
            WRITE(11, *) 'i = ', i, ' | b(i) = ', b(i)
        END DO
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE Debug
    
        
    
    
    !maj à l'Etape 2, 2D
    !subroutine de creation du maillage spatial, contient les coordonnées exactes de chaque pt, en 2D
    SUBROUTINE create_space_grid()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j
    
        
        ALLOCATE(space_grid%x(n_x))
        ALLOCATE(space_grid%y(n_y))
        
        !calcul du pas spatial en x
        dx = l_x / REAL(n_x - 1)

        !calcul du pas spatial en y
        dy = l_y / REAL(n_y - 1)

        !assignation des coordonnées
        DO j = 1, n_y
            DO i = 1, n_x
                space_grid%x(i) = dx * REAL(i-1)
                space_grid%y(j) = dy * REAL(j-1) 
            END DO
        END DO
        
    END SUBROUTINE create_space_grid
    
    
    
    !maj à l'Etape 6, 2D
    !Initialisation aux conditions initiales, elles s'appliquent aux deux composantes u et v
    SUBROUTINE init_solution()
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        ALLOCATE(u(n_x, n_y))
        ALLOCATE(v(n_x, n_y))
        ALLOCATE(p(n_x, n_y))
        
        !Assignation pour chaque position
        u(:, :) = 0_RKind
        v(:, :) = 0_RKind
        p(:, :) = 0_RKind
        
        u(:, n_y) = 1.0_RKind
        
    END SUBROUTINE init_solution
    
    
    
    !maj à l'étape 4, 2D
    !Création et remplissage de la matrice A et allocation des tableaux nécessaires à Jacobi
    SUBROUTINE jacobi_init()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j, k, k_max
        REAL(KIND = RKIND) :: inv_x_2, inv_y_2

        !calcul de la dimension du vecteur solution
        k_max = n_x*n_y

        ALLOCATE(a(k_max, k_max))
        ALLOCATE(a_loc(k_max, k_max))
        ALLOCATE(p_vec(k_max))
        ALLOCATE(p_vec_temp(k_max))
        ALLOCATE(jacobi_r(k_max))
        ALLOCATE(b(k_max))
        


        !Remplissage de la matrice A
        !initialisation
        a(:, :) = 0

        !calcul des coefficients
        inv_x_2 = 1_RKind/(dx**2.0_RKind)
        inv_y_2 = 1_RKind/(dy**2.0_RKind)

        
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                
                !Conditions limites de dérivée nulle
                IF (i == 1) THEN
                    
                    a(k, k) = 1_RKind
                    a(k, k + 1) = -1_RKind
                
                
                ELSE IF (i == n_x) THEN
                    
                    a(k, k) = 1_RKind
                    a(k, k - 1) = -1_RKind
                
                
                ELSE IF (j == 1) THEN
                
                    a(k, k) = 1_RKind
                    a(k, k + n_x) = -1_RKind
                
                
                ELSE IF (j == n_y) THEN
                    
                    a(k, k) = 1_RKind
                    a(k, k - n_x) = -1_RKind
                    
                ELSE 
                    !sinon on applique les coefficients de l'équation
                    a(k, k) = -2_RKind*(inv_x_2 + inv_y_2)
                    a(k, k + 1) = inv_x_2
                    a(k, k - 1) = inv_x_2
                    a(k, k + n_x) = inv_y_2
                    a(k, k - n_x) = inv_y_2
                    

                END IF

            END DO
        END DO
        
    END SUBROUTINE jacobi_init
    
    
    
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
        
        CALL jacobi_init()
        
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
    
    
    
    !remplit le vecteur b (terme de droite de l'equation de poisson)
    SUBROUTINE fill_b ()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: k_max, i, j, k
        
        k_max = n_x*n_y
        
        
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN

                    !application de la CL à b
                    b(k) = 0_RKind
                    
                ELSE 
                    
                    !valeur de b pour tous les autres points
                    b(k) = (u_temp(i+1, j) - u_temp(i-1, j))/(2_RKind*dx) + (v_temp(i, j+1) - v_temp(i, j-1))/(2_RKind*dy)
                    b(k) = b(k)*density/dt
                    !PRINT*, i, j, u_temp(i, j), (u_temp(i+1, j) - u_temp(i-1, j)), density/dt
                END IF
            END DO
        END DO
        
        
    END SUBROUTINE fill_b
    
    
    
    !maj à l'Etape 2, 2D!Calcul la norme 2 d'un vecteur
    SUBROUTINE norm_2(vec, norm)
    
    IMPLICIT NONE
        
        REAL(KIND = Rkind), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vec
        
        INTEGER(KIND = IKIND) :: vec_size, i
        REAL(KIND = RKIND) :: norm
        
        vec_size = SIZE(vec, 1)
        
        norm = SUM(vec(:)**2_RKind)
        
        
        norm = SQRT(norm)
        
    END SUBROUTINE norm_2
    
    
    
    !Methode de Jacobi (resolution iterative de systeme lineaire)
    SUBROUTINE jacobi_method()
    
    IMPLICIT NONE
        
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001     !point d'arret de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000       !Arret forcé de jacobi
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        
        !CALL CPU_TIME(time1)
        
        !Tentative initiale
        
        DO i = 1, n_x
            DO j = 1, n_y
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        p_vec_temp(:) = p_vec(:)
        
        k_max = SIZE(p_vec)
        
        DO i = 1, k_max
            jacobi_r(i) = SUM(a(i, :)*p_vec(:)) - b(i)
        END DO
        
        iteration = 0
        
       lower_norm = 1
       upper_norm = 1 
        DO WHILE (upper_norm/lower_norm > RTol)
            p_vec_temp(:) = p_vec(:)
            DO i = 1, k_max
                
                p_vec(i) = p_vec_temp(i) - jacobi_r(i)/a(i, i)
            END DO
            
            !Integrale de la pression
            ! integral = 0_RKIND
            ! DO i = 1, n_x
            !     DO j = 1, n_y
            !         IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN
            !             integral = integral + p_vec((j-1)*n_x+i)*dx*dy/2
            !         ELSE
            !             integral = integral + p_vec((j-1)*n_x+i)*dx*dy
            !         END IF
            !     END DO
            ! END DO
            ! DO i = 1, n_x
            !     DO j = 1, n_y
            !         IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN
            !             p_vec((j-1)*n_x+i) = p_vec((j-1)*n_x+i) - integral/(dx*dy/2)
            !         ELSE
            !             p_vec((j-1)*n_x+i) = p_vec((j-1)*n_x+i) - integral/(dx*dy)
            !         END IF
            !     END DO
            ! END DO
            integral = SUM(p_vec(:))
            p_vec(:) = p_vec(:) - integral/((n_x)*(n_y))
            
            CALL norm_2(p_vec, lower_norm)
            jacobi_r(:) = (p_vec(:)-p_vec_temp(:))
            CALL norm_2(jacobi_r, upper_norm)
            
            jacobi_r(:) = - b(:)
            DO j = 1, k_max
                jacobi_r(:) = jacobi_r(:) + a(:, j)*p_vec(j)
            END DO
            
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de jacobi atteint (', IterationMax, ' )'
                STOP
            END IF
            
            IF (MOD(iteration, 100) == 0) THEN
                PRINT*, 'iteration = ', iteration, ' | convergence = ', upper_norm/lower_norm
            END IF
            
            !Integrale
            !integral = SUM(p_vec(:))
            !PRINT*, integral
            
            
        END DO
        
        !CALL CPU_TIME(time2)
        
        DO i = 1, n_x
            DO j = 1, n_y
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        integral = 0_RKIND
        DO i = 1, n_x
            DO j = 1, n_y
                IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN
                    integral = integral + p_vec((j-1)*n_x+i)*dx*dy/4
                ELSE
                    integral = integral + p_vec((j-1)*n_x+i)*dx*dy
                END IF
            END DO
        END DO
        
        !PRINT*, 'Jacobi for a grid size of ', n_x, ' : ', time2 - time1, ' seconds (', iteration, ' iterations)'
        PRINT*, 'Jacobi :', iteration, ', iterations | integrale(p) = ', integral
        
    END SUBROUTINE jacobi_method
    
    
    
    !calcul du profil de vitesse pour une itération temporelle
    SUBROUTINE speed_guess()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j !compteur
        INTEGER(KIND = IKIND) :: upwind_x, upwind_y     !permet de déterminer la direction du upwind (1 si backward, 0 sinon)
        
        !Application des conditions limites sur u
        u_temp(1,:) = 0_RKind
        u_temp(n_x,:) = 0_RKind
        u_temp(:,1) = 0_RKind
        u_temp(:,n_y) = 1_RKind
        
        !Application des conditions limites sur v
        v_temp(1,:) = 0_RKind
        v_temp(n_x,:) = 0_RKind
        v_temp(:,1) = 0_RKind
        v_temp(:,n_y) = 0_RKind
        
        
        IF (scheme == 'UR1') THEN

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
                    u_temp(i,j) = u(i,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind) &
                    - u(i,j)*dt/dx*2.0_RKind*(REAL(upwind_x) - 0.5_RKind) - v(i,j)*dt/dy*2.0_RKind*(REAL(upwind_y) - 0.5_RKind)) &
                    + u(i-1,j)*dt*(viscosity/dx**2_RKind + u(i,j)/dx*REAL(upwind_x)) &
                    + u(i,j-1)*dt*(viscosity/dy**2_RKind + v(i,j)/dy*REAL(upwind_y)) &
                    + u(i+1,j)*dt*(viscosity/dx**2_RKind - u(i,j)/dx*REAL(1_RKind - upwind_x)) &
                    + u(i,j+1)*dt*(viscosity/dy**2_RKind - v(i,j)/dy*REAL(1_RKind - upwind_y))
                    
                    !Calcul de v
                    v_temp(i,j) = v(i,j)*(1.0_RKind - 2.0_RKind*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2) &
                    - u(i,j)*dt/dx*2.0_RKIND*(REAL(upwind_x) - 0.5_RKind) - v(i,j)*dt/dy*2.0_RKind*(REAL(upwind_y) - 0.5_RKind)) &
                    + v(i-1,j)*dt*(viscosity/dx**2_RKind + u(i,j)/dx*REAL(upwind_x)) &
                    + v(i,j-1)*dt*(viscosity/dy**2_RKind + v(i,j)/dy*REAL(upwind_y)) &
                    + v(i+1,j)*dt*(viscosity/dx**2_RKind - u(i,j)/dx*REAL(1_RKind - upwind_x)) &
                    + v(i,j+1)*dt*(viscosity/dy**2_RKind - v(i,j)/dy*REAL(1_RKind - upwind_y))
                END DO
            END DO

        ELSE IF (scheme == 'CD4') THEN
        
            

            !calcul du n+1
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    
                    IF ((i == 2) .OR. (i == n_x-1) .OR. (j == 2) .OR. (j == n_y-1)) THEN
                        
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
                        u_temp(i,j) = u(i,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind) &
                        - u(i,j)*dt/dx*2.0_RKind*(REAL(upwind_x) - 0.5_RKind) - v(i,j)*dt/dy*2.0_RKind*(REAL(upwind_y) &
                        - 0.5_RKind)) + u(i-1,j)*dt*(viscosity/dx**2_RKind + u(i,j)/dx*REAL(upwind_x)) &
                        + u(i,j-1)*dt*(viscosity/dy**2_RKind + v(i,j)/dy*REAL(upwind_y)) &
                        + u(i+1,j)*dt*(viscosity/dx**2_RKind - u(i,j)/dx*REAL(1_RKind - upwind_x)) &
                        + u(i,j+1)*dt*(viscosity/dy**2_RKind - v(i,j)/dy*REAL(1_RKind - upwind_y))
                        
                        !Calcul de v
                        v_temp(i,j) = v(i,j)*(1.0_RKind - 2.0_RKind*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2) &
                        - u(i,j)*dt/dx*2.0_RKIND*(REAL(upwind_x) - 0.5_RKind) - v(i,j)*dt/dy*2.0_RKind*(REAL(upwind_y) &
                        - 0.5_RKind)) + v(i-1,j)*dt*(viscosity/dx**2_RKind + u(i,j)/dx*REAL(upwind_x)) &
                        + v(i,j-1)*dt*(viscosity/dy**2_RKind + v(i,j)/dy*REAL(upwind_y)) &
                        + v(i+1,j)*dt*(viscosity/dx**2_RKind - u(i,j)/dx*REAL(1_RKind - upwind_x)) &
                        + v(i,j+1)*dt*(viscosity/dy**2_RKind - v(i,j)/dy*REAL(1_RKind - upwind_y))
                        
                    ELSE
                        
                        !Calcul de u
                        u_temp(i,j) = u(i,j)*(1.0_RKind - dt*(2.0_RKind/(3.0_RKind*dx)*(u(i+1,j) - u(i-1,j)) &
                        - 1.0_RKind/(12_RKIND*dx)*(u(i+2,j) - u(i-2,j)) + 2.0_RKind*viscosity/dx**2.0_RKind + & 
                        2.0_RKind*viscosity/dy**2.0_RKind)) - dt*u(i,j+1)*(2.0_RKind/(3.0_RKind*dy)*v(i,j) - & 
                        viscosity/dy**2.0_RKind) + dt*u(i,j-1)*(2.0_RKind/(3.0_RKind*dy)*v(i,j) + viscosity/dy**2.0_RKind) &
                        + dt*u(i,j+2)*(1.0_RKind/(12_RKIND*dy)*v(i,j)) - dt*u(i,j-2)*(1.0_RKind/(12_RKIND*dy)*v(i,j)) &
                        + dt*viscosity/dx**2.0_RKind*(u(i+1,j) + u(i-1,j))
                    
                        !Calcul de v
                        v_temp(i,j) = v(i,j)*(1.0_RKind - dt*(2.0_RKind/(3.0_RKind*dx)*(v(i+1,j) - v(i-1,j)) &
                        - 1.0_RKind/(12_RKIND*dx)*(v(i+2,j) - v(i-2,j)) + 2.0_RKind*viscosity/dx**2.0_RKind + & 
                        2.0_RKind*viscosity/dy**2.0_RKind)) - dt*v(i,j+1)*(2.0_RKind/(3.0_RKind*dy)*u(i,j) - & 
                        viscosity/dy**2.0_RKind) + dt*v(i,j-1)*(2.0_RKind/(3.0_RKind*dy)*u(i,j) + viscosity/dy**2.0_RKind) &
                        + dt*v(i,j+2)*(1.0_RKind/(12_RKIND*dy)*u(i,j)) - dt*v(i,j-2)*(1.0_RKind/(12_RKIND*dy)*u(i,j)) &
                        + dt*viscosity/dx**2.0_RKind*(v(i+1,j) + v(i-1,j))
                        
                    END IF
                    
                END DO
            END DO
            
            

        END IF
        
    END SUBROUTINE speed_guess


    !Créée Etape 6, résout l'équation de poisson contenant les prédictions de vitesse, pour une itération temporelle
    SUBROUTINE compute_pressure()
    
    IMPLICIT NONE

        !remplissage de b à chaque itération temporelle, là où A est constant
        CALL fill_b()

        !resolution de l'equation de poisson avec méthode de Jacobi
        CALL jacobi_method()

    END SUBROUTINE compute_pressure
    
    
    SUBROUTINE adjust_speed()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        !Application des conditions limites sur u
        u(1,:) = 0_RKind
        u(n_x,:) = 0_RKind
        u(:,1) = 0_RKind
        u(:,n_y) = 1_RKind
        
        !Application des conditions limites sur v
        v(1,:) = 0_RKind
        v(n_x,:) = 0_RKind
        v(:,1) = 0_RKind
        v(:,n_y) = 0_RKind
        
        
        DO i = 2, n_x-1
            DO j = 2, n_y-1
                
                u(i, j) = u_temp(i, j) - ((p(i+1, j) - p(i-1, j))/(2.0_RKind*dx))*dt/density
                v(i, j) = v_temp(i, j) - ((p(i, j+1) - p(i, j-1))/(2.0_RKind*dy))*dt/density
                
            END DO
        END DO
        
    END SUBROUTINE adjust_speed
    
    
    !maj à l'étape 2, 2D
    !Réalise la boucle temporelle
    SUBROUTINE resolution_loop()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), PARAMETER :: NMax = 1000000
        INTEGER(KIND = IKind) :: i
        LOGICAL :: last_iteration
        
        i = 0
        last_iteration = .FALSE.
        
        !Allocation du tableau permettant de stocker les valeurs intermediaires de calcul
        ALLOCATE(u_temp(n_x, n_y))
        ALLOCATE(v_temp(n_x, n_y))
            
        !Boucle temporelle du calcul
        DO WHILE ((last_iteration .EQV. .FALSE.) .AND. (i < NMax))
            
            CALL compute_time_step()
            
            IF (t + dt > t_f) THEN
                !dt = t_f - t
                last_iteration = .TRUE.
                EXIT
            END IF
            
            PRINT*, 'dt = ', dt
            
            !appelle le calcul pour cette itération
            CALL speed_guess()
            
            CALL compute_pressure()
            
            CALL adjust_speed()
            
            i = i + 1
            t = t + dt
            
            !écrit dans un fichier toute les frames
            IF ((MOD(i, frame) == 0) .OR. (last_iteration .EQV. .TRUE.)) THEN
                CALL write_output_file(i)
            END IF
            
            !IF (i>=200 .AND. i <= 210) THEN
            !    CALL debug(i)
            !END IF
            
            PRINT*, 't = ', t
            PRINT*, '-------------------'
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
    

    DEALLOCATE(space_grid%x)
    DEALLOCATE(space_grid%y)
    DEALLOCATE(u)
    DEALLOCATE(v)
    DEALLOCATE(a)
    DEALLOCATE(a_loc)
    DEALLOCATE(b)
    DEALLOCATE(p)
    DEALLOCATE(p_vec)
    DEALLOCATE(p_vec_temp)
    DEALLOCATE(jacobi_r)
    
END PROGRAM main